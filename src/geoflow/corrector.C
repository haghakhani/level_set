/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: get_coef_and_eigen.C,v 1.4 2004/08/11 15:58:46 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define DO_EROSION

#include "../header/hpfem.h"
#include "../header/geoflow.h"
void correct(HashTable* NodeTable, HashTable* El_Table,
	     double dt, MatProps* matprops_ptr, 
	     FluxProps *fluxprops, TimeProps *timeprops,
	     void *EmTemp_in,double *forceint, double *forcebed, 
	     double *eroded, double *deposited)
{
  Element *EmTemp=(Element *) EmTemp_in;
  double *dx=EmTemp->get_dx();
  double dtdx = dt/dx[0];
  double dtdy = dt/dx[1];
  double kactxy[DIMENSION]; 

  double tiny = GEOFLOW_TINY;
  int xp=EmTemp->get_positive_x_side();
  int yp=(xp+1)%4, xm=(xp+2)%4, ym=(xp+3)%4; 

  int ivar,i, j, k;
  double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS]; 
  double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];

  Node* nxp = (Node*) NodeTable->lookup(EmTemp->getNode()+(xp+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxxp[ivar] = nxp->flux[ivar];

  Node* nyp = (Node*) NodeTable->lookup(EmTemp->getNode()+(yp+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxyp[ivar] = nyp->flux[ivar];

  Node* nxm = (Node*) NodeTable->lookup(EmTemp->getNode()+(xm+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxxm[ivar] = nxm->flux[ivar];

  Node* nym = (Node*) NodeTable->lookup(EmTemp->getNode()+(ym+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxym[ivar] = nym->flux[ivar];

  //if(fluxxp[0]!=0||fluxxm[0]!=0||fluxyp[0]!=0||fluxym[0]!=0)
  //printf("fluxxp=...%f, fluxxm=...%f, fluxyp = ...%f,fluxym =... %f \n",fluxxp[0],fluxxm[0],fluxyp[0],fluxym[0]);

#ifdef DO_EROSION
  int do_erosion=1;
#else
  int do_erosion=0;
#endif

#ifdef STOPCRIT_CHANGE_SOURCE
  int IF_STOPPED=EmTemp->get_stoppedflags();
#else
  int IF_STOPPED=!(!EmTemp->get_stoppedflags());
#endif

  double *state_vars=EmTemp->get_state_vars();
  double *prev_state_vars=EmTemp->get_prev_state_vars();
  double *d_state_vars=EmTemp->get_d_state_vars();
  double *gravity=EmTemp->get_gravity();
  double *d_gravity=EmTemp->get_d_gravity();
  double *lap_phi=EmTemp->get_lap_phi();
  double *zeta=EmTemp->get_zeta();
  double *curvature=EmTemp->get_curvature();
  double bedfrict=EmTemp->get_effect_bedfrict();
  double *Influx=EmTemp->get_influx();
  double solid_den=matprops_ptr->den_solid;
  double fluid_den=matprops_ptr->den_fluid;
  double terminal_vel=matprops_ptr->v_terminal;
  double navslip_coef=matprops_ptr->navslip_coef;

  double lscale=matprops_ptr->LENGTH_SCALE;
  double gscale=matprops_ptr->GRAVITY_SCALE;
  double hscale=matprops_ptr->HEIGHT_SCALE;
  double velocity_scale = sqrt(lscale *gscale);
  double momentum_scale = hscale * velocity_scale;

  double Vfluid[DIMENSION], Vsolid[DIMENSION];
  // double volf;

  if ( state_vars[0] < 0 && state_vars[1]>0)
    {
      for (i=0; i<DIMENSION; i++)
	kactxy[i]=*(EmTemp->get_effect_kactxy()+i);

      // fluid velocities
      Vfluid[0]=state_vars[4]/state_vars[1];
      Vfluid[1]=state_vars[5]/state_vars[1];
      Vsolid[0]=state_vars[2]/state_vars[1];
      Vsolid[1]=state_vars[3]/state_vars[1];

      // volume fractions
      //volf = state_vars[1]/state_vars[0];
    }
  else
    {
      for (i=0; i<DIMENSION; i++)
	{
	  kactxy[i]=matprops_ptr->epsilon;
	  Vfluid[i]=0.;
	  Vsolid[i]=0.;
	}
      //volf=1.;
      bedfrict=matprops_ptr->bedfrict[EmTemp->get_material()];
    }

  double V_avg[DIMENSION];
  V_avg[0] = Vsolid[0];//*volf + Vfluid[0]*(1.-volf);
  V_avg[1] = Vsolid[1];//*volf + Vfluid[1]*(1.-volf);
  EmTemp->convect_dryline(V_avg,dt); //this is necessary

  if (state_vars[0] < GEOFLOW_SHORT && state_vars[1]!=0)
    navslip_coef *= state_vars[1];

  double dragforce[2] = {0., 0.};
  int iter=timeprops->iter;
  int keys1=*(EmTemp->pass_key());	
  int keys2=*(EmTemp->pass_key()+1);

  double a,b,c,d;
  a=b=c=d=0;
  a= fluxxp[1]-fluxxm[1];
  b= fluxxp[0]+fluxxm[4];
  c= fluxyp[1]-fluxym[1];
  d= fluxyp[0]+fluxym[4];
//  if(a!=0||b!=0||c!=0||d!=0){
//    printf("\n");
//    printf("The element keys: K1= %d    K2=%d\n",keys1,keys2);
//    printf("the flux result in X: HLL =.... %f , ROE=....%f \n",fluxxp[1]-fluxxm[1],fluxxp[0]+fluxxm[4]); 
//    printf("HLL: fluxxp[1]=....%f    fluxxm[1]=....%f\n", fluxxp[1],fluxxm[1]);
//    printf("ROE: fluxxp[0]=....%f    fluxxm[4]=....%f\n", fluxxp[0],fluxxm[4]);
//    printf("\n"); 
//    printf("the flux result in Y: HLL =.... %f , ROE=....%f \n",fluxyp[1]-fluxym[1],fluxyp[0]+fluxym[4]);
//    printf("HLL: fluxyp[1]=....%f    fluxym[1]=....%f\n", fluxyp[1],fluxym[1]);
//    printf("ROE: fluxyp[0]=....%f    fluxym[4]=....%f\n", fluxyp[0],fluxym[4]);
//    printf(".........................................................................................\n"); 
//  }

  //if (keys1==3563192320 && keys2==0)
  //printf("hello i found you\n");
  correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym,
	   &tiny, &dtdx, &dtdy, &dt, d_state_vars, (d_state_vars+NUM_STATE_VARS), 
	   lap_phi ,&(zeta[0]), &(zeta[1]), curvature,
	   &(matprops_ptr->intfrict), &bedfrict,
	   gravity, kactxy,d_gravity ,&(matprops_ptr->frict_tiny),
	   forceint, forcebed, dragforce, &do_erosion, eroded, Vsolid, Vfluid,
	   &solid_den, &fluid_den, &terminal_vel,
	   &(matprops_ptr->epsilon), &IF_STOPPED, Influx, &navslip_coef);

  //if (state_vars[2]>100) {print_elem_data(EmTemp,matprops_ptr, fluxprops,timeprops); exit(1);}
  EmTemp->put_drag(dragforce);
  *forceint*=dx[0]*dx[1];
  *forcebed*=dx[0]*dx[1];
  *eroded*=dx[0]*dx[1];

  // bool debug=false;
  // if (((*(EmTemp->get_coord()))*(matprops_ptr)->LENGTH_SCALE > 645000 ) && 
  //     ((*(EmTemp->get_coord()))*(matprops_ptr)->LENGTH_SCALE < 645100 ) &&
  //     ((*(EmTemp->get_coord()+1))*(matprops_ptr)->LENGTH_SCALE > 2164800 )&&
  //     ((*(EmTemp->get_coord()+1))*(matprops_ptr)->LENGTH_SCALE > 2164850 )&&
  //     (state_vars[0]<1))
  //   debug=true;
  // if ( debug )
  //   {
  //     double tempg[NUM_STATE_VARS];
  //     for (i=0; i<NUM_STATE_VARS; i++)
  // 	tempg[i] = prev_state_vars[i]-
  // 	  dtdx*(fluxxp[i]-fluxxm[i])-
  // 	  dtdy*(fluxyp[i]-fluxym[i]);
  //     printf("ElemKey: %u\n", *EmTemp->pass_key());
  //     printf("Kactxy = %10.5e, %10.5e\n", kactxy[0], kactxy[1]);
  //     printf("BedFrict: %10.5e: IntFrict: %10.5e\n", bedfrict, matprops_ptr->intfrict);
  //     printf("DtDx=%f , DtDy=%f ,Length scale=%f , Height scale=%f , Gravity scale=%f\n",dtdx,dtdy, lscale , hscale , gscale);

  //     printf("state_vars: \n");
  //     for (i=0; i<NUM_STATE_VARS; i++)
  // 	printf("%10.5e, ", state_vars[i]);
  //     printf("\n");

  //     printf("prev_state_vars: \n");
  //     for (i=0; i<NUM_STATE_VARS; i++)
  // 	printf("%10.5e, ", prev_state_vars[i]);
  //     printf("\n");

  //     printf("Ustore: \n");
  //     for (i=0; i<NUM_STATE_VARS; i++)
  // 	printf("%10.5e, ", tempg[i]);
  //     printf("\n");

  //     printf("fluxes: \n");
  //     for (i=0; i<NUM_STATE_VARS; i++)
  // 	printf("fluxxp:%10.5e, fluxxm:%10.5e, fluxyp:%10.5e, fluxym:%10.5e \n ", fluxxp[i],fluxxm[i],fluxyp[i],fluxym[i]);
  //     printf("\n");
  //     double one=16*2*(2*pow(prev_state_vars[0],2)-3*prev_state_vars[0]+1);
  //     double two=1e-6*(lap_phi[0]+lap_phi[1]);

  // 	// double source =16*((4*pow(prev_state_vars[0],3)-6*pow(prev_state_vars[0],2)+2*prev_state_vars[0])* (lap_phi[0]+lap_phi[1])
  // 	// 			 +(12*pow(prev_state_vars[0],2)-12*prev_state_vars[0]+2)*(pow(d_state_vars[0],2)+pow(d_state_vars[6],2))
  // 	// 			 +prev_state_vars[0]*(d_state_vars[2]+d_state_vars[8]));      
  // 	printf("source term eq0: \n");
  //     printf("source one:%10.5e, source two:%10.5e, lap_phi[0]:%10.5e, lap_phi[1]:%10.5e \n ", one,two,lap_phi[0],lap_phi[1]);

  //     // exit(1);
  //   }
  double ratio=0;
  if (state_vars[0]>5.5){
    ratio = -1;//dabs((state_vars[0]-state_vars[1])/state_vars[1]);
    //printf("the ratio is %f\n",ratio);
  }

  if ( ratio > .1) 
    //*EmTemp->pass_key()==3842346279 && *(EmTemp->pass_key()+1)==2368179492) //((isnan(state_vars[0]))||(state_vars[0]<0))
    {
      printf ("the ratio is %10.5f:\n",ratio);
      double tempU[NUM_STATE_VARS];
      for (i=0; i<NUM_STATE_VARS; i++)
	tempU[i] = prev_state_vars[i]-
	  dtdx*(fluxxp[i]-fluxxm[i])-
	  dtdy*(fluxyp[i]-fluxym[i]);
      printf("ElemKey: %u  ,   %u\n", *EmTemp->pass_key(),*(EmTemp->pass_key()+1));
      printf("Kactxy = %10.5e, %10.5e\n", kactxy[0], kactxy[1]);
      printf("BedFrict: %10.5e: IntFrict: %10.5e\n", bedfrict, matprops_ptr->intfrict);
      printf("DtDx=%f , DtDy=%f ,Length scale=%f , Height scale=%f , Gravity scale=%f\n",dtdx,dtdy, lscale , hscale , gscale);

      printf("state_vars: ");
      for (i=0; i<NUM_STATE_VARS; i++)
	printf("%10.5e, ", state_vars[i]);
      printf("\n");

      printf("prev_state_vars: ");
      for (i=0; i<NUM_STATE_VARS; i++)
	printf("%10.5e, ", prev_state_vars[i]);
      printf("\n");

      printf("Ustore: ");
      for (i=0; i<NUM_STATE_VARS; i++)
	printf("%10.5e, ", tempU[i]);
      printf("\n");

      printf("fluxes: \n");
      for (i=0; i<NUM_STATE_VARS; i++)
	printf("fluxxp:%10.5e, fluxxm:%10.5e, fluxyp:%10.5e, fluxym:%10.5e \n ", fluxxp[i],fluxxm[i],fluxyp[i],fluxym[i]);
      printf("stop from corrector\n");
      		exit(1);
    }

  if(EmTemp->get_stoppedflags()==2) 
    *deposited=state_vars[1]*dx[0]*dx[1];
  else
    *deposited=0.0;

  if(EmTemp->get_stoppedflags()) 
    *eroded=0.0;

  EmTemp->calc_shortspeed(1.0/dt);
  return;
}
