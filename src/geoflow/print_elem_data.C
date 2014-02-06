#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/geoflow.h"

int gogol(HashTable* El_Table ,int b){
double a=El_Table->get_invdxrange();
b=a+b;
return(b);}


void print_elem_data(Element *EmTemp,	MatProps* matprops_ptr, FluxProps* fluxprops_ptr,TimeProps* timeprops_ptr)  
{

	int xp,yp,xm,ym;

	int ivar,i, j, k;
	double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
	double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];

	double *state_vars, *prev_state_vars, *d_state_vars, *gravity, *d_gravity, *dx_ptr, *kactxy;
	double *lap_phi, *zeta, *curvature, bedfrict, *Influx, solid_den;
	double fluid_den, terminal_vel, navslip_coef, lscale, gscale, velocity_scale, momentum_scale, hscale, dt;
       
        timeprops_ptr->incrtime(&dt); 


        double *dx=EmTemp->get_dx();
        double dtdx = dt/dx[0];
        double dtdy = dt/dx[1];


	xp=EmTemp->get_positive_x_side();
	yp=(xp+1)%4; xm=(xp+2)%4; ym=(xp+3)%4;


	state_vars=EmTemp->get_state_vars();
	prev_state_vars=EmTemp->get_prev_state_vars();
	d_state_vars=EmTemp->get_d_state_vars();
	gravity=EmTemp->get_gravity();
	d_gravity=EmTemp->get_d_gravity();
	lap_phi=EmTemp->get_lap_phi();
	zeta=EmTemp->get_zeta();
	curvature=EmTemp->get_curvature();
	bedfrict=EmTemp->get_effect_bedfrict();
	Influx=EmTemp->get_influx();
	solid_den=matprops_ptr->den_solid;
	fluid_den=matprops_ptr->den_fluid;
	terminal_vel=matprops_ptr->v_terminal;
	navslip_coef=matprops_ptr->navslip_coef;

	lscale=matprops_ptr->LENGTH_SCALE;
	gscale=matprops_ptr->GRAVITY_SCALE;
	hscale=matprops_ptr->HEIGHT_SCALE;
	velocity_scale = sqrt(lscale *gscale);
	momentum_scale = hscale * velocity_scale;
//       if ( state_vars[0] > GEOFLOW_TINY)
//                for (i=0; i<DIMENSION; i++)
//                        kactxy[i]=*(EmTemp->get_effect_kactxy()+i);
//        else
//                for (i=0; i<DIMENSION; i++)
//                        kactxy[i]=1;//matprops_ptr->epsilon;




	//printf("Please attention the phi is .......................... %f\n", *(EmTemp->get_state_vars()));
	dx_ptr = EmTemp->get_dx();

	double tempU[NUM_STATE_VARS];
for (i=0; i<NUM_STATE_VARS; i++)
	tempU[i] = prev_state_vars[i]-
	dtdx*(fluxxp[i]-fluxxm[i])-
	dtdy*(fluxyp[i]-fluxym[i]);
	printf("ElemKey: %u      %u\n", *EmTemp->pass_key(),*(EmTemp->pass_key()+1));
//	printf("Kactxy = %10.5e, %10.5e\n", kactxy[0], kactxy[1]);
	printf("BedFrict: %10.5e: IntFrict: %10.5e\n", bedfrict, matprops_ptr->intfrict);
	printf("DtDx=%f , DtDy=%f ,Length scale=%f , Height scale=%f , Gravity scale=%f\n"
			,1,1, lscale , hscale , gscale);

	printf("state_vars: \n");
	for (i=0; i<NUM_STATE_VARS; i++)
		printf("%10.5e, ", state_vars[i]);
	printf("\n");

        printf("d_state_vars: \n");
        for (i=0; i<2*NUM_STATE_VARS; i++)
                {printf("%10.5e, ", d_state_vars[i]);if(i==NUM_STATE_VARS-1) printf("\n");}
printf("\n");

	printf("prev_state_vars: \n");
	for (i=0; i<NUM_STATE_VARS; i++)
		printf("%10.5e, ", prev_state_vars[i]);
	printf("\n");

	printf("Ustore: \n");
	for (i=0; i<NUM_STATE_VARS; i++)
	printf("%10.5e, ", tempU[i]);
	printf("\n");

	printf("fluxes: \n");
	for (i=0; i<NUM_STATE_VARS; i++)
		printf("fluxxp:%10.5e, fluxxm:%10.5e,fluxyp:%10.5e, fluxym:%10.5e \n "
				, fluxxp[i],fluxxm[i],fluxyp[i],fluxym[i]);
        printf("Gravity: \n");
        for (i=0; i<3; i++)
        printf("%10.5e, ", gravity[i]);
        printf("\n");

        printf("d_gravity: \n");
        for (i=0; i<NUM_STATE_VARS; i++)
        printf("%10.5e, ", d_gravity[i]);
        printf("\n");

	
printf("\n");


	return;
}
