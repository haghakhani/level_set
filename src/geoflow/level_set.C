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
 * $Id: step.C 164 2013-08-07 06:15:30Z haghakhani $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define REFINE_THRESHOLD1 /*-.2*GEOFLOW_TINY*/ .2*GEOFLOW_TINY
#define REFINE_THRESHOLD2 /*-.05*GEOFLOW_TINY*/ .05*GEOFLOW_TINY
#define REFINE_THRESHOLD /*-.025*GEOFLOW_TINY*/.025*GEOFLOW_TINY

#include "../header/hpfem.h"
#include "../header/geoflow.h"

double signed_d_f(double p,double eps);
double mini(double a,double b);
double maxi(double a,double b);
double absol (double a);
void initialize_phi(Element *EmTemp, double* norm, double dt,double min,int * elem);
void record_of_phi(HashTable* NodeTable, HashTable* El_Table);
void initialization(HashTable* NodeTable, HashTable* El_Table,
		    double dt, MatProps* matprops_ptr,
		    FluxProps *fluxprops, TimeProps *timeprops)
{
  double norm;
  int i,n,elem=0;
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  Element* Em_Temp;
  double *dx,min,min_dx,max,max_delta,*phi_slope,thresh;

  int flag=0;

  record_of_phi(NodeTable, El_Table);
  min=1;  
  for(i=0; i<El_Table->get_no_of_buckets(); i++) 
    {
      currentPtr = *(buck+i);
      while(currentPtr) 
	{
	  Em_Temp=(Element*)(currentPtr->value);
	       	    
	  if(Em_Temp->get_adapted_flag()>0) {
	    dx=Em_Temp->get_dx();
	    min_dx = mini(dx[0],dx[1]);
	    if(min_dx<min) min=min_dx;
            
	    //Em_Temp->get_slopes(El_Table, NodeTable, matprops_ptr->gamma);
	  }
	  currentPtr=currentPtr->next; 
	}
    }
  int count=0;
thresh =.1*min*min*min;
 
  do{
    norm=0;
    elem=0;
    //elem++;
    max=0;
    for(i=0; i<El_Table->get_no_of_buckets(); i++) 
      {
  	currentPtr = *(buck+i);
  	while(currentPtr) 
  	  {
  	    Em_Temp=(Element*)(currentPtr->value);
	       	    
  	    if(Em_Temp->get_adapted_flag()>0) {
	      Em_Temp->calc_phi_slope(El_Table, NodeTable);
	      //	      phi_slope=Em_Temp->get_phi_slope();
	      //	      max_delta=maxi(maxi(phi_slope[2],phi_slope[3]),maxi(phi_slope[0],phi_slope[1]));
	      //	      if (max_delta>max) max=max_delta;

  	      //Em_Temp->get_slopes(El_Table, NodeTable, matprops_ptr->gamma);
  	    }
	    currentPtr=currentPtr->next; 
  	  }
      }
    // printf("the max is ...... %f  & min is ......%f\n",max,min);
    double time=.1*min;//absol(min/(5*max));

    for(i=0; i<El_Table->get_no_of_buckets(); i++) 
      {
	currentPtr = *(buck+i);
	while(currentPtr) 
	  {
	    Em_Temp=(Element*)(currentPtr->value);
            if(/*(timeprops->iter<2 && */ Em_Temp->get_adapted_flag()>0)/* ||
               ((abs(Em_Temp->get_adapted_flag())==BUFFER)  ||
               (Em_Temp->if_first_buffer_boundary(El_Table,GEOFLOW_TINY     )>0)||
               (Em_Temp->if_first_buffer_boundary(El_Table,REFINE_THRESHOLD1)>0)||
               (Em_Temp->if_first_buffer_boundary(El_Table,REFINE_THRESHOLD2)>0)||
               (Em_Temp->if_first_buffer_boundary(El_Table,REFINE_THRESHOLD) >0)||
	       (Em_Temp->if_next_buffer_boundary(El_Table,NodeTable,REFINE_THRESHOLD)>0) ||
	       (Em_Temp->if_next_buffer_boundary(El_Table,NodeTable,REFINE_THRESHOLD1)>0)||
	       (Em_Temp->if_next_buffer_boundary(El_Table,NodeTable,REFINE_THRESHOLD2)>0)||
               (Em_Temp->if_pile_boundary(El_Table,GEOFLOW_TINY)>0)||
               (Em_Temp->if_pile_boundary(El_Table,REFINE_THRESHOLD1)>0)||
               (Em_Temp->if_pile_boundary(El_Table,REFINE_THRESHOLD2)>0)||
               (Em_Temp->if_pile_boundary(El_Table,REFINE_THRESHOLD)>0))) */

	      {
                if(Em_Temp->if_pile_boundary(El_Table,GEOFLOW_TINY)>1) 
                   flag=1;  
		initialize_phi( Em_Temp,&norm,time,min,&elem);
		*(Em_Temp->get_state_vars()+5)=5;

		//		printf("norm =%f    \n", norm);
	      }
            flag=0;
	    currentPtr=currentPtr->next;      	    
	  }
      }
    
    printf("norm = %e.....dt=%e..... count=%d .....elem=%d..... thresh=%e...... min=%f\n", norm/elem,time,count,elem,thresh,min);
    count++;

    //if(timeprops->iter<2) 
    //thresh=1;
    //else
    //thresh=.5;

  }while((norm/elem)>thresh /* && count<1000*/);//&& count<21)||(sqrt(norm)>thresh && timeprops->iter<100 ));//&& (elem<1);//(sqrt(abs(norm))>.0100);
 
  return;
}

void initialize_phi( Element *EmTemp, double *norm, double dt,double min,int *elem){
  
  double *phi_slope=EmTemp->get_phi_slope();
  double *state_vars=EmTemp->get_state_vars();
  double *prev_state_vars=EmTemp->get_prev_state_vars();
  //double *d_state_vars=EmTemp->get_d_state_vars();
  double delta_p,delta_m;
  double flux_phi,a_p,a_m,b_p,b_m,c_p,c_m,d_p,d_m;  
  //  *norm=0;
  prev_state_vars[0]=state_vars[0];

  if (dabs(state_vars[0])<.5) *elem=*elem+1; 

  a_p=maxi( phi_slope[0] ,0.);
  a_m=mini( phi_slope[0] ,0.);
  b_p=maxi( phi_slope[1] ,0.);
  b_m=mini( phi_slope[1] ,0.);
  c_p=maxi( phi_slope[2] ,0.);
  c_m=mini( phi_slope[2] ,0.);
  d_p=maxi( phi_slope[3] ,0.);
  d_m=mini( phi_slope[3] ,0.);

  if (state_vars[4]>0)
   flux_phi=sqrt(maxi(a_p*a_p , b_m*b_m) + maxi(c_p*c_p , d_m*d_m))-1;
  else if (state_vars[4]<0)
   flux_phi=sqrt(maxi(a_m*a_m , b_p*b_p) + maxi(c_m*c_m , d_p*d_p))-1;
  else
   flux_phi=0;

   //if (state_vars[4]==0) printf("hello I am zero \n");

  state_vars[0]=prev_state_vars[0] - dt*signed_d_f(state_vars[4],min)*flux_phi;

  for(int j=0;j<4;j++)
     if(*(EmTemp->get_neigh_proc()+j) == INIT)   // this is a boundary!
         state_vars[0]=state_vars[4];//1;//prev_state_vars[0];

//  delta_p = sqrt( pow( maxi( phi_slope[0] ,0.) , 2) + pow( mini( phi_slope[1] ,0.) , 2)+
//		  pow( maxi( phi_slope[2] ,0.) , 2) + pow( mini( phi_slope[3] ,0.) , 2));


//  delta_m = sqrt( pow( maxi( phi_slope[1] ,0.) , 2) + pow( mini( phi_slope[0] ,0.) , 2)+
//		  pow( maxi( phi_slope[3] ,0.) , 2) + pow( mini( phi_slope[2] ,0.) , 2));


//  state_vars[0]=prev_state_vars[0] - dt* (maxi(signed_d_f(state_vars[4],min),0.) * delta_p + 
//					  mini(signed_d_f(state_vars[4],min),0.) * delta_m ) + dt * signed_d_f(state_vars[4],min);

  // state_vars[0]=prev_state_vars[0] + dt*( c_sgn(state_vars[4]) -
  // 					  (maxi(c_sgn(state_vars[4]),0.)*delta_p + mini(c_sgn(state_vars[4]),0.)*delta_m ));

  //if((d_state_vars[0]!=0||d_state_vars[6]!=0) && state_vars[4]!=0)
  //printf("this is the place dx=%f dy=%f \n",d_state_vars[0],d_state_vars[6]);
  //  state_vars[0]=prev_state_vars[0] + dt*c_sgn(state_vars[4])*(1-sqrt(pow(d_state_vars[0],2)+pow(d_state_vars[6],2)));

 
  *norm += dabs(state_vars[0]-prev_state_vars[0]);
//    if (dabs(state_vars[0]-prev_state_vars[0])>0) printf("this is the difference .........%f\n", dabs(state_vars[0]-prev_state_vars[0]));
  return;
}

void record_of_phi(HashTable* NodeTable, HashTable* El_Table){

  int i;
  HashEntryPtr      *buck = El_Table->getbucketptr();

  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  *(Curr_El->get_state_vars()+4)= *(Curr_El->get_state_vars());
	  //if (*(Curr_El->get_state_vars()+4)!=*(Curr_El->get_state_vars())) exit(1);
	}
	currentPtr=currentPtr->next;
      }
    }
  return;
}
double mini(double a,double b){
  if (a<b)
    return(a);
  else
    return(b);
}

double maxi(double a,double b){
  if (a>b)
    return(a);
  else
    return(b);
}


double absol (double a){
  if (a<0) return (-a);
  else  return (a);
}

double signed_d_f(double p,double eps){
  return (p/sqrt(p*p+eps*eps));
  //return (.5*(tanh(p/(2*eps))+1));
  //return (tanh(p/(2*eps)));
}
