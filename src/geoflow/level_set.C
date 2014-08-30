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

#define REFINE_THRESHOLD1 .2*GEOFLOW_TINY
#define REFINE_THRESHOLD2 .05*GEOFLOW_TINY
#define REFINE_THRESHOLD .025*GEOFLOW_TINY

#include "../header/hpfem.h"
#include "../header/geoflow.h"

double signed_d_f(double p,double eps);
double mini(double a,double b);
double maxi(double a,double b);
double absol (double a);
void initialize_phi(HashTable* El_Table,Element *EmTemp, double* norm, double dt,double min,int* elem);
void record_of_phi(HashTable* NodeTable, HashTable* El_Table);
int num_nonzero_elem(HashTable *El_Table);
double minidxdy(double x,double y, int* falg);

struct Min_rank {
  double min;
  int rank;
};

//struct Boundry_Point{
//  double xpos;
//  double ypos;
//};

void narrow_bound_layers(HashTable *ElemTable, Element* elem, int num_layer){

  while(num_layer>0){
    // cout<<"num lyer: "<<num_layer<<endl;
    for(int ineigh=0;ineigh<8;ineigh++){
      // cout<<"neighbor:  "<<ineigh<<endl;
      if(*(elem->get_neigh_proc()+ineigh)>=0) //don't check outside map boundary or duplicate neighbor
      {
        Element* ElemNeigh=(Element*) ElemTable->lookup(elem->get_neighbors()+ineigh*KEYLENGTH);
        assert(ElemNeigh);
        *(elem->get_state_vars()+5)=(double) num_layer;
        // cout<<"elem keys are:   "<<*(ElemNeigh->pass_key())<<"  and  "<<*(ElemNeigh->pass_key()+1)<<endl;
        narrow_bound_layers(ElemTable, ElemNeigh,num_layer-1);
      }
    }
    return;
  }
  return;
}

void narrow_bound_marker(HashTable *ElemTable, Element* elem, int num_layer){
  //*(elem->get_state_vars()+5)=0.;//we use here the state_vars[5] as a flag
  if (elem->if_phase_baundary(ElemTable)){
    *(elem->get_state_vars()+5)=(double)num_layer;
    narrow_bound_layers(ElemTable,elem , num_layer-1);
  }
  return;
}


void create_narrow_bound(HashTable* El_Table){
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  int num_layer=4, num=0;

  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
  {
    currentPtr = *(buck+i);
    while(currentPtr)
    {
      Element* Em_Temp=(Element*)(currentPtr->value);

      if(Em_Temp->get_adapted_flag()>0 /* && num<1*/ ){
        *(Em_Temp->get_state_vars()+5)= 0.;
        num++; 
      }

      currentPtr=currentPtr->next;
    }
  }

  cout<<"number of elements are:  "<<num<<endl;
  num=0;

  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
  {
    currentPtr = *(buck+i);
    while(currentPtr)
    {
      Element* Em_Temp=(Element*)(currentPtr->value);

      if(Em_Temp->get_adapted_flag()>0 /* && num<1*/ ){
        narrow_bound_marker(El_Table,Em_Temp,num_layer);
        //if (Em_Temp->if_phase_baundary(El_Table)) cout << "I am working"<<endl;
        //cout<<"now performing:  "<<num++<<endl; 

      }

      currentPtr=currentPtr->next;
    }
  }


  return;
}

void initialization(HashTable* NodeTable, HashTable* El_Table,
    double dt, MatProps* matprops_ptr,
    FluxProps *fluxprops, TimeProps *timeprops, OutLine* outline_ptr,
    int nump, int rank)
{
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  //vector<Bounday_Point> boundary_point;
  //Boundry_Point newBoundaryPoint;

  double *coord;
  int locNumBoundPoint=0;
  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
  {
    currentPtr = *(buck+i);
    while(currentPtr)
    {
      Element* Em_Temp=(Element*)(currentPtr->value);

      if(Em_Temp->get_adapted_flag()>0  && *(Em_Temp->get_state_vars())==0. ) 
        locNumBoundPoint++;

      currentPtr=currentPtr->next;
    }
  }

  double **boundaryPoints=CAllocD2(locNumBoundPoint,2);


  int index=0;
  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
  {
    currentPtr = *(buck+i);
    while(currentPtr)
    {
      Element* Em_Temp=(Element*)(currentPtr->value);
      if(Em_Temp->get_adapted_flag()>0  && *(Em_Temp->get_state_vars())==0. ) {
        //coord=Em_Temp->get_coord();
        boundaryPoints[index][0]=*(Em_Temp->get_coord());
        boundaryPoints[index][1]=*(Em_Temp->get_coord()+1);
        index++;
      }
      currentPtr=currentPtr->next;
    }
  }

  double xrange=(outline_ptr->xminmax[1]-outline_ptr->xminmax[0])/*/matprops_ptr->LENGTH_SCALE*/;
  double yrange=(outline_ptr->yminmax[1]-outline_ptr->yminmax[0])/*/matprops_ptr->LENGTH_SCALE*/;
  double mindist=1000000;//sqrt(xrange*xrange+yrange*yrange);


  for(int i=0; i<El_Table->get_no_of_buckets(); i++)
  {
    currentPtr = *(buck+i);
    while(currentPtr)
    {
      Element* Em_Temp=(Element*)(currentPtr->value);
      if(Em_Temp->get_adapted_flag()>0 && *(Em_Temp->get_state_vars())!=0.) {

        double coef = *(Em_Temp->get_state_vars())>0. ? 1.0: -1.0;//*(Em_Temp->get_state_vars());//*(Em_Temp->get_state_vars())>0. ? 1.0: -1.0;
        double dist=mindist;

        for (int idist=0; idist<locNumBoundPoint; idist++){
          coord=Em_Temp->get_coord();
          double tempdist=sqrt(pow(*(Em_Temp->get_coord())-boundaryPoints[idist][0],2)+pow(*(Em_Temp->get_coord()+1)-boundaryPoints[idist][1],2));
          dist = tempdist < dist ? tempdist : dist;
        }
        dist*=coef;
        *(Em_Temp->get_state_vars())=dist;

      }


      currentPtr=currentPtr->next;
    }
  }

  return;
}
void reinitialization(HashTable* NodeTable, HashTable* El_Table,
    double dt, MatProps* matprops_ptr,
    FluxProps *fluxprops, TimeProps *timeprops, OutLine* outline_ptr,
    int nump, int rank)
{
  double norm;
  int i,n,elem=0;
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  Element* Em_Temp;
  double *dx,min,min_dx,max,max_delta,*phi_slope,thresh=.0001;
  create_narrow_bound(El_Table);
  //cout<<"attention xmin: "<<outline_ptr->xminmax[0]<<" xmax: "<<outline_ptr->xminmax[1]<<" ymin: "<<outline_ptr->yminmax[0]<<" ymax: "<<outline_ptr->yminmax[1]<<endl;
  record_of_phi(NodeTable, El_Table);
  move_data(nump, rank, El_Table, NodeTable,timeprops);

  min=10000;  
  int flagxy=2,flagxymin=2;

  for(i=0; i<El_Table->get_no_of_buckets(); i++) 
  {
    currentPtr = *(buck+i);
    while(currentPtr) 
    {
      Em_Temp=(Element*)(currentPtr->value);

      if(Em_Temp->get_adapted_flag()>0) {
        dx=Em_Temp->get_dx();
        min_dx = minidxdy(dx[0],dx[1],&flagxy);
        if(min_dx<min) {min=min_dx; flagxymin=flagxy;}
      }
      currentPtr=currentPtr->next; 
    }
  }
  Min_rank min_rank,tot_min;
  min_rank.min=min;
  min_rank.rank=rank;

  MPI_Barrier(MPI_COMM_WORLD);
  //we need also to broadcast the flagxymin 
  //Min_rank tot_min;

  MPI_Allreduce(&min_rank,&tot_min,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);

  MPI_Bcast(&flagxymin, 1, MPI_INT,tot_min.rank,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  min=tot_min.min;


  //  MPI_Reduce(&min,&min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  int count=0;

  //if (rank==0){

  double xrange=(outline_ptr->xminmax[1]-outline_ptr->xminmax[0])/*/matprops_ptr->LENGTH_SCALE*/;
  double yrange=(outline_ptr->yminmax[1]-outline_ptr->yminmax[0])/*/matprops_ptr->LENGTH_SCALE*/;

  if (flagxymin==0)
    min/=xrange;
  else if (flagxymin==1)
    min/=yrange;
  else{
    printf("there is an error in levelset.C min computation\n");
    return;
  }
  //}
  //MPI_Bcast(&min, 1, MPI_DOUBLE,0,MPI_COMM_WORLD );
  //MPI_Barrier(MPI_COMM_WORLD);

  double time_inc=.5*min;
  thresh=min*min;
  int flagi=0;
  double total_norm=0.0;
  do{
    norm=0;

    for(i=0; i<El_Table->get_no_of_buckets(); i++) 
    {
      currentPtr = *(buck+i);
      while(currentPtr) 
      {
        Em_Temp=(Element*)(currentPtr->value);
        if(Em_Temp->get_adapted_flag()>0&& (*(Em_Temp->get_state_vars()+5)>0 && *(Em_Temp->get_state_vars()+5)<5 )) {
          Em_Temp->calc_phi_slope(El_Table, NodeTable);
        }
        currentPtr=currentPtr->next; 
      }
    }
    move_data(nump, rank, El_Table, NodeTable,timeprops); 

    for(i=0; i<El_Table->get_no_of_buckets(); i++) 
    {
      currentPtr = *(buck+i);
      while(currentPtr) 
      {
        Em_Temp=(Element*)(currentPtr->value);
        if(Em_Temp->get_adapted_flag()>0&& (*(Em_Temp->get_state_vars()+5)>0 && *(Em_Temp->get_state_vars()+5)<5 ))
          //	       (timeprops->iter==1 && Em_Temp->get_adapted_flag()>0 )|| dabs(*(Em_Temp->get_state_vars()+4))<=.2)
          /*(timeprops->iter<100 &&*/ /*Em_Temp->get_adapted_flag()>0)*/ /*|| 
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
                                                                           (Em_Temp->if_pile_boundary(El_Table,REFINE_THRESHOLD)>0) ))*/
        {

          initialize_phi( El_Table,Em_Temp,&norm,time_inc,/*time_inc*/10*min,&elem);
          //	*(Em_Temp->get_state_vars()+5)=5;

          //		printf("norm =%f    \n", norm);
        }
        currentPtr=currentPtr->next;      	    
      }
    }

    MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    int tot_elem=0;

    if (timeprops->iter>1 )  MPI_Allreduce(&elem,&tot_elem,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    else {
      int nonzelem=num_nonzero_elem(El_Table);

      MPI_Allreduce(&nonzelem,&tot_elem,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    }

    total_norm/=tot_elem;


    if(rank==0) printf("norm=  %e  dt=  %e   count=   %d    thresh= %e  min =  %e \n", total_norm,time_inc,count,thresh,tot_min.min);
    count++;
    if((timeprops->iter>1 && count > 100) || (timeprops->iter==1 && count > 1000) ) 
      //  flagi=1;

      move_data(nump, rank, El_Table, NodeTable,timeprops);//at begining and the end of move_data there is a Barrier, so here I do not need to call MPI_Barrier

  }while((total_norm>thresh) && flagi < 1 );//&& count<21)||(sqrt(norm)>thresh && timeprops->iter<100 ));//&& (elem<1);//(sqrt(abs(norm))>.0100);

  return;
}

void initialize_phi(HashTable *El_Table, Element *EmTemp, double *norm, double dt,double min,int* elem){

  double *phi_slope=EmTemp->get_phi_slope();
  double *state_vars=EmTemp->get_state_vars();
  double *prev_state_vars=EmTemp->get_prev_state_vars();
  //double *d_state_vars=EmTemp->get_d_state_vars();
  double delta_p,delta_m;
  double flux_phi,a_p,a_m,b_p,b_m,c_p,c_m,d_p,d_m; 
  Element  *rev_elem;

  prev_state_vars[0]=state_vars[0];

  if (dabs(state_vars[0])<.5) *elem +=1;

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

  state_vars[0]=prev_state_vars[0] - dt*signed_d_f(state_vars[4],min)*flux_phi;


  //  int flag=1;
  //  for(int j=0;j<4;j++)
  //    if(*(EmTemp->get_neigh_proc()+j) == INIT) {  // this is a boundary!
  //     int rev_side=(j+2)%4;
  //     rev_elem = (Element*)(El_Table->lookup(EmTemp->get_neighbors()+rev_side*KEYLENGTH));//      EmTemp->get_neighbors();
  //     state_vars[0]=*(rev_elem->get_state_vars());//1;//prev_state_vars[0];
  //      flag=0;
  //    }
  //  if (flag)
  *norm += dabs(state_vars[0]-prev_state_vars[0]);
  state_vars[5]=dabs(state_vars[0]-prev_state_vars[0]);
  //  if ((state_vars[0]-prev_state_vars[0])>0) printf("this is the difference .........%e\n", (state_vars[0]-prev_state_vars[0]));
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

double minidxdy(double x,double y, int* flag){
  if (x<y){
    *flag=0;
    return(x);}
  else{
    *flag=1;
    return(y);}
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

//(*(EmTemp->get_dx()+0)<*(EmTemp->get_dx()+1))?*(EmTemp->get_dx()+0):*(EmTemp->get_dx()+1))
int num_nonzero_elem(HashTable *El_Table){
  int               num=0,myid,i;
  HashEntryPtr      currentPtr;
  Element           *Curr_El;
  HashEntryPtr      *buck = El_Table->getbucketptr();
  for(i=0;i<El_Table->get_no_of_buckets();i++)
    if(*(buck+i)){
      currentPtr = *(buck+i);
      while(currentPtr){
        Curr_El=(Element*)(currentPtr->value);
        if(Curr_El->get_adapted_flag()>0){
          num++;
        }
        currentPtr=currentPtr->next;
      }
    }

  return (num);
}
