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
void correct(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
		FluxProps *fluxprops, TimeProps *timeprops, void *EmTemp_in, double *forceint, double *forcebed,
		double *eroded, double *deposited) {
	Element *EmTemp = (Element *) EmTemp_in;
	double *dx = EmTemp->get_dx();
	double dtdx = dt / dx[0];
	double dtdy = dt / dx[1];
	double kactxy[DIMENSION];

	double tiny = GEOFLOW_TINY;
	int xp = EmTemp->get_positive_x_side();
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	int ivar, i, j, k;
	double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
	double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];

	Node* nxp = (Node*) NodeTable->lookup(EmTemp->getNode() + (xp + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxp[ivar] = nxp->flux[ivar];

	Node* nyp = (Node*) NodeTable->lookup(EmTemp->getNode() + (yp + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxyp[ivar] = nyp->flux[ivar];

	Node* nxm = (Node*) NodeTable->lookup(EmTemp->getNode() + (xm + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxm[ivar] = nxm->flux[ivar];

	Node* nym = (Node*) NodeTable->lookup(EmTemp->getNode() + (ym + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxym[ivar] = nym->flux[ivar];

	#ifdef DO_EROSION
	int do_erosion = 1;
#else
	int do_erosion=0;
#endif

#ifdef STOPCRIT_CHANGE_SOURCE
	int IF_STOPPED=EmTemp->get_stoppedflags();
#else
	int IF_STOPPED = !(!EmTemp->get_stoppedflags());
#endif

	double *state_vars = EmTemp->get_state_vars();
	double *prev_state_vars = EmTemp->get_prev_state_vars();
	double *d_state_vars = EmTemp->get_d_state_vars();
	double *gravity = EmTemp->get_gravity();
	double *d_gravity = EmTemp->get_d_gravity();
	double *zeta = EmTemp->get_zeta();
	double *curvature = EmTemp->get_curvature();
	double bedfrict = EmTemp->get_effect_bedfrict();
	double *Influx = EmTemp->get_influx();
	double terminal_vel = matprops_ptr->v_terminal;

	double lscale = matprops_ptr->LENGTH_SCALE;
	double gscale = matprops_ptr->GRAVITY_SCALE;
	double hscale = matprops_ptr->HEIGHT_SCALE;
	double velocity_scale = sqrt(lscale * gscale);
	double momentum_scale = hscale * velocity_scale;

	double Vel[DIMENSION];

	if (state_vars[1] > tiny) {
		for (i = 0; i < DIMENSION; i++)
			kactxy[i] = *(EmTemp->get_effect_kactxy() + i);

		Vel[0] = state_vars[2] / state_vars[1];
		Vel[1] = state_vars[3] / state_vars[1];

	} else {
		for (i = 0; i < DIMENSION; i++) {
			kactxy[i] = matprops_ptr->epsilon;
			Vel[i] = 0.;
		}

		bedfrict = matprops_ptr->bedfrict[EmTemp->get_material()];
	}

	double dragforce[2] = { 0., 0. };

	correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym, &tiny, &dtdx, &dtdy, &dt,
			d_state_vars, (d_state_vars + NUM_STATE_VARS), &(zeta[0]), &(zeta[1]), curvature,
			&(matprops_ptr->intfrict), &bedfrict, gravity, kactxy, d_gravity, &(matprops_ptr->frict_tiny),
			forceint, forcebed, dragforce, &do_erosion, eroded, Vel,
			&terminal_vel, &(matprops_ptr->epsilon), &IF_STOPPED, Influx);

	double ratio = 0;
	if (state_vars[0] > 5.5)
		ratio = -1;


	if (ratio > .1)
	//*EmTemp->pass_key()==3842346279 && *(EmTemp->pass_key()+1)==2368179492) //((isnan(state_vars[0]))||(state_vars[0]<0))
			{
		printf("the ratio is %10.5f:\n", ratio);
		double tempU[NUM_STATE_VARS];
		for (i = 0; i < NUM_STATE_VARS; i++)
			tempU[i] = prev_state_vars[i] - dtdx * (fluxxp[i] - fluxxm[i])
					- dtdy * (fluxyp[i] - fluxym[i]);
		printf("ElemKey: %u  ,   %u\n", *EmTemp->pass_key(), *(EmTemp->pass_key() + 1));
		printf("Kactxy = %10.5e, %10.5e\n", kactxy[0], kactxy[1]);
		printf("BedFrict: %10.5e: IntFrict: %10.5e\n", bedfrict, matprops_ptr->intfrict);
		printf("DtDx=%f , DtDy=%f ,Length scale=%f , Height scale=%f , Gravity scale=%f\n", dtdx, dtdy,
				lscale, hscale, gscale);

		printf("state_vars: ");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", state_vars[i]);
		printf("\n");

		printf("prev_state_vars: ");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", prev_state_vars[i]);
		printf("\n");

		printf("Ustore: ");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", tempU[i]);
		printf("\n");

		printf("fluxes: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("fluxxp:%10.5e, fluxxm:%10.5e, fluxyp:%10.5e, fluxym:%10.5e \n ", fluxxp[i], fluxxm[i],
					fluxyp[i], fluxym[i]);
		printf("stop from corrector\n");
		exit(1);
	}

	return;
}
