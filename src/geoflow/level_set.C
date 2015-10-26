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

#include "../header/hpfem.h"
#include <set>
#include <vector>

typedef double (*Pt2Path)(void* path, double x, double y);
typedef void (*Pt2Grad)(void* path, double x, double y, double* grad);

extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb,
    int *info);

class Edge {
	// this  class holds information of an edge consisting of two elements
	// and also the neighbor number of the second element in the first element
public:
	Edge() {
		elem[0] = NULL;
		elem[1] = NULL;
		neigh_num = 9;
	}
	;

	Edge(Element* elem_l, Element* elem_r, int neigh_n) {
		elem[0] = elem_l;
		elem[1] = elem_r;
		neigh_num = neigh_n;
	}
	;

	bool operator<(const Edge& redge) const {
		if (elem[0] < redge.elem[0] || (elem[0] == redge.elem[0] && elem[1] < redge.elem[1]))
			return true;

		return false;
	}
	;
	Element *elem[2];
	int neigh_num;
};

class Quad {
public:
	Quad() {
		elem[0] = NULL;
		elem[1] = NULL;
		elem[2] = NULL;
		elem[3] = NULL;

	}
	;

	Quad(Element* elem_1, Element* elem_2, Element* elem_3, Element* elem_4) {
		elem[0] = elem_1;
		elem[1] = elem_2;
		elem[2] = elem_3;
		elem[3] = elem_4;
	}
	;

	bool operator<(const Quad& rrect) const {
		if (elem[0] < rrect.elem[0] || (elem[0] == rrect.elem[0] && elem[1] < rrect.elem[1])
		    || (elem[0] == rrect.elem[0] && elem[1] == rrect.elem[1] && elem[2] < rrect.elem[2])
		    || (elem[0] == rrect.elem[0] && elem[1] == rrect.elem[1] && elem[2] == rrect.elem[2]
		        && elem[3] < rrect.elem[3]))
			return true;

		return false;
	}
	;

	Element *elem[4];
};

class Triangle {
public:
	Triangle() {
		elem[0] = NULL;
		elem[1] = NULL;
		elem[2] = NULL;

	}
	;

	Triangle(Element* elem_1, Element* elem_2, Element* elem_3) {
		elem[0] = elem_1;
		elem[1] = elem_2;
		elem[2] = elem_3;
	}
	;

	bool operator<(const Triangle& rrect) const {
		if (elem[0] < rrect.elem[0] || (elem[0] == rrect.elem[0] && elem[1] < rrect.elem[1])
		    || (elem[0] == rrect.elem[0] && elem[1] == rrect.elem[1] && elem[2] < rrect.elem[2]))
			return true;

		return false;
	}
	;

	Element *elem[3];
};

typedef std::set<Edge> EdgeList;
typedef std::set<Quad> QuadList;
typedef std::set<Triangle> TriangList;

class Ellipse {
public:
	Ellipse(double x_center, double y_center, double x_radius, double y_radius, double cosrot,
	    double sinrot) :
			x_center(x_center), y_center(y_center), x_radius(x_radius), y_radius(y_radius), cosrot(
			    cosrot), sinrot(sinrot) {
	}

	const double get_x_center() const {
		return x_center;
	}

	const double get_y_center() const {
		return y_center;
	}

	const double get_x_radius() const {
		return x_radius;
	}

	const double get_y_radius() const {
		return y_radius;
	}

	const double get_cosrot() const {
		return cosrot;
	}

	const double get_sinrot() const {
		return sinrot;
	}

private:
	const double x_center;
	const double y_center;
	const double x_radius;
	const double y_radius;
	const double cosrot;
	const double sinrot;
};

double path_ellipse(void* path, double x, double y) {

	Ellipse* ellipse = (Ellipse*) path;

	const double x_radius = ellipse->get_x_radius();
	const double x_rad_sq = x_radius * x_radius;
	const double y_radius = ellipse->get_y_radius();
	const double y_rad_sq = y_radius * y_radius;
	const double x_center = ellipse->get_x_center();
	const double y_center = ellipse->get_y_center();
	const double cosrot = ellipse->get_cosrot();
	const double sinrot = ellipse->get_sinrot();

	return ((x - x_center) * cosrot + (y - y_center) * sinrot)
	    * ((x - x_center) * cosrot + (y - y_center) * sinrot) / x_rad_sq
	    + ((x - x_center) * sinrot - (y - y_center) * cosrot)
	        * ((x - x_center) * sinrot - (y - y_center) * cosrot) / y_rad_sq - 1.;

}

void grad_path_ellipse(void* path, double x, double y, double* grad) {

	Ellipse* ellipse = (Ellipse*) path;

	const double x_radius = ellipse->get_x_radius();
	const double x_rad_sq = x_radius * x_radius;
	const double y_radius = ellipse->get_y_radius();
	const double y_rad_sq = y_radius * y_radius;
	const double x_center = ellipse->get_x_center();
	const double y_center = ellipse->get_y_center();
	const double cosrot = ellipse->get_cosrot();
	const double sinrot = ellipse->get_sinrot();

	grad[0] = 2. * cosrot * ((x - x_center) * cosrot + (y - y_center) * sinrot) / x_rad_sq
	    + 2. * sinrot * ((x - x_center) * sinrot - (y - y_center) * cosrot) / y_rad_sq;

	grad[1] = 2. * sinrot * ((x - x_center) * cosrot + (y - y_center) * sinrot) / x_rad_sq
	    - 2. * cosrot * ((x - x_center) * sinrot - (y - y_center) * cosrot) / y_rad_sq;

}

void bilinear_interp(Quad quad, double* bilinear_coef) {

//  P=a0+a1*x+a2*y+a3*x*y
//	Ab=x
//	[1,x0,y0,x0y0][a0] [phi0]
//	[1,x1,y1,x1y1][a1] [phi1]
//	[1,x2,y2,x2y2][a2]=[phi2]
//	[1,x3,y3,x3y3][a3] [phi3]

	int dim = 4, one = 1, info, ipiv[dim];

	double A[dim * dim], x[dim], y[dim];

	for (int i = 0; i < dim; ++i) {
		x[i] = *(quad.elem[i]->get_coord());
		y[i] = *(quad.elem[i]->get_coord() + 1);
		bilinear_coef[i] = *(quad.elem[i]->get_state_vars());
		A[i] = 1;
		A[i + 4] = x[i];
		A[i + 8] = y[i];
		A[i + 12] = x[i] * y[i];
	}

//	cout << "Matrix A and vector phi" << endl;
//	for (int i = 0; i < dim; ++i) {
//		cout << "[ " << A[i] << " , " << A[i + 4] << " , " << A[i + 8] << " , " << A[i + 12] << " ]";
//		cout << "[ " << phi[i] << " ]" << endl;
//	}

	dgesv_(&dim, &one, A, &dim, ipiv, bilinear_coef, &dim, &info);

//	cout << "Solution" << endl;
//	for (int i = 0; i < dim; ++i)
//		cout << "[ " << phi[i] << " ]" << endl;

}

double bilinear_surface(void* path, double x, double y) {

//	p=a0+a1x+a2y+a3xy;

	double* bilinear_coef = (double*) path;

	return bilinear_coef[0] + bilinear_coef[1] * x + bilinear_coef[2] * y + bilinear_coef[3] * x * y;
}

void grad_bilinear_surface(void* path, double x, double y, double* grad) {

	//grad[0]=a1+a3y;
	//grad[1]=a2+a3x;

	double* bilinear_coef = (double*) path;

	grad[0] = bilinear_coef[1] + bilinear_coef[3] * y;
	grad[1] = bilinear_coef[2] + bilinear_coef[3] * x;

}

void make_surface(Triangle triangle, double* surface_coef) {
//  P=a0+a1*x+a2*y
//	Ab=x
//	[1,x0,y0][a0] [phi0]
//	[1,x1,y1][a1] [phi1]
//	[1,x2,y2][a2]=[phi2]

	int dim = 3, one = 1, info, ipiv[dim];

	double A[dim * dim], x[dim], y[dim];

	for (int i = 0; i < dim; ++i) {
		x[i] = *(triangle.elem[i]->get_coord());
		y[i] = *(triangle.elem[i]->get_coord() + 1);
		surface_coef[i] = *(triangle.elem[i]->get_state_vars());
		A[i] = 1;
		A[i + 3] = x[i];
		A[i + 6] = y[i];
	}

//	cout << "Matrix A and vector phi" << endl;
//	for (int i = 0; i < dim; ++i) {
//		cout << "[ " << A[i] << " , " << A[i + 4] << " , " << A[i + 8] << " , " << A[i + 12] << " ]";
//		cout << "[ " << phi[i] << " ]" << endl;
//	}

	dgesv_(&dim, &one, A, &dim, ipiv, surface_coef, &dim, &info);

//	cout << "Solution" << endl;
//	for (int i = 0; i < dim; ++i)
//		cout << "[ " << phi[i] << " ]" << endl;

}

double surface(void* path, double x, double y) {

//	p=a0+a1*x+a2*y;

	double* surface_coef = (double*) path;

	return surface_coef[0] + surface_coef[1] * x + surface_coef[2] * y;
}

void grad_surface(void* path, double x, double y, double* grad) {

	//grad[0]=a1;
	//grad[1]=a2;

	double* surface_coef = (double*) path;

	grad[0] = surface_coef[1];
	grad[1] = surface_coef[2];

}

void initialize_distance(Element* elem, Pt2Path& p_value_func, Pt2Grad& p_grad_func, void* ctx,
    double min_dx) {

	// this part is from reinitialization proposed in Chopp's, SOME IMPROVEMENTS OF THE FAST MARCHING METHOD
	// the difference is that instead of bicubic interpolation, we used bilinear interpolation to find the distance
	// near the interface. For the first time step we have the initial shape of pile and we can compute the
	// distance to the interface exactly, and we do not need to bilinear interpolation

	double d1[2] = { 0., 0. }, d2[2] = { 0., 0. }, grad[2], pvalue, grad_dot, x_new, y_new, x_old,
	    y_old, x_half, y_half, epslon = .01;

// initializing the solution
	x_new = x_old = *(elem->get_coord());
	y_new = y_old = *(elem->get_coord() + 1);

	double& phi_old = *(elem->get_state_vars());

	int iter = 0, max_iter = 20;// Chopp's paper says it normally has to converge after 4 or 5 iterations
	double sgn = 1.;
	const double treshold = 1e-3 * min_dx * min_dx;

	do {
		iter++;

		p_grad_func(ctx, x_new, y_new, grad);
		pvalue = p_value_func(ctx, x_new, y_new);

		grad_dot = grad[0] * grad[0] + grad[1] * grad[1];

		d1[0] = -pvalue * grad[0] / grad_dot;
		d1[1] = -pvalue * grad[1] / grad_dot;
		x_half = x_new + d1[0];
		y_half = y_new + d1[1];

		d2[0] = (x_old - x_new) - (x_old - x_new) * grad[0] / grad_dot * grad[0];
		d2[1] = (y_old - y_new) - (y_old - y_new) * grad[1] / grad_dot * grad[1];

		x_new = x_half + d2[0];
		y_new = y_half + d2[1];

	} while (sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d2[0] * d2[0] + d2[1] * d2[1]) > treshold
	    && iter < max_iter);

	if (phi_old < 0.)
		sgn = -1.;

	double phi_new = sqrt((x_new - x_old) * (x_new - x_old) + (y_new - y_old) * (y_new - y_old));

	// this is for the case that a point is computed twice, and we take the smaller value
	if (phi_new < fabs(phi_old))
		phi_old = sgn * phi_new;
//	std::cout << point.get_phi() << std::endl;
}

void adjacent_to_interface(HashTable* El_Table, EdgeList& accepted) {

	HashEntryPtr* buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)

		if (*(buck + i)) {

			HashEntryPtr currentPtr = *(buck + i);

			while (currentPtr) {

				Element* Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					double *phi = (Curr_El->get_state_vars());

					for (int ineigh = 0; ineigh < 8; ineigh++)

						if (*(Curr_El->get_neigh_proc() + ineigh) >= 0) {
							// this condition is to avoid elements on the boundary, or duplicated neighbors

							unsigned* neigh_key = Curr_El->get_neighbors();

							Element* ElemNeigh = (Element*) El_Table->lookup(neigh_key + ineigh * KEYLENGTH);

							double neighb_phi = *(ElemNeigh->get_state_vars());

							if (neighb_phi * phi[0] < 0.) {

								int yp, xm, ym, xp = Curr_El->get_positive_x_side();
								yp = (xp + 1) % 4;
								xm = (xp + 2) % 4;
								ym = (xp + 3) % 4;

								int orientation = abs(ineigh - xp);

								if (orientation % 4 == 0 || orientation % 4 == 1)
									accepted.insert(Edge(Curr_El, ElemNeigh, ineigh));
								else
									accepted.insert(
									    Edge(ElemNeigh, Curr_El, ElemNeigh->which_neighbor(Curr_El->pass_key())));
							}

						}
				}
				currentPtr = currentPtr->next;
			}

		}
}

void calc_phi_slope(HashTable* El_Table, HashTable* NodeTable) {

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		HashEntryPtr currentPtr = *(buck + i);
		while (currentPtr) {
			Element* Em_Temp = (Element*) (currentPtr->value);
			if (Em_Temp->get_adapted_flag() > 0)
				Em_Temp->calc_phi_slope(El_Table, NodeTable);

			currentPtr = currentPtr->next;
		}
	}
}

void test_nbflag(HashTable* El_Table) {

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		HashEntryPtr currentPtr = *(buck + i);
		while (currentPtr) {
			Element* Em_Temp = (Element*) (currentPtr->value);
			if (Em_Temp->get_adapted_flag() > 0 && *(Em_Temp->get_nbflag()))
				cout << "Error this should not happen" << endl;
			currentPtr = currentPtr->next;
		}
	}
}

void reset_nbflag(HashTable* El_Table) {

	HashEntryPtr *buck = El_Table->getbucketptr();
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		HashEntryPtr currentPtr = *(buck + i);
		while (currentPtr) {
			Element* Em_Temp = (Element*) (currentPtr->value);
			if (Em_Temp->get_adapted_flag() > 0)
				*(Em_Temp->get_nbflag()) = 0;
			currentPtr = currentPtr->next;
		}
	}
}

void update_phi(HashTable* El_Table, double min_dx, double* norm, int* elem) {

	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Element* Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0 && *(Curr_El->get_nbflag()) != 1
				// I did note get good result from the following condition
				    && fabs(*(Curr_El->get_state_vars())) <= 10 * min_dx) {
					// with the last condition we narrow the range of update
					// to maximum of 10 cell-width far the the interface

					(*elem)++;

					double *phi_slope = Curr_El->get_phi_slope();
					double *state_vars = Curr_El->get_state_vars();
					double *prev_state_vars = Curr_El->get_prev_state_vars();
					double delta_p, delta_m;
					double flux_phi, a_p, a_m, b_p, b_m, c_p, c_m, d_p, d_m;
					const double CFL = .5;
					const double dt = CFL * min_dx;

					// based on Sussman Smerka 1994

					prev_state_vars[0] = state_vars[0];

					a_p = max(phi_slope[0], 0.);
					a_m = min(phi_slope[0], 0.);
					b_p = max(phi_slope[1], 0.);
					b_m = min(phi_slope[1], 0.);
					c_p = max(phi_slope[2], 0.);
					c_m = min(phi_slope[2], 0.);
					d_p = max(phi_slope[3], 0.);
					d_m = min(phi_slope[3], 0.);

					if (state_vars[4] > 0.)
						flux_phi = sqrt(max(a_p * a_p, b_m * b_m) + max(c_p * c_p, d_m * d_m)) - 1;
					else if (state_vars[4] < 0.)
						flux_phi = sqrt(max(a_m * a_m, b_p * b_p) + max(c_m * c_m, d_p * d_p)) - 1;
					else
						flux_phi = 0.;

					state_vars[0] = state_vars[0] - dt * sign(state_vars[4]) * flux_phi;

					*norm += (state_vars[0] - prev_state_vars[0]) * (state_vars[0] - prev_state_vars[0]);

					// this is just to see which elements are updating
					state_vars[5] = 1.;

				}

				currentPtr = currentPtr->next;
			}
		}
}

void record_of_phi(HashTable* El_Table) {

	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Element* Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					*(Curr_El->get_state_vars() + 4) = *(Curr_El->get_state_vars());

				currentPtr = currentPtr->next;
			}
		}
}

void make_quad_trangle(HashTable* El_Table, EdgeList& accepted, QuadList& quadlist,
    TriangList& trianglist) {

	Element* elem[4];
	int neigh_num;
	int debug = 0;

	for (EdgeList::iterator it = accepted.begin(); it != accepted.end(); ++it) {
		debug++;
//		cout<<debug<<endl;

		if (it->elem[0]->get_adapted_flag() > 0) {
			elem[0] = it->elem[0];
			elem[1] = it->elem[1];
			neigh_num = it->neigh_num;
		} else {
			elem[0] = it->elem[1];
			elem[1] = it->elem[0];
			neigh_num = elem[0]->which_neighbor(elem[1]->pass_key());					//= (it->neigh_num + 2) % 8;
		}

		// given the neighbor number, we can build 2 Quads, but we do not care.
		// what we need to do is to find 4 points for making approximated surface
		// from the bilinear map, and then computing the distance to the interface

		unsigned* elem_key = elem[0]->get_neighbors();
		elem[2] = (Element*) El_Table->lookup(elem_key + ((neigh_num + 1) % 8) * KEYLENGTH);
		unsigned* neigh_key;

		if (elem[2] && elem[2]->get_adapted_flag() > 0) {
			// first condition is to avoid the case that elem[2] is not available

			neigh_key = elem[2]->get_neighbors();
			elem[3] = (Element*) El_Table->lookup(neigh_key + neigh_num * KEYLENGTH);
		} else {
			// this means that this proc has not access to this direction, so we have to change the direction
			elem[2] = (Element*) El_Table->lookup(elem_key + ((neigh_num - 1 + 8) % 8) * KEYLENGTH);

			if (elem[2]) {
				// the condition is to avoid the case that elem[2] is not available
				neigh_key = elem[2]->get_neighbors();
				elem[3] = (Element*) El_Table->lookup(neigh_key + neigh_num * KEYLENGTH);
			}
		}

		if (elem[3]) {
			quadlist.insert(Quad(elem[0], elem[1], elem[2], elem[3]));
		} else {
			// this is for the case that the element is in the corner, and three of its sides are in
			// in the other processors, in this case we just make a triangle
			trianglist.insert(Triangle(elem[0], elem[1], elem[2]));
		}

//		if (!(elem[0] && elem[1] && elem[2] && elem[3])) {
//			cout << "adaption flags are: " << elem[0]->get_adapted_flag() << " , "
//					<< elem[1]->get_adapted_flag() << " , " << elem[2]->get_adapted_flag() << endl;
//			cout << "x of quad: " << *(elem[0]->get_coord()) << " " << *(elem[1]->get_coord()) << " "
//					<< *(elem[2]->get_coord()) << endl;
//			cout << "y of quad: " << *(elem[0]->get_coord() + 1) << " " << *(elem[1]->get_coord() + 1)
//					<< " " << *(elem[2]->get_coord() + 1) << endl;
//			cout << endl;
//		}

	}

}

void pde_reinitialization(HashTable* El_Table, HashTable* NodeTable, TimeProps* timeprops,
    double min_dx, int nump, int rank) {

	const double threshold = .5 * min_dx * min_dx * min_dx;
	double norm, total_norm, normalized_norm;
	int elem, tot_elem;

	record_of_phi(El_Table);

	int iter = 0;

	do {
		iter++;
		elem = 0;
		tot_elem = 0;
		norm = 0.;
		normalized_norm = 0.;
		calc_phi_slope(El_Table, NodeTable);
		update_phi(El_Table, min_dx, &norm, &elem);

		// after updating this data have to be transmitted into others
		move_data(nump, rank, El_Table, NodeTable, timeprops);
		MPI_Allreduce(&norm, &total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&elem, &tot_elem, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		normalized_norm = sqrt(total_norm) / tot_elem;

		// CFL number is 0.5 and we use 12 iteration to make sure at least 6 neighbor element has been updated
	} while (normalized_norm > threshold && iter < 13);

	if (rank == 0)
		cout << "norm: " << normalized_norm << " threshold is: " << threshold << endl;

}

void reinitialization(HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr,
    TimeProps *timeprops, PileProps *pileprops_ptr, int nump, int rank) {

	reset_nbflag(El_Table);

	EdgeList accepted;
	adjacent_to_interface(El_Table, accepted);

	QuadList rectangles;
	TriangList trianglist;

// with each Edge we can build 2 rectangles, but this procedure can be repeated
// for an edge, so make a new list of edges to prevent making these again
	make_quad_trangle(El_Table, accepted, rectangles, trianglist);

	double min_dx = 0.;
	find_min_dx(El_Table, &min_dx);

	Pt2Path p2path = bilinear_surface;
	Pt2Grad p2grad = grad_bilinear_surface;

	double bilinear_coef[4];

	for (QuadList::iterator it = rectangles.begin(); it != rectangles.end(); ++it) {

		bilinear_interp(*it, bilinear_coef);

		for (int j = 0; j < 4; ++j) {
			initialize_distance(it->elem[j], p2path, p2grad, (void*) bilinear_coef, min_dx);
			*((it->elem[j])->get_nbflag()) = 1;
		}
	}

	p2path = surface;
	p2grad = grad_surface;

	double surface_coef[4];

	for (TriangList::iterator it = trianglist.begin(); it != trianglist.end(); ++it) {

		make_surface(*it, bilinear_coef);

		for (int j = 0; j < 3; ++j) {
			initialize_distance(it->elem[j], p2path, p2grad, (void*) surface_coef, min_dx);
			*((it->elem[j])->get_nbflag()) = 1;
		}
	}

	pde_reinitialization(El_Table, NodeTable, timeprops, min_dx, nump, rank);
}

void initialization(HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr,
    TimeProps *timeprops, PileProps *pileprops_ptr, int nump, int rank) {

// data in pileprops_ptr are scaled
	Ellipse ellipse(pileprops_ptr->xCen[0], pileprops_ptr->yCen[0], pileprops_ptr->majorrad[0],
	    pileprops_ptr->minorrad[0], pileprops_ptr->cosrot[0], pileprops_ptr->sinrot[0]);

	reset_nbflag(El_Table);
	EdgeList accepted;
	adjacent_to_interface(El_Table, accepted);

	double min_dx = 0.;
	find_min_dx(El_Table, &min_dx);

	Pt2Path p2path = path_ellipse;
	Pt2Grad p2grad = grad_path_ellipse;

	for (EdgeList::iterator it = accepted.begin(); it != accepted.end(); ++it)
		for (int j = 0; j < 2; ++j)
			if (it->elem[j]->get_adapted_flag() > 0) {		// we do not want to update ghost element
				initialize_distance(it->elem[j], p2path, p2grad, (void*) &ellipse, min_dx);
				*((it->elem[j])->get_nbflag()) = 1;
			}

	pde_reinitialization(El_Table, NodeTable, timeprops, min_dx, nump, rank);

}

int num_nonzero_elem(HashTable *El_Table) {
	int num = 0, myid;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	//printf("the numebr of buckets are: %d", El_Table->get_no_of_buckets());

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {
					num++;
				}
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}
