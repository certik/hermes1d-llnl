// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "mesh.h"
#include "iterator.h"
#include "adapt.h"

Element::Element() 
{
  x1 = x2 = 0;
  p = 0; 
  dof = NULL;
  sons[0] = sons[1] = NULL; 
  active = 1;
  level = 0;
  id = -1;
  dof_size = 0;
}

Element::Element(double x_left, double x_right, int deg, int n_eq) 
{
  x1 = x_left;
  x2 = x_right;
  p = deg; 
  dof_size = n_eq;
  this->dof_alloc();
  if (dof == NULL) error("Not enough memory in Element().");
  sons[0] = sons[1] = NULL; 
  active = 1;
  level = 0;
  id = -1;
}

unsigned Element::is_active() 
{
  return this->active;
}

void Element::refine(int3 cand) 
{
  int hp_refinement = cand[0];
  if(!hp_refinement) {
    this->p = cand[1];
  }
  else {
    double x1 = this->x1;
    double x2 = this->x2;
    double midpoint = (x1 + x2)/2.; 
    this->sons[0] = new Element(x1, midpoint, cand[1], dof_size);
    this->sons[1] = new Element(midpoint, x2, cand[2], dof_size);
    this->sons[0]->level = this->level + 1; 
    this->sons[1]->level = this->level + 1; 
    // copying negative dof to sons if any
    for(int c=0; c<dof_size; c++) {
      if (this->dof[c][0] < 0) this->sons[0]->dof[c][0] = this->dof[c][0];
      if (this->dof[c][1] < 0) this->sons[1]->dof[c][1] = this->dof[c][1];
    }
    this->active = 0;
  }
}

void Element::init(double x1, double x2, int p_init, int n_eq)
{
  this->x1 = x1;
  this->x2 = x2;
  this->p = p_init;
  this->dof_size = n_eq;
  // allocate element dof arrays for all solution components 
  // and length MAX_POLYORDER
  this->dof_alloc();
}

void Element::dof_alloc() 
{
  this->dof = (int**)malloc(dof_size*sizeof(int*));
  if(this->dof == NULL) error("Element dof_alloc() failed.");
  // c is solution component
  for(int c=0; c<dof_size; c++) {
    this->dof[c] = new int[MAX_POLYORDER + 1];
    // important for th etreatment of boundary conditions
    for(int i=0; i<MAX_POLYORDER + 1; i++) this->dof[c][i] = 0;
  }
}

// return coefficients for all shape functions on the element m,
// for all solution components
void Element::get_coeffs(double *y_prev, 
                         double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM],
                         double bc_left_dir_values[MAX_EQN_NUM],
                         double bc_right_dir_values[MAX_EQN_NUM])
{
  if (!this->is_active()) error("Internal in calculate_elem_coeffs().");
  int dof_size = this->dof_size;
  for(int c=0; c<dof_size; c++) {
    // coeff of the left vertex function
    if (this->dof[c][0] == -1) coeffs[c][0] = bc_left_dir_values[c];
    else coeffs[c][0] = y_prev[this->dof[c][0]];
    // coeff of the right vertex function
    if (this->dof[c][1] == -1) coeffs[c][1] = bc_right_dir_values[c];
    else coeffs[c][1] = y_prev[this->dof[c][1]];
    //completing coeffs of bubble functions
    for (int j=2; j<=this->p; j++) {
        coeffs[c][j] = y_prev[this->dof[c][j]];
    }
  }
}

// transforms point 'x_phys' from element (x1, x2) to (-1, 1)
double inverse_map(double x1, double x2, double x_phys) 
{
  double c = x1 - x2;
  double a = -2/c;
  double b = (x1 + x2)/c; 
  return a*x_phys + b;
}

// Evaluate solution and its derivatives in the points 'x_phys' 
// in the element (coeffs[][] array provided).
void Element::get_solution(double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                           int pts_num, double x_phys[MAX_PTS_NUM], 
                           double val_phys[MAX_EQN_NUM][MAX_PTS_NUM], 
                           double der_phys[MAX_EQN_NUM][MAX_PTS_NUM])
{
  double x1 = this->x1;
  double x2 = this->x2;
  double jac = (x2-x1)/2.; // Jacobian of reference map
  int dof_size = this->dof_size;
  int p = this->p;
  double x_ref[MAX_PTS_NUM];
  // transforming points to (-1, 1)
  for (int i=0 ; i<pts_num; i++) x_ref[i] = inverse_map(x1, x2, x_phys[i]);
  // filling the values and derivatives
  for(int c=0; c<dof_size; c++) { 
    for (int i=0 ; i<pts_num; i++) {
      der_phys[c][i] = val_phys[c][i] = 0;
      for(int j=0; j<=p; j++) {
        val_phys[c][i] += coeff[c][j]*lobatto_fn_tab_1d[j](x_ref[i]);
        der_phys[c][i] += coeff[c][j]*lobatto_der_tab_1d[j](x_ref[i]);
      }
      der_phys[c][i] /= jac;
    }
  }
} 

// Evaluate solution and its derivatives in the points 'x_phys' 
// in the element (coeffs[][] array not provided).
void Element::get_solution(double x_phys[MAX_PTS_NUM], int pts_num,  
                           double val_phys[MAX_EQN_NUM][MAX_PTS_NUM], 
                           double der_phys[MAX_EQN_NUM][MAX_PTS_NUM],
                           double *y_prev, 
                           double bc_left_dir_values[MAX_EQN_NUM],
                           double bc_right_dir_values[MAX_EQN_NUM])
{
  double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->get_coeffs(y_prev, coeff, bc_left_dir_values, bc_right_dir_values);
  this->get_solution(coeff, pts_num, x_phys, val_phys, der_phys);
} 

// evaluate solution and its derivative 
// at the point x_phys (coeffs[][] array provided)
void Element::get_solution_point(double x_phys,
			         double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                                 double val[MAX_EQN_NUM], double der[MAX_EQN_NUM])
{
  double x1 = this->x1;
  double x2 = this->x2;
  double jac = (x2-x1)/2.;
  int dof_size = this->dof_size; 
  int p = this->p;
  // transforming point x_phys to (-1, 1)
  double x_ref = inverse_map(x1, x2, x_phys);
  for(int c=0; c<dof_size; c++) {
    der[c] = val[c] = 0;
    for(int j=0; j<=p; j++) {
      val[c] += coeff[c][j]*lobatto_fn_tab_1d[j](x_ref);
      der[c] += coeff[c][j]*lobatto_der_tab_1d[j](x_ref);
    }
    der[c] /= jac;
  }
} 

// evaluate solution and its derivative 
// at the point x_phys (coeffs[][] array not provided)
void Element::get_solution_point(double x_phys,
				 double val_phys[MAX_EQN_NUM], double der_phys[MAX_EQN_NUM],
                                 double *y_prev, double *bc_left_dir_values,
                                 double *bc_right_dir_values)
{
  double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->get_coeffs(y_prev, coeff, bc_left_dir_values, bc_right_dir_values); 
  this->get_solution_point(x_phys, coeff, val_phys, der_phys);
} 

// copying elements including their refinement trees 
// inactive Dirichlet DOF are replicated
// FIXME - the recursive version is slow, improve it!
void Element::copy_sons_recursively(Element *e_trg) {
  // if element has been refined
  if(this->sons[0] != NULL) {
    int p_left = this->sons[0]->p;
    int p_right = this->sons[1]->p;
    int3 cand = {1, p_left, p_right};
    e_trg->refine(cand);
    // left son
    this->sons[0]->copy_sons_recursively(e_trg->sons[0]);
    // right son
    this->sons[1]->copy_sons_recursively(e_trg->sons[1]);
  }
}

// gets physical coordinate of a reference point x_ref \in [-1,1]
double Element::get_x_phys(double x_ref) 
{
  double a = this->x1;
  double b = this->x2;
  return (a+b)/2. + x_ref*(b-a)/2.;
}

Mesh::Mesh() {
  n_eq = 0;
  n_base_elem = 0;
  n_active_elem = 0;
  n_dof = 0;
  base_elems = NULL;
}

// creates equidistant mesh with uniform polynomial degree of elements
Mesh::Mesh(double a, double b, int n_base_elem, int p_init, int n_eq)
{
  // domain end points
  left_endpoint = a;
  right_endpoint = b;
  // number of equations
  this->n_eq = n_eq;
  // print the banner (only once)
  static int n_calls = 0;
  n_calls++;
  if (n_calls == 1) intro();
  // check maximum number of equations
  if(n_eq > MAX_EQN_NUM) 
  error("Maximum number of equations exceeded (set in common.h)");
  // arrays for boundary conditions
  for (int i=0; i<MAX_EQN_NUM; i++) {
    this->bc_left_dir_values[i] = 0;
    this->bc_right_dir_values[i] = 0;
  }
  // number of base elements
  this->n_base_elem = n_base_elem;
  // number of active elements
  this->n_active_elem = n_base_elem;
  // allocate element array
  this->base_elems = new Element[this->n_base_elem];     
  if (base_elems == NULL) error("Not enough memory in Mesh::create().");
  if (p_init > MAX_POLYORDER) 
    error("Max element order exceeded (set in common.h).");
  // element length
  double h = (b - a)/this->n_base_elem;          
  // fill initial element array
  for(int i=0; i<this->n_base_elem; i++) {         
    this->base_elems[i].init(a + i*h, a + i*h + h, p_init, n_eq);
  }
  this->assign_elem_ids();
}

// creates mesh using a given array of n_base_elem+1 points (pts_array)
// and an array of n_base_elem polynomial degrees (p_array)
Mesh::Mesh(int n_base_elem, double *pts_array, int *p_array, int n_eq)
{
  // domain end points
  left_endpoint = pts_array[0];
  right_endpoint = pts_array[n_base_elem];
  // number of equations
  this->n_eq = n_eq;
  // print the banner (only once)
  static int n_calls = 0;
  n_calls++;
  if (n_calls == 1) intro();
  // check maximum number of equations
  if(n_eq > MAX_EQN_NUM) 
  error("Maximum number of equations exceeded (set in common.h)");
  // arrays for boundary conditions
  for (int i=0; i<MAX_EQN_NUM; i++) {
    this->bc_left_dir_values[i] = 0;
    this->bc_right_dir_values[i] = 0;
  }
  // number of base elements
  this->n_base_elem = n_base_elem;
  // number of active elements
  this->n_active_elem = n_base_elem;
  // allocate element array
  this->base_elems = new Element[this->n_base_elem];     
  if (base_elems == NULL) error("Not enough memory in Mesh::create().");
  // fill initial element array
  for(int i=0; i<this->n_base_elem; i++) {
    if (p_array[i] > MAX_POLYORDER) 
      error("Max element order exceeded (set in common.h).");
    // polynomial degree
    this->base_elems[i].p = p_array[i];
    // allocate element dof arrays for all solution components 
    // and length MAX_POLYORDER
    this->base_elems[i].dof_alloc();
    // define element end points
    this->base_elems[i].x1 = pts_array[i];
    this->base_elems[i].x2 = pts_array[i+1];
  }
  this->assign_elem_ids();
}

// caution - this is expensive (traverses the entire tree 
// from the beginning until the element is found)
void Mesh::refine_single_elem(int id, int3 cand)
{
    Iterator I(this);
    Element *e;
    while ((e = I.next_active_element()) != NULL) {
        printf("%d %d\n", e->id, id);
        if (e->id == id) {
            e->refine(cand);
            if (cand[0] == 1) this->n_active_elem++; // hp-refinement
            return;
        }
    }
    error("refine_single_elem: Element not found.");
}

// performs mesh refinement using a list of elements to be 
// refined and a list of corresponding polynomial degree 
// pairs for the sons
void Mesh::refine_elems(int elem_num, int *id_array, int3 *cand_array)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id == id_array[count]) {
            if (count >= elem_num)
                error("refine_multi_elems: not enough elems specified");
            e->refine(cand_array[count]);
            this->n_active_elem++;
            count++;
        }
    }
}

// splits the indicated elements and 
// increases poly degree in sons by one
void Mesh::reference_refinement(int start_elem_id, int elem_num)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id >= start_elem_id && e->id < start_elem_id + elem_num) {
	    if (count >= elem_num) return;
            int3 cand = {1, e->p + 1, e->p + 1};
            e->refine(cand);
            if (cand[0] == 1) this->n_active_elem++; // if hp-refinement
            count++;
        }
    }
}

void Mesh::set_bc_left_dirichlet(int eq_n, double val)
{
  this->bc_left_dir_values[eq_n] = val;
  // deactivate the corresponding dof for the left-most
  // element and all his descendants adjacent to the 
  // left boundary
  Element *e = this->base_elems + 0;
  do {
    e->dof[eq_n][0] = -1;
    e = e->sons[0];
  } while (e != NULL);
}

void Mesh::set_bc_right_dirichlet(int eq_n, double val)
{
  this->bc_right_dir_values[eq_n] = val;
  // deactivate the corresponding dof for the right-most
  // element and all his descendants adjacent to the 
  // right boundary
  Element *e = this->base_elems + this->n_base_elem - 1;
  do {
    e->dof[eq_n][1] = -1;
    e = e->sons[1];
  } while (e != NULL);

}

// define element connectivities (global dof)
int Mesh::assign_dofs()
{
  Iterator *I = new Iterator(this);
  // (1) enumerate vertex dofs
  int count_dof = 0;
  // loop over solution components
  for(int c=0; c<this->n_eq; c++) {    
    Element *e;
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
      if (e->dof[c][0] != -1) e->dof[c][0] = count_dof++; 
      if (e->dof[c][1] != -1) e->dof[c][1] = count_dof; 
      else count_dof--;
    }
    count_dof++;
    // (2) enumerate bubble dofs
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
      for(int j=2; j<=e->p; j++) {
        e->dof[c][j] = count_dof;
        count_dof++;
      }
    }
  }
  this->n_dof = count_dof;

  // enumerate elements
  this->assign_elem_ids();

  // print element connectivities
  if(DEBUG_ELEM_DOF) {
    printf("Printing element DOF arrays:\n");
    printf("Elements = %d\n", this->n_base_elem);
    printf("DOF = %d", this->n_dof);
    for(int c = 0; c<this->n_eq; c++) {
      I->reset();
      Element *e;
      while ((e = I->next_active_element()) != NULL) {
        printf("\nElement (%g, %g), id = %d, p = %d\n ", 
               e->x1, e->x2, e->id, e->p); 
        for(int j = 0; j<e->p + 1; j++) {
          printf("dof[%d][%d] = %d\n ", c, j, e->dof[c][j]);
        }
      }
    }
    printf("\n"); 
  }

  delete I;
  return this->n_dof;
}

int Mesh::assign_elem_ids()
{
    Iterator *I = new Iterator(this);
    int count_id = 0;
    Element *e;
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
        e->id = count_id++;
    }
    delete I;
}

Element* Mesh::first_active_element()
{
  Element *e = base_elems;
  while(!e->is_active()) {
    e = e->sons[0];
  }
  return e;
}

Element* Mesh::last_active_element()
{
  Element *e = base_elems + n_base_elem - 1;
  while(!e->is_active()) {
    e = e->sons[1];
  }
  return e;
}

int Element::create_cand_list(int p_ref_left, int p_ref_right, int3 *cand_list) 
{
    int counter = 0;
    // p->p+1
    if (this->p + 1 <= max(p_ref_left, p_ref_right)) {
      cand_list[counter][0] = 0;
      cand_list[counter][1] = this->p + 1;
      cand_list[counter][2] = -1;
      counter++;
    }
    // p->p+2
    if (this->p + 2 <= max(p_ref_left, p_ref_right)) {
      cand_list[counter][0] = 0;
      cand_list[counter][1] = this->p + 2;
      cand_list[counter][2] = -1;
      counter++;
    }
    // so that the if statements below are simplified:
    if (p_ref_left == -1) p_ref_left = 100000;
    if (p_ref_right == -1) p_ref_right = 100000;
    // p -> (p/2, p/2) 
    int base_p = this->p / 2; 
    if (base_p < 1) base_p = 1;
    if (base_p <= p_ref_left && base_p <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p;
        cand_list[counter][2] = base_p;
        counter++;
    }
    // p -> (p/2+1, p/2) 
    if (base_p + 1 <= p_ref_left && base_p <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p+1;
        cand_list[counter][2] = base_p;
        counter++;
    }
    // p -> (p/2, p/2+1) 
    if (base_p <= p_ref_left && base_p + 1 <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p;
        cand_list[counter][2] = base_p+1;
        counter++;
    }
    // p -> (p/2+1, p/2+1) 
    if (base_p + 1 <= p_ref_left && base_p + 1 <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p+1;
        cand_list[counter][2] = base_p+1;
        counter++;
    }
    // p -> (p/2+2, p/2) 
    if (base_p + 2 <= p_ref_left && base_p <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p+2;
        cand_list[counter][2] = base_p;
        counter++;
    }
    // p -> (p/2, p/2+2) 
    if (base_p <= p_ref_left && base_p + 2 <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p;
        cand_list[counter][2] = base_p+2;
        counter++;
    }
    // p -> (p/2+1, p/2+2) 
    if (base_p + 1 <= p_ref_left && base_p + 2 <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p+1;
        cand_list[counter][2] = base_p+2;
        counter++;
    }
    // p -> (p/2+2, p/2+1) 
    if (base_p + 2 <= p_ref_left && base_p + 1 <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p+2;
        cand_list[counter][2] = base_p+1;
        counter++;
    }
    // p -> (p/2+2, p/2+2) 
    if (base_p + 2 <= p_ref_left && base_p + 2 <= p_ref_right) {
        cand_list[counter][0] = 1;
        cand_list[counter][1] = base_p+2;
        cand_list[counter][2] = base_p+2;
        counter++;
    }

    return counter;
}

void Element::print_cand_list(int num_cand, int3 *cand_list) 
{
  printf("Element (%g, %g): refinement candidates:\n", this->x1, this->x2);
  for (int i=0; i < num_cand; i++) { 
    printf("%d %d %d\n", cand_list[i][0], cand_list[i][1], cand_list[i][2]);
  }
}

// transformation of k-th shape function defined on Gauss points 
// corresponding to 'order' to physical interval (a,b)
void element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  int pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0 ; i<pts_num; i++) {
    // change function values and derivatives to interval (a, b)
    val[i] = lobatto_fn_tab_1d[k](ref_tab[i][0]);
    double jac = (b-a)/2.; 
    der[i] = lobatto_der_tab_1d[k](ref_tab[i][0]) / jac; 
  }
};

// transformation of k-th shape function at the reference 
// point x_ref to physical interval (a,b).
void element_shapefn_point(double x_ref, double a, double b, 
		           int k, double *val, double *der) {
    // change function values and derivatives to interval (a, b)
    *val = lobatto_fn_tab_1d[k](x_ref);
    double jac = (b-a)/2.; 
    *der = lobatto_der_tab_1d[k](x_ref) / jac; 
}

// Replicate mesh
Mesh *Mesh::replicate()
{
  // copy base mesh element array, use dummy poly degree (p_init
  // cannot be used since p-refinements on base mesh may have taken 
  // place)
  int p_dummy = -1; 
  Mesh *mesh_ref = new Mesh(this->left_endpoint, this->right_endpoint, 
			    this->n_base_elem, p_dummy, this->n_eq);

  // copy element degrees on base mesh
  for(int i=0; i<this->n_base_elem; i++) {
    mesh_ref->get_base_elems()[i].p = this->get_base_elems()[i].p;
  }

  // copy number of base elements
  mesh_ref->set_n_base_elem(this->get_n_base_elem());
    
  // copy arrays of Dirichlet boundary conditions
  for(int c = 0; c<this->n_eq; c++) {
    mesh_ref->bc_left_dir_values[c] = this->bc_left_dir_values[c]; 
    mesh_ref->bc_right_dir_values[c] = this->bc_right_dir_values[c]; 
  }

  // copy inactive Dirichlet DOF on first and last element of base mesh
  for(int c=0; c<this->n_eq; c++) {
    if(this->base_elems[0].dof[c][0] == -1) {
      Element *e_left = mesh_ref->get_base_elems();
      e_left->dof[c][0] = -1;
    }
    if(this->base_elems[this->n_base_elem-1].dof[c][1] == -1) {
      Element *e_right = mesh_ref->get_base_elems()+this->n_base_elem-1;
      e_right->dof[c][1] = -1;
    }
  }

  // replicate tree structure on all base mesh elements,
  // this includes replication of inactive Dirichlet DOF
  for(int i=0; i<this->n_base_elem; i++) {
    Element *e_src = this->base_elems + i;
    Element *e_trg = mesh_ref->get_base_elems() + i;
    e_src->copy_sons_recursively(e_trg);
  }

  // enumerate elements in reference mesh
  mesh_ref->assign_elem_ids();

  return mesh_ref;
}

void Mesh::plot(const char* filename) 
{
    FILE *f = fopen(filename, "wb");
    if(f == NULL) error("problem opening file in Mesh::plot().");
    Iterator I(this);
    Element *e;
    while ((e = I.next_active_element()) != NULL) {
      fprintf(f, "%g %d\n", e->x1, 0);
      fprintf(f, "%g %d\n", e->x1, e->p);
      fprintf(f, "%g %d\n", e->x2, e->p);
      fprintf(f, "%g %d\n\n", e->x2, 0);
    }
    fclose(f);
    printf("Mesh written to %s.\n", filename);
}

// Plots the error between the reference and coarse mesh solutions
// if the reference refinement was p-refinement
void Mesh::plot_element_error_p(FILE *f[MAX_EQN_NUM], Element *e, Element *e_ref, 
                                double *y_prev, double *y_prev_ref,
                                int subdivision)
{
  int n_eq = this->get_n_eq();
  int pts_num = subdivision + 1;
  double x1 = e->x1;
  double x2 = e->x2;

  // calculate point array
  double x_phys[MAX_PTS_NUM];  
  double h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx[MAX_EQN_NUM][MAX_PTS_NUM];
  e->get_solution(x_phys, pts_num, phys_u, phys_dudx, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  // get fine mesh solution values and derivatives
  double phys_u_ref[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx_ref[MAX_EQN_NUM][MAX_PTS_NUM];
  e->get_solution(x_phys, pts_num, phys_u_ref, phys_dudx_ref, y_prev_ref, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int c=0; c < n_eq; c++) {
    for (int i=0; i < pts_num; i++) {
      fprintf(f[c], "%g %g\n", x_phys[i], phys_u_ref[c][i] - phys_u[c][i]);
    }
    fprintf(f[c], "\n");
  }
}

// Plots the error between the reference and coarse mesh solutions
// if the reference refinement was hp-refinement
void Mesh::plot_element_error_hp(FILE *f[MAX_EQN_NUM], Element *e, Element *e_ref_left, 
                                 Element *e_ref_right, 
                                 double *y_prev, double *y_prev_ref,
                                 int subdivision)
{
  int n_eq = this->get_n_eq();
  // we will be using two intervals of 50% length
  subdivision /= 2;
  int pts_num = subdivision + 1;

  // First: left half (-1, 0)
  double x1 = e_ref_left->x1;
  double x2 = e_ref_left->x2;

  // calculate point array
  double x_phys_left[MAX_PTS_NUM];  
  double h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys_left[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u_left[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx_left[MAX_EQN_NUM][MAX_PTS_NUM];
  e->get_solution(x_phys_left, pts_num, phys_u_left, phys_dudx_left, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  // get fine mesh solution values and derivatives
  double phys_u_ref_left[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx_ref_left[MAX_EQN_NUM][MAX_PTS_NUM];
  e_ref_left->get_solution(x_phys_left, pts_num, phys_u_ref_left, phys_dudx_ref_left, y_prev_ref, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int c=0; c < n_eq; c++) {
    for (int i=0; i < pts_num; i++) {
      fprintf(f[c], "%g %g\n", x_phys_left[i], phys_u_ref_left[c][i] - phys_u_left[c][i]);
    }
  }

  // Second: right half (0, 1)
  x1 = e_ref_right->x1;
  x2 = e_ref_right->x2;

  // calculate point array
  double x_phys_right[MAX_PTS_NUM];  
  h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys_right[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u_right[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx_right[MAX_EQN_NUM][MAX_PTS_NUM];
  e->get_solution(x_phys_right, pts_num, phys_u_right, phys_dudx_right, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  // get fine mesh solution values and derivatives
  double phys_u_ref_right[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx_ref_right[MAX_EQN_NUM][MAX_PTS_NUM];
  e_ref_right->get_solution(x_phys_right, pts_num, phys_u_ref_right, phys_dudx_ref_right, 
                            y_prev_ref, this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int c=0; c < n_eq; c++) {
    for (int i=0; i < pts_num; i++) {
      fprintf(f[c], "%g %g\n", x_phys_right[i], phys_u_ref_right[c][i] - phys_u_right[c][i]);
    }
    fprintf(f[c], "\n");
  }

}

// Plots the error wrt. the exact solution (if available)
void Mesh::plot_element_error_exact(FILE *f[MAX_EQN_NUM], Element *e, 
                                    double *y_prev, exact_sol_type exact_sol, int subdivision)
{
  int n_eq = this->get_n_eq();
  int pts_num = subdivision + 1;
  double x1 = e->x1;
  double x2 = e->x2;

  // calculate point array
  double x_phys[MAX_PTS_NUM];  
  double h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u[MAX_EQN_NUM][MAX_PTS_NUM], phys_dudx[MAX_EQN_NUM][MAX_PTS_NUM];
  e->get_solution(x_phys, pts_num, phys_u, phys_dudx, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int i=0; i < pts_num; i++) {
    double exact_sol_point[MAX_EQN_NUM];
    double exact_der_point[MAX_EQN_NUM];
    exact_sol(x_phys[i], exact_sol_point, exact_der_point);
    for (int c=0; c < n_eq; c++) {
      fprintf(f[c], "%g %g\n", x_phys[i], exact_sol_point[c] - phys_u[c][i]);
    }
  }
  for (int c=0; c < n_eq; c++) fprintf(f[c], "\n");
}

// Plots the error between the reference and coarse mesh solutions
void Mesh::plot_error_est(const char *filename, Mesh* mesh_ref, 
		      double* y_prev, double* y_prev_ref, int subdivision)
{
  int n_eq = this->get_n_eq();

  FILE *f_array[MAX_EQN_NUM];
  char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
  for(int c=0; c<n_eq; c++) {
    if(n_eq == 1)
        sprintf(final_filename[c], "%s", filename);
    else
        sprintf(final_filename[c], "%s_%d", filename, c);
    f_array[c] = fopen(final_filename[c], "wb");
    if(f_array[c] == NULL) error("problem opening file in plot_error_est().");
  }

  // simultaneous traversal of 'this' and 'mesh_ref'
  Element *e;
  Iterator *I = new Iterator(this);
  Iterator *I_ref = new Iterator(mesh_ref);
  while ((e = I->next_active_element()) != NULL) {
    Element *e_ref = I_ref->next_active_element();
    if (e->level == e_ref->level) { // element 'e' was not refined in space
                                    // for reference solution
      plot_element_error_p(f_array, e, e_ref, y_prev, y_prev_ref, subdivision);
    }
    else { // element 'e' was refined in space for reference solution
      Element* e_ref_left = e_ref;
      Element* e_ref_right = I_ref->next_active_element();
      plot_element_error_hp(f_array, e, e_ref_left, e_ref_right, 
                            y_prev, y_prev_ref, subdivision);
    }
  }

  for(int c=0; c<n_eq; c++) {
    fclose(f_array[c]);
    printf("Error function written to %s.\n", final_filename[c]);
  }
}

// Plots the error wrt. the exact solution (if available)
void Mesh::plot_error_exact(const char *filename,  
		            double* y_prev, exact_sol_type exact_sol, int subdivision)
{
  int n_eq = this->get_n_eq();

  FILE *f_array[MAX_EQN_NUM];
  char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
  for(int c=0; c<n_eq; c++) {
    if(n_eq == 1)
        sprintf(final_filename[c], "%s", filename);
    else
        sprintf(final_filename[c], "%s_%d", filename, c);
    f_array[c] = fopen(final_filename[c], "wb");
    if(f_array[c] == NULL) error("problem opening file in plot_error_exact().");
  }

  // traversal of 'this'
  Element *e;
  Iterator *I = new Iterator(this);
  while ((e = I->next_active_element()) != NULL) {
    plot_element_error_exact(f_array, e, y_prev, exact_sol, subdivision);
  }

  for(int c=0; c<n_eq; c++) {
    fclose(f_array[c]);
    printf("Exact solution error written to %s.\n", final_filename[c]);
  }
}

// Refine all elements in the id_array list whose id_array >= 0
void Mesh::adapt(double threshold, Mesh *mesh_ref, double *y_prev, double *y_prev_ref, 
                 double *err_squared_array) 
{
  // Find element with largest error
  double max_elem_error = 0;
  for(int i=0; i<this->get_n_active_elem(); i++) {
    double elem_error = sqrt(err_squared_array[i]);
    if (elem_error > max_elem_error) {
      max_elem_error = elem_error;
    }
  }

  // Create auxiliary array of element indices
  int id_array[MAX_ELEM_NUM];
  for(int i=0; i<this->get_n_active_elem(); i++) {
   if(sqrt(err_squared_array[i]) < threshold*max_elem_error) id_array[i] = -1; 
   else id_array[i] = i;
  }

  // Print elements to be refined
  printf("Elements to be refined:\n");
  for (int i=0; i<this->get_n_active_elem(); i++) {
    if (id_array[i] >= 0) printf("Elem[%d], error = %g\n", id_array[i], 
                                 sqrt(err_squared_array[i]));
  }

  int adapt_list[MAX_ELEM_NUM];
  int num_to_adapt = 0;

  // Create list of elements to be refined, in increasing order
  for (int i=0; i<this->get_n_active_elem(); i++) printf("id_array[%d] = %d\n", i, id_array[i]);
  for (int i=0; i<this->get_n_active_elem(); i++) {
    if (id_array[i] >= 0) {
      adapt_list[num_to_adapt] = id_array[i];
      num_to_adapt++;
    }
  }
 
  // Debug: Printing list of elements to be refined
  printf("refine_elements(): Elements to be refined:\n");
  for (int i=0; i<num_to_adapt; i++) printf("Elem[%d]\n", adapt_list[i]);

  Iterator *I = new Iterator(this);
  Iterator *I_ref = new Iterator(mesh_ref);

  // Simultaneous traversal of 'this' and 'mesh_ref'.
  // For each element, create a list of refinement 
  // candidates and have it checked. 
  Element *e = I->next_active_element();
  Element *e_ref = I_ref->next_active_element();
  int counter = 0;
  while (counter != num_to_adapt) {
    if (e->id == adapt_list[counter]) {
      counter++;
      int choice = 0;
      int3 cand_list[MAX_CAND_NUM];    // Every refinement candidates consists of three
                                      // numbers: 1/0 whether it is a p-candidate or not,
                                      // and then either one or two polynomial degrees
      // debug:
      //e->print_cand_list(num_cand, cand_list);
      if (e->level == e_ref->level) { // element 'e' was not refined in space
                                      // for reference solution
        int num_cand = e->create_cand_list(e_ref->p, -1, cand_list);
        choice = select_hp_refinement_ref_p(num_cand, cand_list, e, e_ref, y_prev_ref, 
                                            this->bc_left_dir_values,
	  		                    this->bc_right_dir_values);
      }
      else { // element 'e' was refined in space for reference solution
        Element* e_ref_left = e_ref;
        Element* e_ref_right = I_ref->next_active_element();
        int num_cand = e->create_cand_list(e_ref_left->p, e_ref_right->p, cand_list);
        choice = select_hp_refinement_ref_hp(num_cand, cand_list, e, e_ref_left, 
                                             e_ref_right, y_prev_ref, 
                                             this->bc_left_dir_values,
			                     this->bc_right_dir_values);
      }
      Element *e_last = e;
      e = I->next_active_element();
      e_ref = I_ref->next_active_element();
      // perform the refinement
      printf("Refining element (%g, %g), cand = (%d %d %d)\n", e_last->x1, e_last->x2, 
             cand_list[choice][0], cand_list[choice][1], cand_list[choice][2]);
      e_last->refine(cand_list[choice]);
      if(cand_list[choice][0] == 1) this->n_active_elem++; 
    }
    else {
      e = I->next_active_element();
      e_ref = I_ref->next_active_element();
      if (e->level != e_ref->level) { 
        e_ref = I_ref->next_active_element();
      }    
    }
  }
}

