#pragma once

void sample(double &pm, double &um, double &s, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr, double &d, double &u, double &p);
void prefun(double &f, double &fd, double &p, double &dk, double &pk, double &ck);
double guessp(double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr);
void starpu(double &p, double &u, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr);

// Functions

double initial_density(double x);
double initial_pressure(double x);
double initial_velocity(double x);
void linear(double dl, double ul, double pl, double dr, double ur, double pr, double &d, double &u, double &p);
void linear_check(double dl, double ul, double pl, double dr, double ur, double pr, int &left, int &middle, int &right, int numb);
void runge(double *massiv);
double gyugonio(double p1, double ro1, double p2/*�� ������� ������*/);
double sw_speed2(double ro1, double u1, double p1, double ro2 /*�� ������� ������*/, double p2 /*�� ������� ������*/);
double sw_speed(double ro1, double ro2, double u1, double u2);
double* finite_difference(int numb, double *mas);

void analitical_RW(FILE* file_name, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r, double numb);
void analitical_SW(int numcells, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r, double *res_p, double *res_u, double* res_d, double timer);
void analitical_riemann(int numcells, double p1, double ro1, double u1, double p2, double ro2, double u2, double *sol_p, double *sol_u);
void analitical_riemann_modeling(int numcells, double ro1, double u1, double p1, double ro2, double u2, double p2, double timer, /*output*/double *all_d, double *all_u, double *all_p);
void analitical_writing_into_file(int numcells, double* R_D, double*R_U, double*R_P, double timer);
void difference_analitical_riemann_Linf(int numb, double* R, double*U, double*P, double* R_D, double*R_U, double*R_P, double &delta_ro, double &delta_u, double &delta_p);
void difference_analitical_riemann_L1(int numb, double* R, double*U, double*P, double* R_D, double*R_U, double*R_P, double &sum_ro, double &sum_u, double &sum_p);
double RW_prop(int digit, double x, double numb, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r);
void outline_integral_riemann(int numcells, double timer, double tau, const double tt1, const double tt2, double xx1, double xx2, double* R, double*U, double*P, double*RE, double*S, /*output*/ double sum[4][4]);
void gnuplot_one_iteration(int numcells);
void gnuplot_RW_DIFF(int numcells);
void gnuplot_RW_NUM_ANALITIC(int numcells);
void gnuplot_P_PLUS_PG(int numcells);
void gnuplot_all_iterations_NC(int numb);
void gnuplot_one_it_NC();
void gnuplot_conservative(int numb);
void gnuplot_five_t_steps(int numb);
void gnuplot_n_smooth_NC(int numb);
void gnuplot_n_smooth_NC2(int numcells, int* n_r, int* n_u, int* n_p);
void gnuplot_analitical_riemann(int numcells, double* R, double*U, double*P, double* R_D, double*R_U, double*R_P);
void gnuplot_n_smooth(int numb);
void gnuplot_n_smooth2(int numcells, int* sw1_r, int* sw1_u, int* sw1_p, int* sw2_r, int* sw2_u, int* sw2_p, int* sw3_r, int* sw3_u, int* sw3_p);
void gnuplot_n_smooth3(int numcells, int* sw1_r, int* sw1_u, int* sw1_p, int* sw2_r, int* sw2_u, int* sw2_p, int* sw3_r, int* sw3_u, int* sw3_p);
void gnuplot_analitical_riemann2(int numcells, int* n_r, int* n_u, int* n_p);
void nonlinear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss);
void linear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss, int last);
void iteration(int numb);
void iteration_bound(int numb);
void first_step_validation(int numcells, double c_c, double timer, double *R, double *U, double *P, double *UFLUX);
void boundary_conditions(int numcells, double *dss, double *uss, double *pss, double *R, double *U, double *P, double *FR, double *FRU, double *FRE);