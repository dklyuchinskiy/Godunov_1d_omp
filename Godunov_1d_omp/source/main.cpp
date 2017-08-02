#include "definitions.h"
#include "support.h"


/* Initial condition */
double time_max;
double explosion;


/********************************************/

double inflections[NUM_ITER] = { 0 };

int *sw1_num_p, *sw2_num_p, *sw3_num_p;
int *sw1_num_u, *sw2_num_u, *sw3_num_u;
int *sw1_num_r, *sw2_num_r, *sw3_num_r;

double F_ro_I[NUM_ITER] = { 0 }, F_ro_S[NUM_ITER] = { 0 }, F_ro_M[NUM_ITER] = { 0 }, F_ro_E[NUM_ITER] = { 0 };
double ITER_TIME[NUM_ITER] = { 0 };
double delta_D[NUM_ITER] = { 0 }, delta_U[NUM_ITER] = { 0 }, delta_P[NUM_ITER] = { 0 };
double l1_D[NUM_ITER] = { 0 }, l1_U[NUM_ITER] = { 0 }, l1_P[NUM_ITER] = { 0 };

double LOOP_TIME[LOOPS][OMP_CORES] = { 0 };

void iteration(int numb)
{
	double *R, *R1, *R2,    // density
		*P, *P1, *P2,		// pressure
		*U, *U1, *U2,		// velocity
		*RS, *S, *S_diff,   // entropy
		*RU, *RU1, *RU2,    // moment of impulse
		*RE, *RE1, *RE2;	// total energy

	double *FR,		// density flux
	       *FRU,	// moment flux
		   *FRP,	// pressure gradient
		   *FRE,	// energy flux
		   *FRUS;   // entropy flux
	double *UFLUX;	// velocity flux

	int i, numcells, start_print, jump_print, left_index, right_index, imesh, k;
	double	timer, tau, dx, len, x, x_NC, ds = 0, us = 0, ps = 0, es, cs = 0, ss, term, D_num, rp, pg, hh, *uss, *pss, *dss;

	double *x_layer, *x_layer_NC, *E;

	double *all_exact_P, *all_exact_U, *all_exact_R, *all_exact_RE, *all_exact_S;
	double *diff_riem_P, *diff_riem_U, *diff_riem_R, *diff_riem_RE, *diff_riem_S;

	double delta_ro, delta_u, delta_p;

	double* right_parts;

	int* w_num_p, *w_num_r, * w_num_u;

	char FileName[255], FileName2[255], FileName3[255], FileName4[255];

	FILE *fout, *fmesh;
	FILE *fout2, *fout3, *file;
	FILE* fout_NC;
	FILE* fout4;

	double sum_m[4][4] = { 0 };

	int check1 = 0, check2 = 0, check3 = 0, check5[5] = { 0 };

	double riemann_U, riemann_P, riemann_D[2];

	double *diff; //память под массив выделяется внутри функции

	double dtdx; // dt / dx

	/* Set number of cells */
	numcells = nmesh[numb];	//	N = 100 * 3^K = 100 * 3^(NUM_MESH)
	start_print = nprnt[numb];	/*	(3^K-1)/2
								K = 0	-> 0
								K = 1	-> 1
								K = 2	-> 4
								K = 3	-> 13
								K = 4	-> 40
								K = 5	-> 121
								*/
	int numcells3 = numcells / 3;
	jump_print = numcells / 100;		//	N / 100

	/* Boundary */
	if (PROBLEM == 3)	// periodic
	{
		left_index = numcells - 2;
		right_index = 1;
	}
	else				// no periodic
	{
		left_index = 1;
		right_index = numcells - 2;
	}

	/* Domain */
	len = LENGTH;
	dx = len / double(numcells);	// step 
	x_layer_NC = new double[numcells];

	/* Create arrays */
	R = (double*)_mm_malloc(numcells*sizeof(double),32);
	P = (double*)_mm_malloc(numcells*sizeof(double), 32);
	U = (double*)_mm_malloc(numcells*sizeof(double), 32);
	S = (double*)_mm_malloc(numcells*sizeof(double), 32);
	E = (double*)_mm_malloc(numcells*sizeof(double), 32);
	RU = (double*)_mm_malloc(numcells*sizeof(double), 32);
	RE = (double*)_mm_malloc(numcells*sizeof(double), 32);
	RS = (double*)_mm_malloc(numcells*sizeof(double), 32);
	S_diff = (double*)_mm_malloc(numcells*sizeof(double), 32);
	
	all_exact_U = new double[numcells];
	all_exact_R = new double[numcells];
	all_exact_P = new double[numcells];
	all_exact_RE = new double[numcells];
	all_exact_S = new double[numcells];
	right_parts = new double[numcells];

	diff_riem_P = new double[numcells];
	diff_riem_R = new double[numcells];
	diff_riem_U = new double[numcells];
	diff_riem_RE = new double[numcells];
	diff_riem_S = new double[numcells];

	w_num_p = new int[N_smooth]; w_num_r = new int[N_smooth]; w_num_u = new int[N_smooth];


	sw1_num_p = new int[N_smooth]; sw2_num_p = new int[N_smooth]; sw3_num_p = new int[N_smooth];
	sw1_num_r = new int[N_smooth]; sw2_num_r = new int[N_smooth]; sw3_num_r = new int[N_smooth];
	sw1_num_u = new int[N_smooth]; sw2_num_u = new int[N_smooth]; sw3_num_u = new int[N_smooth];

	FR = (double*)_mm_malloc((numcells+1)*sizeof(double), 32);
	FRU = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	FRP = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	FRE = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	UFLUX = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	dss = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	uss = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	pss = (double*)_mm_malloc((numcells + 1)*sizeof(double), 32);
	x_layer = (double*)_mm_malloc((numcells)*sizeof(double), 32);

	double max_x2 = 0;
	double pg_max = 0;

	double es_diff;

	double D_analit = 0;

	double l1, r1, l2, r2;

	int omp_chunk = numcells / OMP_CORES;

	w_num_p[0:N_smooth] = 0; w_num_r[0:N_smooth] = 0;  w_num_u[0:N_smooth] = 0;

	int count2 = 0;

	w_num_p[0:N_smooth] = 0;
	w_num_r[0:N_smooth] = 0;
	w_num_u[0:N_smooth] = 0;

	sw1_num_p[0:N_smooth] = 0;
	sw2_num_p[0:N_smooth] = 0;
	sw3_num_p[0:N_smooth] = 0;

	sw1_num_u[0:N_smooth] = 0;
	sw2_num_u[0:N_smooth] = 0;
	sw3_num_u[0:N_smooth] = 0;

	sw1_num_r[0:N_smooth] = 0;
	sw2_num_r[0:N_smooth] = 0;
	sw3_num_r[0:N_smooth] = 0;

	int c_c = 0;

	double wtime;

	double u1 = 0, u2 = 0, u3 = 0, u_loc = 0, u_max = 0;

	double loop_full[LOOPS] = { 0 };

	/*********************************THE PROGRAM BELOW******************************/
	sprintf(FileName2,"first step analysis_%c.dat", TYPE);
	file = fopen(FileName2,"w");

	printf("\nIteration %d\n", numb + 1);

	clock_t start_t, end_t;
	clock_t loop_time;
	start_t = clock();

	/* Initial conditions */

	printf("Loop 0\n");
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
	{
		wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static)
		for (int i = 0; i < numcells; i++)
		{
			x_layer[i] = i*dx + 0.5*dx;          // that are middles of cells
			R[i] = initial_density(x_layer[i]);
			P[i] = initial_pressure(x_layer[i]);
			U[i] = initial_velocity(x_layer[i]);
		}
		wtime = omp_get_wtime()- wtime;
		LOOP_TIME[0][omp_get_thread_num()] += wtime;
		printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[0][omp_get_thread_num()]);
	}


	/*********************** Boundary conditions ******************************/

	boundary_conditions(numcells, dss, uss, pss, R, U, P, FR, FRU, FRE);
	
	/*********************** Boundary conditions ******************************/
	printf("Loop 1\n");
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
	{
		wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static)
		for (int i = 0; i < numcells; i++)
		{
			RU[i] = R[i] * U[i];
			RE[i] = P[i] / (GAMMA - 1.0) + 0.5 * R[i] * U[i] * U[i]; // full mechanic energy
		}
		wtime = omp_get_wtime() - wtime;
		LOOP_TIME[1][omp_get_thread_num()] += wtime;
		printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[1][omp_get_thread_num()]);
	}


#if PROBLEM==0
	D_analit = (R[numcells - 1] * U[numcells - 1] - R[0] * U[0]) / (R[numcells - 1] - R[0]);
	printf("Analitical speed: %10.8lf\n\n", D_analit);
#elif PROBLEM == 12
	double uc_left, uc_right;
	
//	uc_left = U[0] +0.43 - sqrt(GAMMA*P[0] / R[0]);
	//uc_right = U[numcells] + 0.43  - sqrt(GAMMA*P[numcells] / R[numcells]);
    uc_left = U[0] - sqrt(GAMMA*P[0] / R[0]);
	uc_right = U[numcells-5]  - sqrt(GAMMA*P[numcells-5] / R[numcells-5]);
	printf("U-C left: %8.6lf\nU-C right: %8.6lf\nmiddle: %8.6lf\n",uc_left,uc_right, (uc_left+uc_right)/2);
	system("pause");

#elif PROBLEM==2
	double ro_right = 1.271413930046081;
	double ro_left = 1.0;
	double u_right = 0.292868067614595;
	double u_left = 0;

	D_analit = (ro_right * u_right - ro_left * u_left) / (ro_right - ro_left);
//	printf("Analitical speed: %10.8lf\n\n", D_analit);
//	system("pause");

	
#elif (PROBLEM==8 || PROBLEM==4)
/*	printf("-----------OLD WAY----------\n");
	printf("\n--------first sw-------\n");
	printf("dens: %lf\n", gyugonio(st_P3, st_R3, st_P2));
	printf("sw speed: %lf\n", sw_speed(st_R3, st_R2, st_U3, st_U2));
	printf("u after sw: %lf\n", st_U2);
	printf("\n---------second sw---------\n");
	printf("dens: %lf\n", gyugonio(st_P2, st_R2, st_P1));
	printf("sw speed: %lf\n", sw_speed(st_R2, st_R1, st_U2, st_U1));
	

	printf("\n---------summa---------\n");
	printf("sw speed: %lf\n", sw_speed(st_R3, st_R1, st_U3, st_U1));
	D_analit = sw_speed(st_R3, st_R1, st_U3, st_U1);*/
	//D_analit = 3.48;

	printf("---------------NEW WAY-----------\n");
	printf("\n--------first sw-------\n");
	//printf("dens after sw: %lf\n", gyugonio(st_P3, st_R3, st_P2));
	double D_sw1 = sw_speed2(st_R3, st_U3, st_P3, gyugonio(st_P3, st_R3, st_P2), st_P2);
	printf("sw speed1 without u3: %lf\n", D_sw1);
	printf("rup after sw1: %lf, %lf, %lf\n", st_R2, st_U2, st_P2);
	system("pause");
	printf("\n--------second sw-------\n");
	//printf("dens after sw: %lf\n", st_R1);
	double D_sw2 = sw_speed2(st_R2, st_U2, st_P2, st_R1, st_P1);
	printf("sw speed2 without u3: %lf\n", D_sw2);
	printf("rup after sw2: %lf, %lf, %lf\n", st_R1, st_U1, st_P1);
	system("pause");
#endif

	timer = 0.0;

	time_max = time_max_array[PROBLEM];


	/***** Вычислительный цикл метода Годунова *****/

	int last = 0;

#ifdef FLUX_COUNT
	FILE* array_flux[N_bound];
	for (int i = 0; i < N_bound; i++)
	{
		sprintf(FileName4, "FLUXES_%c_%d_P%d_N%d.dat", TYPE, i, int(PROBLEM), int(numcells));
		array_flux[i] = fopen(FileName4, "w");
	}
#endif

	
	while (timer < time_max)
		/************************************* The beggining of iteration*****************************/
	{
		c_c++;

#pragma omp parallel firstprivate(u1,u2,u3,u_loc) shared(u_max) num_threads(OMP_CORES)
		{
#pragma omp for schedule(static, omp_chunk) nowait
			for (int i = 0; i < numcells; i++)// not parallel and vectorized: dependences
			{
				u1 = U[i] + sqrt(GAMMA*P[i] / R[i]);
				u2 = U[i] - sqrt(GAMMA*P[i] / R[i]);
				u3 = U[i];

				if (u1 > u_loc) u_loc = u1;
				if (u2 > u_loc) u_loc = u2;
				if (u3 > u_loc) u_loc = u3;
			}
#pragma omp critical
			if (u_loc > u_max) u_max = u_loc;
		}

#ifdef CFL_SWITCH
		if (timer < time_max / 2)
		{
			tau = CFL08*dx / u_max;
		}
		else
		{
			tau = CFL04*dx / u_max;
		}
#else
		tau = CFL04*dx / u_max;
#endif

		if (timer + tau > time_max)
		{
			tau = time_max - timer; // if the last time's step is bigger than distance to "time_max"
			last = 1;
		}

		dtdx = tau / dx;


		if (timer < CROSS_POINT)
		{
			//printf("NONLINEAR\n");
			nonlinear_solver(numcells, R, U, P, dss, uss, pss);
		}
		else
		{
			//printf("LINEAR\n");
			if (last) printf("Loop 2\n");
			loop_time = clock();

			linear_solver(numcells, R, U, P, dss, uss, pss, last);

			loop_full[2] += (double)(clock() - loop_time) / CLOCKS_PER_SEC;
			if (last) printf("Full time loop: %lf\n", loop_full[2]);
		}

		// computation of conservative variables

		if (last) printf("Loop 3\n");
		loop_time = clock();
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
		{
			wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static) 
			for (int i = 1; i < numcells; i++)
			{
				UFLUX[i] = uss[i];
				FR[i] = dss[i] * uss[i];
				FRU[i] = dss[i] * uss[i] * uss[i] + pss[i];
				FRE[i] = (pss[i] / (GAMMA - 1.0) + 0.5*dss[i] * uss[i] * uss[i])*uss[i] + pss[i] * uss[i];
			}
			wtime = omp_get_wtime() - wtime;
			LOOP_TIME[3][omp_get_thread_num()] += wtime;
			if (last) printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[3][omp_get_thread_num()]);
		}
		loop_full[3] += (double)(clock() - loop_time) / CLOCKS_PER_SEC;
		if (last) printf("Full time loop: %lf\n", loop_full[3]);


		/****************************************************************
		IMPORTANT: Solving of this problem is in conservative variables!
		They are mass (R), impulse (RU), energy (RE)
		****************************************************************/

		/* Computations of conservations laws*/

		if (last) printf("Loop 4\n");
		loop_time = clock();

#pragma omp parallel private(wtime) num_threads(OMP_CORES)
		{
			wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static) // vectorized and auto_parallelized
			for (int i = 0; i < numcells; i++)
			{
				R[i] = R[i] - dtdx * (FR[i + 1] - FR[i]);
				RU[i] = RU[i] - dtdx * (FRU[i + 1] - FRU[i]);
				RE[i] = RE[i] - dtdx * (FRE[i + 1] - FRE[i]);
			}
			wtime = omp_get_wtime() - wtime;
			LOOP_TIME[4][omp_get_thread_num()] += wtime;
			if (last) printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[4][omp_get_thread_num()]);
		}
		loop_full[4] += (double)(clock() - loop_time)/ CLOCKS_PER_SEC;
		if (last) printf("Full time loop: %lf\n", loop_full[4]);

		if (last) printf("Loop 5\n");
		loop_time = clock();
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
		{
			wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static)
			for (int i = 0; i < numcells; i++)    // over-computation of velocity and pressure for the next step ( n --> n+1 ) - as R, RU, RE (n+1 step)
			{
				U[i] = RU[i] / R[i]; //используется дальше в программе!!!!
				P[i] = (GAMMA - 1.0) * (RE[i] - 0.5 * RU[i] * U[i]);
				S[i] = P[i] / pow(R[i], GAMMA);
			}
			wtime = omp_get_wtime() - wtime;
			LOOP_TIME[5][omp_get_thread_num()] += wtime;
			if (last) printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[5][omp_get_thread_num()]);
		}
		loop_full[5] += (double)(clock() - loop_time) / CLOCKS_PER_SEC;
		if (last) printf("Full time loop: %lf\n", loop_full[5]);

#if PROBLEM==0

		/*******************difference analit and numeric solutions************/
		/**********************shock wave*************************************/

		/*	analitical_SW(numcells, initial_pressure(0.05), initial_density(0.05), initial_velocity(0.05), initial_pressure(0.2), initial_density(0.2), initial_velocity(0.2), shw_analit_p, shw_analit_u, shw_analit_d, timer);
		for (int j = 0; j < numcells; j++)
		{
		shw_diff_p[j] = P[j] - shw_analit_p[j];
		shw_diff_u[j] = U[j] - shw_analit_u[j];
		shw_diff_d[j] = R[j] - shw_analit_d[j];
		}*/
#endif

		//************ CURRENT DOMAIN **********************//

	
#if (PROBLEM==17)
		const double xx1 = 0.0, xx2 = 10.0, tt1 = 2.8;  //0.05. 0.25. 0.4 - норм
#else
		const double xx1 = 0.4, xx2 = 0.95, tt1 = 0.35;  //0.05. 0.25. 0.4 - норм
#endif

#if (PROBLEM==0 || PROBLEM==1)
		const double tt2 = 0.495;  // берем не самую верхнюю границу (0.5), так как она не считается!!! она равна нулю по моему алгоритму реализации
#elif (PROBLEM==17)
		const double tt2 = 3.4;
#else
		const double tt2 = 0.295;
#endif

		//************ CURRENT DOMAIN **********************//

		//analitical_riemann_modeling(numcells, initial_density(0.1), initial_velocity(0.1), initial_pressure(0.1), initial_density(0.9), initial_velocity(0.9), initial_pressure(0.9), timer, all_exact_R, all_exact_U, all_exact_P);

#ifdef DIFF_ANALIT_RIEMANN
		for (int i = 0; i < numcells; i++)
		{
			all_exact_RE[i] = (all_exact_P[i] / (GAMMA - 1.0) + 0.5*all_exact_R[i] * all_exact_U[i] * all_exact_U[i])*all_exact_U[i] + all_exact_P[i] * all_exact_U[i];
			all_exact_S[i] = all_exact_P[i] / pow(all_exact_R[i], GAMMA);
		}

		diff_riem_R[0:numcells] = R[0:numcells] - all_exact_R[0:numcells];
		diff_riem_U[0:numcells] = U[0:numcells] - all_exact_U[0:numcells];
		diff_riem_P[0:numcells] = P[0:numcells] - all_exact_P[0:numcells];
		diff_riem_RE[0:numcells] = RE[0:numcells] - all_exact_RE[0:numcells];
		diff_riem_S[0:numcells] = S[0:numcells] - all_exact_S[0:numcells];
#endif
		timer += tau;

#ifdef DEBUG
		if (c_c <= 2)
		{
			first_step_validation(file, numcells, c_c, timer, R, U, P, dss, uss, pss);
		}
#endif

		/************** cчет интегралов по контуру ****************/
#ifdef INTEGRAL
#ifndef DIFF_ANALIT_RIEMANN
		outline_integral_riemann(numcells, timer, tau, tt1, tt2, xx1, xx2, R, U, P, RE, S, sum_m); //numerical solution on timer with new P U R 
#else
		outline_integral_riemann(numcells, timer, tau, tt1, tt2, xx1, xx2, diff_riem_R, diff_riem_U, diff_riem_P, diff_riem_RE, diff_riem_S, sum_m); //difference num and exact
#endif
#endif

		/*************** cчет интегралов по контуру ****************/

#ifdef OUTPUT_N_SMOOTH
		gnuplot_n_smooth_steps(numcells, timer, tau, R, U, P, S_diff);
#endif

#ifdef FLUX_COUNT
		/************** расчет движения потоков ********************/
		double t[N_bound] = { 0 };
		int t_ind[N_bound] = { 0 };
		int numcells_flux;
		numcells_flux = numcells;

		for (int i = 0; i < N_bound; i++)
		{
			t_ind[i] =  i * numcells_flux / N_bound;
			t[i] = (t_ind[i] + 0.5)*dx + UFLUX[t_ind[i]] * timer;
			fprintf(array_flux[i], "%lf %lf %lf\n", t[i], timer, UFLUX[t_ind[i]]);
		}
#endif

	} /******************************************* The end of iteration**************************/

#ifdef FLUX_COUNT
	for (int i = 0; i < N_bound; i++)
		fclose(array_flux[i]);

#endif

#ifdef OUTPUT_N_SMOOTH
#ifdef PRINT
#if(PROBLEM==2 || PROBLEM==9)
	gnuplot_analitical_riemann2(numcells, w_num_r, w_num_u, w_num_p);
#else
#ifdef NC2
	gnuplot_n_smooth_NC2(numcells, w_num_r, w_num_u, w_num_p);
#endif
	gnuplot_n_smooth2(numcells, sw1_num_r, sw1_num_u, sw1_num_p, sw2_num_r, sw2_num_u, sw2_num_p, sw3_num_r, sw3_num_u, sw3_num_p);
#if (PROBLEM == 10)
	gnuplot_n_smooth3(numcells, sw1_num_r, sw1_num_u, sw1_num_p, sw2_num_r, sw2_num_u, sw2_num_p, sw3_num_r, sw3_num_u, sw3_num_p);
#endif
#endif
#endif
#endif

#if PROBLEM==2
#ifdef PRINT
	gnuplot_n_smooth2(numcells, sw1_num_r, sw1_num_u, sw1_num_p, sw2_num_r, sw2_num_u, sw2_num_p, sw3_num_r, sw3_num_u, sw3_num_p);
#endif
#endif

#ifdef DIFF_ANALIT_RIEMANN
	free(diff_riem_P);
	free(diff_riem_R);
	free(diff_riem_U);
	free(diff_riem_RE);
	free(diff_riem_S);
#endif

	
	_mm_free(FR);
	_mm_free(FRU);
	_mm_free(FRP);
	_mm_free(FRE);
	_mm_free(UFLUX);
	_mm_free(dss);
	_mm_free(uss);
	_mm_free(pss);



/*#if(PROBLEM==2)
	analitical_riemann_modeling(numcells, initial_density(0.0), initial_velocity(0.0), initial_pressure(0.0), initial_density(0.9), initial_velocity(0.9), initial_pressure(0.9), time_max_array[PROBLEM], all_exact_R, all_exact_U, all_exact_P);
#ifdef PRINT	
	gnuplot_analitical_riemann(numcells, R, U, P, all_exact_R, all_exact_U, all_exact_P);
#endif
#endif*/

	/*********Итог счета интегралов по контуру**********/
#ifdef DIFF_ANALIT_RIEMANN
#ifdef L1-NORM
	difference_analitical_riemann_Linf(numb, R, U, P, all_exact_R, all_exact_U, all_exact_P, delta_ro, delta_u, delta_p);
	delta_D[numb] = delta_ro;
	delta_U[numb] = delta_u;
	delta_P[numb] = delta_p;
	
	difference_analitical_riemann_L1(numb, R, U, P, all_exact_R, all_exact_U, all_exact_P, delta_ro, delta_u, delta_p);
	l1_D[numb] = delta_ro;
	l1_U[numb] = delta_u;
	l1_P[numb] = delta_p;

	free(all_exact_U);
	free(all_exact_R);
	free(all_exact_P);
	free(all_exact_RE);
	free(all_exact_S);
#endif
#endif

	/*********Итог счета интегралов по контуру**********/
#ifdef INTEGRAL
	F_ro_I[numb] = sum_m[0][1] - sum_m[1][1] + sum_m[2][1] - sum_m[3][1];
	printf("sum_l_I: %30.28lf\nsum_t_I: %30.28lf\nsum_r_I: %30.28lf\nsum_b_I: %30.28lf\n", sum_m[3][1],sum_m[0][1],sum_m[2][1], sum_m[1][1]);
	F_ro_M[numb] = sum_m[0][0] - sum_m[1][0] + sum_m[2][0] - sum_m[3][0];
	printf("\nsum_l_M: %30.28lf\nsum_t_M: %30.28lf\nsum_r_M: %30.28lf\nsum_b_M: %30.28lf\n", sum_m[3][0] , sum_m[0][0] , sum_m[2][0] , sum_m[1][0]);
	F_ro_S[numb] = sum_m[0][2] - sum_m[1][2] + sum_m[2][2] - sum_m[3][2];
	printf("\nsum_l_S: %30.28lf\nsum_t_S: %30.28lf\nsum_r_S: %30.28lf\nsum_b_S: %30.28lf\n", sum_m[3][2] , sum_m[0][2] , sum_m[2][2] , sum_m[1][2]);
	F_ro_E[numb] = sum_m[0][3] - sum_m[1][3] + sum_m[2][3] - sum_m[3][3];
#endif
	/*********Итог счета интегралов по контуру**********/

	end_t = clock();
	double duration = (double)(end_t - start_t) / CLOCKS_PER_SEC;

	printf("Iteration %d time: %lf\n", numb, duration);
	ITER_TIME[numb] = duration;

	//  ВЫВОД!!!
#ifndef OUTPUT_N_SMOOTH
#ifndef RP
	sprintf(FileName, "N%04d_P%1d_SLV%1d_TERM%.0lf.dat", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);
	fout = fopen(FileName, "w");
#else
	sprintf(FileName, "RP_N%04d_P%1d_SLV%1d_TERM%.0lf.dat", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);
	fout = fopen(FileName, "w");
#endif


	//system("pause");
	/*******************OUTPUT***************************************/



#ifdef FIRST
	for (i = start_print; i < numcells; i += jump_print)  // вывод всегда по 100 точек ( с первой итерации которые )
#else
	for (i = 0; i < numcells; i++)  // вывод всегда по 100 точек ( с первой итерации которые )
#endif
	{
		x = i*dx + 0.5*dx;
#if (PROBLEM==0)
		x_NC = (i*dx + 0.5*dx) - D_analit*time_max_array[0];  //NC
#elif (PROBLEM==1)
		x_NC = (i*dx + 0.5*dx - 0.1) / (time_max_array[0] * (numb + 1));  //NC
#endif

		/********************************************************************************************************************************
		| 1 |    2    |   3   |     4    |          5        |         6       |      7      |          8        |      9       |   10  |
		| x | density | speed | pressure | velocity of sound |       entropy   |      term    |	                 |   pg/pg_max  | h_pow |
		********************************************************************************************************************************/
		ds = R[i];
		us = U[i];
		ps = P[i];
		cs = sqrt(GAMMA*P[i] / R[i]);
		ss = P[i] / pow(R[i], GAMMA);


#ifndef NC
		fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf \n", x, ds, us, ps, cs, ss);
#else 
		//fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf %30.24lf %30.24lf %30.24lf \n", x_NC, ds, us, ps, cs, es, rp, pg, hh);
		fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf\n", x_layer_NC[i], ds, us, ps, cs, ss);
#endif

	}

#endif
	/************************** Analysis after count on one iteration *********************/

#if (PROBLEM==0)
#ifdef SW_FINITE_DIFF
	printf("************Pressure*****************\n");
	diff = finite_difference(numb, P);

	printf("************Right parts*****************\n");
	diff = finite_difference(numb, pressure_gradient);
#endif
#endif

#if (PROBLEM==1)
#ifdef RW_NUM_ANALITICAL
	   rw_diff_num_analit(numb, numcells, R, U, P);
#endif
#endif

	_mm_free(R);
	_mm_free(P);
	_mm_free(U);
	_mm_free(S_diff);
	_mm_free(S);
	_mm_free(RU);
	_mm_free(RE);
	_mm_free(RS);

	fclose(file);

#ifdef PRINT	
#ifndef OUTPUT_N_SMOOTH
	fclose(fout);
#endif
#endif

}



/************************************************************************************/
/**********************************  MAIN  ******************************************/
/************************************************************************************/

int main()
{
	clock_t start;

	int i = 0;
	int run[NUM_ITER];
#ifdef INTEGRAL
	run[0:NUM_ITER] = 1;
#else
	run[0] = 0;
	run[1] = 1;
	run[2] = 0;
	run[3] = 0;
	run[4] = 0;
	run[5] = 0;
	run[6] = 0;
#endif

#if 1
	
#ifdef NC
	NC;
#endif

#ifdef FIRST
	FIRST;
#else
	SECOND;
#endif

#ifdef PRINT
	PRINT;
#endif

	start = clock(); // start of program's work
	for (i = 0; i < NUM_ITER; i++) // Iterations   //there is dependence between iterations!!! its impossible to start new iteration before last ends
	{
#ifdef OTKOL
		if (run[i] == 1) iteration_bound(i);
#else
		if (run[i] == 1) iteration(i);
#endif
		printf("Printing...\n");

		/******************Different settings and tasks***************************/

#ifdef PRINT
#ifdef SIMPLE
		gnuplot_one_iteration(nmesh[i]);
#endif
#ifdef OUTPUT_N_SMOOTH
#ifndef NC
	//	gnuplot_n_smooth(i);
#else
	//	gnuplot_n_smooth_NC(i);
#endif
#endif

#ifdef RW_NUM_ANALITICAL
#if(PROBLEM==1)
		//	gnuplot_RW_DIFF(nmesh[i]);
		gnuplot_RW_NUM_ANALITIC(nmesh[i]);
#endif


#endif

#if (PROBLEM==2)
	//	gnuplot_conservative(i);

#endif
#if (PROBLEM==3)
#ifdef FIVE_T_STEPS
		gnuplot_five_t_steps(i);
#endif
#endif
#endif
	}

	double duration;
	duration = (double)(clock() - start) / CLOCKS_PER_SEC;

#ifdef PRINT
#ifdef NC
	gnuplot_one_it_NC();  // i - current number of iterations
#endif
#endif

	/******************Different settings and tasks***************************/

#ifdef INTEGRAL
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_M: %30.28lf\n", F_ro_M[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_I: %30.28lf\n", F_ro_I[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_S: %30.28f\n", F_ro_S[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_E: %30.28lf\n", F_ro_E[i]);
	}
	printf("\n");
	printf("Massa\n");
	runge(F_ro_M);
	printf("Impulse\n");
	runge(F_ro_I);
	printf("Energy\n");
	runge(F_ro_E);
	printf("Entropy\n");
	runge(F_ro_S);
#endif
	printf("\n**************************\n");
	for (int i=0; i<NUM_ITER;i++)
		printf("Iter %d time: %lf\n",i + 1,ITER_TIME[i]);

#ifdef DIFF_ANALIT_RIEMANN
	printf("----------------------\nL infinity norm\n----------------------\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("delta_ro: %30.28lf\n", delta_D[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("delta_u: %30.28lf\n", delta_U[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("delta_p: %30.28f\n", delta_P[i]);
	}
	printf("\n");

	printf("\n");
	printf("Density in rarify wave\n");
	runge(delta_D);
	printf("Velocity in rarify wave\n");
	runge(delta_U);
	printf("Pressure in rarify wave\n");
	runge(delta_P);

	printf("----------------------\nL1 norm\n----------------------\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("l1_ro: %30.28lf\n", l1_D[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("l1_u: %30.28lf\n", l1_U[i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("l1_p: %30.28f\n", l1_P[i]);
	}
	printf("\n");

	printf("\n");
	printf("Density in rarify wave\n");
	runge(l1_D);
	printf("Velocity in rarify wave\n");
	runge(l1_U);
	printf("Pressure in rarify wave\n");
	runge(l1_P);
#endif

	printf("\nElapsed time: %lf sec\n", duration);
#endif

	system("pause");
	return 0;
}
