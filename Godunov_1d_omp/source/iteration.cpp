#include "definitions.h"
#include "support.h"

void iteration(int numb, double F_ro[], double ITER_TIME[])
{
	double *R,				   // density
		   *P,				   // pressure
	       *U,  			   // velocity
     	   *RU,				   // moment of impulse
    	   *RE,			       // total energy
	       *S,				   // entropy
		   *S_diff,
		   *S_prev;  

	double *FR,		           // density flux
		   *FRU,	           // moment flux
		   *FRE,	           // energy flux
	       *UFLUX;	           // velocity flux

	double *uss, *pss, *dss;                                         // boundaries
	double *exact_R, *exact_U, *exact_P, *exact_RE, *exact_S;		 // exact values
	double *diff_R, *diff_U, *diff_P, *diff_RE, *diff_S;             // diff of analyt and numerical functions
	double *diff;									    			 // for this array the memory is allocated inside function 
	double *x_init, *x_layer;                                        // coordinates

	int numcells, start_print, jump_print, left_index, right_index, imesh;
	int iter = 0, count = 0, last = 0;
	int check1 = 0, check2 = 0, check3 = 0, check5[5] = { 0 };

	double timer, time_max, tau, dx, dtdx, len, x, x_NC, wtime, CFL;
	double ds = 0, us = 0, ps = 0, es = 0, es_diff = 0, cs = 0;
	double u1 = 0, u2 = 0, u3 = 0, u_loc = 0, u_max = 0;
	double l1, r1, l2, r2;	
	double D_analit = 0;
	double delta_ro, delta_u, delta_p;
	double cfl_time = 0;

	char FileName[255], FileName2[255], FileName3[255], FileName4[255];
	double sum_m[4][4] = { 0 };
	
	double loop_full[LOOPS] = { 0 };
	double LOOP_TIME[LOOPS][OMP_CORES] = { 0 };

	FILE *fout, *fout_NC, *fmesh, *file;

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

	int omp_chunk = numcells / OMP_CORES;

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

	/* Create arrays */
	mem_alloc(numcells, &R, 32);
	mem_alloc(numcells, &P, 32);
	mem_alloc(numcells, &U, 32);
	mem_alloc(numcells, &S, 32);
	mem_alloc(numcells, &RU, 32);
	mem_alloc(numcells, &RE, 32);
	mem_alloc(numcells, &S_diff, 32);
	mem_alloc(numcells, &S_prev, 32);
	mem_alloc(numcells, &S_diff, 32);

	mem_alloc(numcells, &x_layer, 32);
	mem_alloc(numcells, &x_init, 32);

	mem_alloc(numcells, &exact_R, 32);
	mem_alloc(numcells, &exact_U, 32);
	mem_alloc(numcells, &exact_P, 32);
	mem_alloc(numcells, &exact_S, 32);
	mem_alloc(numcells, &exact_RE, 32);

	mem_alloc(numcells, &diff_R, 32);
	mem_alloc(numcells, &diff_U, 32);
	mem_alloc(numcells, &diff_P, 32);
	mem_alloc(numcells, &diff_S, 32);
	mem_alloc(numcells, &diff_RE, 32);

	mem_alloc(numcells + 1, &FR, 32);
	mem_alloc(numcells + 1, &FRU, 32);
	mem_alloc(numcells + 1, &FRE, 32);
	mem_alloc(numcells + 1, &UFLUX, 32);
	mem_alloc(numcells + 1, &dss, 32);
	mem_alloc(numcells + 1, &uss, 32);
	mem_alloc(numcells + 1, &pss, 32);

	int w_num_p[N_smooth] = { 0 };
	int w_num_r[N_smooth] = { 0 };
	int w_num_u[N_smooth] = { 0 };

	int sw1_num[3][N_smooth] = { 0 };
	int sw2_num[3][N_smooth] = { 0 };
	int sw3_num[3][N_smooth] = { 0 };

	double t_flux[N_bound] = { 0 };

	/*********************************************************************/

#ifdef DEBUG
	sprintf(FileName2, "first step analysis_%c.dat", TYPE);
	file = fopen(FileName2, "w");
#endif
	
	printf("\nIteration %d\n", numb + 1);

	clock_t start_t, end_t;
	clock_t loop_time;
	start_t = clock();

	/* Mesh */
#pragma omp for simd schedule(simd:static)
	for (int i = 0; i < numcells; i++)
	{
		x_init[i] = i*dx + 0.5*dx;          // that are middles of cells
		x_layer[i] = x_init[i];
	}

	/* Initial conditions */
	printf("Loop 0\n");
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
	{
		wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static)
		for (int i = 0; i < numcells; i++)
		{
			R[i] = initial_density(x_init[i]);
			P[i] = initial_pressure(x_init[i]);
			U[i] = initial_velocity(x_init[i]);
			S_prev[i] = log(P[i] / pow(R[i], GAMMA));
		}
		wtime = omp_get_wtime() - wtime;
		LOOP_TIME[0][omp_get_thread_num()] += wtime;
		printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[0][omp_get_thread_num()]);
	}

	/* Computation of RU and RE */
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

	/*** To get some information before the start ***/
	inf_before_start(numcells, R, U, P, D_analit);


	/***** Computational loop *****/

#ifdef FLUX_COUNT
	FILE* array_flux[N_bound];
	for (int i = 0; i < N_bound; i++)
	{
		sprintf(FileName4, "FLUXES/FLUXES_%c_%d_P%d_N%d.dat", TYPE, i, int(PROBLEM), int(numcells));
		array_flux[i] = fopen(FileName4, "w+");
	}
#endif

	/* Main loop of computational algrithm */
	iter = 0;
	timer = 0.0;
	time_max = time_max_array[PROBLEM];

	while (timer < time_max)
	{
		iter++;
		
		loop_time = clock();
#pragma omp parallel firstprivate(u1,u2,u3,u_loc) shared(u_max) num_threads(OMP_CORES)
		{
			/* CFL condition */
#pragma omp for schedule(static, omp_chunk) nowait
			for (int i = 0; i < numcells; i++)
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
		cfl_time += (double)(clock()-loop_time) / CLOCKS_PER_SEC;
		
#ifdef CFL_SWITCH
		CFL = timer < time_max / 2 ? CFL08 : CFL04;
		tau = CFL*dx / u_max;
#else
		tau = CFL04*dx / u_max;
#endif

		if (timer + tau > time_max)
		{
			tau = time_max - timer; // if the last time's step is bigger than distance to "time_max"
			last = 1;
		}
		if (last) printf("CFL loop: %lf\n", cfl_time);
		dtdx = tau / dx;

		/* Boundary conditions */
		boundary_conditions(numcells, dss, uss, pss, R, U, P);

		/* Nonlinear or linear Godunov scheme */
		if (timer < CROSS_POINT)
		{
			nonlinear_solver(numcells, R, U, P, dss, uss, pss);
		}
		else
		{
			if (last) printf("Loop 2\n");
			loop_time = clock();

			linear_solver(numcells, R, U, P, dss, uss, pss, last);

			loop_full[2] += (double)(clock() - loop_time) / CLOCKS_PER_SEC;
			if (last) printf("Full time loop: %lf\n", loop_full[2]);
		}

		/* Computation of flux variables */
		if (last) printf("Loop 3\n");
		loop_time = clock();
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
		{
			wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static) 
			for (int i = 0; i <= numcells; i++)
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


		/* Computation of conservations laws (in conservative variables!) */
		if (last) printf("Loop 4\n");
		loop_time = clock();
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
		{
			wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static)
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
		loop_full[4] += (double)(clock() - loop_time) / CLOCKS_PER_SEC;
		if (last) printf("Full time loop: %lf\n", loop_full[4]);

		/* Over-computation of velocity, pressure and entropy */
		if (last) printf("Loop 5\n");
		loop_time = clock();
#pragma omp parallel private(wtime) num_threads(OMP_CORES)
		{
			wtime = omp_get_wtime();
#pragma omp for simd schedule(simd:static)
			for (int i = 0; i < numcells; i++) 
			{
				U[i] = RU[i] / R[i]; 
				P[i] = (GAMMA - 1.0) * (RE[i] - 0.5 * RU[i] * U[i]);
				S[i] = log(P[i] / pow(R[i], GAMMA));
				S_diff[i] = S[i] - S_prev[i];
				S_prev[i] = S[i];
			}
			wtime = omp_get_wtime() - wtime;
			LOOP_TIME[5][omp_get_thread_num()] += wtime;
			if (last) printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[5][omp_get_thread_num()]);
		}
		loop_full[5] += (double)(clock() - loop_time) / CLOCKS_PER_SEC;
		if (last) printf("Full time loop: %lf\n", loop_full[5]);

#if (PROBLEM==0)
#ifdef DIFF_ANALYT
		difference_SW(inumcells, timer, R, U, P, diff_R, diff_U, diff_P, exact_R, exact_U, exact_P);
#endif
#endif

		timer += tau;

#ifdef DEBUG
		if (iter <= 2) first_step_validation(file, numcells, iter, timer, R, U, P, dss, uss, pss);
#endif

		/* Checking of accuracy of the solution */
#ifdef INTEGRAL
#ifndef DIFF_ANALIT_RIEMANN
		outline_integral_riemann(numcells, timer, tau, (double)T1, (double)T2, (double)X1, (double)X2, x_init, R, U, P, RE, S, sum_m); //numerical solution on timer with new P U R 
#else
		outline_integral_riemann(numcells, timer, tau, (double)T1, (double)T2, (double)X1, (double)X2, x_init, diff_R, diff_U, diff_P, diff_RE, diff_S, sum_m); //difference num and exact
#endif
#endif
		/* Output to file during computations */
#if defined(OUTPUT_N_SMOOTH) && defined(PRINT)
#if (PROBLEM == 18)
		file_n_smooth_steps(numcells, timer, tau, x_layer, R, U, P, RE, S, S_diff, UFLUX);
#else
		file_n_smooth_steps(numcells, timer, tau, x_init, R, U, P, RE, S, S_diff, UFLUX);
#endif
#endif
		/* Euler coordinates */
		for (int i = 0; i < numcells; i++)
			x_layer[i] = x_layer[i] + tau * UFLUX[i];

		/* Moving of boundaries */
#ifdef FLUX_COUNT
		flux_count(array_flux, iter, numcells, timer, tau, t_flux, UFLUX);
#endif

	} /******************************************* The end of iteration**************************/

#ifdef FLUX_COUNT
	for (int i = 0; i < N_bound; i++)
		fclose(array_flux[i]);
#endif

	/* Analytical computations */
#if defined(DIFF_ANALIT_RIEMANN) && defined(L1_NORM)
	difference_analitical_riemann_Linf(numb, R, U, P, exact_R, exact_U, exact_P, delta_ro, delta_u, delta_p);
	delta_RUP[0][numb] = delta_ro;
	delta_RUP[1][numb] = delta_u;
	delta_RUP[2][numb] = delta_p;

	difference_analitical_riemann_L1(numb, R, U, P, exact_R, exact_U, exact_P, delta_ro, delta_u, delta_p);
	l1_RUP[0][numb] = delta_ro;
	l1_RUP[1][numb] = delta_u;
	l1_PUP[2][numb] = delta_p;
#endif

	/* Checking of accuracy of the solution */
	int ldf = 4;

	// Matrix [4 x numb]. Let's use column major.
#ifdef INTEGRAL
	F_ro[1 + ldf * numb] = sum_m[0][1] - sum_m[1][1] + sum_m[2][1] - sum_m[3][1];
	printf("sum_l_I: %30.28lf\nsum_t_I: %30.28lf\nsum_r_I: %30.28lf\nsum_b_I: %30.28lf\n", sum_m[3][1], sum_m[0][1], sum_m[2][1], sum_m[1][1]);
	F_ro[0 + ldf * numb] = sum_m[0][0] - sum_m[1][0] + sum_m[2][0] - sum_m[3][0];
	printf("\nsum_l_M: %30.28lf\nsum_t_M: %30.28lf\nsum_r_M: %30.28lf\nsum_b_M: %30.28lf\n", sum_m[3][0], sum_m[0][0], sum_m[2][0], sum_m[1][0]);
	F_ro[3 + ldf * numb] = sum_m[0][2] - sum_m[1][2] + sum_m[2][2] - sum_m[3][2];
	printf("\nsum_l_S: %30.28lf\nsum_t_S: %30.28lf\nsum_r_S: %30.28lf\nsum_b_S: %30.28lf\n", sum_m[3][2], sum_m[0][2], sum_m[2][2], sum_m[1][2]);
	F_ro[2 + ldf * numb] = sum_m[0][3] - sum_m[1][3] + sum_m[2][3] - sum_m[3][3];
#endif

	end_t = clock();
	double duration = (double)(end_t - start_t) / CLOCKS_PER_SEC;

	printf("Iteration %d time: %lf\n", numb + 1, duration);
	ITER_TIME[numb] = duration;


	/* Analysis after count on one iteration */
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

	/* The final output */
#ifdef PRINT
#ifdef OUTPUT_N_SMOOTH
	gnuplot_n_smooth2(numcells, sw1_num, sw2_num, sw3_num);
#else
	gnuplot_last_step(numcells, dx, D_analit, R, U, P);
#if (PROBLEM==2 || PROBLEM==9)
	gnuplot_analitical_riemann2(numcells, w_num_r, w_num_u, w_num_p);
#endif
#endif
#endif

	mem_free(&R);
	mem_free(&P);
	mem_free(&U);
	mem_free(&S_diff);
	mem_free(&S);
	mem_free(&RU);
	mem_free(&RE);

	mem_free(&FR);
	mem_free(&FRU);
	mem_free(&FRE);
	mem_free(&UFLUX);
	mem_free(&dss);
	mem_free(&uss);
	mem_free(&pss);

	mem_free(&exact_U);
	mem_free(&exact_R);
	mem_free(&exact_P);
	mem_free(&exact_RE);
	mem_free(&exact_S);

	mem_free(&diff_U);
	mem_free(&diff_R);
	mem_free(&diff_P);
	mem_free(&diff_RE);
	mem_free(&diff_S);

	mem_free(&x_layer);
	mem_free(&x_init);

#ifdef DEBUG
	fclose(file);
#endif
}