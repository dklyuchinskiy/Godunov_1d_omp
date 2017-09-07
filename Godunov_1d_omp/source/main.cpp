#include "definitions.h"
#include "support.h"

int main()
{
	clock_t start;

	int i = 0;
	int run[NUM_ITER];

	// M, I, E, S
	double F_ro[4*NUM_ITER] = { 0 };
	double ITER_TIME[NUM_ITER] = { 0 };

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

	start = clock(); // start of program's work
	for (i = 0; i < NUM_ITER; i++) // Iterations   //there is dependence between iterations!!! its impossible to start new iteration before last ends
	{
#ifdef OTKOL
		if (run[i] == 1) iteration_bound(i);
#else
		if (run[i] == 1) iteration(i, F_ro, ITER_TIME);
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
	printf("good\n");
	/******************Different settings and tasks***************************/
	int ldf = 4;
#ifdef INTEGRAL
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_M: %30.28lf\n", F_ro[0 + ldf * i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_I: %30.28lf\n", F_ro[1 + ldf * i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_E: %30.28lf\n", F_ro[2 + ldf * i]);
	}
	printf("\n");
	for (i = 0; i < NUM_ITER; i++)
	{
		printf("F_ro_S: %30.28lf\n", F_ro[3 + ldf * i]);
	}
	printf("\n");
	printf("Massa\n");
	runge(F_ro, ldf, 0);
	printf("Impulse\n");
	runge(F_ro, ldf, 1);
	printf("Energy\n");
	runge(F_ro, ldf, 2);
	printf("Entropy\n");
	runge(F_ro, ldf, 3);
#endif
	printf("\n**************************\n");
	for (int i = 0; i < NUM_ITER; i++)
		printf("Iter %d time: %lf\n", i + 1, ITER_TIME[i]);

	printf("\nElapsed time: %lf sec\n", duration);
#endif

	system("pause");
	return 0;
}
