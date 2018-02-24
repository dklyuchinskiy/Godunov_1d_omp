#include "definitions.h"
#include "support.h"

int main()
{
	double start;
	double duration;

	int i = 0;
	int run[NUM_ITER];

	// M, I, E, S
	double F_ro[4*NUM_ITER] = { 0 };
	double ITER_TIME[NUM_ITER] = { 0 };

#ifdef INTEGRAL
	for (int i = 0; i < NUM_ITER; i++)
		run[i] = 1;
#else
	run[0] = 0;
	run[1] = 0;
	run[2] = 0;
	run[3] = 0;
	run[4] = 0;
	run[5] = 1;
	run[6] = 0;
#endif

#if 1

	start = omp_get_wtime();
	for (i = 0; i < NUM_ITER; i++) // Iterations   //there is dependence between iterations!!! its impossible to start new iteration before last ends
	{
		if (run[i] == 1) iteration(i, F_ro, ITER_TIME);

#if (defined(PRINT) && defined(SIMPLE))
		gnuplot_one_iteration(nmesh[i]);
#endif
	}

	// P_0 : timing 0.35
	double time = 0.3450;
	// P_1 : timing 0.3500
	time = 0.0375;
	// P_2 : timing 0.0075, 0.015, 0.0225, 0.03, 0.075, 0.2250
	//time = 0.2250;

	// P_7: timing 0.1500
	time = 0.1500;
	// P_4: timing 
	time = 0.04;
	// gnuplot_all_iter_one_time(run, 4, time);

	duration = omp_get_wtime() - start;

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

	printf("\n************************************\n");
	for (int i = 0; i < NUM_ITER; i++)
		printf("Iter %d time: %lf\n", i + 1, ITER_TIME[i]);

	printf("\nElapsed time: %lf sec\n", duration);
#endif

	system("pause");
	return 0;
}
