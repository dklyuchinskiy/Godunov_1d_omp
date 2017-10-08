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
	for (int i = 0; i < NUM_ITER; i++)
		run[i] = 1;
#else
	run[0] = 1;
	run[1] = 1;
	run[2] = 1;
	run[3] = 1;
	run[4] = 1;
	run[5] = 0;
	run[6] = 0;
#endif

#if 1

	start = clock();
	for (i = 0; i < NUM_ITER; i++) // Iterations   //there is dependence between iterations!!! its impossible to start new iteration before last ends
	{
		if (run[i] == 1) iteration(i, F_ro, ITER_TIME);

#if (defined(PRINT) && defined(SIMPLE))
		gnuplot_one_iteration(nmesh[i]);
#endif
	}

	double duration;
	duration = (double)(clock() - start) / CLOCKS_PER_SEC;

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
