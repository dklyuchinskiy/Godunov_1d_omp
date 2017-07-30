#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <omp.h>
#include <xmmintrin.h>

/************** MACRO ****************/

//#define NC printf("Coordinates of shockwave is setted up\n");  // It influences on output in .dat file
//#define NC2 prog // (x-x0-Dt)/h


//#define FIRST printf("The first output: only first 100 dots\n"); // For 1 task: NC, SECOND 
#define SECOND printf("The second output: printing all dots\n");
#define PRINT printf("Printing with GNUPLOT is set up\n");
//#define P_PLUS_PG printf("Printing new task about pressure and pressure gradient");
//#define SW_FINITE_DIFF  printf("Its a work with shock wave and right parts. We are determinating the wavelength\n")
//#define RW_NUM_ANALITICAL  printf("We are conculating the difference between numeric and analitical solvers\n")   // если установлен этот макрос, то print не нуженropy 
//#define ENTROPY_RP  printf("Right parts are in entropy equation. We are calculating through entropy conservation law.");
//#define SIMPLE printf("The simple print")
//#define MY_OFFLOAD  printf("GFX offload")
//#define INTEGRAL printf("Integrals are computed")
//#define FIVE_T_STEPS
#define OUTPUT_N_SMOOTH
//#define DIFF_ANALIT_RIEMANN
//#define PRINT_TAU
//#define BOOST_VEC
//#define BOOST
#define OMP
//#define RP
//#define L1-NORM
//#define ENTROPY_CHECK
//#define SW_POINTS_PRINT
//#define CFL_SWITCH


//#define BOUND_COND

//#define NEW_VERSION  // без перекидки значений с 1 на 0 элемент массива


#ifdef INTEGRAL
#define RUNGE
#endif

#define NUM_ITER 7

#define GRID 3

#define N_smooth 40

#define OMP_CORES 4
#define LOOPS 6

#define PROBLEM		12
/*
0 - shock wave
1 - rarify wave
2 - shock tube
3 - periodic continious
4 - задача с отклом
5 - two shock tube (распад разрыва) в x=4 и x=6
6 - волна разрежени€ x=6 догон€ет ударную волну в x=7; x=20
7 - две волны разрежени€ в разные стороны
8 - одна ударна€ догон€ет другую; расчет количества точек на волнах; strong sw
9 - распад разрыва; тест —ода
10 - эксперименты с гиперболическим синусом
11 - ударна€ волна встречает область разрежени€
12 - тест “оро 1
13 - тест “оро 2
14 - тест “оро 3
15 - тест “оро 4
16 - тест “оро 5
17 - shock tube [0:10], x = 5 - discontinuity
18 - shock wave with bound conditions
19 - shock wave with U=u(t)
20 - shock wave intersect another shock wave. Bound cond
*/

/*************************************/

#define EXACT_DISC

#ifdef EXACT_DISC
#define CROSS_POINT time_max_array[PROBLEM]
#define TYPE 'E'
#else
//#define CROSS_POINT 1.875  problem 1
#define CROSS_POINT 0.000
#define TYPE 'L'
#endif

/*************************************/

#if (PROBLEM == 8 )
#define DELTA 0.005
#else
#define DELTA 0.005
#endif


#define CFL	  0.4		// cfl number       изменили CFL с 0.4 до 0.8 - пересчитать пор€дки точности!
#define CFL04 CFL
#define CFL08 0.8

// вз€ть меньшие  уранты!
#define PI			3.1415926535897932

#define A_TERM		1.0			// term = - A * \vec{u} h^{2K-1} (dp/dx)^{2K}
#define K_TERM		2.0			

#define GAMMA 1.4
#define R0 8.93
#define C0 3.97



#define C1 1.0


/**************************************/

#define X_G(x) (pow((x), GAMMA))
#define SQ_2(x) ((x)*(x)/2)

/**************************************/

#define RUNGE_KUTTA	0	/*	1 - Runge-Kutta method
0 - Classic method
*/

/*		* Timers *
0.5  - shock wave
0.5  - rarify wave
0.25 - shock tube
0.25 - periodic continious
*/

/***************************************/
#if PROBLEM==4
//#define P4_ONE_WAVE
#endif

// SW 1 // SW 2 // SW 3

#define st_R3 1.0
#define st_P3 1.0
#define st_U3 0.0


#define st_P2 1.4017
#define st_R2 gyugonio(st_P3, st_R3, st_P2)
#define st_U2 ((sw_speed2(st_R3, st_U3, st_P3, st_R2, st_P2)*(1.0 - (st_R3)/(st_R2))) + ((st_R3)*(st_U3)/(st_R2)))

#ifdef P4_ONE_WAVE
#define st_P1 5.0
#else
#define st_P1 3.288303 //5
#endif 
#define st_R1 gyugonio(st_P2, st_R2, st_P1)
#define st_U1 ((sw_speed2(st_R2, st_U2, st_P2, st_R1, st_P1)*(1.0 - (st_R2)/(st_R1))) + ((st_R2)*(st_U2)/(st_R1)))

/**************************************************/

#define st_th_P2 1.0
#define st_th_R2 1.0
#define st_th_U2 0.0

#define st_th_P1 3.5
#define st_th_R1 gyugonio(st_th_P2, st_th_R2, st_th_P1)
#define st_th_U1 ((sw_speed2(st_th_R2, st_th_U2, st_th_P2, st_th_R1, st_th_P1)*(1.0 - (st_th_R2)/(st_th_R1))) + ((st_th_R2)*(st_th_U2)/(st_th_R1)))



// если в точном распаде есть характеристика с нулевой скоростью то линейный распад непременим! не имеем право примен€ть линейный распад!
// в этом случае толькой нелинейный распад!
// линеаризировать можно, только знать когда!
// но это не снимает вопросов которые поставлены, схема применима не всегда

#if (PROBLEM==0 || PROBLEM==1)
#define DISC_POINT 0.1         // 0.3 - rarification wave test
#elif (PROBLEM == 2 || PROBLEM == 9)
#define DISC_POINT 0.5
#elif (PROBLEM == 4)
#define DISC_POINT 0.0
#elif (PROBLEM == 5)
#define DISC_POINT 0
#elif (PROBLEM == 6)
#define DISC_POINT 1.0
#elif (PROBLEM == 7)
#define DISC_POINT 0.5
#elif (PROBLEM == 8) || (PROBLEM==20)
#define DISC_POINT 0.2
#elif (PROBLEM == 10)
#define DISC_POINT 0.2
#elif (PROBLEM == 11)
#define DISC_POINT 0.2
#elif (PROBLEM == 12)
#define DISC_POINT 0.3
#elif (PROBLEM == 13)
#define DISC_POINT 0.5
#elif (PROBLEM == 14)
#define DISC_POINT 0.5
#elif (PROBLEM == 15)
#define DISC_POINT 0.4
#elif (PROBLEM == 16)
#define DISC_POINT 0.8
#elif (PROBLEM == 17)
#define DISC_POINT 5.0
#elif (PROBLEM == 18)
#define DISC_POINT 0.0
#elif (PROBLEM == 19)
#define DISC_POINT 0.0
#else
#define DISC_POINT 0.0
#endif

#if (PROBLEM == 5 || PROBLEM == 17)
#define LENGTH 10.0
#elif (PROBLEM == 19)
#define LENGTH 20.0
#else
#define LENGTH 1.0
#endif

#define NUM_MESH	4			// number of mesh

extern double time_max_array[];

//				 0    1    2    3     4     5      6      7       8
extern int nmesh[];
extern int nprnt[];

extern double inflections[];

extern double boundary[][3];

extern double LOOP_TIME[][OMP_CORES];

extern double F_ro_I[], F_ro_S[], F_ro_M[], F_ro_E[];
extern double delta_D[], delta_U[], delta_P[];
extern double l1_D[], l1_U[], l1_P[];
