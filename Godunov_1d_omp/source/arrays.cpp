#include "definitions.h"
#include "support.h"

/**********************************************
This source file contains general constants
and arrays used to set some known conditions
(for example, the maximal time of computations,
number of grid points for each iteration, 
left and right points of output domain in plot file
************************************************/

/***** Constants *****/
const double g1 = (GAMMA - 1.0) / (2.0*GAMMA);
const double g2 = (GAMMA + 1.0) / (2.0*GAMMA);
const double g3 = 2.0*GAMMA / (GAMMA - 1.0);
const double g4 = 2.0 / (GAMMA - 1.0);
const double g5 = 2.0 / (GAMMA + 1.0);
const double g6 = (GAMMA - 1.0) / (GAMMA + 1.0);
const double g7 = (GAMMA - 1.0) / 2.0;
const double g8 = GAMMA - 1.0;

// density, velocity, pressure, energy, entropy, difference entropy, integral difference entropy
char prop[7] = { 'R', 'U', 'P', 'E', 'S', 'D', 'G' };

/***** Maximal time of computations *****/
#ifdef SW_POINTS_PRINT
double time_max_array[21] = { 0.5, 0.5, 0.5, 1.00, 0.5, 4.0, 0.6, 8.0, 0.4, 0.2, 0.3, 0.3, 0.2, 0.15, 0.012, 0.035, 0.012, 4.0, 0.5. 8, 0.5 };
#else
double time_max_array[21] = { 0.5, 0.5, 0.3, 1.00, 0.1, 4.0, 0.6, 0.2, 0.5, 0.2, 0.8, 0.6, 0.3, 0.15, 0.012, 0.035, 0.012, 4.0, 0.7, 8, 0.5 };
#endif

// Grid points for each iteration
#if (GRID == 3)
int nmesh[10] = { 100, 300, 900, 2700, 8100, 24300, 32768, 65536, 656100 };
//		           0    1    2    3     4     5      6      7       8   
#elif (GRID == 2)
int nmesh[10] = { 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200 };
#endif
int nprnt[10] = { 0, 1, 4, 13, 40, 121, 364, 1093, 3280 };


/******** Boundaries for plot file **************/

/**** R, U, P, -, S, diff S, diff S integral ****/

float left_SW[7] = { 0.9f, -0.1f, 0.9f, 1.1f, -0.001f, -0.0001f, -0.0001f};
float right_SW[7] = { 1.35f, 0.35f, 1.5f, 1.3f, 0.008f, 0.0002f, 0.0002f };

float left_19[7] = { 0.0f, -1.0f, 0.0f, 1.1f, 0.998f, -0.0001f, -0.0001f };
float right_19[7] = { 15.f, 5.0f, 15.f, 1.3f, 1.005f, 0.0001f, 0.0001f };

float left_6[7] = { 1.0f, -0.8f, 0.0f, 1.1f, -0.88f, -0.005f, -0.005f };
float right_6[7] = { 5.0f, 0.6f, 4.5f, 1.3f, -0.7f, 0.005f, 0.005f };

float left_7[7] = { 0.2f, -1.4f, 0.2f, 1.1f, 0.0f, 0.0f, 0.0f };
float right_7[7] = { 1.0f, 1.4f, 1.0f, 1.3f, 0.25f, 0.001f, 0.001f };

float left_10[7] = { 1.0f, -0.2f, 1.0f, 1.1f, 0.0f, -0.00001f, -0.00001f };
float right_10[7] = { 1.5f, 0.35f, 1.7f, 1.3f, 0.011f, 0.0001f, 0.0001f };

float left_11[7] = { 1.4f, 0.0f, 1.4f, 1.1f, -0.28f, -0.0001f, -0.0001f };
float right_11[7] = { 2.8f, 0.7f, 3.2f, 1.3f, -0.262f, 0.00025f, 0.00025f };

float left_SW_cs[7] = { 0.9f, -1.1f, 0.9f, 1.1f, 0.998f, 0.0f, 0.0f };
float right_SW_cs[7] = { 1.35f, 1.1f, 1.5f, 1.3f, 1.005f, 0.0000005f, 0.0000005f };

float left_ST[7] = { 0.9f, -0.1f, 0.9f, 1.1f, -0.3f, -0.0003f, -0.0003f};
float right_ST[7] = { 2.1f, 0.4f, 2.1f, 1.3f, 0.1f, 0.0003f, 0.0003f};

float left_RW[7] = { 1.5f, -0.35f, 1.3f, 1.05f, -0.270f, -0.0001f, -0.0001f };
float right_RW[7] = { 2.2f, 0.1f, 2.2f, 1.25f, -0.278f, 0.0001f, 0.0001f};

float left_sodd[7] = { 0.0f, -0.1f, 0.0f, 1.05f, 0.75f, 0.0f, 0.0f };
float right_sodd[7] = { 1.1f, 1.0f, 1.1f, 1.25f, 0.76f, 0.0000005f, 0.0000005f };

float left_2RR[7] = { 0.f, -1.0f, 0.0f, 1.05f, 0.75f, 0.0f, 0.0f };
float right_2RR[7] = { 10.f, 1.f, 12.0f, 1.25f, 0.76f, 0.0000005f, 0.0000005f };

float left_RW_SW[7] = { 0.f, -1.0f, 0.0f, 1.05f, 0.75f, 0.0f, 0.0f };
float right_RW_SW[7] = { 8.f, 1.f, 6.0f, 1.25f, 0.76f, 0.0000005f, 0.0000005f };

float left_RW_RW[7] = { 0.0f, -0.5f, 0.0f, 1.05f, 0.98f, -0.0005f, -0.0005f };
float right_RW_RW[7] = { 4.f, 2.f, 4.0f, 1.25f, 1.1f, 0.0005f, 0.0005f };

float left_RW1_RW1[7] = { 0.f, -2.0f, 0.0f, 1.05f, 0.75f, 0.0f, 0.0f };
float right_RW1_RW1[7] = { 10.f, 2.f, 10.0f, 1.25f, 0.76f, 0.0000005f, 0.0000005f };

float left_RW2_RW2[7] = { -0.5f, -1.5f, -0.5f, 1.05f, 0.98f, 0.98f, 0.98f };
float right_RW2_RW2[7] = { 8.f, 6.f, 8.0f, 1.25f, 1.1f, 1.1f, 1.1f };

float left_RW3_RW3[7] = { -0.5f, -1.5f, -0.5f, 1.05f, 0.0f, 0.0f, 0.0f };
float right_RW3_RW3[7] = { 8.f, 6.f, 8.0f, 1.25f, 2.0f, 2.0f, 2.0f };

float left_SW_SW[7] = { 1.0f, 2.0f, 0.0f, 1.1f, 0.0f, 0.0f, 0.0f };
float right_SW_SW[7] = { 7.0f, -2.0f, 6.0f, 1.3f, 0.008f, 0.0002f, 0.0002f };

float left_RP[3] = { 0, 0, 0 };
float right_RP[3] = { 4, 4, 4 };

char dip[3] = { 'D', 'I', 'P' };
float percents[NUM_ITER];

// iteration.cpp

double delta_RUP[3][NUM_ITER] = { 0 };
double l1_RUP[3][NUM_ITER] = { 0 };
