#include "definitions.h"
#include "support.h"

/***** ��������� �� ���������� �������� *****/
double
g1 = (GAMMA - 1.0) / (2.0*GAMMA),
g2 = (GAMMA + 1.0) / (2.0*GAMMA),
g3 = 2.0*GAMMA / (GAMMA - 1.0),
g4 = 2.0 / (GAMMA - 1.0),
g5 = 2.0 / (GAMMA + 1.0),
g6 = (GAMMA - 1.0) / (GAMMA + 1.0),
g7 = (GAMMA - 1.0) / 2.0,
g8 = GAMMA - 1.0;


/**************************************************/

char prop[6] = { 'R', 'U', 'P', 'C', 'S', 'D' };

#ifdef SW_POINTS_PRINT
double time_max_array[21] = { 0.5, 0.5, 0.5, 1.00, 0.5, 4.0, 0.6, 8.0, 0.4, 0.2, 0.3, 0.3, 0.2, 0.15, 0.012, 0.035, 0.012, 4.0, 0.5. 8, 0.5 };
#else
double time_max_array[21] = { 0.6, 0.5, 0.3, 1.00, 0.5, 4.0, 0.6, 0.2, 0.5, 0.2, 0.3, 0.3, 0.3, 0.15, 0.012, 0.035, 0.012, 4.0, 0.7, 4, 0.5 };
#endif

//				  0    1    2    3     4     5      6      7       8   9
#if (GRID == 3)
int nmesh[10] = { 100, 300, 900, 2700, 8100, 24400, 32768, 218700, 656100 };
#elif (GRID == 2)
int nmesh[10] = { 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200 };
#endif
int nprnt[10] = { 0, 1, 4, 13, 40, 121, 364, 1093, 3280 };

float left_SW[6] = { 0.9f, -0.1f, 0.9f, 1.1f, -0.001f, -0.0001f };
float right_SW[6] = { 1.35f, 0.35f, 1.5f, 1.3f, 0.008f, 0.0002f };

float left_19[6] = { 0.9f, -0.1f, 0.9f, 1.1f, 0.998f, -0.0001f };
float right_19[6] = { 15.f, 5.0f, 15.f, 1.3f, 1.005f, 0.0001f };

//float left_7[6] = { 0.3f, -1.0f, 0.2f, 1.1f, 0.0f, 0.0f };
//float right_7[6] = { 1.0f, 1.0f, 1.0f, 1.3f, 0.25f, 0.001f };

float left_7[6] = { 0.2f, -1.0f, 0.2f, 1.1f, 0.0f, 0.0f };
float right_7[6] = { 2.1f, 1.0f, 1.0f, 1.3f, 0.25f, 0.001f };

float left_SW_cs[6] = { 0.9f, -1.1f, 0.9f, 1.1f, 0.998f, 0.0f };
float right_SW_cs[6] = { 1.35f, 1.1f, 1.5f, 1.3f, 1.005f, 0.0000005f };

float left_ST[6] = { 0.9f, -0.1f, 0.9f, 1.1f, -0.3f, -0.0002f };
float right_ST[6] = { 2.1f, 0.35f, 2.1f, 1.3f, 0.05f, 0.0003f };

float left_RW[6] = { 1.5f, -0.35f, 1.3f, 1.05f, 0.95f, 0.0f };
float right_RW[6] = { 2.1f, 0.1f, 2.1f, 1.25f, 1.05f, 0.0000005f };

float left_sodd[6] = { 0.0f, -0.1f, 0.0f, 1.05f, 0.75f, 0.0f };
float right_sodd[6] = { 1.1f, 1.0f, 1.1f, 1.25f, 0.76f, 0.0000005f };

float left_2RR[6] = { 0.f, -1.0f, 0.0f, 1.05f, 0.75f, 0.0f };
float right_2RR[6] = { 10.f, 1.f, 12.0f, 1.25f, 0.76f, 0.0000005f };

float left_RW_SW[6] = { 0.f, -1.0f, 0.0f, 1.05f, 0.75f, 0.0f };
float right_RW_SW[6] = { 8.f, 1.f, 6.0f, 1.25f, 0.76f, 0.0000005f };

float left_RW_RW[6] = { 0.0f, -0.5f, 0.0f, 1.05f, 0.98f, -0.0005f };
float right_RW_RW[6] = { 4.f, 2.f, 4.0f, 1.25f, 1.1f, 0.0005f };

float left_RW1_RW1[6] = { 0.f, -2.0f, 0.0f, 1.05f, 0.75f, 0.0f };
float right_RW1_RW1[6] = { 10.f, 2.f, 10.0f, 1.25f, 0.76f, 0.0000005f };

float left_RW2_RW2[6] = { -0.5f, -1.5f, -0.5f, 1.05f, 0.98f, 0.98f };
float right_RW2_RW2[6] = { 8.f, 6.f, 8.0f, 1.25f, 1.1f, 1.1f };

float left_RW3_RW3[6] = { -0.5f, -1.5f, -0.5f, 1.05f, 0.0f, 0.0f };
float right_RW3_RW3[6] = { 8.f, 6.f, 8.0f, 1.25f, 2.0f, 2.0f };

//float left_RP[3] = { 1, 0, 1};
//float right_RP[3] = { 2, 1, 2};

float left_RP[3] = { 0, 0, 0 };
float right_RP[3] = { 4, 4, 4 };

char dip[3] = { 'D', 'I', 'P' };

float percents[NUM_ITER];

double boundary[NUM_ITER][3] = { 0 };
