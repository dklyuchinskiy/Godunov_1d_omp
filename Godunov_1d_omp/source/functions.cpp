#include "definitions.h"
#include "support.h"

/**************************************************
Получение плотности, давления и скорости в точке
**************************************************/
void sample(double &pm, double &um, double &s, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr,
	double &d, double &u, double &p)
{

	double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

	if (s <= um)
	{
		// точка слева от контактного разрыва
		if (pm <= pl)
		{
			// левая волна разрежения
			shl = ul - cl;

			if (s <= shl)
			{
				// точка слева
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				cml = cl * pow(pm / pl, g1);
				stl = um - cml;

				if (s > stl)
				{
					// точка слева от контактного разрыва
					d = dl * pow(pm / pl, 1.0 / GAMMA);
					u = um;
					p = pm;
				}
				else
				{
					// точка в волне разрежения
					u = g5 * (cl + g7*ul + s);
					c = g5 * (cl + g7*(ul - s));
					d = dl*pow(c / cl, g4);
					p = pl*pow(c / cl, g3);
				}
			}
		}
		else
		{
			// левая ударная волна
			pml = pm / pl;
			sl = ul - cl * sqrt(g2*pml + g1);

			if (s <= sl)
			{
				// точка слева от ударной волны
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				// точка слева от контактного разрыва
				d = dl * (pml + g6) / (pml*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
	}
	else
	{
		// точка справа от контактного разрыва
		if (pm > pr)
		{
			// правая ударная волна
			pmr = pm / pr;
			sr = ur + cr * sqrt(g2*pmr + g1);

			if (s >= sr)
			{
				// точка справа от ударной волны
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				// точка справа от контактного разрыва
				d = dr * (pmr + g6) / (pmr*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
		else
		{
			// правая волна разрежения
			shr = ur + cr;

			if (s >= shr)
			{
				// точка справа от волны разрежения
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				cmr = cr * pow(pm / pr, g1);
				str = um + cmr;

				if (s <= str)
				{
					// точка справа от контактного разрыва
					d = dr * pow(pm / pr, 1.0 / GAMMA);
					u = um;
					p = pm;
				}
				else
				{
					// точка в левой волне разрежения
					u = g5 * (-cr + g7*ur + s);
					c = g5 * (cr - g7*(ur - s));
					d = dr * pow(c / cr, g4);
					p = pr * pow(c / cr, g3);
				}
			}
		}
	}

}

/**************************************************
Решение для одной из частей
**************************************************/
void prefun(double &f, double &fd, double &p,
	double &dk, double &pk, double &ck)
{
	double ak, bk, pratio, qrt;

	if (p <= pk)
	{
		// волна разрежения
		pratio = p / pk;
		f = g4*ck*(pow(pratio, g1) - 1.0);
		fd = (1.0 / (dk*ck))*pow(pratio, -g2);

	}
	else
	{
		// ударная волна
		ak = g5 / dk;
		bk = g6*pk;
		qrt = sqrt(ak / (bk + p));
		f = (p - pk)*qrt;
		fd = (1.0 - 0.5*(p - pk) / (bk + p)) * qrt;
	}

}


/**************************************************
Вычисление начального приближения
**************************************************/
double guessp(double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr)
{
	double	cup, gel, ger,
		pmax, pmin, ppv, pq, pm,
		ptl, ptr,
		qmax, quser, um;

	/*** Вычисление приближения давления из PVRS решения Римана ***/
	quser = 2.0;
	cup = 0.25 * (dl + dr) * (cl + cr);
	ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
	ppv = ppv > 0.0 ? ppv : 0.0;
	pmin = pl > pr ? pr : pl;
	pmax = pl > pr ? pl : pr;
	qmax = pmax / pmin;

	if (qmax <= quser && pmin <= ppv && ppv <= pmax)
	{
		pm = ppv;
	}
	else
	{
		if (ppv < pmin)
		{
			// две волны разрежения
			pq = pow(pl / pr, g1);
			um = (pq*ul / cl + ur / cr + g4*(pq - 1.0)) / (pq / cl + 1.0 / cr);
			ptl = 1.0 + g7*(ul - um) / cl;
			ptr = 1.0 + g7*(um - ur) / cr;
			pm = 0.5*(pl*pow(ptl, g3) + pr*pow(ptr, g3));
		}
		else
		{
			// две ударных волны
			gel = sqrt((g5 / dl) / (g6*pl + ppv));
			ger = sqrt((g5 / dr) / (g6*pr + ppv));
			pm = (gel*pl + ger*pr - (ur - ul)) / (gel + ger);
		}

	}
	return pm;
}


/**************************************************
Определение давления и скорости контактного разрыва
**************************************************/
void starpu(double &p, double &u, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr)
{
	int i;

	double pstart,	// начальное приближение давления 
		pold,	// предыдущее приближение давления
		udiff,	// разность скоростей
		change;	// разность давлений

	double	fl, fld, // функции давления
		fr, frd;

	/*** Вычисление начального приближения ***/
	pstart = guessp(dl, ul, pl, cl, dr, ur, pr, cr);

	/*** Предыдущее приближение ***/
	pold = pstart;

	/*** Разность скоростей ***/
	udiff = ur - ul;

	/*** Метод Ньютона для определения давления ***/

	// разность между разными приближениями давлений
	change = 10.0 * EPS;

	for (i = 0; i<MAX_ITER && change>EPS; i++)
	{
		// решение для левой части
		prefun(fl, fld, pold, dl, pl, cl);

		// решение для правой части
		prefun(fr, frd, pold, dr, pr, cr);

		// очередное приближение давления
		p = pold - (fl + fr + udiff) / (fld + frd);

		// разность между разными приближениями давлений
		change = 2.0 * fabs((p - pold) / (p + pold));

		// если давление отрицательное, до обнуляем его
		if (p < 0.0) p = 0.0;

		pold = p;
	}

	// определение скорости
	u = 0.5 * (ul + ur + fr - fl);
}

void nonlinear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss)
{

	double um = 0, pm = 0;
	double s_char;

	/***** Параметры ударной трубы для соседних ячеек *****/
	// Параметры слева
	double	dl,  // плотность
		ul,  // скорость
		pl,  // давление
		cl;  // скорость звука

			 // Параметры справа
	double	dr,  // плотность
		ur,  // скорость
		pr,  // давление
		cr;  // скорость звука

	//====================================== EXACT RIEMANN PROBLEM =========================
	for (int i = 1; i < numcells; i++)
	{
		// левая сторона
		dl = R[i - 1];
		pl = P[i - 1];
		ul = U[i - 1];
		cl = sqrt(GAMMA*pl / dl);
		// правая сторона
		dr = R[i];
		pr = P[i];
		ur = U[i];
		cr = sqrt(GAMMA*pr / dr);

		starpu(pm, um, dl, ul, pl, cl, dr, ur, pr, cr);

		// Решение задачи распада разрыва
		sample(pm, um, s_char, dl, ul, pl, cl, dr, ur, pr, cr, dss[i], uss[i], pss[i]);
	}
}

void boundary_conditions(int numcells, double *dss, double *uss, double *pss, double *R, double *U, double *P)
{
	double c0, s0, cn, sn;

#if (PROBLEM==18 || PROBLEM==20)
	// set pressure
	pss[0] = initial_pressure(0.0);
	//	pss[numcells] = initial_pressure(LENGTH);

	c0 = sqrt(GAMMA*P[0] / R[0]);
	//	cn = sqrt(GAMMA*P[numcells] / R[numcells]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	//double rn_const = U[numcells] + P[numcells] / (R[numcells] * cn);

	uss[0] = l0_const + pss[0] / (R[0] * c0);
	//	uss[numcells] = rn_const - pss[numcells] / (R[numcells] * cn);

	s0 = log(P[0] / pow(R[0], GAMMA));
	//	sn = log(P[numcells] / pow(R[numcells], GAMMA));

	//dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);
	//	dss[numcells] = pow(pss[numcells] / sn, 1.0 / GAMMA);

	//R3 = dl - dl / cl * (bigU - ul);
	dss[0] = R[0] + R[0] / c0 * (uss[0] - U[0]);

	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);

#elif (PROBLEM == 19)

	uss[0] = timer;
	c0 = sqrt(GAMMA*P[0] / R[0]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	double rn_const = U[numcells] + P[numcells] / (R[numcells] * Cn);

	pss[0] = (uss[0] - l0_const)*(R[0] * c0);

	s0 = log(P[0] / pow(R[0], GAMMA));
	dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);

#elif (PROBLEM == 4)
	// set pressure
	pss[0] = 4.0 / exp(3.0*timer / time_max);
	//	pss[numcells] = initial_pressure(LENGTH);

	c0 = sqrt(GAMMA*P[0] / R[0]);
	//	cn = sqrt(GAMMA*P[numcells] / R[numcells]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	double rn_const = U[numcells] + P[numcells] / (R[numcells] * cn);

	uss[0] = l0_const + pss[0] / (R[0] * c0);
	//	uss[numcells] = rn_const - pss[numcells] / (R[numcells] * cn);

	s0 = log(P[0] / pow(R[0], GAMMA));
	//	sn = log(P[numcells] / pow(R[numcells], GAMMA));

	dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);
	//	dss[numcells] = pow(pss[numcells] / S[numcells], 1.0 / GAMMA);

#elif (PROBLEM == 3)
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[numcells], uss[numcells], pss[numcells]);
#else
	linear(R[0], U[0], P[0], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);
#endif
}

void flux_count(FILE* *array_flux, int iter, int numcells, double timer, double tau, double *t, double *UFLUX)
{
	int t_ind[N_bound] = { 0 };
	int numcells_flux;
	numcells_flux = numcells;

	double dx = LENGTH / double(numcells);

	for (int i = 0; i < N_bound; i++)
	{
		t_ind[i] = i * numcells_flux / N_bound;
	}

	if (iter == 1)
	{
		for (int i = 0; i < N_bound; i++)
			t[i] = (t_ind[i] + 0.5)*dx;
	}
	else
	{
		for (int i = 0; i < N_bound; i++)
		{
			//t[i] = t[i] + UFLUX[t_ind[i]] * tau;
			t[i] = t[i];

			fprintf(array_flux[i], "%lf %lf %lf\n", t[i], timer, UFLUX[t_ind[i]]);
		}
	}
		//t[i] = (t_ind[i] + 0.5)*dx - UFLUX[t_ind[i]] * timer;
}

void linear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss, int last)
{
	// Параметры слева
	double	dl,  // плотность
		ul,  // скорость
		pl,  // давление
		cl;  // скорость звука

			 // Параметры справа
	double	dr,  // плотность
		ur,  // скорость
		pr,  // давление
		cr;  // скорость звука

	double wtime = 0;
	double bigU, bigP, bigS, bigR, help, hl, hr, R3, R4;


#pragma omp parallel private(ul,pl,dl,ur,pr,dr,cl,cr,hl,hr,bigU,bigP,bigR,wtime) num_threads(OMP_CORES)
	{
		wtime = omp_get_wtime();
#pragma omp for schedule(dynamic,64) 
		for (int i = 1; i < numcells; i++)
		{
			ul = U[i - 1];
			pl = P[i - 1];
			dl = R[i - 1];

			ur = U[i];
			pr = P[i];
			dr = R[i];

			cl = sqrt(GAMMA*pl / dl);
			cr = sqrt(GAMMA*pr / dr);
			hl = 1.0 / (dl*cl);
			hr = 1.0 / (dr*cr);

			if (ul > cl)
			{
				bigP = pl;
				bigU = ul;
				bigR = dl;

			}
			else if (ur < -cr)
			{
				bigP = pr;
				bigU = ur;
				bigR = dr;
			}
			else
			{
				bigP = (ul - ur + pl / (dl*cl) + pr / (dr*cr)) / (hl + hr);
				bigU = (dl*cl*ul + dr*cr*ur + pl - pr) / (dl*cl + dr*cr);

				R3 = dl - dl / cl * (bigU - ul);
				R4 = dr + dr / cr * (bigU - ur);
				if (bigU > 0) bigR = R3;
				else bigR = R4;
				
			}

			uss[i] = bigU;
			pss[i] = bigP;
			dss[i] = bigR;

		}
		wtime = omp_get_wtime() - wtime;
		LOOP_TIME[2][omp_get_thread_num()] += wtime;
		if (last) printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[2][omp_get_thread_num()]);
	}
}

/* Linear solver */
void linear(double dl, double ul, double pl, double dr, double ur, double pr, double &d, double &u, double &p) // передача параметров по ссылке
{
	double bigU, bigP, bigS, bigR, cl, cr, help, hl, hr, R3, R4;

	cl = sqrt(GAMMA*pl / dl);
	cr = sqrt(GAMMA*pr / dr);
	hl = 1.0 / (dl*cl);
	hr = 1.0 / (dr*cr);


	
	if (ul > cl)
	{
		bigP = pl;
		bigU = ul;
		bigR = dl;

	}
	else if (ur < -cr)
	{
		bigP = pr;
		bigU = ur;
		bigR = dr;
	}
	else
	{
		bigP = (ul - ur + pl / (dl*cl) + pr / (dr*cr)) / (hl + hr);
		bigU = (dl*cl*ul + dr*cr*ur + pl - pr) / (dl*cl + dr*cr);
	/*	if (bigU >= 0) bigS = pl / pow(dl, GAMMA);
		else bigS = pr / pow(dr, GAMMA);
		help = bigP / bigS;
		bigR = pow(help, 1.0 / GAMMA);*/
		R3 = dl - dl / cl * (bigU - ul);
		R4 = dr + dr / cr * (bigU - ur);
		if (bigU > 0) bigR = R3;
		else bigR = R4;
	}
	
	u = bigU;
	p = bigP;
	d = bigR;

}

void first_step_validation(FILE *file3, int numcells, int c_c, double timer, double *R, double *U, double *P, double *dss, double *uss, double *pss)
{
	double len = LENGTH;
	double dx = len / double(numcells);
	double x;


	fprintf(file3, "cc: %d timer: %lf\n\n", c_c, timer);
	
	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;
		fprintf(file3, "%lf %lf %lf %lf %lf %lf %lf\n", x, R[i], U[i], P[i], dss[i], uss[i], pss[i]);
	}
	fprintf(file3,"*********************************************************\n\n");
	
}


void linear_check(double dl, double ul, double pl, double dr, double ur, double pr, int &left, int &middle, int &right, int numb)
{
	double bigU, bigP, bigS, bigR, cl, cr, help, hl, hr;
	cl = sqrt(GAMMA*pl / dl);
	cr = sqrt(GAMMA*pr / dr);
	hl = 1.0 / (dl*cl);
	hr = 1.0 / (dr*cr);

	if (numb == 1)
	{
		if (ul > cl || ur < -cr)
		{
			if (ul > cl)
			{
				left++;
			}
			if (ur < -cr)
			{
				right++;
			}
		}
		else
		{
			middle++;
		}
	}

	if (numb == 2)
	{
		if (ul + cl < 0 || ur - cr > 0)
		{
			if (ul + cl < 0)
			{
				left++;
			}
			if (ur - cr > 0)
			{
				right++;
			}
		}
		else
		{
			middle++;
		}
	}


}

double gyugonio(double p1, double ro1, double p2/*за ударной волной*/)
{
	return ro1*((GAMMA + 1.0)*p2 + (GAMMA - 1.0)*p1) / ((GAMMA - 1.0)*p2 + (GAMMA + 1.0)*p1);
	//  R2 = ro2*((GAMMA + 1.0)*P + (GAMMA - 1.0)*p2) / ((GAMMA - 1.0)*P + (GAMMA + 1.0)*p2);
}

double sw_speed(double ro1, double ro2, double u1, double u2)
{
	return (ro1*u1 - ro2*u2) / (ro1 - ro2);

}

double sw_speed2(double ro1, double u1, double p1, double ro2 /*за ударной волной*/, double p2 /*за ударной волной*/)
{
	return u1 + sqrt(ro2*(p2 - p1) / ro1 / (ro2 - ro1));
}

double law(double t)
{
	return sqrt(t);
}

#ifdef BOOST_VEC
#pragma omp declare simd
#endif
double initial_density(double x)
{
	switch (PROBLEM)
	{
	case 19:
	case 18:
	case 0:	if (x <= DISC_POINT)                  // разрыв значений в начальный момент времени в точке x=0.2 ( для задачи с УВ )
		return 1.271413930046081;
			else
				return 1.0;
		/*case 1:	if (x <= DISC_POINT)
		return 1.551608179649565;
		else
		return 2.0;*/
	case 1:	if (x <= DISC_POINT)
		return 1.0;
			else
				return 0.585507;
	case 2:	if (x <= 0.5)
		return 2.0;
			else
				return 1.0;
	case 3:	return 1.0;
		/*	case 4: if (x <= DISC_POINT) return 1.271413930046081;
		else if (x >= 1.0 - DISC_POINT) return 1.271413930046081;
		else return 1.0;*/
	case 4: return R0;

	case 5: if (x <= 3.0) return 6;
			else if (x >= 7.0) return 6;
			else return 4.0;
	case 6: if (x <= 0.2) return 1.36205;
			else if (x <= 0.3) return 4.875;
			else return 3;
	case 7: if (x < DISC_POINT) return 1.0;
			else return 2.0;
	case 8: if (x <= 0.05) return st_R1;
			else if (x <= 0.3) return st_R2;
			else return st_R3;
	case 9: if (x <= DISC_POINT) return 1.0;
			else return 0.125;
	case 10: if (x <= DISC_POINT) return st_th_R1;
			 else return st_th_R2;
			 //case 10: return -0.55*tanh(2 * PI*(x - 0.5)) + 1.546; // weak
			 //case 10: return -0.9*tanh(2 * PI*(x - 0.5)) + 1.9; // strong
	case 11: if (x <= 0.1) return st_R2;
			 else if (x <= 0.6) return st_R3;
			 else 0.2;
			 // TORO tests
			 // 1
	case 12: if (x <= DISC_POINT) return 1.0;
			 else return 0.125;
			 // 2
	case 13: return 1.0;
		// 3
	case 14: return 1.0;
		// 4
	case 15: if (x <= DISC_POINT) return 5.99924;
			 else return 5.99242;
			 // 5
	case 16: return 1.0;
	case 17:	if (x <= 5.0) return 2.0;
				else return 1.0;
	case 20: if (x <= 0.00) return st_R1;
			else if (x <= DISC_POINT) return st_R2;
			else return st_R3;
	}
	return 0;
}
#ifdef BOOST_VEC
#pragma omp declare simd
#endif
double initial_pressure(double x)
{
	switch (PROBLEM)
	{
	case 19:
	case 18:
	case 0:	if (x <= DISC_POINT)
		return 1.401789770179879;
			else
				return 1.0;
		/*case 1:	if (x <= DISC_POINT)
		return 1.401789770179879;
		else
		return 2.0;*/
	case 1:	if (x <= DISC_POINT)
		return 1.0;
			else
				return 0.466727;

	case 2:	if (x <= 0.5)
		return 2.0;
			else return 1.0;
	case 3:	return 1.0;
		/*case 4: if (x <= DISC_POINT) return 1.401789770179879;
		else if (x >= 1.0 - DISC_POINT) return 1.401789770179879;
		else return 1.0;*/
	case 4: return 0.0;

	case 5: if (x <= 3.0) return 7;
			else if (x >= 7.0) return 7;
			else return 2.0;
	case 6: if (x <= 0.2) return 0.67099;
			else if (x <= 0.3) return 4;
			else return 2;
	case 7: return 1.0;
	case 8: if (x <= 0.05) return st_P1;
			else if (x <= 0.3) return st_P2;
			else return st_P3;
	case 9: if (x <= DISC_POINT) return 1.0;
			else return 0.1;
			//case 10: return -tanh(2 * PI*(x - 0.5)) + 2; //weak

			//case 10: return -2 * tanh(2 * PI*(x - 0.5)) + 3; //strong


	case 10: if (x <= DISC_POINT) return st_th_P1;
			 else return st_th_P2;
	case 11: if (x <= 0.1) return st_P2;
			 else return st_P3;
			 // TORO tests
			 // 1
	case 12: if (x <= DISC_POINT) return 1.0;
			 else return 0.1;
			 // 2
	case 13: return 0.4;
		// 3
	case 14: if (x <= DISC_POINT) return 1000.0;
			 else return 0.01;
			 // 4
	case 15: if (x <= DISC_POINT) return 460.894;
			 else return 46.095;
			 // 5
	case 16: if (x <= DISC_POINT) return 1000.0;
			 else return 0.01;
	case 17: if (x <= 5.0) return 2.0;
			 else return 1.0;
	case 20: if (x <= 0.00) return st_P1;
			else if (x <= DISC_POINT) return st_P2;
			else return st_P3;
	}
	return 0;
}
#ifdef BOOST_VEC
#pragma omp declare simd
#endif
double initial_velocity(double x)
{
	switch (PROBLEM)
	{
	case 19:
	case 18:
	case 0:	if (x <= DISC_POINT)
		return 0.292868067614595;
			else
				return 0.0;
		/*case 1:	if (x <= DISC_POINT)
		return -0.292868067614595;
		else
		return 0.0;*/

	case 1:	if (x <= DISC_POINT)
		return 0.75 + 0.4533;// -0.43;
			else
				return 1.386 + 0.4533;// -0.43;

	case 2:	return 0.0;
		//if (x <= 0.5)  2 ударных волны
		//	return 1.0;
		//else return 0.0;
	case 3:	return sin(2.0*PI*x);
		/*case 4: if (x <= DISC_POINT) return 0.292868067614595;
		else if (x >= 1.0 - DISC_POINT) return 0.292868067614595;
		else return 0.0;*/
	case 4: return 0;

	case 5: if (x <= 3.0) return 0;
			else if (x >= 7.0) return 0.0;
			else return 0.0;
	case 6: if (x <= 0.2) return -0.7;
			else if (x <= 0.3) return 0.50637;
			else return 0;
	case 7: if (x <= DISC_POINT) return -1.0;
			else return 1.0;
	case 8: if (x <= 0.05) return st_U1;
			else if (x <= 0.3) return st_U2;
			else return st_U3;
	case 9: return 0.0;
	case 10: if (x <= DISC_POINT) return st_th_U1;
			 else return st_th_U2;
			 //	case 10: return -0.51*tanh(2 * PI*(x - 0.5)) + 0.51; // weak

			 //case 10: return -0.8*tanh(2 * PI*(x - 0.5)) + 0.8; // strong
	case 11: if (x <= 0.1) return st_U2;
			 else return st_U3;
			 // TORO tests
			 // 1
	case 12: if (x <= DISC_POINT) return 0.75;// +0.2;// +0.4533;
			 else return 0.0;// +0.2;// 0.4533;
			 // 2
	case 13: if (x <= DISC_POINT) return -2.0;
			 else return 2.0;
			 // 3
	case 14: return 0.0;
		// 4
	case 15: if (x <= DISC_POINT) return 19.5975;
			 else return -6.19633;
			 // 5
	case 16: return -19.59745;
	case 17: return 0.0;
	case 20: if (x <= 0.00) return st_U1;
			else if (x <= DISC_POINT) return st_U2;
			else return st_U3;
	}
	return 0;
}

void difference_SW(int numcells, double timer, double *R, double *U, double *P, 
	                          double *shw_diff_d, double *shw_diff_u, double *shw_diff_p, 
	                          double *shw_analit_d, double *shw_analit_u, double *shw_analit_p)
{
	/*******************difference analit and numeric solutions************/
	/**********************shock wave*************************************/

	analitical_SW(numcells, initial_pressure(0.05), initial_density(0.05), initial_velocity(0.05),
		                    initial_pressure(0.2), initial_density(0.2), initial_velocity(0.2),
		                    shw_analit_p, shw_analit_u, shw_analit_d, timer);
	for (int j = 0; j < numcells; j++)
	{
	shw_diff_p[j] = P[j] - shw_analit_p[j];
	shw_diff_u[j] = U[j] - shw_analit_u[j];
	shw_diff_d[j] = R[j] - shw_analit_d[j];
	}
}

void analitical_RW(FILE* file_name, double ip_l, double id_l, double iu_l,
	double ip_r, double id_r, double iu_r, double numb)
{
	FILE* out_ul, *out_pl, *out_dl;
	FILE* out_ur, *out_pr, *out_dr;
	out_pl = fopen("data_pl.txt", "w");
	out_ul = fopen("data_ul.txt", "w");
	out_dl = fopen("data_dl.txt", "w");

	out_pr = fopen("data_pr.txt", "w");
	out_ur = fopen("data_ur.txt", "w");
	out_dr = fopen("data_dr.txt", "w");

	double xl, xr, c, l0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	//tt = 1 / (u2 + c2 - D);
	tt = numb*0.1;   // only this one have NUMB count since 1  !!!

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	fprintf(file_name, "u(x)=%lf*(x-0.1)/%lf+(%lf)\n", q3, tt, q11);

	xl = (iu_l - q11)*tt / q3 + 0.1;
	xr = (iu_r - q11)*tt / q3 + 0.1;
	fprintf(file_name, "xu_l=%lf\nxu_r=%lf\n", xl, xr);
	fprintf(out_ul, "0 %lf\n%lf %lf", iu_l, xl, iu_l);
	fprintf(out_ur, "%lf %lf\n1 %lf", xr, iu_r, iu_r);

	// для давления  p=(q4*x-q5)^q6

	q4 = q1 / (sqrt(GAMMA)*pow(A, 1 / (2 * GAMMA)));
	q5 = q4*l0;
	q6 = q2*GAMMA;
	fprintf(file_name, "p(x)=(%lf*(x-0.1)/%lf-(%lf))**%lf\n", q4, tt, q5, q6);  // почему, если разрвы в 0.1 мы должны отнять

	xl = (pow(ip_l, (1 / q6)) + q5)*tt / q4 + 0.1;
	xr = (pow(ip_r, (1 / q6)) + q5)*tt / q4 + 0.1;
	fprintf(file_name, "xp_l=%lf\nxp_r=%lf\n", xl, xr);

	fprintf(out_pl, "0 %lf\n%lf %lf", ip_l, xl, ip_l);
	fprintf(out_pr, "%lf %lf\n1 %lf", xr, ip_r, ip_r);

	// для плотности ro=(q7*x-q8)^q2

	q7 = q1 / sqrt(GAMMA*A);
	q8 = q7*l0;
	fprintf(file_name, "ro(x)=(%lf*(x-0.1)/%lf-(%lf))**%lf\n", q7, tt, q8, q2);

	xl = (pow(id_l, (1 / q2)) + q8)*tt / q7 + 0.1;
	xr = (pow(id_r, (1 / q2)) + q8)*tt / q7 + 0.1;
	fprintf(file_name, "xd_l=%lf\nxd_r=%lf\n\n", xl, xr);

	fprintf(out_dl, "0 %lf\n%lf %lf", id_l, xl, id_l);
	fprintf(out_dr, "%lf %lf\n1 %lf", xr, id_r, id_r);

	fprintf(file_name, "id_l(x)=%lf\n", id_l);
	fprintf(file_name, "ip_l(x)=%lf\n", ip_l);
	fprintf(file_name, "iu_l(x)=%lf\n\n", iu_l);

	fprintf(file_name, "id_r(x)=%lf\n", id_r);
	fprintf(file_name, "ip_r(x)=%lf\n", ip_r);
	fprintf(file_name, "iu_r(x)=%lf\n", iu_r);

	fclose(out_dl);
	fclose(out_dr);
	fclose(out_ul);
	fclose(out_ur);
	fclose(out_pl);
	fclose(out_pr);

	return;
}

void analitical_SW(int numcells, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r, double *res_p, double *res_u, double* res_d, double timer)
{
	double D_analit;
	double x;
	double *x_lay;
	x_lay = new double[numcells];
	double dx;
	dx = LENGTH / double(numcells);
	for (int i = 0; i < numcells; i++)
		x_lay[i] = 0;

	D_analit = (id_r * iu_r - id_l * iu_l) / (id_r - id_l);

	x = D_analit * timer + 0.1; // discontinuity point

	for (int i = 0; i < numcells; i++)
	{
		x_lay[i] = i*dx + 0.5*dx;
		if (x_lay[i] <= x)
		{
			res_p[i] = ip_l;
			res_d[i] = id_l;
			res_u[i] = iu_l;
		}
		else
		{
			res_p[i] = ip_r;
			res_d[i] = id_r;
			res_u[i] = iu_r;
		}
	}

}

/*modeling of exact solution of riemann problem*/
void analitical_riemann(int numcells, double p1, double ro1, double u1, double p2, double ro2, double u2, double *sol_p, double *sol_u)
{
	//double P = 0, U = 0;  //искомые давление и скорость (постоянные значения слева и справа от контактного разрыва)
	double R1, E1; //искомые значения плотности и внутренней энергии СЛЕВА от контактного разрыва
	double R2, E2; //искомые значения плотности и внутренней энергии СПРАВА от контактного разрыва

	// начальные данные

	double c1; // давление, плотность, скорость СЛЕВА
	double c2; // давление, плотность, скорость СПРАВА



	c1 = sqrt(GAMMA*p1 / ro1);
	c2 = sqrt(GAMMA*p2 / ro2);

	double a1; // массовая скорость слева
	double a2; // массовая скорость справа


	double P_prev = 0;
	double P_now = 0;

	double delta = 0;
	double drob = 0;
	double h1 = 0, h2 = 0;
	int l = 0;  //счетчик

	P_prev = 20; // P_prev > p1 и P_prev > p2;    ударные волны СЛЕВА и СПРАВА


	/////////////// ПЕРЕРАСЧЕТ ПО ВИДОИЗМЕНЕННЫМ ФОРМУЛАМ///////////////////////////
	do
	{
		if (P_prev > p1) a1 = sqrt(ro1*((GAMMA + 1)*P_prev / 2.0 + (GAMMA - 1)*p1 / 2.0));
		else
		{
			h2 = (GAMMA - 1.0) / (2.0 * GAMMA);
			h1 = P_prev / p1;
			drob = (1.0 - P_prev / p1) / (1.0 - pow(h1, h2));
			a1 = (GAMMA - 1.0) / (2.0 * GAMMA)*ro1*c1*drob;
		}
		drob = 0;

		if (P_prev>p2) a2 = sqrt(ro2*((GAMMA + 1)*P_prev / 2.0 + (GAMMA - 1)*p2 / 2.0));
		else
		{
			h2 = (GAMMA - 1) / (2.0 * GAMMA);
			h1 = P_prev / p2;
			drob = (1.0 - P_prev / p2) / (1.0 - pow(h1, h2));
			a2 = (GAMMA - 1) / (2.0 * GAMMA)*ro2*c2*drob;
		}
		drob = 0;

		double fi_P1 = 0;
		double alfa = 0;
		double help = 0, h3 = 0, h4 = 0;
		double z = 0;


		z = P_prev / (p1 + p2);
		h3 = (GAMMA + 1) / 2.0 / GAMMA;
		h4 = (GAMMA - 1) / 2.0 / GAMMA;
		help = ((GAMMA - 1) / 3.0 / GAMMA)*((1.0 - z) / (pow(z, h3)*(1.0 - pow(z, h4)))) - 1;
		if (help>0) alfa = help;
		else alfa = 0;
		fi_P1 = (a2*p1 + a1*p2 + a1*a2*(u1 - u2)) / (a1 + a2);
		P_now = (alfa*P_prev + fi_P1) / (1 + alfa);

		delta = P_now - P_prev;

		// переприсваивание

		P_prev = P_now;
		P_now = 0;

		l++;
	} while (fabs(delta)>0.0000001);

	*sol_p = P_prev;
	*sol_u = (a1*u1 + a2*u2 + p1 - p2) / (a1 + a2);


	//	printf("Pressure: %lf\nSpeed: %lf\nIterations: %i\n", *sol_p, *sol_u, l);
}

void rw_diff_num_analit(int numb, int numcells, double *R, double *U, double *P)
{
	/* Rarify wave - numeric vs analitic */
	double c, l0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;
	double iu_l, id_l, ip_l, iu_r, id_r, ip_r;
	double x, x_NC, xl, xr;
	int check1 = 0, check2 = 0;

	double dx = LENGTH / double(numcells);

	iu_l = initial_velocity(0.05);
	id_l = initial_density(0.05);
	ip_l = initial_pressure(0.05);

	iu_r = initial_velocity(0.2);
	id_r = initial_density(0.2);
	ip_r = initial_pressure(0.2);

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	tt = (numb + 1)*0.1;

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	xl = (iu_l - q11)*tt / q3 + DISC_POINT; // счет по старому, по итерациям - какая итерация, такое и время. 
	xr = (iu_r - q11)*tt / q3 + DISC_POINT; // в нашем случае ЭТО НЕВЕРНО, так как теперь итерации отвечают только за ШАГ СЕТКИ, время должно быть ФИКСИРОВАНО
									 // здесь ошибка при вычислении точного решения!!!


									 // 0-U, 1-P, 2-R, 3-RU, 4-RE
	double** difference_RW;
	difference_RW = new double*[5];

	for (int j = 0; j < 5; j++)
		difference_RW[j] = new double[numcells];

	double* RW_R, *RW_P, *RW_U;
	RW_R = new double[numcells];
	RW_P = new double[numcells];
	RW_U = new double[numcells];

	int counter = 0, counter_all = 0, counter2 = 0;

	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;

		if (x < xl)
		{
			RW_U[i] = iu_l; // FOR RARIFY WAVE
			RW_P[i] = ip_l;
			RW_R[i] = id_l;
			difference_RW[0][i] = U[i] - RW_U[i];
			difference_RW[1][i] = P[i] - RW_P[i];
			difference_RW[2][i] = R[i] - RW_R[i];
			difference_RW[3][i] = 0;
			difference_RW[4][i] = 0;
			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter2++;
			}

		}
		if (x >= xl && x <= xr)
		{
			counter_all++;
			RW_U[i] = RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[0][i] = U[i] - RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			RW_P[i] = RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[1][i] = P[i] - RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter++;
				counter2++;
			}

			RW_R[i] = RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[2][i] = R[i] - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			difference_RW[3][i] = R[i] * U[i] - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)*RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[4][i] = R[i] * (P[i] / pow(R[i], GAMMA) + SQ_2(U[i])) - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)*(RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r) / pow(RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r), GAMMA) + SQ_2(RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)));
		}
		if (x > xr)
		{
			RW_U[i] = iu_r;
			RW_P[i] = ip_r;
			RW_R[i] = id_r;

			difference_RW[0][i] = U[i] - RW_U[i];
			difference_RW[1][i] = P[i] - RW_P[i];
			difference_RW[2][i] = R[i] - RW_R[i];
			difference_RW[3][i] = 0;
			difference_RW[4][i] = 0;
			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter2++;
			}

		}

	}

	check1 = 0;
	check2 = 0;

	int *i_helper;
	i_helper = new int[10];

	double *x_helper;
	x_helper = new double[10];

	double xl_num, xr_num;

	/****************Boundary of numerical rarify wave******************/

	int i_mem_left = 0, i_mem_right = 0;

	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;
		if (x >= xl && check1 == 0)
		{
			xl_num = x;   // по координате x
			i_mem_left = i;  // по счетчику i
			check1 = 1;
		}
		if (x >= xr && check2 == 0)
		{
			xr_num = x;
			i_mem_right = i;
			check2 = 1;
		}
	}
	printf("%lf %lf", xl_num, xr_num);
	/****************Boundary of numerical rarify wave******************/

	if (numb > 0)
	{
		int helper = counter_all / 10;

		for (int j = 0; j < 10; j++)
		{
			i_helper[j] = i_mem_left + j*helper;
		}
		printf("\n");
		for (int i = 0; i < numcells; i++)
		{
			x = i*dx + 0.5*dx;
			for (int j = 0; j < 10; j++)
			{
				if (i == i_helper[j])
				{
					x_helper[j] = x - (U[i] + sqrt(GAMMA*P[i] / R[i]))*time_max_array[PROBLEM];
					//		printf("%lf\n", x_helper[j]);
				}
			}
		}
	}

	printf("\nxl: %lf, xr: %lf\n", xl, xr);
	printf("points in rarify wave %d %d\n", counter, counter_all);
	printf("points2 %d %d\n", counter2, numcells);
	percents[numb] = float(counter) / float(counter_all) * 100.0f;  //when in int i counter is devided by counter all, its a деление нацело, so the result of 1/4=0;
	printf("percents of the middle: %f\n", percents[numb]);
	percents[numb] = float(counter2) / float(numcells) * 100.0f;  //when in int i counter is devided by counter all, its a деление нацело, so the result of 1/4=0;
	printf("percents of all stream [0:1]: %f\n", percents[numb]);


	FILE *out4;
	FILE *out5;
	char *FileName2, *FileName3;
	FileName2 = new char[64];
	FileName3 = new char[64];
	sprintf(FileName2, "N%04d_RW_difference.dat", numcells);
	sprintf(FileName3, "N%04d_RW_NUM_ANALITIC.dat", numcells);
	out4 = fopen(FileName2, "w");
	out5 = fopen(FileName3, "w");
	check1 = 0;
	check2 = 0;
	for (int i = 0; i < numcells; i++)
	{
#ifndef NC
		x = i*dx + 0.5*dx;
		if (x >= xl && check1 == 0)
		{
			fprintf(out5, "left b: %lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
			check1 = 1;
			continue;
		}
		if (x >= xr && check2 == 0)
		{
			fprintf(out5, "right b: %lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
			check2 = 1;
			continue;
		}
		fprintf(out4, "%lf %lf %lf %lf\n", x, difference_RW[0][i], difference_RW[1][i], difference_RW[2][i]);
		fprintf(out5, "%lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
#else

		// Выведем в новых координатах
		x_NC = (i*dx + 0.5*dx - 0.1) / (time_max_array[0] * (numb + 1));
		fprintf(out4, "%lf %lf %lf %lf\n", x_NC, difference_RW[0][i], difference_RW[1][i], difference_RW[2][i]);
		fprintf(out5, "%lf %lf %lf %lf %lf %lf %lf\n", x_NC, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i]);
#endif
	}
	fclose(out4);
	fclose(out5);

}

void analitical_riemann_modeling(int numcells, double ro1, double u1, double p1, double ro2, double u2, double p2, double timer,
	/*output*/double *all_d, double *all_u, double *all_p)
{
	double static P, U; //solution of Riemann problem
	double R1, R2;
	int numcells2 = numcells / 2;
	double c1, C, prop, r0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;
	int static ex_sol = 0;

	double v1, V;

	double *xx;
	xx = new double[numcells];

	double dx;
	dx = LENGTH / double(numcells);

	if (ex_sol == 0)
	{
		analitical_riemann(numcells, p1, ro1, u1, p2, ro2, u2, &P, &U); //для всех шагов сетки и шагов по времени нам нужно получить точное решение всего 1 раз
	}

	c1 = sqrt(GAMMA*p1 / ro1);
	r0 = u1 + 2 * c1 / (GAMMA - 1); // для левой волны разрежения; r0 инвариант - постоянен
	A = p1 / (pow(ro1, GAMMA));
	C = (r0 - U)*(GAMMA - 1) / 2.0;

	V = U - C;
	v1 = u1 - c1;
	//находим плотность R1 за волной разрежения через cвойство постоянности инварианта римана r0
	R1 = GAMMA*P / (C*C);
	//находим плотность sol_d[0] за ударной волной c помощью адиабаты Гюгонио
	R2 = ro2*((GAMMA + 1.0)*P + (GAMMA - 1.0)*p2) / ((GAMMA - 1.0)*P + (GAMMA + 1.0)*p2);

	/*крнтактный разрыв*/
	double x_KR;
	x_KR = U*timer + 0.5;
	//	printf("cont razr %lf \n", x_KR);

	/*Строим решение справа от x=0.5 - ударная волна*/
	double D_analit;
	double x;

	D_analit = (ro2* u2 - R2 * U) / (ro2 - R2);
	x = D_analit * timer + 0.5; // discontinuity point
	//	printf("shock wave %lf \n", x);

	double c1_rw, C_rw, l0, A_rw;

	c1_rw = sqrt(GAMMA*p1 / ro1);
	l0 = u1 - 2 * c1_rw / (GAMMA - 1); // для правой волны разрежения; l0 инвариант - постоянен
	A_rw = p1 / (pow(ro1, GAMMA));
	C_rw = (U - l0)*(GAMMA - 1) / 2.0;

	if (ex_sol == 0)
	{
		printf("u1-c1: %lf\nU-C: %lf, D_analit: %lf, D_cont razr: %lf\n", v1, V, D_analit, U);
		printf("t=0.2: x1: %lf, x2: %lf, x3: %lf, x4: %lf\n", v1*0.2 + 0.5, V*0.2 + 0.5, D_analit*0.2 + 0.5, U*0.2 + 0.5);
		printf("rariry wave t=0.25: x1: %lf, x2: %lf\n", (u1 + c1_rw)*0.25 + DISC_POINT, (U + C)*0.25 + DISC_POINT);
		ex_sol = 1;
	}

	for (int i = 0; i < numcells; i++)
	{
		xx[i] = i*dx + 0.5*dx;
	}

	for (int i = numcells2; i < numcells; i++)
	{
		if (xx[i] <= x)
		{
			all_p[i] = P;
			all_u[i] = U;
		}
		else
		{
			all_p[i] = p2;
			all_d[i] = ro2;
			all_u[i] = u2;
		}
	}

	/*отдельно для плотности расчет, учитываем скорость контактного разрыва*/
	for (int i = 0; i < numcells; i++)
	{
		if (xx[i] >= x_KR && xx[i] <= x) all_d[i] = R2;
	}


	/*Строим решение слева от x=0.5 - волна разрежения*/
	double xl, xr;

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	q4 = q1 / (sqrt(GAMMA)*pow(A, 1 / (2 * GAMMA)));
	q6 = q2*GAMMA;
	q7 = q1 / sqrt(GAMMA*A);

	xl = (u1 - q1*r0)*timer / q3 + 0.5;
	xr = (U - q1*r0)*timer / q3 + 0.5;
	//	printf("%lf %lf\n", xl, xr);

	for (int i = 0; i < numcells2; i++)
	{
		if (xx[i] < xl)
		{
			all_p[i] = p1;
			all_d[i] = ro1;
			all_u[i] = u1;
		}
		else if (xx[i] >= xl && xx[i] <= xr)
		{
			all_p[i] = pow(q4*(r0 - (xx[i] - 0.5) / timer), q6);
			all_d[i] = pow(q7*(r0 - (xx[i] - 0.5) / timer), q2);
			all_u[i] = q3*(xx[i] - 0.5) / timer + q1*r0;
		}
		else
		{
			all_p[i] = P;
			all_u[i] = U;
		}
	}

	for (int i = 0; i < numcells; i++)
	{
		if (xx[i] >= xr && xx[i] <= x_KR) all_d[i] = R1;
	}

	free(xx);
}

void difference_analitical_riemann_Linf(int numb, double* R, double*U, double*P, double* R_D, double*R_U, double*R_P, double &delta_ro, double &delta_u, double &delta_p)
{
	double *x_lay;
	x_lay = new double[nmesh[numb]];
	x_lay[0:nmesh[numb]] = 0;

	double dx;
	dx = LENGTH / double(nmesh[numb]);
	double* difference_D;
	double* difference_U;
	double* difference_P;

	delta_ro = 0;
	delta_u = 0;
	delta_p = 0;

	int new_numcells = int(0.41*nmesh[numb]);

	difference_P = new double[new_numcells];
	difference_U = new double[new_numcells];
	difference_D = new double[new_numcells];

	for (int i = 0; i < new_numcells; i++)
	{
		x_lay[i] = i*dx + 0.5*dx;
		difference_D[i] = fabs(R[i] - R_D[i]);
		difference_U[i] = fabs(U[i] - R_U[i]);
		difference_P[i] = fabs(P[i] - R_P[i]);
	}

	for (int i = 0; i < new_numcells; i++)
	{
		if (difference_D[i]>delta_ro) delta_ro = difference_D[i];
		if (difference_U[i]>delta_u) delta_u = difference_U[i];
		if (difference_P[i]>delta_p) delta_p = difference_P[i];
	}

	delta_ro = delta_ro*dx;
	delta_u = delta_u*dx;
	delta_p = delta_p*dx;
}

void difference_analitical_riemann_L1(int numb, double* R, double*U, double*P, double* R_D, double*R_U, double*R_P, double &sum_ro, double &sum_u, double &sum_p)
{
	double *x_lay;
	x_lay = new double[nmesh[numb]];
	x_lay[0:nmesh[numb]] = 0;

	double dx;
	dx = LENGTH / double(nmesh[numb]);
	double* difference_D;
	double* difference_U;
	double* difference_P;

	sum_ro = 0;
	sum_u = 0;
	sum_p = 0;

	int new_numcells = int(0.41*nmesh[numb]);

	difference_P = new double[new_numcells];
	difference_U = new double[new_numcells];
	difference_D = new double[new_numcells];

	for (int i = 0; i < new_numcells; i++)
	{
		x_lay[i] = i*dx + 0.5*dx;
		difference_D[i] = fabs(R[i] - R_D[i]);
		difference_U[i] = fabs(U[i] - R_U[i]);
		difference_P[i] = fabs(P[i] - R_P[i]);
		//if (numb == 0 || numb==1) printf("%lf\n", difference_D[i]);
	}



	for (int i = 0; i < new_numcells; i++)
	{
		sum_ro += difference_D[i];
		sum_u += difference_U[i];
		sum_p += difference_P[i];
	}

	sum_ro = sum_ro*dx;
	sum_u = sum_u*dx;
	sum_p = sum_p*dx;
}

/* 0 - U, 1 - P, 2 - D */
double RW_prop(int digit, double x, double numb, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r)
{
	double c, prop, l0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	tt = (numb + 1)*0.1; 	//tt = 1 / (u2 + c2 - D);

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	if (digit == 0)
	{
		return prop = q3*(x - 0.1) / tt + q11;
	}

	// для давления  p=(q4*x-q5)^q6

	q4 = q1 / (sqrt(GAMMA)*pow(A, 1 / (2 * GAMMA)));
	q5 = q4*l0;
	q6 = q2*GAMMA;
	if (digit == 1)
	{
		return prop = pow((q4*(x - 0.1) / tt - q5), q6);

	}

	// для плотности ro=(q7*x-q8)^q2

	q7 = q1 / sqrt(GAMMA*A);
	q8 = q7*l0;

	if (digit == 2)
	{
		return prop = pow((q7*(x - 0.1) / tt - q8), q2);
	}

	return 0;
}

void outline_integral_riemann(int numcells, double timer, double tau, double tt1, double tt2, double xx1, double xx2, double* R, double*U, double*P, double*RE, double*S,
	/*output*/ double sum[4][4])
{
	int static check1 = 0;
	int static check2 = 0;

	double static sum_l_M, sum_l_I, sum_l_S, sum_l_E;
	double static sum_r_M, sum_r_I, sum_r_S, sum_r_E;
	double static sum_t_M, sum_t_I, sum_t_S, sum_t_E;
	double static sum_b_M, sum_b_I, sum_b_S, sum_b_E;

	int static numcells_check;

	if (numcells == 100) numcells_check = 100;

	/*проверка на новый пространственный шаг*/
	if (numcells_check != numcells)
	{
		numcells_check = numcells;
		sum_l_M = 0, sum_l_I = 0, sum_l_S = 0, sum_l_E = 0;
		sum_r_M = 0, sum_r_I = 0, sum_r_S = 0, sum_r_E = 0;
		sum_t_M = 0, sum_t_I = 0, sum_t_S = 0, sum_t_E = 0;
		sum_b_M = 0, sum_b_I = 0, sum_b_S = 0, sum_b_E = 0;
		check1 = 0;
		check2 = 0;
	}

	double dx = LENGTH / double(numcells);
	double *xx = new double[numcells];
	int omp_chuck = numcells / OMP_CORES;

	if (timer >= tt1 && timer <= tt2)
	{
		int l_bound = int((xx1 - 0.5*dx) / dx);
		int r_bound = int((xx2 - 0.5*dx) / dx);
		//printf("atata: %d %d\n", l_bound, r_bound);

		/*************massa*****************/

		sum_l_M += R[l_bound] * U[l_bound] * tau;
		sum_r_M += R[r_bound] * U[r_bound] * tau;

		/***********impulse****************/

		sum_l_I += (P[l_bound] + R[l_bound] * U[l_bound] * U[l_bound]) * tau;
		sum_r_I += (P[r_bound] + R[r_bound] * U[r_bound] * U[r_bound]) * tau;

		/***********entropy****************/

		sum_l_S += R[l_bound] * S[l_bound] * U[l_bound] * tau;
		sum_r_S += R[r_bound] * S[r_bound] * U[r_bound] * tau;

		/***********energy****************/

		sum_l_E += (RE[l_bound] + P[l_bound]) * U[l_bound] * tau;
		sum_r_E += (RE[r_bound] + P[r_bound]) * U[r_bound] * tau;
	}

#pragma omp parallel for num_threads(OMP_CORES) schedule(guided)
	for (int i = 0; i < numcells; i++)
	{
		xx[i] = i*dx + 0.5*dx;
	}

	if (timer >= tt1 && check1 == 0)
	{

#pragma omp parallel for simd reduction(+:sum_b_M,sum_b_I,sum_b_S,sum_b_E) schedule(guided) num_threads(OMP_CORES)
		for (int i = 0; i < numcells; i++)
		{
			if (xx[i] >= xx1 && xx[i] <= xx2)
			{
				sum_b_M += R[i] * dx;
				sum_b_I += R[i] * U[i] * dx;
				sum_b_S += R[i] * S[i] * dx;
				sum_b_E += RE[i] * dx;
			}
		}

		check1 = 1;
	}


	if (timer >= tt2 && check2 == 0)
	{
//#pragma omp parallel for simd reduction(+:sum_t_M,sum_t_I,sum_t_S,sum_t_E) schedule(guided) num_threads(OMP_CORES)
		for (int i = 0; i < numcells; i++)
		{
			if (xx[i] >= xx1 && xx[i] <= xx2)
			{
				sum_t_M += R[i] * dx;
				sum_t_I += R[i] * U[i] * dx;
				sum_t_S += R[i] * S[i] * dx;
				sum_t_E += RE[i] * dx;
			}
		}
		check2 = 1;
	}

	sum[0][0] = sum_t_M;
	sum[0][1] = sum_t_I;
	sum[0][2] = sum_t_S;
	sum[0][3] = sum_t_E;

	sum[1][0] = sum_b_M;
	sum[1][1] = sum_b_I;
	sum[1][2] = sum_b_S;
	sum[1][3] = sum_b_E;

	sum[2][0] = sum_r_M;
	sum[2][1] = sum_r_I;
	sum[2][2] = sum_r_S;
	sum[2][3] = sum_r_E;

	sum[3][0] = sum_l_M;
	sum[3][1] = sum_l_I;
	sum[3][2] = sum_l_S;
	sum[3][3] = sum_l_E;

	free(xx);
}

void inf_before_start(int numcells, double *R, double *U, double *P, double &D_analit)
{ 
#if PROBLEM==0
	D_analit = (R[numcells - 1] * U[numcells - 1] - R[0] * U[0]) / (R[numcells - 1] - R[0]);
	printf("Analitical speed: %10.8lf\n\n", D_analit);
#elif PROBLEM == 12
	double uc_left, uc_right;

	//	uc_left = U[0] +0.43 - sqrt(GAMMA*P[0] / R[0]);
	//uc_right = U[numcells] + 0.43  - sqrt(GAMMA*P[numcells] / R[numcells]);
	uc_left = U[0] - sqrt(GAMMA*P[0] / R[0]);
	uc_right = U[numcells - 5] - sqrt(GAMMA*P[numcells - 5] / R[numcells - 5]);
	printf("U-C left: %8.6lf\nU-C right: %8.6lf\nmiddle: %8.6lf\n", uc_left, uc_right, (uc_left + uc_right) / 2);
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
}


/**************************************************/

double* finite_difference(int numb, double *mas)
{
	int val;
	val = nmesh[numb];
	double *dif;
	dif = new double[val - 2];
	int check = 0;
	long int p = 0;
	long int m = 0;
	long int zero = 0;


	for (int i = 0; i < val - 2; i++)
	{
		dif[i] = mas[i + 2] - 2 * mas[i + 1] + mas[i];
		if (i == 0)
		{
			if (dif[i] > 0) { check = 1; p++; continue; }
			if (dif[i] < 0) { check = 0; m++; continue; }
			if (dif[i] == 0) { check = 2; zero++; continue; }
		}

		if (dif[i] > 0)
		{
			if (check == 1) { p++; check = 1; }
			else {
				p++;
				if (check == 0) { if (m > 30) { printf("-: %d\n", m); } m = 0; }
				if (check == 2) { if (zero > 30) { printf("0: %d\n", zero); } zero = 0; }
				check = 1;
			}

		}
		if (dif[i] < 0)
		{
			if (check == 0) { m++; check = 0; }
			else {
				m++;
				if (check == 1) { if (p > 30) { printf("+: %d\n", p); } p = 0; }
				if (check == 2) { if (zero > 30) { printf("0: %d\n", zero); } zero = 0; }
				check = 0;
			}
		}

		if (dif[i] == 0)
		{
			if (check == 2) { zero++; check = 2; }
			else {
				zero++;
				if (check == 0) { if (m > 30) { printf("-: %d\n", m); } m = 0; }
				if (check == 1) { if (p > 30) { printf("+: %d\n", p); } p = 0; }
				check = 2;
			}
		}

		if (i == val - 3)
		{
			if (check == 0) printf("-: %d\n", m);
			if (check == 1) printf("+: %d\n", p);
			if (check == 2) printf("0: %d\n", zero);
		}

	}

	return dif; // возвращаем указатель на массив!!! его имя!
}

void file_exact_diff(int numcells, double *exact_R, double *exact_U, double *exact_P, double *exact_RE, double *exact_S, double *diff_R, double *diff_U, double *diff_P, double time)
{
}

void runge(double *massiv, int lda, int numb)
{
	double value[NUM_ITER - 2]; // without boundary points
	for (int i = 0; i < NUM_ITER - 2; i++)
	{
		//value[i] = log(fabs((fabs(massiv[i]) - fabs(massiv[i + 1])) / (fabs(massiv[i + 1]) - fabs(massiv[i + 2])))) / log(double(GRID));
		value[i] = log(fabs((massiv[numb + lda * i] - massiv[numb + lda * (i + 1)]) / (massiv[numb + lda * (i + 1)] - massiv[numb + lda * (i + 2)]))) / log(double(GRID));

		printf("Precision order for iteration double %d: %lf\n", i + 1, value[i]);
	}
}

void analitical_writing_into_file(int numcells, double* R_D, double*R_U, double*R_P, double timer)
{
	FILE* fout;
	FILE* plot;
	double x, x_NC;
	double dx = LENGTH / double(numcells);
	char name1[255], name2[255];

	double D_analit;
	double ro_right, ro_left, u_right, u_left;

	ro_right = 1.271413930046081;
	ro_left = 1.0;
	u_right = 0.292868067614595;
	u_left = 0;

	D_analit = (ro_right * u_right - ro_left * u_left) / (ro_right - ro_left);

#ifndef NC
	sprintf(name1, "workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_analit_%6.4lf.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, timer);
#else
	sprintf(name1, "workspace/%03d/NC_N%03d_P%1d_SLV%1d_TERM%.0lf_analit_%6.4lf.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, timer);
#endif

	fout = fopen(name1, "w");

	for (int i = 0; i < numcells; i++)
	{

#ifndef NC
		x = i*dx + 0.5*dx;
#else
#ifdef NC2
		x_NC = (i*dx + 0.5*dx - DISC_POINT - D_analit*timer) / dx;
#else
		x_NC = i*dx + 0.5*dx - DISC_POINT - D_analit*timer;
#endif

#endif


		/**********************************
		| 1 |    2    |   3   |     4    |
		| x | density | speed | pressure |
		***********************************/
#ifndef NC
		fprintf(fout, "%9.6lf %lf %lf %lf \n", x, R_D[i], R_U[i], R_P[i]);
#else
		fprintf(fout, "%9.6lf %lf %lf %lf \n", x_NC, R_D[i], R_U[i], R_P[i]);
#endif
	}
	fclose(fout);
}



//#define BOUND_COND

/************************************************
**************************************************
***************************************************/
#ifdef OTKOL
void iteration_bound(int numb)
{
	double *R, *R1, *R2,		// density
		*P, *P1, *P2,		// pressure
		*U, *U1, *U2,		// velocity
		*RS, *S_a, *S_b, *S, *S_diff,
		*RU, *RU1, *RU2,		// moment of impulse
		*RE, *RE1, *RE2,		// total energy
	*UBOUND, *C;

	double *FR,		// density flux
		*FRU,	// moment flux
		*FRP,	// pressure gradient
		*FRE,	// energy flux
		*FRUS;     //entropy flux
	double *UFLUX;	// velocity flux
	int i, numcells, start_print, jump_print, left_index, right_index, imesh, k;
	double	timer, tau, umax, u1, u2, u3, dx, len, x, x_NC, ds = 0, us = 0, ps = 0, es, cs = 0, ss, term, D_num, rp, pg, hh, *uss, *pss, *dss;

	double *x_layer, *x_layer_NC, *E;

	double *all_exact_P, *all_exact_U, *all_exact_R, *all_exact_RE, *all_exact_S;
	double *diff_riem_P, *diff_riem_U, *diff_riem_R, *diff_riem_RE, *diff_riem_S;

	double delta_ro, delta_u, delta_p;

	double* right_parts;

	int* w_num_p, *w_num_r, *w_num_u;

	char FileName[255];
	char FileName2[255];
	char FileName3[255];
	char FileName4[255];

	FILE *fout, *fmesh;
	FILE *fout2, *fout3;
	FILE* fout_NC;
	FILE* f_flux, *f_flux2;

	//	double moving_mesh[11], moving_mesh1[11], moving_mesh2[11];

	double sum_m[4][4] = { 0 };

	int check1 = 0, check2 = 0, check3 = 0, check5[5] = { 0 };
	int proverka[N_smooth] = { 0 };

	double riemann_U, riemann_P, riemann_D[2];

	double *diff; //память под массив выделяется внутри функции

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
	R = new double[numcells];
	P = new double[numcells];
	U = new double[numcells];
	C = new double[numcells];
	S_b = new double[numcells];
	S_a = new double[numcells];
	S_diff = new double[numcells];
	S = new double[numcells];
	E = new double[numcells];
	RU = new double[numcells];
	RE = new double[numcells];
	RS = new double[numcells];

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
	int *sw1_num_p, *sw2_num_p, *sw3_num_p;
	int *sw1_num_u, *sw2_num_u, *sw3_num_u;
	int *sw1_num_r, *sw2_num_r, *sw3_num_r;
	sw1_num_p = new int[N_smooth]; sw2_num_p = new int[N_smooth]; sw3_num_p = new int[N_smooth];
	sw1_num_r = new int[N_smooth]; sw2_num_r = new int[N_smooth]; sw3_num_r = new int[N_smooth];
	sw1_num_u = new int[N_smooth]; sw2_num_u = new int[N_smooth]; sw3_num_u = new int[N_smooth];



#ifdef NEW_VERSION
	UBOUND = new double[numcells + 1];
	FR = new double[numcells + 1];
	FRU = new double[numcells + 1];
	FRP = new double[numcells + 1];
	FRE = new double[numcells + 1];
	FRUS = new double[numcells + 1];
	UFLUX = new double[numcells + 1];
	dss = new double[numcells + 1];
	uss = new double[numcells + 1];
	pss = new double[numcells + 1];
#else
	FR = new double[numcells - 1];
	FRU = new double[numcells - 1];
	FRP = new double[numcells - 1];
	FRE = new double[numcells - 1];
	FRUS = new double[numcells - 1];
	UFLUX = new double[numcells - 1];
	dss = new double[numcells - 1];
	uss = new double[numcells - 1];
	pss = new double[numcells - 1];
#endif
	x_layer = new double[numcells];

	double max_x2 = 0;
	double pg_max = 0;

	double D_analit = 0;

	double es_av, es_diff;
	double l1, r1, l2, r2;

	w_num_p[0:N_smooth] = 0; w_num_r[0:N_smooth] = 0;  w_num_u[0:N_smooth] = 0;

	int count2 = 0;

	/*********************************THE PROGRAM BELOW******************************/

	FILE* entropy, *file3;
	char buf[255];


	/* Initial conditions */
#ifdef BOOST
#pragma loop count min(64)
#pragma unroll(5)
#endif
	// vectorized and auto_parallized
#ifdef OMP
#pragma omp parallel for simd schedule(guided)
#endif
	for (i = 0; i < numcells; i++)
	{
		x_layer[i] = i*dx + 0.5*dx;          // that are middles of cells
		R[i] = initial_density(x_layer[i]);
		P[i] = initial_pressure(x_layer[i]);
		U[i] = initial_velocity(x_layer[i]);
	}


#ifdef BOOST
#pragma parallel always
#endif
#ifdef OMP
#pragma omp parallel for simd schedule(guided) 	// vectorized and auto_parallized
#endif
	for (i = 0; i < numcells; i++)
	{
		RU[i] = R[i] * U[i];
		RE[i] = (P[i]-(R[i]-R0) * C0 * C0) / (GAMMA - 1.0) + 0.5 * R[i] * U[i] * U[i]; // full mechanic energy
	}

	printf("\nIteration %d\n", numb + 1);

	timer = 0.0;
	double time_max;
	double explosion;
	time_max = time_max_array[PROBLEM];

#ifdef OUTPUT_N_SMOOTH

	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
	{
		time_control[i] = (i + 1)*k_step;
		//		printf("%lf\n", time_control[i]);
	}
#endif

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
	/***** Вычислительный цикл метода Годунова *****/

	// определяем характеристику
	double s_char = 0.0;
	double pm = 0, um = 0;

	FILE* array_flux[N_smooth];
	for (int i = 0; i < N_smooth; i++)
	{
		sprintf(FileName4, "FLUXES_%d_P%d_N%d.dat", i, int(PROBLEM), int(numcells));
		array_flux[i] = fopen(FileName4, "w");
	}

	while (timer < time_max)
		/************************************* The beggining of iteration*****************************/
	{
		c_c++;
		//printf("c_c: %d %lf %lf %lf\n", c_c, timer, time_max, dx);

		// cfl condition  (условие Курранта)
		umax = 0.0;

		//dependence! not vectorized and parallelized!
		for (int i = 0; i < numcells; i++)// not parallel and vectorized: dependences
		{

			u1 = fabs(U[i] + sqrt(GAMMA*P[i] / R[i]));
			u2 = fabs(U[i] - sqrt(GAMMA*P[i] / R[i]));
			u3 = fabs(U[i]);
			if (u1 > umax) umax = u1;
			if (u2 > umax) umax = u2;
			if (u3 > umax) umax = u3;
		}
#ifdef CFL_SWITCH
		if (timer < time_max / 2)
		{
			tau = CFL08*dx / umax;
		}
		else
		{
			tau = CFL04*dx / umax;
		}
#else
		tau = CFL04*dx / umax;
#endif

#ifdef PRINT_TAU
		printf("tau: %4.3d %4.3d %4.3d %10.9lf %5.4lf\n", new_timer, count_razmaz, schet, tau, tau / dx);
#endif
		if (timer + tau > time_max)
		{
			tau = time_max - timer; // if the last time's step is bigger than distance to "time_max"
		}


		/* Boundary conditions */
#ifdef NEW_VERSION
#ifdef BOUND_COND
#if (PROBLEM==18 || PROBLEM==20)
		// set pressure
		pss[0] = initial_pressure(0.0);
		//	pss[numcells] = initial_pressure(LENGTH);

		C[0] = sqrt(GAMMA*P[0] / R[0]);
		//	C[numcells] = sqrt(GAMMA*P[numcells] / R[numcells]);

		double l0_const = U[0] - P[0] / (R[0] * C[0]);
		double rn_const = U[numcells] + P[numcells] / (R[numcells] * C[numcells]);

		uss[0] = l0_const + pss[0] / (R[0] * C[0]);
		//	uss[numcells] = rn_const - pss[numcells] / (R[numcells] * C[numcells]);

		S[0] = log(P[0] / pow(R[0], GAMMA));
		//	S[numcells] = log(P[numcells] / pow(R[numcells], GAMMA));

		dss[0] = pow(pss[0] / S[0], 1.0 / GAMMA);
		//	dss[numcells] = pow(pss[numcells] / S[numcells], 1.0 / GAMMA);

#elif (PROBLEM == 19)

		uss[0] = timer;
		C[0] = sqrt(GAMMA*P[0] / R[0]);

		double l0_const = U[0] - P[0] / (R[0] * C[0]);
		double rn_const = U[numcells] + P[numcells] / (R[numcells] * C[numcells]);

		pss[0] = (uss[0] - l0_const)*(R[0] * C[0]);

		S[0] = log(P[0] / pow(R[0], GAMMA));
		dss[0] = pow(pss[0] / S[0], 1.0 / GAMMA);
#elif (PROBLEM == 4)
		// set pressure
		pss[0] = 4.0 / exp(3.0*timer / time_max);
		//	pss[numcells] = initial_pressure(LENGTH);

		C[0] = sqrt(GAMMA*P[0] / R[0]);
		//	C[numcells] = sqrt(GAMMA*P[numcells] / R[numcells]);

		double l0_const = U[0] - P[0] / (R[0] * C[0]);
		double rn_const = U[numcells] + P[numcells] / (R[numcells] * C[numcells]);

		uss[0] = l0_const + pss[0] / (R[0] * C[0]);
		//	uss[numcells] = rn_const - pss[numcells] / (R[numcells] * C[numcells]);

		S[0] = log(P[0] / pow(R[0], GAMMA));
		//	S[numcells] = P[numcells] / pow(R[numcells], GAMMA);

		dss[0] = pow(pss[0] / S[0], 1.0 / GAMMA);
		//	dss[numcells] = pow(pss[numcells] / S[numcells], 1.0 / GAMMA);

#endif

	/*	pss[numcells] = initial_pressure(LENGTH);
		uss[numcells] = initial_velocity(LENGTH);
		dss[numcells] = initial_density(LENGTH);*/

		linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);

		FR[0] = dss[0] * uss[0];
		FRU[0] = dss[0] * uss[0] * uss[0] + pss[0];
		FRE[0] = ((pss[0] - (dss[0] - R0) * C0 * C0) / (GAMMA - 1.0) + 0.5*dss[0] * uss[0] * uss[0])*uss[0] + pss[0] * uss[0];
		FR[numcells] = dss[numcells] * uss[numcells];
		FRU[numcells] = dss[numcells] * uss[numcells] * uss[numcells] + pss[numcells];
		FRE[numcells] = ((pss[numcells] - (dss[numcells] - R0) * C0 * C0) / (GAMMA - 1.0) + 0.5*dss[numcells] * uss[numcells] * uss[numcells])*uss[numcells] + pss[numcells] * uss[numcells];

		UFLUX[0] = uss[0];
		UFLUX[numcells] = uss[numcells];	

#endif
#endif
		/* Boundary conditions */
		// shock wave configuration
#ifdef BOOST_VEC
#pragma distribute point
#endif
#ifdef NEW_VERSION
		if (timer < CROSS_POINT)
		{
			//====================================== EXACT RIEMANN PROBLEM =========================
			for (int i = 1; i < numcells; i++)
			{
				// левая сторона
				dl = R[i - 1];
				pl = P[i - 1];
				ul = U[i - 1];
				cl = sqrt(GAMMA*pl / dl);
				// правая сторона
				dr = R[i];
				pr = P[i];
				ur = U[i];
				cr = sqrt(GAMMA*pr / dr);

				starpu(pm, um);

				// Решение задачи распада разрыва
				sample(pm, um, s_char, ds, us, ps);
				UFLUX[i] = us;
				FR[i] = ds * us;
				FRU[i] = ds * us * us + ps;
				FRE[i] = ((ps - (ds - R0) * C0 * C0) / (GAMMA - 1.0) + 0.5*ds * us * us)*us + ps * us;
				//		ss = ps / pow(ds, GAMMA);
				FRP[i] = ps;
			}
		}
		else
		{
#ifdef OMP
#pragma distribute point
#pragma omp parallel for simd schedule(guided) // vectorized and auto_parallelized
#endif
			for (int i = 1; i < numcells; i++)
			{
				linear(R[i - 1], U[i - 1], P[i - 1], R[i], U[i], P[i], dss[i], uss[i], pss[i]);      // computation of the Riemann problem in each two cells
		
				UFLUX[i] = uss[i];
				FR[i] = dss[i] * uss[i];
				FRU[i] = dss[i] * uss[i] * uss[i] + pss[i];
				FRE[i] = ((pss[i] - (dss[i] - R0) * C0 * C0) / (GAMMA - 1.0) + 0.5*dss[i] * uss[i] * uss[i])*uss[i] + pss[i] * uss[i];
				FRP[i] = pss[i];
			}

		}
#else
#ifdef EXACT_SOLUTION
		for (int i = 1; i < numcells; i++)
		{
			// левая сторона
			dl = R[i - 1];
			pl = P[i - 1];
			ul = U[i - 1];
			cl = sqrt(gamma*pl / dl);
			// правая сторона
			dr = R[i];
			pr = P[i];
			ur = U[i];
			cr = sqrt(gamma*pr / dr);
			/*			if (c_c == 1 && i == 1)
			{
			printf("dl,pl,ul,cl: %lf %lf %lf %lf\ndlr,pr,ur,cr: %lf %lf %lf %lf\n", dl, pl, ul, cl, dr, pr, ur, cr);
			printf("one %lf %lf %lf %lf %lf %lf\n", pm, um, s_char, ds, us, ps);
			}*/
			starpu(pm, um);
			// Решение задачи распада разрыва
			//			if (c_c == 1 && i == 1) printf("two %lf %lf %lf %lf %lf %lf\n", pm, um, s_char, ds, us, ps);
			sample(pm, um, s_char, ds, us, ps);
			/*			if (c_c == 1 && i == 1)
			{
			printf("thr %lf %lf %lf %lf %lf %lf\n", pm, um, s_char, ds, us, ps);
			system("pause");
			}*/

			UFLUX[i - 1] = us;
			FR[i - 1] = ds * us;
			FRU[i - 1] = ds * us * us + ps;
			FRP[i - 1] = ps;
			FRE[i - 1] = (ps / (GAMMA - 1.0) + 0.5*ds * us * us)*us + ps * us;
#else
#ifdef OMP
#pragma distribute point
#pragma omp parallel for simd schedule(guided) // vectorized and auto_parallelized
#endif
		for (int i = 1; i < numcells; i++)
		{
			linear(R[i - 1], U[i - 1], P[i - 1], R[i], U[i], P[i], dss[i - 1], uss[i - 1], pss[i - 1]);      // computation of the Riemann problem in each two cells

			// finding: Pressure, Density, Velocity on the boundary of each sell
			//			printf("i: %d r: %lf  u: %lf  p: %lf\n", i, dss[i - 1], uss[i - 1], pss[i - 1]);
			UFLUX[i - 1] = uss[i - 1];
			FR[i - 1] = dss[i - 1] * uss[i - 1];
			FRU[i - 1] = dss[i - 1] * uss[i - 1] * uss[i - 1] + pss[i - 1];
			FRP[i - 1] = pss[i - 1];
			FRE[i - 1] = (pss[i - 1] / (GAMMA - 1.0) + 0.5*dss[i - 1] * uss[i - 1] * uss[i - 1])*uss[i - 1] + pss[i - 1] * uss[i - 1];
#endif
		}
#endif


#ifdef NEW_VERSION		
#ifndef BOUND_COND
#if PROBLEM!=3
		linear(R[0], U[0], P[0], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
		linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);
		FR[0] = dss[0] * uss[0];
		FRU[0] = dss[0] * uss[0] * uss[0] + pss[0];
		FRE[0] = ((pss[0] - (dss[0] - R0) * C0 * C0) / (GAMMA - 1.0) + 0.5*dss[0] * uss[0] * uss[0])*uss[0] + pss[0] * uss[0];
		FR[numcells] = dss[numcells] * uss[numcells];
		FRU[numcells] = dss[numcells] * uss[numcells] * uss[numcells] + pss[numcells];
		FRE[numcells] = ((pss[numcells] - (dss[numcells] - R0) * C0 * C0) / (GAMMA - 1.0) + 0.5*dss[numcells] * uss[numcells] * uss[numcells])*uss[numcells] + pss[numcells] * uss[numcells];

#else
		linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
		linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[numcells], uss[numcells], pss[numcells]);
		FR[0] = dss[0] * uss[0];
		FRU[0] = dss[0] * uss[0] * uss[0] + pss[0];
		FRE[0] = (pss[0] / (GAMMA - 1.0) + 0.5*dss[0] * uss[0] * uss[0])*uss[0] + pss[0] * uss[0];
		FR[numcells] = dss[numcells] * uss[numcells];
		FRU[numcells] = dss[numcells] * uss[numcells] * uss[numcells] + pss[numcells];
		FRE[numcells] = (pss[numcells] / (GAMMA - 1.0) + 0.5*dss[numcells] * uss[numcells] * uss[numcells])*uss[numcells] + pss[numcells] * uss[numcells];
#endif
#endif
#endif

		//Comparing u and bigU
/*		if (c_c == 10)
		{
			file3 = fopen("comparing_U.dat", "w+");
			for (int i = 0; i < numcells; i++)
			{
				x_layer[i] = i*dx + 0.5*dx;
				fprintf(file3, "%lf %lf %lf\n", x_layer[i], U[i], UFLUX[i]);
			}
			fclose(file3);
		}*/
		//system("pause");

		/****************************************************************
		IMPORTANT: Solving of this problem is in conservative variables!
		(R), impulse (RU), energy (RE)
		****************************************************************/

		/* Computations of conservations laws*/
#ifdef BOOST
#pragma loop count min(32)
#pragma parallel always
#endif	
#ifdef OMP
		//#pragma distribute point
#pragma omp parallel for simd schedule(guided) // vectorized and auto_parallelized
#endif
#ifdef NEW_VERSION
		for (int i = 0; i < numcells; i++)
		{
			R[i] = R[i] - tau * (FR[i + 1] - FR[i]) / dx;
#ifdef RP
			term = -U[i] * pow(dx, 2.0*K_TERM - 1.0)*pow((FRP[i + 1] - FRP[i]) / dx, 2.0*K_TERM);  // it's definition of the right part of equations
			right_parts[i] = tau*term;
			RU[i] = RU[i] - tau * (FRU[i + 1] - FRU[i]) / dx + right_parts[i];  // new, modificated laws (with right part)
#else
			RU[i] = RU[i] - tau * (FRU[i + 1] - FRU[i]) / dx;  // without right parts
#endif
			RE[i] = RE[i] - tau * (FRE[i + 1] - FRE[i]) / dx;
		}
#else
		//		printf("\n\n");
		for (int i = 1; i < numcells - 1; i++)
		{
			R[i] = R[i] - tau * (FR[i] - FR[i - 1]) / dx;
#ifdef RP
			term = -U[i] * pow(dx, 2.0*K_TERM - 1.0)*pow((FRP[i] - FRP[i - 1]) / dx, 2.0*K_TERM);  // it's definition of the right part of equations
			right_parts[i] = tau*term;
			RU[i] = RU[i] - tau * (FRU[i] - FRU[i - 1]) / dx + right_parts[i];  // new, modificated laws (with right part)
#else
			RU[i] = RU[i] - tau * (FRU[i] - FRU[i - 1]) / dx;  // without right parts

#endif
			RE[i] = RE[i] - tau * (FRE[i] - FRE[i - 1]) / dx;
			//		printf("i: %d r: %lf  ru: %lf  re: %lf\n", i, R[i], RU[i], RE[i]);
		}
#endif
	
#ifndef NEW_VERSION
		// boundary has already setted up
		R[0] = R[left_index];			  // R[0] - left boundary                   lenght(R)=numcells
		R[numcells - 1] = R[right_index];   // R[numcells-1] - right boundary
		RU[0] = RU[left_index];
		RU[numcells - 1] = RU[right_index];
		RE[0] = RE[left_index];
		RE[numcells - 1] = RE[right_index];
#endif


		// moment -> velocity
#ifdef BOOST
#pragma loop count min(16)
#pragma parallel always
#endif
		// vectorized and auto_parallelized
#ifdef OMP
#pragma omp parallel for simd schedule(guided)
#endif
		for (int i = 0; i < numcells; i++)    // over-computation of velocity and pressure for the next step ( n --> n+1 ) - as R, RU, RE (n+1 step)
		{
			U[i] = RU[i] / R[i]; //используется дальше в программе!!!!
			P[i] = (GAMMA - 1.0) * (RE[i] - 0.5 * RU[i] * U[i]) + (R[i] - R0) * C0 * C0;
			S[i] = log(P[i] / pow(R[i], GAMMA));
		}


		//************ CURRENT DOMAIN **********************//

		// распад разрыва
		/* 1)   x=[0;0.46]
		2)   x=[0.46;0.61]
		3)   x=[0.61;1]    */
		//const double xx1 = 0.05, xx2 = 0.46, tt1 = 0.14; 0.295
		//const double xx1 = 0.05, xx2 = 0.56, tt1 = 0.14; 0.295

		// new calculations until t=0.3; no boundaries; no right parts

		/* Rarification wave

		1) x=[0;0.36]
		2) x=[0.36;0.42]
		3) x=[0:51;0:64] t=0.1
		x=[0:56;0.78] t=0.2
		*/

		/* 2SW: 0.3-0.5 left, 0.65-0.9 right */
		/* for RW: 0-1 distance */

		/* очень хороший случай: 0.2-1, 0.245-0.5*/
		//const double xx1 = 8, xx2 = 17, tt1 = 2.0;
#if (PROBLEM==17)
		const double xx1 = 0.0, xx2 = 10.0, tt1 = 2.8;  // 0.05. 0.25. 0.4 - норм
#else
		const double xx1 = 0.00, xx2 = 0.95, tt1 = 0.35;   // 0.05. 0.25. 0.4 - норм
#endif


		//const double xx1 = 0.05, xx2 = 0.36, tt1 = 0.148; -- удиничные порядки

		// контур вокруг столкновения двух волн: 8, 17, 2.0, 5.0
		// столкновение сильных волн (P_8): x[0.3, 0.8], y[0.16, 0.3]
		// одна новая волна (P_8): x[0.7,0.95], y[0.2875,0.3625]

		// 6 задача: 13, 18, 4.0, 7.0
		// 8 задача: 8, 18, 0, 3

#if (PROBLEM==0 || PROBLEM==1)
		const double tt2 = 0.495;  // берем не самую верхнюю границу (0.5), так как она не считается!!! она равна нулю по моему алгоритму реализации
#elif (PROBLEM==17)
		const double tt2 = 3.4;
#else
		const double tt2 = 0.495;
#endif

		//************ CURRENT DOMAIN **********************//

		//	analitical_riemann_modeling(numcells, initial_density(0.1), initial_velocity(0.1), initial_pressure(0.1), initial_density(0.9), initial_velocity(0.9), initial_pressure(0.9), timer, all_exact_R, all_exact_U, all_exact_P);

		timer += tau;

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

#pragma noparallel
		for (int k = 0; k < N_smooth; k++)
		{

			// точности: 0.005 - for problem 1
			// 0.0005 - for problem 2
			//	if (fabs(timer - time_control[k]) <= 0.005 && proverka[k] == 0)
			if (fabs(timer - time_control[k]) <= tau && proverka[k] == 0)
			{

				sprintf(FileName, "workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4lf.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (char)TYPE, time_control[k]);
				sprintf(FileName2, "workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4lf_NC.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (char)TYPE, time_control[k]);

#ifdef TIME				
				if (numb == TIME && (fabs(time_control[k] - 0.0015) < k_step || fabs(time_control[k] - 0.0035) < k_step || fabs(time_control[k] - 0.015) < k_step || fabs(time_control[k] - 0.08) < k_step))
				{
#endif
					fout = fopen(FileName, "w");
					fprintf(fout, "Timer: %lf\n", timer);

#ifdef NC
					fout_NC = fopen(FileName2, "w");
					fprintf(fout_NC, "Timer: %lf\n", timer);
#endif
#pragma noparallel
					for (int i = 0; i < numcells; i++)  // вывод всегда по 100 точек ( с первой итерации которые )
					{
#ifndef NC
						//x_layer[i] = i*dx + 0.5*dx-0.4533*timer;
						//x_layer[i] = i*dx + 0.5*dx;

						//x_layer[i] = i*dx + 0.5*dx+UFLUX[i]*timer;
#if PROBLEM==19
						x_layer[i] = i*dx + 0.5*dx+timer*timer;
						
#elif PROBLEM==18
						//x_layer[i] = i*dx + 0.5*dx+initial_velocity(0.0)*timer;
						x_layer[i] = i*dx + 0.5*dx;
#else
						x_layer[i] = i*dx + 0.5*dx;
#endif

#else
						//x_layer_NC[i] = i*dx + 0.5*dx - D_analit*timer - DISC_POINT;      //без деления на t^alpha
#ifdef alpha
						x_layer_NC[i] = (i*dx + 0.5*dx - D_analit*timer - DISC_POINT) / (C1*pow(timer, alpha));
#endif

#ifdef NC2
						x_layer_NC[i] = (i*dx + 0.5*dx - DISC_POINT - D_analit*timer) / dx;
#else
						x_layer_NC[i] = (i*dx + 0.5*dx - DISC_POINT - D_analit*timer);
						//	x_layer_NC[i] = (i*dx + 0.5*dx - 0.43*timer);
#endif
#endif

						//	x_layer_NC[i] = (i*dx + 0.5*dx - DISC_POINT) / (D_analit*pow(timer, alpha)); //еще один случай

						ds = R[i];
#if PROBLEM == 12
						us = U[i] - 0.4533;
#else
						us = U[i];
#endif
						ps = P[i];
						cs = sqrt(GAMMA*P[i] / R[i]);
						es = log(P[i] / pow(R[i], GAMMA));
						if (i != 0 && i != numcells - 1) es_av = (P[i - 1] / pow(R[i - 1], GAMMA) + P[i] / pow(R[i], GAMMA) + P[i + 1] / pow(R[i + 1], GAMMA)) / 3.0;
						else if (i == 0) es_av = (P[i] / pow(R[i], GAMMA) + P[i + 1] / pow(R[i + 1], GAMMA)) / 2.0;
						else es_av = (P[i - 1] / pow(R[i - 1], GAMMA) + P[i] / pow(R[i], GAMMA)) / 2.0;
						es_diff = S_diff[i];
						//rp = us + cs;

						//проверка толщины ударной волны
#ifdef TIME
						if (numb == TIME)
						{
							if (P[i] < (initial_pressure(0.05) - DELTA) && P[i] > (initial_pressure(0.95) + DELTA)) w_num_p[k]++;
							if (R[i] < (initial_density(0.05) - DELTA) && R[i] > (initial_density(0.95) + DELTA)) w_num_r[k]++;
							if (U[i] < (initial_velocity(0.05) - DELTA) && U[i] > (initial_velocity(0.95) + DELTA)) w_num_u[k]++;
						}
#endif
#ifdef SW_POINTS_PRINT
						//-----------------------------


#if (PROBLEM == 0 || PROBLEM == 1)
						if (P[i] < (initial_pressure(0.05) - DELTA) && P[i] > (initial_pressure(0.95) + DELTA)) w_num_p[k]++;
						if (R[i] < (initial_density(0.05) - DELTA) && R[i] > (initial_density(0.95) + DELTA)) w_num_r[k]++;
						if (U[i] < (initial_velocity(0.05) - DELTA) && U[i] > (initial_velocity(0.95) + DELTA)) w_num_u[k]++;
#elif (PROBLEM == 2)
						// shock wave checking
						if (P[i] < (1.401789770179879 - DELTA) && P[i] > (1.0 + DELTA)) w_num_p[k]++;
						if (R[i] < (1.271413930046081 - DELTA) && R[i] > (1.0 + DELTA)) w_num_r[k]++;
						if (U[i] < (0.292868067614595 - DELTA) && U[i] > (0.0 + DELTA) && i>numcells / 2) w_num_u[k]++;
#elif (PROBLEM == 8)
						if (P[i] < (st_P1 - DELTA - 0.005) && P[i] > (st_P2 + DELTA)) sw1_num_p[k]++;
						if (P[i] < (st_P2 - DELTA - 0.005) && P[i] > (st_P3 + DELTA)) sw2_num_p[k]++;
						if (P[i] < (3.2627) && P[i] > (st_P3 + DELTA)) sw3_num_p[k]++;

						if (R[i] < (st_R1 - DELTA) && R[i] > (st_R2 + DELTA)) sw1_num_r[k]++;
						if (R[i] < (st_R2 - DELTA) && R[i] > (st_R3 + DELTA)) sw2_num_r[k]++;
						if (R[i] < (st_R1 - 0.03) && R[i] > (st_R3 + DELTA)) sw3_num_r[k]++;

						if (U[i] < (st_U1 - DELTA) && U[i] > (st_U2 + DELTA)) sw1_num_u[k]++;
						if (U[i] < (st_U2 - DELTA) && U[i] > (st_U3 + DELTA)) sw2_num_u[k]++;
						if (U[i] < (st_U1 - 0.03) && U[i] > (st_U3 + DELTA)) sw3_num_u[k]++;
#elif (PROBLEM == 4)
						if (P[i] < (st_P1 - DELTA - 0.005) && P[i] > (st_P2 + DELTA)) sw1_num_p[k]++;
						if (P[i] < (st_P2 - DELTA - 0.005) && P[i] > (st_P3 + DELTA)) sw2_num_p[k]++;
						if (P[i] < (st_P1 - DELTA - 0.02) && P[i] > (st_P3 + DELTA)) sw3_num_p[k]++;

						if (R[i] < (st_R1 - DELTA) && R[i] > (st_R2 + DELTA)) sw1_num_r[k]++;
						if (R[i] < (st_R2 - DELTA) && R[i] > (st_R3 + DELTA)) sw2_num_r[k]++;
						if (R[i] < (st_R1 - 0.03) && R[i] > (st_R3 + DELTA)) sw3_num_r[k]++;

						if (U[i] < (st_U1 - DELTA) && U[i] > (st_U2 + DELTA)) sw1_num_u[k]++;
						if (U[i] < (st_U2 - DELTA) && U[i] > (st_U3 + DELTA)) sw2_num_u[k]++;
						if (U[i] < (st_U1 - 0.03) && U[i] > (st_U3 + DELTA)) sw3_num_u[k]++;


#endif
#endif

						//-----------------------------

						fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf %lf\n", x_layer[i], ds, us, ps, cs, es, es_diff);
#ifdef NC
						fprintf(fout_NC, "%9.6lf %lf %lf %lf %lf %lf \n", x_layer_NC[i], ds, us, ps, cs, es);
#endif

					}

					fclose(fout);
#ifdef NC
					fclose(fout_NC);

#ifdef TIME
				}
#endif
#endif
#if(PROBLEM==2 || PROBLEM==9)
				analitical_riemann_modeling(numcells, initial_density(0.05), initial_velocity(0.05), initial_pressure(0.05), initial_density(0.95), initial_velocity(0.95), initial_pressure(0.95), timer, all_exact_R, all_exact_U, all_exact_P);
				analitical_writing_into_file(numcells, all_exact_R, all_exact_U, all_exact_P, time_control[k]);
#endif
				proverka[k] = 1;

			}
		}

#endif

		/************** расчет движения потоков ********************/
		double t[N_smooth] = { 0 };
		int t_ind[N_smooth] = { 0 };
		int numcells_flux;
		if (PROBLEM == 19) numcells_flux = numcells;
		else numcells_flux = numcells;

		for (int i = 0; i < N_smooth; i++)
		{
			t_ind[i] =  i * numcells_flux / N_smooth;
			t[i] = (t_ind[i] + 0.5)*dx + UFLUX[t_ind[i]] * timer;
			fprintf(array_flux[i], "%lf %lf %lf\n", t[i], timer, UFLUX[t_ind[i]]);
		}


	} /******************************************* The end of iteration**************************/

	for (int i = 0; i < N_smooth; i++)
	fclose(array_flux[i]);

#ifdef OUTPUT_N_SMOOTH
#ifdef PRINT
#if(PROBLEM==2 || PROBLEM==9)
	//gnuplot_analitical_riemann2(numcells, w_num_r, w_num_u, w_num_p);
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

	free(diff_riem_P);
	free(diff_riem_R);
	free(diff_riem_U);
	free(diff_riem_RE);
	free(diff_riem_S);

	free(FR);
	free(FRU);
	free(FRP);
	free(FRE);
	free(FRUS);
	free(UFLUX);
	free(dss);
	free(uss);
	free(pss);

#ifdef ENTROPY_CHECK
	if (numcells == 100 || numcells == 300)
	{
		//	fclose(entropy);
	}
#endif

#ifdef OUTPUT_N_SMOOTH
	FILE *sw_wide;
#ifdef TIME
	if (numb == TIME)
	{
		sprintf(FileName2, "sw_width_N%04d_delta%4.3lf_Nsmooth%05d.dat", numcells, DELTA, N_smooth);
		sw_wide = fopen(FileName2, "w");
		for (int k = 0; k < N_smooth; k++)
		{
			fprintf(sw_wide, "%lf %d %d %d\n", time_control[k], rw_num_r[k], rw_num_u[k], rw_num_p[k]);
		}
		fclose(sw_wide);
	}
#endif
#endif

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
#endif
#endif
	/*********Итог счета интегралов по контуру**********/

	free(all_exact_U);
	free(all_exact_R);
	free(all_exact_P);
	free(all_exact_RE);
	free(all_exact_S);

	/*********Итог счета интегралов по контуру**********/
#ifdef INTEGRAL
	F_ro_I[numb] = sum_m[0][1] - sum_m[1][1] + sum_m[2][1] - sum_m[3][1];
	printf("sum_l_I: %30.28lf\nsum_t_I: %30.28lf\nsum_r_I: %30.28lf\nsum_b_I: %30.28lf\n", sum_m[3][1], sum_m[0][1], sum_m[2][1], sum_m[1][1]);
	F_ro_M[numb] = sum_m[0][0] - sum_m[1][0] + sum_m[2][0] - sum_m[3][0];
	printf("\nsum_l_M: %30.28lf\nsum_t_M: %30.28lf\nsum_r_M: %30.28lf\nsum_b_M: %30.28lf\n", sum_m[3][0], sum_m[0][0], sum_m[2][0], sum_m[1][0]);
	F_ro_S[numb] = sum_m[0][2] - sum_m[1][2] + sum_m[2][2] - sum_m[3][2];
	printf("\nsum_l_S: %30.28lf\nsum_t_S: %30.28lf\nsum_r_S: %30.28lf\nsum_b_S: %30.28lf\n", sum_m[3][2], sum_m[0][2], sum_m[2][2], sum_m[1][2]);
	F_ro_E[numb] = sum_m[0][3] - sum_m[1][3] + sum_m[2][3] - sum_m[3][3];
#endif
	/*********Итог счета интегралов по контуру**********/


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
		ss = log(P[i] / pow(R[i], GAMMA));
		rp = right_parts[i];


#ifndef NC
		fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf %40.38lf \n", x, ds, us, ps, cs, ss, rp);
#else 
		//fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf %30.24lf %30.24lf %30.24lf \n", x_NC, ds, us, ps, cs, es, rp, pg, hh);
		fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf\n", x_layer_NC[i], ds, us, ps, cs, ss);
#endif

	}

#endif
	/************************** Analysis after count on one iteration *********************/

#ifdef P_PLUS_PG
	sprintf(FileName, "N%04d_P%1d_SLV%1d_TERM%.0lf_P_PLUS_PG_NC.dat", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);
	fout3 = fopen(FileName, "w");
	for (i = 0; i < numcells; i++)  // вывод всегда по 100 точек ( с первой итерации которые )
	{
		x_NC = (i*dx + 0.5*dx) - D_analit*time_max_array[0] * (numb + 1);  //NC
		ps = P[i];
		pg = pressure_gradient[i];

		fprintf(fout3, "%9.6lf %30.24lf %30.24lf \n", x_NC, ps, pg / pg_max);
	}
#endif

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
	/* Rarify wave - numeric vs analitic */
	double c, l0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;
	double iu_l, id_l, ip_l, iu_r, id_r, ip_r;
	double xl, xr;

	iu_l = initial_velocity(0.05);
	id_l = initial_density(0.05);
	ip_l = initial_pressure(0.05);

	iu_r = initial_velocity(0.2);
	id_r = initial_density(0.2);
	ip_r = initial_pressure(0.2);

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	tt = (numb + 1)*0.1;

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	xl = (iu_l - q11)*tt / q3 + 0.1; // счет по старому, по итерациям - какая итерация, такое и время. 
	xr = (iu_r - q11)*tt / q3 + 0.1; // в нашем случае ЭТО НЕВЕРНО, так как теперь итерации отвечают только за ШАГ СЕТКИ, время должно быть ФИКСИРОВАНО
	// здесь ошибка при вычислении точного решения!!!


	// 0-U, 1-P, 2-R, 3-RU, 4-RE
	double** difference_RW;
	difference_RW = new double*[5];

	for (int j = 0; j < 5; j++)
		difference_RW[j] = new double[numcells];

	double* RW_R, *RW_P, *RW_U;
	RW_R = new double[numcells];
	RW_P = new double[numcells];
	RW_U = new double[numcells];

	int counter = 0, counter_all = 0, counter2 = 0;


	//	cilk::reducer_opadd<int>counter_all(0);
	//	cilk::reducer_opadd<int>counter(0);
	for (i = 0; i<numcells; i++)
	{
		x = i*dx + 0.5*dx;

		if (x < xl)
		{
			RW_U[i] = iu_l; // FOR RARIFY WAVE
			RW_P[i] = ip_l;
			RW_R[i] = id_l;
			difference_RW[0][i] = U[i] - RW_U[i];
			difference_RW[1][i] = P[i] - RW_P[i];
			difference_RW[2][i] = R[i] - RW_R[i];
			difference_RW[3][i] = 0;
			difference_RW[4][i] = 0;
			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter2++;
			}

		}
		if (x >= xl && x <= xr)
		{
			counter_all++;
			RW_U[i] = RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[0][i] = U[i] - RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			RW_P[i] = RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[1][i] = P[i] - RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter++;
				counter2++;
			}

			RW_R[i] = RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[2][i] = R[i] - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			difference_RW[3][i] = R[i] * U[i] - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)*RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[4][i] = R[i] * (P[i] / pow(R[i], GAMMA) + SQ_2(U[i])) - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)*(RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r) / pow(RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r), GAMMA) + SQ_2(RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)));
		}
		if (x>xr)
		{
			RW_U[i] = iu_r;
			RW_P[i] = ip_r;
			RW_R[i] = id_r;

			difference_RW[0][i] = U[i] - RW_U[i];
			difference_RW[1][i] = P[i] - RW_P[i];
			difference_RW[2][i] = R[i] - RW_R[i];
			difference_RW[3][i] = 0;
			difference_RW[4][i] = 0;
			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter2++;
			}

		}

	}

	check1 = 0;
	check2 = 0;

	int *i_helper;
	i_helper = new int[10];

	double *x_helper;
	x_helper = new double[10];

	double xl_num, xr_num;

	/****************Boundary of numerical rarify wave******************/

	int i_mem_left = 0, i_mem_right = 0;

	for (i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;
		if (x >= xl && check1 == 0)
		{
			xl_num = x;   // по координате x
			i_mem_left = i;  // по счетчику i
			check1 = 1;
		}
		if (x >= xr && check2 == 0)
		{
			xr_num = x;
			i_mem_right = i;
			check2 = 1;
		}
	}
	printf("%lf %lf", xl_num, xr_num);
	/****************Boundary of numerical rarify wave******************/

	if (numb > 0)
	{
		int helper = counter_all / 10;
		//	printf("h: %i\n", helper);

		for (int j = 0; j < 10; j++)
		{
			i_helper[j] = i_mem_left + j*helper;
			//		printf("%d\n", i_helper[j]);
		}
		printf("\n");
		for (i = 0; i < numcells; i++)
		{
			x = i*dx + 0.5*dx;
			for (int j = 0; j < 10; j++)
			{
				if (i == i_helper[j])
				{
					x_helper[j] = x - (U[i] + sqrt(GAMMA*P[i] / R[i]))*time_max_array[PROBLEM];
					//		printf("%lf\n", x_helper[j]);
				}
			}
		}
	}

	printf("\nxl: %lf, xr: %lf\n", xl, xr);
	printf("points in rarify wave %d %d\n", counter, counter_all);
	printf("points2 %d %d\n", counter2, numcells);
	percents[numb] = float(counter) / float(counter_all) * 100.0f;  //when in int i counter is devided by counter all, its a деление нацело, so the result of 1/4=0;
	printf("percents of the middle: %f\n", percents[numb]);
	percents[numb] = float(counter2) / float(numcells) * 100.0f;  //when in int i counter is devided by counter all, its a деление нацело, so the result of 1/4=0;
	printf("percents of all stream [0:1]: %f\n", percents[numb]);
	//printf("%lf %lf\n", xl, xr);

	/**************** счет интегралов по контуру (верх - низ)********/

	printf("DXDXDXDXD: %lf\n", dx);
#ifdef INTEGRAL
	boundary[numb][0] = integral_x(numcells, dx, difference_RW[2]);
	boundary[numb][1] = integral_x(numcells, dx, difference_RW[3]);
	boundary[numb][2] = integral_x(numcells, dx, difference_RW[4]);
#endif

	/**************** счет интегралов по контуру (верх - низ)********/

	FILE *out4;
	FILE *out5;
	char *FileName2, *FileName3;
	FileName2 = new char[64];
	FileName3 = new char[64];
	sprintf(FileName2, "N%04d_RW_difference.dat", numcells);
	sprintf(FileName3, "N%04d_RW_NUM_ANALITIC.dat", numcells);
	out4 = fopen(FileName2, "w");
	out5 = fopen(FileName3, "w");
	check1 = 0;
	check2 = 0;
	for (i = 0; i < numcells; i++)
	{
#ifndef NC
		x = i*dx + 0.5*dx;
		if (x >= xl && check1 == 0)
		{
			fprintf(out5, "left b: %lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
			check1 = 1;
			continue;
		}
		if (x >= xr && check2 == 0)
		{
			fprintf(out5, "right b: %lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
			check2 = 1;
			continue;
		}
		fprintf(out4, "%lf %lf %lf %lf\n", x, difference_RW[0][i], difference_RW[1][i], difference_RW[2][i]);
		fprintf(out5, "%lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
#else

		// Выведем в новых координатах
		x_NC = (i*dx + 0.5*dx - 0.1) / (time_max_array[0] * (numb + 1));
		fprintf(out4, "%lf %lf %lf %lf\n", x_NC, difference_RW[0][i], difference_RW[1][i], difference_RW[2][i]);
		fprintf(out5, "%lf %lf %lf %lf %lf %lf %lf\n", x_NC, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i]);
#endif
	}
	fclose(out4);
	fclose(out5);

#endif
#endif

	free(R);
	free(P);
	free(U);
	free(RU);
	free(RE);
	free(RS);
	free(S_a);
	free(S_b);

	free(right_parts);



	//system("pause");
#ifdef PRINT	
#ifndef OUTPUT_N_SMOOTH
	fclose(fout);
#endif
#endif

}
#endif