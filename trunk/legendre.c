#include "makespl.h"
#include "piv_ge_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA:	liczbę używanych f. bazowych można ustawić przez wartość
			zmiennej środowiskowej APPROX_BASE_SIZE
			np. export APPROX_BASE_SIZE=3
*/
static double
P(int i, double x)
{
	switch (i)
	{
		case 0:
			return 1;
		case 1:
			return x;
		case 2:
			return ((3*pow(x,2) - 1) / 2);
		case 3:
			return ((5*pow(x,3) - 3*x) / 2);
		case 4:
			return ((35*pow(x,4) - 30*pow(x,2) + 3) / 8);
		case 5:
			return ((63*pow(x,5) - 70*pow(x,3) + 15*x) / 8);
		case 6:
			return ((231*pow(x,6) - 315*pow(x,4) + 105*pow(x,2) - 5) / 16);
		case 7:
			return ((429*pow(x,7) - 693*pow(x,5) + 315*pow(x,3) - 35*x ) / 16);
		case 8:
			return ((6435*pow(x,8) - 12012*pow(x,6) + 6930*pow(x,4) - 1260*pow(x,2) + 35) / 128);
		case 9:
			return ((12155*pow(x,9) - 25740*pow(x,7) + 18018*pow(x,5) - 4620 * pow(x,3) + 315*x) / 128);
		case 10:
			return ((46189*pow(x,10) - 109395*pow(x,8) + 90090*pow(x,6) - 30030*pow(x,4) + 3465*pow(x,2) - 63) / 256);
		case 11:
			return ((88179*pow(x,11) - 230945*pow(x,9) + 218790*pow(x,7) - 90090*pow(x,5) + 15015*pow(x,3) - 693*x) / 256);
		default:
			return 0;
	}
}

/* Pierwsza pochodna P */
static double
dP(int i, double x)
{
	switch (i)
	{
		case 0:
			return 0;
		case 1:
			return 1;
		case 2:
			return (3*x);
		case 3:
			return ((5*pow(x,2) - 1) * 3 / 2);
		case 4:
			return ((7*pow(x,2) - 3) * x * 5 / 2);
		case 5:
			return ((21*pow(x,4) - 14*pow(x,2) + 1) * 15 / 8);
		case 6:
			return ((33*pow(x,4) - 30*pow(x,2) + 5) * x * 21 / 8 );
		case 7: 
			return ((429*pow(x,6) - 495*pow(x,4) + 135*pow(x,2) - 5) * 7 / 16);
		case 8:
			return ((715*pow(x,6) - 1001*pow(x,4) + 385*pow(x,2) - 35) * x * 9 / 16);
		case 9:
			return ((2431*pow(x,8) - 4004*pow(x,6) + 2002*pow(x,4) - 308*pow(x,2) + 7) * 45 / 128);
		case 10:
			return ((4199*pow(x,8) - 7956*pow(x,6) + 4914*pow(x,4) - 1092*pow(x,2) + 63) * x * 55 / 128);
		case 11:
			return ((29393*pow(x,10) - 62985*pow(x,8) + 46410*pow(x,6) - 13650*pow(x,4) + 1365*pow(x,2) - 21) * 33 / 256);
		default:
			return 0;
	}
}

/* Druga pochodna P */
static double
d2P(int i, double x)
{
	switch (i)
	{
		case 0:
			return 0;
		case 1:
			return 0;
		case 2:
			return 3;
		case 3:
			return (15*x);
		case 4:
			return ((7*pow(x,2) - 1) * 15 / 2);
		case 5:
			return ((3*pow(x,2) - 1) * x * 105 / 2);
		case 6:
			return ((33*pow(x,4) - 18*pow(x,2) + 1) * 105 / 8);
		case 7:
			return ((143*pow(x,4) - 110*pow(x,2) + 15) * x * 63 / 8);
		case 8:
			return ((143*pow(x,6) - 143*pow(x,4) + 33*pow(x,2) - 1) * 315 / 16);
		case 9:
			return ((221*pow(x,6) - 273*pow(x,4) + 91*pow(x,2) - 7) * x * 495 / 16);
		case 10:
			return ((4199*pow(x,8) - 6188*pow(x,6) + 2730*pow(x,4) - 364*pow(x,2) + 7) * 495 / 256);
		case 11:
			return ((2261*pow(x,8) - 3876*pow(x,6) + 2142*pow(x,4) - 420*pow(x,2) + 21) * x * 2145 / 128); 
		default:
			return 0;
	}
}

/* Trzecia pochodna P */
static double
d3P(int i, double x)
{
	switch (i)
	{	
		case 0:
			return 0;
		case 1:
			return 0;
		case 2:
			return 0;
		case 3:
			return 15;
		case 4:
			return (105*x);
		case 5:
			return ((9*pow(x,2) - 1) * 105 / 2);
		case 6:
			return ((11*pow(x,2) - 3) * x * 315 / 2);
		case 7:
			return ((143*pow(x,4) - 66*pow(x,2) + 3) * 315 / 8);
		case 8:
			return ((39*pow(x,4) - 26*pow(x,2) + 3) * x * 3465 / 8 );
		case 9:
			return ((221*pow(x,6) - 195*pow(x,4) + 39*pow(x,2) - 1) * 3465 / 16);
		case 10:
			return ((323*pow(x,6) - 357*pow(x,4) + 105*pow(x,2) - 7) * x * 6435 / 16);
		case 11:
			return ((969*pow(x,8) - 1292*pow(x,6) + 510*pow(x,4) - 60*pow(x,2) + 1) * 45045 / 128);
		default:
			return 0;
	}
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, P(i, x[k]) * P(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * P(j, x[k]));
	}

	if (piv_ge_solver(eqs)) {
		if (eqs != NULL) free_matrix(eqs);
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, 1) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] =x[0];
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * P  (k, xx);
				spl->f1[i] += ck * dP (k, xx);
				spl->f2[i] += ck * d2P(k, xx);
				spl->f3[i] += ck * d3P(k, xx);
			}
		}
	}

free_matrix(eqs);
return;
}
