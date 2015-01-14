#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

double P0 ()
{
	return 1;
}

double P1 (double x)
{
	return x;
}

double P2 (double x)
{
	return ( 3*x*x - 1 ) /2;
}

double P3 (double x)
{
	return ( 5*x*x*x - 3*x ) /2;
}

void
make_spl (points_t * pts, spline_t * spl)
{
	matrix_t *eqs = NULL;
	double *x = pts->x;
	double *y = pts->y;
	int n = pts->n;
	int i, j, k;
	int nb = 4; 
	double suma = 0;

	eqs = make_matrix (nb, nb + 1);


	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
		{
			suma = 0;
			for (k = 0; k < pts->n; k++)
			{
				if (i==0)
				{	
					if (j==0)
						suma += P0();
					else if (j==1)
						suma += P1(x[k]);
					else if (j==2)
						suma += P2(x[k]);
					else if (j==3)
						suma += P3(x[k]);
				}
				else if (i==1)
				{
					if (j==0)
						suma += P1(x[k]);
					if (j==1)
						suma += 2 * P1(x[k]) * P1(x[k]);
					if (j==2)
						suma += 2 * P1(x[k]) * P2(x[k]);
					if (j==3)
						suma += 2 * P1(x[k]) * P3(x[k]);
				}
				else if (i==2)
				{
					if (j==0)
						suma += P2(x[k]);
					if (j==1)
						suma += 2 * P2(x[k]) * P1(x[k]);
					if (j==2)
						suma += 2 * P2(x[k]) * P2(x[k]);
					if (j==3)
						suma += 2 * P2(x[k]) * P3(x[k]);
				}
				else if (i==3)
				{
					if (j==0)
						suma += P3(x[k]);
					if (j==1)
						suma += 2 * P3(x[k]) * P1(x[k]);
					if (j==2)
						suma += 2 * P3(x[k]) * P2(x[k]);
					if (j==3)
						suma += 2 * P3(x[k]) * P3(x[k]);
				}
				
			}
			add_to_entry_matrix(eqs, j, i, suma);
		}	
		suma = 0;
		for (k = 0; k < n; k++)
		{
			if (j==0)
				suma += y[k];
			if (j==1)
				suma += 2 * y[k] * P1(x[k]);
			if (j==2)
				suma += 2 * y[k] * P2(x[k]);
			if (j==3)
				suma += 2 * y[k] * P3(x[k]);
		}
		add_to_entry_matrix(eqs, j, nb, suma);
	}


	if (piv_ge_solver (eqs))
	{
		spl->n = 0;
		return;
	}

	if (alloc_spl (spl, 1) == 0)
	{
		double a0 = get_entry_matrix(eqs, 0, nb);
		double a1 = get_entry_matrix(eqs, 1, nb);
		double a2 = get_entry_matrix(eqs, 2, nb);
		double a3 = get_entry_matrix(eqs, 3, nb);
		double xx = spl->x[0] = pts->x[0];
		spl->f[0] = 2*a0 + a1*P1(xx) + a2*P2(xx) + a3*P3(xx);
		spl->f1[0] = 3 * ( xx * a2 + a3 * (5 * xx * xx - 1) / 2 );
		spl->f2[0] = 3 * ( a2 + 5 * a3 * xx );
		spl->f3[0] = 15 * a3;
	}
}
