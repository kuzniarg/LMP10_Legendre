#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

static double P (int n, double x)
{
	if (n==0) return 1;
	else if (n==1) return x;
	else if (n==2) return ( 3*x*x - 1 ) /2;
	else if (n==3) return ( 5*x*x*x - 3*x ) /2;
}

static double wartosc (int i, int j, double x)
{
	if ( i==0 || j==0 )
		return (P (i, x) * P (j, x));
	else
		return (2 * P (i, x) * P (j, x));
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
	double a0, a1, a2, a3, xx;
	
	eqs = make_matrix (nb, nb + 1);


	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
		{
			suma = 0;
			for (k = 0; k < pts->n; k++)
			{
				suma += wartosc (i, j, x[k]);
			}
			add_to_entry_matrix(eqs, j, i, suma);
		}	
		suma = 0;
		for (k = 0; k < n; k++)
		{
			if (j==0)
				suma += y[k];
			else 
				suma += 2 * y[k] * P ( j, x[k] );
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
		a0 = get_entry_matrix(eqs, 0, nb);
		a1 = get_entry_matrix(eqs, 1, nb);
		a2 = get_entry_matrix(eqs, 2, nb);
		a3 = get_entry_matrix(eqs, 3, nb);
		
		xx = spl->x[0] = pts->x[0];
		
		spl->f[0] = 2*a0 + a1*P(1, xx) + a2*P(2, xx) + a3*P(3, xx);
		spl->f1[0] = 3 * ( xx * a2 + a3 * (5 * xx * xx - 1) / 2 );
		spl->f2[0] = 3 * ( a2 + 5 * a3 * xx );
		spl->f3[0] = 15 * a3;
	}

	free (eqs->e);
	free (eqs);
}
