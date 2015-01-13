#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE*/

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
  double a = x[0];
  double b = x[pts->n - 1];
  int i, j, k;
  int nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv = getenv ("APPROX_BASE_SIZE");

  if (nbEnv != NULL && atoi (nbEnv) > 0)
    nb = atoi (nbEnv);

  eqs = make_matrix (nb, nb + 1);


  for (j = 0; j < nb; j++)
    {
      for (i = 0; i < nb; i++)
	for (k = 0; k < pts->n; k++)
	  add_to_entry_matrix (eqs, j, i,
			       fi (a, b, nb, i, x[k]) * fi (a, b, nb, j,
							    x[k]));

      for (k = 0; k < pts->n; k++)
	add_to_entry_matrix (eqs, j, nb, y[k] * fi (a, b, nb, j, x[k]));
    }


  if (piv_ge_solver (eqs))
    {
      spl->n = 0;
      return;
    }

  if (alloc_spl (spl, nb) == 0)
    {
      for (i = 0; i < spl->n; i++)
	{
	  double xx = spl->x[i] = a + i * (b - a) / (spl->n - 1);
	  xx += 10.0 * DBL_EPSILON;	// zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
	  spl->f[i] = 0;
	  spl->f1[i] = 0;
	  spl->f2[i] = 0;
	  spl->f3[i] = 0;
	  for (k = 0; k < nb; k++)
	    {
	      double ck = get_entry_matrix (eqs, k, nb);
	      spl->f[i] += ck * fi (a, b, nb, k, xx);
	      spl->f1[i] += ck * dfi (a, b, nb, k, xx);
	      spl->f2[i] += ck * d2fi (a, b, nb, k, xx);
	      spl->f3[i] += ck * d3fi (a, b, nb, k, xx);
	    }
	}
    }


}
