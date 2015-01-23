#include "makespl.h"
#include "piv_ge_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
zmiennej środowiskowej APPROX_BASE_SIZE
*/
/*
* Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
* - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
*/
double
P (int i, double x)
{
  switch (i)
    {
    case 0:
      return 1;
    case 1:
      return x;
    case 2:
      return ((3 * pow (x, 2) - 1) / 2);
    case 3:
      return ((5 * pow (x, 3) - 3 * x) / 2);
    case 4:
      return ((35 * pow (x, 4) - 30 * pow (x, 2) + 3) / 8);
    case 5:
      return ((63 * pow (x, 5) - 70 * pow (x, 3) + 15 * x) / 8);
    case 6:
      return ((231 * pow (x, 6) - 315 * pow (x, 4) + 105 * pow (x, 2) -
	       5) / 16);
    case 7:
      return ((429 * pow (x, 7) - 693 * pow (x, 5) + 315 * pow (x, 3) -
	       35 * x) / 16);
    case 8:
      return ((6435 * pow (x, 8) - 12012 * pow (x, 6) + 6930 * pow (x, 4) -
	       1260 * pow (x, 2) + 35) / 128);
    case 9:
      return ((12155 * pow (x, 9) - 25740 * pow (x, 7) + 18018 * pow (x, 5) -
	       4620 * pow (x, 3) + 315 * x) / 128);
    case 10:
      return ((46189 * pow (x, 10) - 109395 * pow (x, 8) +
	       90090 * pow (x, 6) - 30030 * pow (x, 4) + 3465 * pow (x,
								     2) -
	       63) / 256);
    case 11:
      return ((88179 * pow (x, 11) - 230945 * pow (x, 9) +
	       218790 * pow (x, 7) - 90090 * pow (x, 5) + 15015 * pow (x,
								       3) -
	       693 * x) / 256);
    default:
      return 0;
    }
}

/* Pierwsza pochodna fi */
double
dP (int i, double x)
{
  switch (i)
    {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return (3 * x);
    case 3:
      return ((5 * pow (x, 2) - 1) * 3 / 2);
    case 4:
      return ((7 * pow (x, 2) - 3) * x * 5 / 2);
    case 5:
      return ((21 * pow (x, 4) - 14 * pow (x, 2) + 1) * 15 / 8);
    case 6:
      return ((33 * pow (x, 4) - 30 * pow (x, 2) + 5) * x * 21 / 8);
    case 7:
      return ((429 * pow (x, 6) - 495 * pow (x, 4) + 135 * pow (x, 2) -
	       5) * 7 / 16);
    case 8:
      return ((715 * pow (x, 6) - 1001 * pow (x, 4) + 385 * pow (x, 2) -
	       35) * x * 9 / 16);
    case 9:
      return ((2431 * pow (x, 8) - 4004 * pow (x, 6) + 2002 * pow (x, 4) -
	       308 * pow (x, 2) + 7) * 45 / 128);
    case 10:
      return ((4199 * pow (x, 8) - 7956 * pow (x, 6) + 4914 * pow (x, 4) -
	       1092 * pow (x, 2) + 63) * x * 55 / 128);
    case 11:
      return ((29393 * pow (x, 10) - 62985 * pow (x, 8) + 46410 * pow (x, 6) -
	       13650 * pow (x, 4) + 1365 * pow (x, 2) - 21) * 33 / 256);
    default:
      return 0;
    }
}

/* Druga pochodna fi */
double
d2P (int i, double x)
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
      return (15 * x);
    case 4:
      return ((7 * pow (x, 2) - 1) * 15 / 2);
    case 5:
      return ((3 * pow (x, 2) - 1) * x * 105 / 2);
    case 6:
      return ((33 * pow (x, 4) - 18 * pow (x, 2) + 1) * 105 / 8);
    case 7:
      return ((143 * pow (x, 4) - 110 * pow (x, 2) + 15) * x * 63 / 8);
    case 8:
      return ((143 * pow (x, 6) - 143 * pow (x, 4) + 33 * pow (x, 2) -
	       1) * 315 / 16);
    case 9:
      return ((221 * pow (x, 6) - 273 * pow (x, 4) + 91 * pow (x, 2) -
	       7) * x * 495 / 16);
    case 10:
      return ((4199 * pow (x, 8) - 6188 * pow (x, 6) + 2730 * pow (x, 4) -
	       364 * pow (x, 2) + 7) * 495 / 256);
    case 11:
      return ((2261 * pow (x, 8) - 3876 * pow (x, 6) + 2142 * pow (x, 4) -
	       420 * pow (x, 2) + 21) * x * 2145 / 128);
    default:
      return 0;
    }
}

/* Trzecia pochodna fi */
double
d3P (int i, double x)
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
      return (105 * x);
    case 5:
      return ((9 * pow (x, 2) - 1) * 105 / 2);
    case 6:
      return ((11 * pow (x, 2) - 3) * x * 315 / 2);
    case 7:
      return ((143 * pow (x, 4) - 66 * pow (x, 2) + 3) * 315 / 8);
    case 8:
      return ((39 * pow (x, 4) - 26 * pow (x, 2) + 3) * x * 3465 / 8);
    case 9:
      return ((221 * pow (x, 6) - 195 * pow (x, 4) + 39 * pow (x, 2) -
	       1) * 3465 / 16);
    case 10:
      return ((323 * pow (x, 6) - 357 * pow (x, 4) + 105 * pow (x, 2) -
	       7) * x * 6435 / 16);
    case 11:
      return ((969 * pow (x, 8) - 1292 * pow (x, 6) + 510 * pow (x, 4) -
	       60 * pow (x, 2) + 1) * 45045 / 128);
    default:
      return 0;
    }
}

double
fi (double a, double b, int n, int i, double x)
{
  return P (i, x);
}

/* Pierwsza pochodna fi */
double
dfi (double a, double b, int n, int i, double x)
{
  return dP (i, x);
}

/* Druga pochodna fi */
double
d2fi (double a, double b, int n, int i, double x)
{
  return d2P (i, x);
}

/* Trzecia pochodna fi */
double
d3fi (double a, double b, int n, int i, double x)
{
  return d3P (i, x);
}

/* Pomocnicza f. do rysowania bazy */
double
xfi (double a, double b, int n, int i, FILE * out)
{
  double h = (b - a) / (n - 1);
  double h3 = h * h * h;
  int hi[5] = { i - 2, i - 1, i, i + 1, i + 2 };
  double hx[5];
  int j;
  for (j = 0; j < 5; j++)
    hx[j] = a + h * hi[j];
  fprintf (out, "# nb=%d, i=%d: hi=[", n, i);
  for (j = 0; j < 5; j++)
    fprintf (out, " %d", hi[j]);
  fprintf (out, "] hx=[");
  for (j = 0; j < 5; j++)
    fprintf (out, " %g", hx[j]);
  fprintf (out, "]\n");
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
#ifdef DEBUG
#define TESTBASE 500
  {
    FILE *tst = fopen ("debug_base_plot.txt", "w");
    double dx = (b - a) / (TESTBASE - 1);
    for (j = 0; j < nb; j++)
      xfi (a, b, nb, j, tst);
    for (i = 0; i < TESTBASE; i++)
      {
	fprintf (tst, "%g", a + i * dx);
	for (j = 0; j < nb; j++)
	  {
	    fprintf (tst, " %g", fi (a, b, nb, j, a + i * dx));
	    fprintf (tst, " %g", dfi (a, b, nb, j, a + i * dx));
	    fprintf (tst, " %g", d2fi (a, b, nb, j, a + i * dx));
	    fprintf (tst, " %g", d3fi (a, b, nb, j, a + i * dx));
	  }
	fprintf (tst, "\n");
      }
    fclose (tst);
  }
#endif
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
#ifdef DEBUG
  write_matrix (eqs, stdout);
#endif
  if (piv_ge_solver (eqs))
    {
      spl->n = 0;
      return;
    }
#ifdef DEBUG
  write_matrix (eqs, stdout);
#endif
  if (alloc_spl (spl, 1) == 0)
    {
      for (i = 0; i < spl->n; i++)
	{
	  double xx = spl->x[i] = (a + b) / 2;
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
#ifdef DEBUG
  {
    FILE *tst = fopen ("debug_spline_plot.txt", "w");
    double dx = (b - a) / (TESTBASE - 1);
    for (i = 0; i < TESTBASE; i++)
      {
	double yi = 0;
	double dyi = 0;
	double d2yi = 0;
	double d3yi = 0;
	double xi = a + i * dx;
	for (k = 0; k < nb; k++)
	  {
	    yi += get_entry_matrix (eqs, k, nb) * fi (a, b, nb, k, xi);
	    dyi += get_entry_matrix (eqs, k, nb) * dfi (a, b, nb, k, xi);
	    d2yi += get_entry_matrix (eqs, k, nb) * d2fi (a, b, nb, k, xi);
	    d3yi += get_entry_matrix (eqs, k, nb) * d3fi (a, b, nb, k, xi);
	  }
	fprintf (tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi);
      }
    fclose (tst);
  }
#endif
  free_matrix (eqs);
  return;
}
