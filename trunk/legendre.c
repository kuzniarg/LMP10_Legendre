#include "makespl.h"
#include <stdio.h>

double f (double x)
{
	return ((3 * x * x - 1)/2);
}
void
make_spl (points_t * pts, spline_t * spl)
{
	alloc_spl (spl, 1);
	int i, n = pts->n;
	double ck, roznica;
	roznica = (pts->x[n - 1] - pts->x[0]) / (n - 1);
	for (i = 0; i < 1; i++)
	{
		ck = (pts->y[n - 1] * f(pts->x[n-1])) / (f(pts->x[0]) * f(pts->x[0]));

		spl->x[i] = pts->x[0];
		spl->f[i] = ck * (3 * spl->x[i] * spl->x[i] - 1)/6;
		spl->f1[i] = ck * 3 * spl->x[i];
		spl->f2[i] = ck * 3;
		spl->f3[i] = 0; 
	}
	return;



}
