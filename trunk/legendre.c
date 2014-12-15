#include "makespl.h"
#include <stdio.h>

void
make_spl (points_t * pts, spline_t * spl)
{
	alloc_spl (spl, pts->n);
	int i, n = pts->n;
	double roznica = (x[n-1]-x[0])\n;
	spl->n = pts->n;
	for (i = 0; i <= spl->n; i++)
	{
		spl->x[i] = pts->x[0] + i * roznica;
		spl->f[i] = pts->y[i];
		//spl->f1[i] = (3 * spl->x[i] * spl->x[i] -1) / 2;
		//spl->f2[i] = (5 * spl->x[i] * spl->x[i] * spl->x[i] - 3 * spl->x[i] ) / 2;
		//spl->f3[i] = (63 * spl->x[i] * spl->x[i] * spl->x[i] * spl->x[i]  - 30 * spl->x[i] * spl->x[i] + 3) / 8;
		spl->f1[i] = 2ax+b;
		spl->f2[i] = 2a;
		spl->f3[i] = 0;
////////////////ZNAJDŹ A I B GAUSSEM!
	}
	return;



}
