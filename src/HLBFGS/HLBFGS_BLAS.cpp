#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "HLBFGS_BLAS.h"
#include <cmath>

double HLBFGS_DDOT(const size_t n, const double *x, const double *y)
{
	double result = 0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:result)
#endif
	for (ptrdiff_t i = 0; i < (ptrdiff_t)n; i++)
	{
		result += x[i] * y[i];
	}
	return result;
}

void HLBFGS_DAXPY(const size_t n, const double alpha, const double *x,
		double *y)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (ptrdiff_t i = 0; i < (ptrdiff_t)n; i++)
	{
		y[i] += alpha * x[i];
	}
}

double HLBFGS_DNRM2(const size_t n, const double *x)
{
	double result = 0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:result)
#endif
	for (ptrdiff_t i = 0; i < (ptrdiff_t)n; i++)
	{
		result += x[i] * x[i];
	}
	return std::sqrt(result);
}

void HLBFGS_DSCAL(const size_t n, const double a, double *x)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (ptrdiff_t i = 0; i < (ptrdiff_t)n; i++)
	{
		x[i] *= a;
	}
}

