#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <random>
constexpr int NMAX = 5;
constexpr double l = 0.6;


double A = 1.0, k = 0;
double p[NMAX], r[NMAX];	// sirurile de pozitii si de raze

double v(int k, double loco_p[]);
double fn1(const gsl_vector q[], void* params);

int main(void)
{
	int j;

	FILE *fp;
	char numefis[200];
	sprintf(numefis, "E:\\Stoleriu\\C\\special\\3d\\res\\2019\\Elastic\\Minimization\\minim_%d.dat", NMAX);
	fp = fopen(numefis, "w");

	// varianta mea favorita de generat numere aleatorii conform C++11 e cu std::random_device
	// (varianta C-clasic, cu functia rand() e invechita si n-ar mai trebui folosita)
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(-0.05, 0.05);
	// la fiecare apel de tipul dist(mt) vom primi un double intre -0.05 si +0.05 (distrib. uniforma)


	p[0] = 0.4 + dist(mt);
	r[0] = 0.2;
	printf("0 @ %lf\n", p[0]);
	fprintf(fp, "%lf %lf ", p[0], v(0, p));
	for (j = 1; j < NMAX; j++)
	{
		r[j] = 0.2;
		p[j] = p[j - 1] + r[j - 1] + l + r[j] + dist(mt);
		printf("%d @ %lf\n", j, p[j]);
		fprintf(fp, "%lf %lf ", p[j], v(j, p));
	}
	fprintf(fp, "\n");

	const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2rand;
	gsl_multimin_fminimizer* s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc(NMAX);
	for (int i = 0; i < NMAX; i++)
	{
		gsl_vector_set(x, i, p[i]);
	}

	/* Set initial step sizes */
	ss = gsl_vector_alloc(NMAX);
	gsl_vector_set_all(ss, 0.01);


	double* par = NULL;
	/* Initialize method and iterate */
	minex_func.n = NMAX;
	minex_func.f = fn1;
	minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc(T, NMAX);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-3);

		if (status == GSL_SUCCESS)
		{
			printf("converged to minimum at\n");
		}

		printf("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), s->fval, size);

		for (j = 0; j < NMAX; j++)
		{
			fprintf(fp, "%lf %lf ", gsl_vector_get(s->x, j), -A * sin(5 * gsl_vector_get(s->x, j) * M_PI));
		}
		fprintf(fp, "\n");

	} while (status == GSL_CONTINUE && iter < 100);

	for (int i = 0; i < NMAX; i++)
	{
		printf("%d -> %6.4lf\n", i, gsl_vector_get(s->x, i));
	}

	fclose(fp);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	return status;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

double v(int k, double loco_p[])
{
	return -A * sin(5 * loco_p[k] * M_PI);
}

double fn1(const gsl_vector q[], void* params)
{
	(void)(params); /* avoid unused parameter warning */
	double sv = 0, sx = 0, f;

	double nou_p[NMAX];
	for (int i = 0; i < NMAX; i++)
	{
		nou_p[i] = gsl_vector_get(q, i);
	}

	for (int i = 0; i < NMAX - 1; i++)
	{
		sv = sv + v(i, nou_p);
		sx = sx + (nou_p[i] - nou_p[i - 1] - r[i] - r[i - 1] - l) * (nou_p[i] - nou_p[i - 1] - r[i] - r[i - 1] - l);
	}
	sv += v(NMAX - 1, nou_p);

	f = sv + (k / 2.0) * sx;
	return f;
}

