#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <random>
#include <windows.h>

#define metode_fdf 1
//#undef metode_fdf

constexpr int Npart = 300;			// numar de particule
constexpr int Npasi = 1000;			// numar de pasi (de pozitii de echilibru)
constexpr double l0 = 0.6;			// lungimea resortului
constexpr double r_mic = 0.2;		// raza mica
constexpr double R_mare = 0.25;		// raza mare
constexpr double A = 1.0;			// "adancimea" gropii de potential
constexpr double perioada = 5.0;	// perioada potentialului de suprafata
constexpr double k_el = 0.0;		// constanta elastica

double p[Npart], r[Npart];	// sirurile de pozitii si de raze

double v(double loco_p);
double fn1(const gsl_vector q[], void *params);
void dfn1(const gsl_vector q[], void *params, gsl_vector *df);
void fdfn1(const gsl_vector q[], void *params, double *f, gsl_vector *df);

int main(void)
{
	int j;
	DWORD starttime, elapsedtime;

	FILE *fp;
	char numefis[200];
	sprintf(numefis, "E:\\Stoleriu\\C\\special\\3d\\res\\2019\\Elastic\\Minimization\\minim_%d.dat", Npart);
	fp = fopen(numefis, "w");
	fclose(fp);

	// varianta mea favorita de generat numere aleatorii conform C++11 e cu std::random_device
	// (varianta C-clasic, cu functia rand() e invechita si n-ar mai trebui folosita)
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(-0.01, 0.01);
	// la fiecare apel de tipul dist(mt) vom primi un double intre -0.01 si +0.01 (distrib. uniforma)


	p[0] = 0.4 /*+ dist(mt)*/;
	r[0] = r_mic;
	for (j = 1; j < Npart; j++)
	{
		p[j] = p[j - 1] + r[j - 1] + l0 + r[j] /*+ dist(mt)*/;
		r[j] = r_mic;
	}

#ifdef metode_fdf
	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
	gsl_multimin_fdfminimizer* s;
	gsl_vector* x;
	gsl_multimin_function_fdf minex_func;
#else
 	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
 	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
#endif

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc(Npart);
	for (int i = 0; i < Npart; i++)
	{
		gsl_vector_set(x, i, p[i]);
	}

#ifndef metode_fdf 
	/* Set initial step sizes */
	ss = gsl_vector_alloc(Npart);		// pentru metode bazate pe gradient nu trebuie asa ceva
	gsl_vector_set_all(ss, /*0.01*/r_mic/100.0);
#endif

	double *par = NULL;
	/* Initialize method and iterate */
	minex_func.n = Npart;
	minex_func.f = fn1;
#ifdef metode_fdf
	minex_func.df = dfn1;
	minex_func.fdf = fdfn1;
#endif
	minex_func.params = par;

#ifdef metode_fdf
	s = gsl_multimin_fdfminimizer_alloc(T, Npart);
	gsl_multimin_fdfminimizer_set(s, &minex_func, x, 0.01, 1.0e-4);
#else
 	s = gsl_multimin_fminimizer_alloc(T, Npart);
 	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
#endif


	starttime = timeGetTime();
	//echilibru initial
	do
	{
		iter++;
#ifdef metode_fdf
		status = gsl_multimin_fdfminimizer_iterate(s);
#else
		status = gsl_multimin_fminimizer_iterate(s);
#endif		

		if (status)
			break;

#ifdef metode_fdf
		status = gsl_multimin_test_gradient(s->gradient, 1.0e-3);
		printf("%5zd %10.3e %10.3e f() = %7.3f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), s->f);
#else
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-3);
		if (!(iter % 10000))
		{
			printf("%5zd %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), s->fval, size);
		}
#endif


	} while (status == GSL_CONTINUE /*&& iter < 100*/);		//eliminare supapa de siguranta iter<100
	
	elapsedtime = timeGetTime() - starttime;
	printf("\n DONE IN %ld milliseconds\n", elapsedtime);


	fp = fopen(numefis, "w");
	for (int i = 0; i < Npart; i++)
	{
		fprintf(fp, "%d %20.16lf\n", i, gsl_vector_get(s->x, i));
	}
	fclose(fp);

	gsl_vector_free(x);
#ifdef metode_fdf
	gsl_multimin_fdfminimizer_free(s);
#else
 	gsl_vector_free(ss);
 	gsl_multimin_fminimizer_free(s);
#endif

	return status;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

double v(double loco_p)
{
	return -A * sin(perioada * loco_p * M_PI);
}

double fn1(const gsl_vector q[], void *params)
{
	(void)(params); /* avoid unused parameter warning */
	double sv = 0, sx = 0, f;
	double xi, xip1;
	
	for (int i = 0; i < (Npart - 1); i++)
	{
		xi = gsl_vector_get(q, i);
		xip1 = gsl_vector_get(q, (i + 1));
 		sv = sv + v(xi);
  		sx = sx + (xip1 - xi - r[i] - r[i+1] - l0) * (xip1 - xi - r[i] - r[i+1] - l0);
	}

	sv += v(gsl_vector_get(q, (Npart-1)));

	f = sv + (k_el / 2.0) * sx;
	return f;
}

void dfn1(const gsl_vector q[], void *params, gsl_vector *df)
{
	(void)(params); /* avoid unused parameter warning */
	double xim1, xi, xip1, deriv_xi;

	xi   = gsl_vector_get(q, 0);
	xip1 = gsl_vector_get(q, 1);
	deriv_xi = -perioada * A * M_PI * cos(perioada * xi * M_PI) - k_el * (xip1 - xi - r[0] - r[1] - l0);
	gsl_vector_set(df, 0, deriv_xi);

	for (int i = 1; i < (Npart - 1); i++)
	{
		xim1 = xi;
		xi   = xip1;
		xip1 = gsl_vector_get(q, (i + 1));
		deriv_xi = -perioada * A * M_PI * cos(perioada * xi * M_PI) - k_el * (r[i-1] - r[i+1] + xim1 - 2.0 * xi + xip1);
		gsl_vector_set(df, i, deriv_xi);
	}

	xim1 = xi;
	xi = xip1;
	deriv_xi = -perioada * A * M_PI * cos(perioada * xi * M_PI) + k_el * (xi - xim1 - r[Npart-2] - r[Npart-1] - l0);
	gsl_vector_set(df, (Npart-1), deriv_xi);
}

void fdfn1(const gsl_vector q[], void *params, double *f, gsl_vector *df)
{
	*f = fn1(q, params);
	dfn1(q, params, df);
}