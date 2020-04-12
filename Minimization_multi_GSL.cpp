#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <random>
constexpr int Npart = 500;			// numar de particule
constexpr int Npasi = 1000;			// numar de pasi (de pozitii de echilibru)
constexpr double l0 = 0.6;			// lungimea resortului
constexpr double r_mic = 0.2;		// raza mica
constexpr double R_mare = 0.25;		// raza mare
constexpr double A = 1.0;			// "adancimea" gropii de potential
constexpr double k_el = 0.0;		// constanta elastica


double p[Npart], r[Npart];	// sirurile de pozitii si de raze

//double v(int k, double loco_p[]);
double v(double loco_p);
double fn1(const gsl_vector q[], void* params);

int main(void)
{
	int j;

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

	const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer* s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc(Npart);
	for (int i = 0; i < Npart; i++)
	{
		gsl_vector_set(x, i, p[i]);
	}

	/* Set initial step sizes */
	ss = gsl_vector_alloc(Npart);
	gsl_vector_set_all(ss, /*0.01*/r_mic/100.0);


	double *par = NULL;
	/* Initialize method and iterate */
	minex_func.n = Npart;
	minex_func.f = fn1;
	minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc(T, Npart);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	//echilibru initial
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-3);

		if (!(iter % 10000))
		{
			printf("%5zd %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), s->fval, size);
		}

	} while (status == GSL_CONTINUE /*&& iter < 100*/);		//eliminare supapa de siguranta iter<100
	


	fp = fopen(numefis, "w");
	for (int i = 0; i < Npart; i++)
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

//double v(int k, double loco_p[])
double v(double loco_p)
{
	//return -A * sin(5 * loco_p[k] * M_PI);
	return -A * sin(5 * loco_p * M_PI);
}

double fn1(const gsl_vector q[], void* params)
{
	(void)(params); /* avoid unused parameter warning */
	double sv = 0, sx = 0, f;

// 	double nou_p[Npart];				//mai rapid daca nu facem o copie locala a gsl_vector q[] la fiecare apel al functiei?
// 	for (int i = 0; i < Npart; i++)
// 	{
// 		nou_p[i] = gsl_vector_get(q, i);
// 	}

	for (int i = 0; i < Npart - 1; i++)
	{
// 		sv = sv + v(i, nou_p);
// 		sx = sx + (nou_p[i] - nou_p[i - 1] - r[i] - r[i - 1] - l0) * (nou_p[i] - nou_p[i - 1] - r[i] - r[i - 1] - l0);
		sv = sv + v(gsl_vector_get(q, i));
 		sx = sx + (gsl_vector_get(q, i) - gsl_vector_get(q, (i+1)) - r[i] - r[i+1] - l0) * 
			      (gsl_vector_get(q, i) - gsl_vector_get(q, (i+1)) - r[i] - r[i+1] - l0);
	}
//	sv += v(Npart - 1, nou_p);
	sv += v(gsl_vector_get(q, (Npart-1)));

	f = sv + (k_el / 2.0) * sx;
	return f;
}

// void dfn1(const gsl_vector q[], void* params, gsl_vector* df)
// {
// 	(void)(params); /* avoid unused parameter warning */
// 
// 	for (int i = 0; i < Npart - 1; i++)
// 	{
// 
// 		sv = sv + v(gsl_vector_get(q, i));
// 		sx = sx + (gsl_vector_get(q, i) - gsl_vector_get(q, (i - 1)) - r[i] - r[i - 1] - l0) *
// 			(gsl_vector_get(q, i) - gsl_vector_get(q, (i - 1)) - r[i] - r[i - 1] - l0);
// 	}
// 	//	sv += v(Npart - 1, nou_p);
// 	sv += v(gsl_vector_get(q, (Npart - 1)));
// 
// 
// }
