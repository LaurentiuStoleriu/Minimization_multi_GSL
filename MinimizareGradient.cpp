#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#define NMAX 50

FILE* fout;
const char output[] = "iesire.txt";

double energie(const gsl_vector* v, void* params);
void gradient_energie(const gsl_vector* v, void* params, gsl_vector* df);
void init_functii(const gsl_vector* x, void* params, double* f, gsl_vector* df);

double A = 1.0, k = 0.0, w = 5.0 * M_PI, L = 1.0;
double pozInit[NMAX], razeInit[NMAX];

int main()
{
	int status, i;

	const gsl_multimin_fdfminimizer_type* T;
	gsl_multimin_fdfminimizer* s;

	double* par = NULL;

	gsl_vector* x;
	gsl_multimin_function_fdf functie;

	functie.n = NMAX;
	functie.f = energie;
	functie.df = gradient_energie;
	functie.fdf = init_functii;
	functie.params = par;

	fout = fopen(output, "wt");
	if (!fout)
	{
		printf("EROARE! Nu pot deschide fisierul %s\nTermin programul.", output);
		return 1;
	}

	//Initializari
	fprintf(fout, "POZITII INITIALE\n");
	pozInit[0] = 0.4;
	razeInit[0] = 0.2;
	fprintf(fout, "0 %lf\n", pozInit[0]);
	for (i = 1; i < NMAX; i++)
	{
		razeInit[i] = 0.2;
		pozInit[i] = pozInit[i - 1] + L + 2 * razeInit[i];
		fprintf(fout, "%d %lf\n", i, pozInit[i]);
	}

	x = gsl_vector_alloc(NMAX);
	for (i = 0; i < NMAX; i++) gsl_vector_set(x, i, pozInit[i]);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, NMAX);

	gsl_multimin_fdfminimizer_set(s, &functie, x, 0.01, 1e-4);

	do
	{
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status) break;

		status = gsl_multimin_test_gradient(s->gradient, 1e-3);
	} while (status == GSL_CONTINUE);

	fprintf(fout, "\nPOZITII FINALE\n");
	for (i = 0; i < NMAX; i++) fprintf(fout, "%d %lf\n", i, gsl_vector_get(s->x, i));

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	return 0;
}


//Functie
double energie(const gsl_vector* v, void* params)
{
	(void)(params);

	double energieGrPot = 0, energieElast = 0;
	double poz[NMAX];
	int i;

	for (i = 0; i < NMAX; i++) poz[i] = gsl_vector_get(v, i);

	for (i = 0; i < NMAX - 1; i++)
	{
		energieGrPot += -A * sin(w * poz[i]);
		energieElast += (k / 2.0) * (poz[i + 1] - poz[i] - razeInit[i + 1] - razeInit[i] - L) * (poz[i + 1] - poz[i] - razeInit[i + 1] - razeInit[i] - L);
	}

	return energieGrPot + energieElast;
}


//Gradientul functiei
void gradient_energie(const gsl_vector* v, void* params, gsl_vector* df)
{
	(void)(params);

	double poz[NMAX];
	int i;

	for (i = 0; i < NMAX; i++) poz[i] = gsl_vector_get(v, i);

	gsl_vector_set(df, 0, -w * A * cos(w * poz[0]) - k * (poz[1] - poz[0] - razeInit[1] - razeInit[0] - L));
	for (i = 1; i < NMAX - 1; i++) gsl_vector_set(df, i, -w * A * cos(w * poz[i]) + k * (poz[i + 1] - poz[i] - razeInit[i + 1] - razeInit[i] - L));
	gsl_vector_set(df, i, -w * A * cos(w * poz[i]) + k * (poz[i] - poz[i - 1] - razeInit[i] - razeInit[i - 1] - L));
}


//Construiesc gradientul si functia
void init_functii(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	*f = energie(x, params);
	gradient_energie(x, params, df);
}