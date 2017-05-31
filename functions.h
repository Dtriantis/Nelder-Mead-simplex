#include <math.h>
#include <stdlib.h>
#include <time.h>

/* Our test functions
 * Arguments:
 * 	- num: Function id with integer values
 * 	-   x: Current point, or in the specific case of nelder-mead algorithm, the current simplex
*/ 

typedef struct Simplex
{
	double **vertices;
	double *fvals;
}Simplex;

typedef struct Fstats
{
	double mean;
	double std;
}Fstats;

typedef struct Counters
{
	int shrink_num; 
	int reflection_num;
	int contraction_num;
	int expansion_num;
	int num_of_iters;
	int num_of_evals;
}Counters;


double functions(int func_num, double *point, int len);
double *initial_point(int min,int max, int N);
Simplex *simplex_allocation(int dim);
void simplex_fvals(Simplex *sim, int dim, int func_num);
Simplex *initial_simplex(double *init_p, int dim, int func_num);
void sort_vertices(Simplex *sim, int dim);
Fstats stats(double *fvals, int dim);
double *mean(Simplex *sim, int dim);
double *point_transformation(double *mean_point, double *worst, float t, int dim);
void nelder_mead_step(Simplex *sim, double *mean_point,int dim,int func_num, Counters *counts, int k);
void printtofile(int func_num, int dim, double minimum, double *minimizer, int shrink_num, int reflection_num, int expansion_num, int contraction_num, int num_of_iters,int num_of_evals);
