#include <stdio.h>
#include "functions.h"
#define KMAX 1000000
#define eps 0.000000001
void main(int argc, char *argv[])
{
	int i, j , dim, func_num, min_res_pos=0, num_of_restarts,k,l;
	float min,max;
	double *mean_point=NULL,min_res=0,norm_sum=0,max_f=0,max_nrm=0;
	Fstats fs = {0, 0};
	FILE *fp;
	FILE *fp1;
	fp = fopen("/home/means.txt ","w");
	fp1 = fopen("/home/vertices.txt ","w");
	srand(time(NULL));

	if(argc < 6)
	{
		printf("Missing arguments...\nExiting!\n");
		exit(0);
	}
	else if(argc > 6)
	{
		printf("Too much arguments...\nExiting!\n");
		exit(0);
	}
	dim = atoi(argv[1]);
	func_num = atoi(argv[2]);
	min = atoi(argv[3]);
	max = atoi(argv[4]);
	num_of_restarts = atoi(argv[5]);
	
	Counters *counts = malloc(num_of_restarts*sizeof(Counters));
	for(i=0;i<num_of_restarts;i++)
	{
		counts[i].num_of_iters=0;
		counts[i].num_of_evals=0;
		counts[i].shrink_num=0;
		counts[i].reflection_num=0;
		counts[i].contraction_num=0;
		counts[i].expansion_num=0;
	}
	
	double *fvals_restarts = malloc((num_of_restarts)*sizeof(double));
	double **minimizer_restarts = malloc((num_of_restarts)*sizeof(double*));
	for (i=0;i<num_of_restarts;i++){
		minimizer_restarts[i] = malloc(dim*sizeof(double));
	}
	
	double *init=NULL;
	Simplex *sim=NULL;

	//Restarts for keeping the best and avoid stagnation
	for (i=0;i<num_of_restarts;i++){

		init = initial_point(min,max,dim);
		sim = initial_simplex(init, dim, func_num);
		counts[i].num_of_evals+=dim+1;
		sort_vertices(sim, dim);
		mean_point = mean(sim, dim);
		fs = stats(sim->fvals,dim);
		
		fprintf(fp,"Number of Restart %d\n", i);
		fprintf(fp,"\n");


	//Nelder-Mead iterations
		
		//First termination criterion
// 		while(fs.std > eps && counts[i].num_of_iters <= KMAX)
// 		{
// 			counts[i].num_of_iters++;
// 			nelder_mead_step(sim, mean_point, dim, func_num, counts,i);
// 			simplex_fvals(sim, dim, func_num);
// 			counts[i].num_of_evals+=dim+1;
// 			sort_vertices(sim, dim);
// 			fs = stats(sim->fvals,dim);
// 			
// 			fprintf(fp, "%f", fs.mean);
// 			fprintf(fp,"\n");
// 
// 			free(mean_point);
// 			mean_point = mean(sim, dim);
// 		
// 		}
		
		//second termination Critereon
		
		for(j=0;j<dim;j++)
			norm_sum += pow(sim->vertices[dim][j]-sim->vertices[0][j],2);
			
		max_f = sim->fvals[dim] - sim->fvals[0];
		max_nrm = sqrt(norm_sum);
		
		while((max_f > eps && max_nrm > eps) && (counts[i].num_of_iters <= KMAX))
		{
			
			counts[i].num_of_iters++;
			nelder_mead_step(sim, mean_point, dim, func_num, counts,i);
			simplex_fvals(sim, dim, func_num);
			counts[i].num_of_evals+=dim+1;
			sort_vertices(sim, dim);
			
			for(k=0;k<dim+1;k++)
			{
				for(l=0;l<dim;l++){
					
					fprintf(fp1, "%f ", sim->vertices[k][l]);
				}
			}
			fprintf(fp1, "\n");
			fs = stats(sim->fvals,dim);
			
			norm_sum = 0;
			for(j=0;j<dim;j++)
				norm_sum += pow(sim->vertices[dim][j]-sim->vertices[0][j],2);
			
			max_f = sim->fvals[dim] - sim->fvals[0];
			max_nrm = sqrt(norm_sum);
			
			fprintf(fp, "%f", fs.mean);
			fprintf(fp,"\n");

			free(mean_point);
			mean_point = mean(sim, dim);
			
		}
		
		fvals_restarts[i] = sim->fvals[0];
		for(j=0;j<dim;j++){
			minimizer_restarts[i][j] = sim->vertices[0][j];
		}
	
		free(init);	
		free(sim);
		init=NULL;
		sim=NULL;
	}
	
	
	//Find minimum
	min_res = fvals_restarts[0];
	for(j=1;j<num_of_restarts;j++)
	{
		if(fvals_restarts[j] < min_res)
		{
			min_res = fvals_restarts[j];
			min_res_pos = j;
		}	
		
		
	}
	fprintf(fp,"Best Min was found in %d\n", min_res_pos);

	printtofile(func_num, dim, min_res, minimizer_restarts[min_res_pos], counts[min_res_pos].shrink_num, counts[min_res_pos].reflection_num, counts[min_res_pos].expansion_num, counts[min_res_pos].contraction_num, counts[min_res_pos].num_of_iters, counts[min_res_pos].num_of_evals);	
	fclose(fp);	
	fclose(fp1);
}

