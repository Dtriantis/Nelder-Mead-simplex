#include "functions.h"
#include <stdio.h>
/* Our test functions
 * Arguments:
 * 	- num: Function id with integer values
 * 	-   x: Current point, or in the specific case of nelder-mead algorithm, the current simplex
*/ 


double functions(int num, double *vec, int dim)
{
	int i;
	double sum=0, theta, F;
	// Simple testing function: Σxi^2
	if (num==1)
	{
		for(i=0;i<dim;i++)
			sum+=(vec[i])*(vec[i]);
		return sum;
	}
	//Freudenstein and Roth
	if (num==2)
	{
		for(i=0;i<dim;i+=2){
			sum += pow((-13 + vec[i] + ((5 - vec[i+1])*vec[i+1] - 2)*vec[i+1]),2)
				+ pow((-29 +  vec[i] + ((vec[i+1] + 1)*vec[i+1] -14)*vec[i+1]),2);
		}
		return sum;
	}
	// Powell
	if (num==3)
	{
		for(i=0;i<dim;i+=2){
			sum += pow((10000*vec[i]*vec[i+1] - 1),2)
				+ pow((exp(-vec[i]) + (exp(-vec[i+1])) - 1.0001),2);
		}
		return sum;
	}
	// Beale
	if (num==4)
	{
		for(i=0;i<dim;i+=2){
			sum += pow((1.5 - vec[i]*(1 - vec[i+1])),2)
				+ pow((2.25 - vec[i]*(1 - pow(vec[i+1],2))),2)
				+ pow((2.625 - vec[i]*(1 - pow(vec[i+1],3))),2);
		}
		return sum;
	}
	//Helical Valley
	if (num==5)
	{
		if(vec[0]>0){
			theta = (double)(1/2*M_PI) * atan(vec[1]/vec[0]);
		}else{
			theta = (double)(1/2*M_PI) * atan(vec[1]/vec[0]) + 0.5;
		}
		
		sum = pow((10*(vec[2] - 10*theta)),2)
			+pow((10*((sqrt(pow(vec[0],2)+pow(vec[1],2)))-1)),2)
			+pow(vec[2],2);
			
		return sum;
	}
	
	//Wood function
	
	if(num==6)
	{
		for(i=0;i<dim;i+=4)
		{
			sum += 100*pow(vec[i+1]-vec[i]*vec[i],2)+pow(1-vec[i],2)+90*pow(vec[i+3]-vec[i+2]*vec[i+2],2)
			+ pow((1-vec[i+2]),2) + 10.1*(pow(vec[i+1]-1,2)+pow(vec[i+3]-1,2)) + 19.8*(vec[i+1]-1)*(vec[i+3]-1);
		}
		return sum;
	}
	
	//Biggs EXP6 function
	if(num==7)
	{
		float t;
		double y;
		double fi;
		for(i=0;i<dim;i++)
		{
			t = 0.1*(i+1);
			y = exp(-t) - 5*exp(-10*t) + 3*exp(-4*t);
			fi = vec[2]*exp(-t*vec[0]) - vec[3]*exp(-t*vec[1]) + vec[5]*exp(-t*vec[4]) - y;
			sum += pow(fi,2);
		}
		return sum; 
	}
	
	//Rosenbrock function
	if(num==8)
	{
		for(i=0;i<dim;i+=2)
		{
			sum += 100*pow(vec[i+1]-vec[i]*vec[i],2) + pow(1-vec[i],2);
		}
		return sum;
	}
	
	// Powell Singular function
	if(num==9)
	{
		for(i=0;i<dim;i+=4)
		{
			sum += pow(vec[i]+10*vec[i+1],2) + 5*pow(vec[i+2]-vec[i+3],2) + pow(vec[i+1]-2*vec[i+2],4) + 10*pow(vec[i]-vec[i+3],4);
		}
		return sum;
	}
	
	//Himmelblau
	if (num==10)
	{
		for(i=0;i<dim;i+=2){
			
			sum+= pow((pow(vec[i],2)+vec[i+1]-11),2)
				+ pow((vec[i] + pow(vec[i+1],2) - 7),2);			
		}
		return sum;
	} 
	
	//Grienwank
	if(num==11)
	{
		double prod=1;
		for(i=0;i<dim;i++)
		{
			sum += pow(vec[i],2);
			prod *= cos(vec[i]/sqrt(i+1));
		}
		F = (double)1/4000*sum-prod +1;

		return F;
	}
	
	//Ackley's function
	if(num==12)
	{
		double sum2=0;
		for(i=0;i<dim;i++)
		{
			sum+=vec[i]*vec[i];
			sum2+=cos(2*M_PI*vec[i]);
		}
		sum=(float)1/dim*sum;
		sum2=(float)1/dim*sum2;
		F = -20*exp((float)-1/5*sqrt(sum))-exp(sum2)+20+M_E;
		return F;
	}
	
	//Rastrigins
	if(num==13)
	{
		for(i=0;i<dim;i++){
			sum += pow(vec[i],2)-10*cos(2*M_PI*vec[i]);		
		}
		
		F=10*dim + sum;

		return F;
	}
	//Schwefel
	if (num==14)
	{
		for(i=0;i<dim;i++){
			sum+= vec[i]*sin(sqrt(abs(vec[i])));
		}
		return -sum;
	}
	
}		

/*
 * Initial random point uniformly in [min,max]
*/
double *initial_point(int min, int max, int N)
{
	int i;
	double *x = malloc(N*sizeof(double));
	
	for (i=0;i<N;i++)
	{
		x[i] = (double)rand()/(RAND_MAX)*(max-min) + min;
	}
	return x;
}

Simplex *simplex_allocation(int dim)
{
	int i,j;
	Simplex *sim = malloc(sizeof(Simplex));
	sim->vertices = malloc((dim+1)*sizeof(double*));
	sim->fvals = malloc((dim+1)*sizeof(double));
	for(i=0;i<dim+1;i++)
	{
		sim->vertices[i]=malloc(dim*sizeof(double));
	}
	/*
	for(i=0;i<dim+1;i++)
	{
		sim->fvals[i]=0;
		for(j=0;j<dim;j++)
		{
			sim->vertices[i][j]=0;
		}
	}
	*/
	return sim;
}

// Function evaluations of simplex vertices
void simplex_fvals(Simplex *sim, int dim, int func_num)
{
	int i;
	for(i=0;i<dim+1;i++)
	{
		sim->fvals[i]=functions(func_num, sim->vertices[i], dim);
	}	
	
}

/*
 * Initial simplex based on the initial random point: every point is perturbed along a specific coordinate (most used).
*/ 

Simplex *initial_simplex(double *init_p, int dim, int func_num)
{
	int i,j;
	double **unit = malloc((dim)*sizeof(double*));
	Simplex *sim = simplex_allocation(dim);
	float h=0.05;
	for(i=0;i<dim;i++)
	{
		unit[i]=malloc(dim*sizeof(double));
	}
	for(i=0;i<dim;i++)
		sim->vertices[0][i]=init_p[i];
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			if(i==j)
				unit[i][j]=1;
			else
				unit[i][j]=0;
		}		
	}			
	for(i=0;i<dim;i++)
	{
		if(sim->vertices[0][i]==0)
			h=0.00025;
		
		for(j=0;j<dim;j++)
		{	

			sim->vertices[i+1][j]=sim->vertices[0][j]+h*unit[i][j];
		}		
	}
	simplex_fvals(sim, dim, func_num);
	for(i=0;i<dim;i++)
		free(unit[i]);
	free(unit);
	return sim;
}

//Sort (BubbleSort) simplex function vertices by their function values
void sort_vertices(Simplex *sim, int dim)
{
	int swapped=1,j=0,i,k;
	double temp_f=0;
	double *temp_x = malloc(dim*sizeof(double));
	
	while(swapped)
	{
		swapped=0;
		j++;
		for(i=0;i<dim+1-j;i++)
		{
			if(sim->fvals[i]>sim->fvals[i+1])
			{
				temp_f=sim->fvals[i+1];
				sim->fvals[i+1]=sim->fvals[i];
				sim->fvals[i]=temp_f;
				
				for(k=0;k<dim;k++)
				{
					temp_x[k] = sim->vertices[i+1][k];
					sim->vertices[i+1][k]=sim->vertices[i][k];
					sim->vertices[i][k]=temp_x[k];
				}
				swapped=1;
			}
		}
	}
	
	free(temp_x);
}

Fstats stats(double *fvals, int dim)
{
	double sum=0, diff=0;
	int i;
	Fstats fs;
	for(i=0;i<dim+1;i++)
		sum+=fvals[i];
	fs.mean=sum/(dim+1);
	sum=0;
	for(i=0;i<dim+1;i++)
	{
		diff = fvals[i]-fs.mean;
		sum+=pow(diff,2);
	}
	fs.std = sum/dim;
	return fs;
}

// The mean of the simplex n best vertices
double *mean(Simplex *sim, int dim)
{
	int i,j;
	double *mean_point = malloc(dim*sizeof(double));
	for(i=0;i<dim;i++)
		mean_point[i]=0;
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			mean_point[i]+=sim->vertices[j][i];
		} 
	} 
	for(i=0;i<dim;i++)
		mean_point[i]/=dim;
	return mean_point;
}

/*
 * The simplex point transformation:
 * x' = xavg + t(xworst - xavg);
 * x': the new point
 * xworst: the point with the worst function values
 * xavg: The simplex mean
 * t: parameter that changes depending on the operation (reflection, expansion, contraction, shrinkage)
*/
double *point_transformation(double *mean_point, double *worst, float t, int dim)
{
	int i;
	double *new_p = malloc(dim*sizeof(double));
	for(i=0;i<dim;i++)
	{
		new_p[i] = mean_point[i] + t*(worst[i]-mean_point[i]);
	}
	
	return new_p;
}

// One full step of nelder-mead algorithm
void nelder_mead_step(Simplex *sim, double *mean_point,int dim,int func_num, Counters *counts, int k)
{
	int i,j,flag=0;
	double fref,fexpa,fcontrin,fcontrout;
	// Reflection point
	double *xref = point_transformation(mean_point,sim->vertices[dim],-1,dim);
	fref = functions(func_num, xref, dim);
	counts[k].num_of_evals++;
	//Reflection
	if((fref >= sim->fvals[0]) && (fref < sim->fvals[dim-1]))
	{
		counts[k].reflection_num++;
		for(i=0;i<dim;i++)
			sim->vertices[dim][i] = xref[i];
	}
	//Expansion
	else if(fref < sim->fvals[0])
	{
		double *xexpa = point_transformation(mean_point,sim->vertices[dim],-2,dim);
		fexpa = functions(func_num, xexpa, dim);
		counts[k].num_of_evals++;
		if(fexpa < fref)
		{
			counts[k].expansion_num++;
			for(i=0;i<dim;i++)
				sim->vertices[dim][i] = xexpa[i];
		}
		else
			counts[k].reflection_num++;
			for(i=0;i<dim;i++)
				sim->vertices[dim][i] = xref[i];
		free(xexpa);
	}
	else if(fref >= sim->fvals[dim-1])
	{
		//Outside contraction
		if((fref >= sim->fvals[dim-1]) && (fref < sim->fvals[dim]))
		{
			double *xcontrout = point_transformation(mean_point,sim->vertices[dim],(float)-1/2,dim);
			fcontrout = functions(func_num, xcontrout, dim);
			counts[k].num_of_evals++;
			if(fcontrout <= fref)
			{
				counts[k].contraction_num++;
				for(i=0;i<dim;i++)
					sim->vertices[dim][i] = xcontrout[i];
			}
			else
				flag = 1;
			free(xcontrout);
		}
		//Inside contraction 
		else
		{
			double *xcontrin = point_transformation(mean_point,sim->vertices[dim],(float)1/2,dim);
			fcontrin = functions(func_num, xcontrin, dim);
			counts[k].num_of_evals++;
			if(fcontrin < sim->fvals[dim])
			{
				counts[k].contraction_num++;
				for(i=0;i<dim;i++)
					sim->vertices[dim][i] = xcontrin[i];
			}
			else
				flag = 1;
			free(xcontrin);
		}	
		
		//Neither outside nor inside contraction was performed, so shrink simplex
		
		if(flag)
		{
			counts[k].shrink_num++;
			for(i=1;i<dim+1;i++)
			{
				for(j=0;j<dim;j++)
					sim->vertices[i][j] = (sim->vertices[0][j] + sim->vertices[i][j])/2;
			}
		}
	}
	free(xref);
}

void printtofile(int func_num, int dim, double minimum, double *minimizer, int shrink_num, int reflection_num, int expansion_num, int contraction_num,int num_of_iters,int num_of_evals)
{
	int i;
	FILE *fp;
	
	fp = fopen("/home/vaggos/Desktop/Master/4o/Optimization/Project/experiments.txt","a");
	
	switch (func_num)	{
		case 1:
			fprintf(fp, "Function : Sum of squares\n");
			printf("Function : Sum of squares\n");
			break;

		case 2:
			fprintf(fp, "Function : Freudenstein and Roth");
			printf("Function : Freudenstein and Roth");
			break;

		case 3:
			fprintf(fp, "Function : Powell");
			printf("Function : Powell");
			break;
		case 4:
			fprintf(fp, "Function : Beale");
			printf("Function : Beale");
			break;
		case 5:
			fprintf(fp, "Function : Helical Valley");
			printf("Function : Helical Valley");
			break;
		case 6:
			fprintf(fp, "Function : Wood");
			printf("Function : Wood");
			break;
		case 7:
			fprintf(fp, "Function : BIGGS EXP6");
			printf("Function : BIGGS EXP6");
			break;
		case 8:
			fprintf(fp, "Function : Rosenbrock");
			printf("Function : Rosenbrock");
			break;
		case 9:
			fprintf(fp, "Function : Powell Singular Function");
			printf("Function : Powell Singular Function");
			break;
		case 10:
			fprintf(fp, "Function : Himmelblau");
			printf("Function : Himmelblau");
			break;
		case 11:
			fprintf(fp, "Function : Grienwank");
			printf("Function : Grienwank");
			break;
		case 12:
			fprintf(fp, "Function : Ackley's");
			printf("Function : Ackley's");
			break;
		case 13:
			fprintf(fp, "Function : Rastrigins");
			printf("Function : Rastrigins");
			break;
		case 14:
			fprintf(fp, "Function : Schwefel");
			printf("Function : Schwefel");
			break;
	}
	
	fprintf(fp,"Dimension: %d\n", dim);
	printf("Dimension: %d\n", dim);
	fprintf(fp,"Function evaluation: %f\n", minimum);
	printf("Function evaluation: %f\n", minimum);
	
	for(i=0;i<dim;i++){
		fprintf(fp,"Minimizer: %f ", minimizer[i]);
		printf("Minimizer: %f ", minimizer[i]);
	}
	fprintf(fp,"\n");
	printf("\n");
	
	fprintf(fp,"Number of Shrinks: %d\n", shrink_num);
	fprintf(fp,"Number of reflections: %d\n", reflection_num);
	fprintf(fp,"Number of expansion: %d\n", expansion_num);
	fprintf(fp,"Number of contractions: %d\n", contraction_num);
	fprintf(fp,"Reflections/(Expansions+Contractions): %f\n", (float)reflection_num/(expansion_num+contraction_num+reflection_num));
	
	printf("Number of Shrinks: %d\n", shrink_num);
	printf("Number of reflections: %d\n", reflection_num);
	printf("Number of expansion: %d\n", expansion_num);
	printf("Number of contractions: %d\n", contraction_num);
	printf("Reflections/(Expansions+Contractions): %f\n", (float)reflection_num/(expansion_num+contraction_num));
	
	fprintf(fp,"Number of Iterations: %d\n", num_of_iters);
	fprintf(fp,"Number of function Evaluations: %d\n", num_of_evals);
	
	printf("Number of Iterations: %d\n", num_of_iters);
	printf("Number of function Evaluations: %d\n", num_of_evals);
	fprintf(fp,"\n");
	printf("\n");
	fclose(fp);
}