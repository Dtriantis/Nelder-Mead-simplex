#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void main(int argc, char *argv[])
{
	int i,num,dim;
	float min,max;
	double *vec,dom;
	num = atoi(argv[1]);
	dim = atoi(argv[2]);
	
	vec = malloc(dim*sizeof(double));
	
	for(i=0;i<dim;i++){
		printf("Dwse init\n");
		scanf("%lf",&dom);
		vec[i]=dom;
	}
	double sum=0, theta, F;
	// Simple testing function: Σxi^2
	if (num==1)
	{
		for(i=0;i<dim;i++)
			sum+=(vec[i])*(vec[i]);
		printf("Fval: %f\n",sum);
	}
	//Freudenstein and Roth
	if (num==2)
	{
		for(i=0;i<dim;i+=2){
			sum += pow((-13 + vec[i] + ((5 - vec[i+1])*vec[i+1] - 2)*vec[i+1]),2)
				+ pow((-29 +  vec[i] + ((vec[i+1] + 1)*vec[i+1] -14)*vec[i+1]),2);
		}
		printf("Fval: %f\n",sum);
	}
	// Powell
	if (num==3)
	{
		for(i=0;i<dim;i+=2){
			sum += pow((10000*vec[i]*vec[i+1] - 1),2)
				+ pow((exp(-vec[i]) + (exp(-vec[i+1])) - 1.0001),2);
		}
		printf("Fval: %f\n",sum);
	}
	// Beale
	if (num==4)
	{
		for(i=0;i<dim;i+=2){
			sum += pow((1.5 - vec[i]*(1 - vec[i+1])),2)
				+ pow((2.25 - vec[i]*(1 - pow(vec[i+1],2))),2)
				+ pow((2.625 - vec[i]*(1 - pow(vec[i+1],3))),2);
		}
		printf("Fval: %f\n",sum);
	}
	//Helical Valley
	if (num==5)
	{
		if(vec[0]>0){
			theta = ((double)1/2*M_PI) * atan2(vec[1],vec[0]);
		}else{
			theta = ((double)1/2*M_PI) * atan2(vec[1],vec[0]) + 0.5;
		}
		
		sum = pow((10*(vec[2] - 10*theta)),2)
			+pow((10*((sqrt(pow(vec[0],2)+pow(vec[1],2)))-1)),2)
			+pow(vec[2],2);
			
		printf("Fval: %f\n",sum);
	}
	
	//Wood function
	
	if(num==6)
	{
		for(i=0;i<dim;i+=4)
		{
			sum += 100*pow(vec[i+1]-vec[i]*vec[i],2)+pow(1-vec[i],2)+90*pow(vec[i+3]-vec[i+2]*vec[i+2],2)
			+ pow((1-vec[i+2]),2) + 10.1*(pow(vec[i+1]-1,2)+pow(vec[i+3]-1,2)) + 19.8*(vec[i+1]-1)*(vec[i+3]-1);
		}
		printf("Fval: %f\n",sum);
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
		printf("Fval: %f\n",sum);
	}
	
	//Rosenbrock function
	if(num==8)
	{
		for(i=0;i<dim;i+=2)
		{
			sum += 100*pow(vec[i+1]-vec[i]*vec[i],2) + pow(1-vec[i],2);
		}
		printf("Fval: %f\n",sum);
	}
	
	// Powell Singular function
	if(num==9)
	{
		for(i=0;i<dim;i+=4)
		{
			sum += pow(vec[i]+10*vec[i+1],2) + 5*pow(vec[i+2]-vec[i+3],2) + pow(vec[i+1]-2*vec[i+2],4) + 10*pow(vec[i]-vec[i+3],4);
		}
		printf("Fval: %f\n",sum);
	}
	
	//Himmelblau
	if (num==10)
	{
		for(i=0;i<dim;i+=2){
			
			sum+= pow((pow(vec[i],2)+vec[i+1]-11),2)
				+ pow((vec[i] + pow(vec[i+1],2) - 7),2);			
		}
		printf("Fval: %f\n",sum);
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

		printf("Fval: %f\n",F);
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
		sum=(double)1/dim*sum;
		sum2=(double)1/dim*sum2;
		F = -20*exp((double)-1/5*sqrt(sum))-exp(sum2)+20+M_E;
		printf("Fval: %f\n",F);
	}
	
	//Rastrigins
	if(num==13)
	{
		for(i=0;i<dim;i++){
			sum += pow(vec[i],2)-10*cos(2*M_PI*vec[i]);		
		}
		
		F=10*dim + sum;

		printf("Fval: %f\n",F);
	}
	//Schwefel
	if (num==14)
	{
		for(i=0;i<dim;i++){
			sum+= vec[i]*sin(sqrt(abs(vec[i])));
		}
		printf("Fval: %f\n",-sum);;
	}
	
}	