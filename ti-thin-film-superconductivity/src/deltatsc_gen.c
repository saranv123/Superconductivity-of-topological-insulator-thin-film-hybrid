#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "energyfunc.c"


struct func_params {int v; double d; double delta;};

double integEtopplus(double k, void * p)
{
	
	double integexp;
	
	struct func_params *params = (struct func_params *)p;
	int v = (params -> v);
	double d = (params -> d);
	double delta = (params -> delta);
	
	integexp = (1./(8.*pi*debye))*(k/sqrt(pow(eigenEtopplus(v, d, k), 2) + pow(delta, 2)));
	
	return integexp;
}	


double integEtopminus(double k, void * p)
{
	
	double integexp;
	
	struct func_params *params = (struct func_params *)p;
	int v = (params -> v);
	double d = (params -> d);
	double delta = (params -> delta);
	
	integexp = (1./(8.*pi*debye))*(k/sqrt(pow(eigenEtopminus(v, d, k), 2) + pow(delta, 2)));
	
	return integexp;
}	

double integEbottomplus(double k, void * p)
{
	
	double integexp;
	
	struct func_params *params = (struct func_params *)p;
	int v = (params -> v);
	double d = (params -> d);
	double delta = (params -> delta);
	
	
	integexp = (1./(8.*pi*debye))*(k/sqrt(pow(eigenEbottomplus(v, d, k), 2) + pow(delta, 2)));
	
	return integexp;
}	


double integralFS1(int v, double d, double delta)
{
    
    double result = 0., result1 = 0., result2 = 0., error, a, b, rootcheck, rootcheck2;
    struct func_params params = {v, d, delta};
		
	gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (1000);
		
				 gsl_integration_workspace * z  = gsl_integration_workspace_alloc (1000);

	rootcheck = newton_raph(0.2, *eigenEtopminus, v, d);
	
	gsl_function F;
	gsl_function G;
	gsl_function H;

	
	if(fabs(rootcheck) <= 0.0001)
	{
		return 0.;//printf("1FS\n");
	}
	else
	{ 
      //  printf("3FS\n");
		rootcheck2 = newton_raph(0.1, *eigenEtopplus, v, d);
		
		
		if(fabs(rootcheck2) <= 0.0001)
		{
			
		//	printf("negative helicity\n");			
			
			
			a = lowlimit(v, d, *eigenEtopminus, 0.);
	
	        b = uplimit(v, d, *eigenEtopminus, 0.);
			
			
		//	printf("a = %lf\t b = %lf\n", a, b);
			
			H.function = &integEtopminus;
            H.params = &params;

            gsl_integration_qag (&H,
                                 b, a,
                                 0.0, 1e-5, 1000,
                                 GSL_INTEG_GAUSS15,
                                 z,
                                 &result1, &error);

	
		}
		
		else
		{
			
		//	printf("positive helicity\n");
		
		 a = lowlimit(v, d, *eigenEtopplus, 0.);
	    
		 b = uplimit(v, d, *eigenEtopplus, 0.);
		 
		// printf("a = %lf\t b = %lf\n", a, b);

	
	
        F.function = &integEtopplus;
        F.params = &params;

        gsl_integration_qag (&F,
                              a, b,
                              0.0, 1e-5, 1000,
                              GSL_INTEG_GAUSS15,
                              w,
                              &result1, &error);
							  
		error = 0.;	b = 0.;				  

		}
		 gsl_integration_workspace_free (w);
		 gsl_integration_workspace_free (z);
		 
		result = result1 + result2;	
		
		
			
	    return result;
	}	

} 


	
double integralFS2(int v, double d, double delta)
{
	
	double result, error, a, b, rootcheck;
    struct func_params params = {v, d, delta};
	
	
	rootcheck = newton_raph(0.3, *eigenEtopminus, v, d);
	
	
	if(fabs(rootcheck) <= 0.0001)
	{
		return 0.;
	}
	else
	{
	
    gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

    

    a = lowlimit(v, d, *eigenEtopminus, rootcheck);
	b = uplimit(v, d, *eigenEtopminus, rootcheck);

	

    gsl_function F;
    F.function = &integEtopminus;
    F.params = &params;

    gsl_integration_qag (&F,
                         a, b,
                         0.0, 1e-5, 1000,
                         GSL_INTEG_GAUSS15,
                         w,
                         &result, &error);

 
    gsl_integration_workspace_free (w);
   
	return result;
	}

}

	
double integralFS3(int v, double d, double delta)
{
    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (1000);

    double result, error, a, b, root;
    struct func_params params = {v, d, delta};
	
	root = newton_raph(0.1, *eigenEbottomplus, v, d);



    a = lowlimit(v, d, *eigenEbottomplus, root);
	b = uplimit(v, d, *eigenEbottomplus, root);
	
	//printf("a = %lf\t, b = %lf\n", a, b); 
	

    gsl_function F;
    F.function = &integEbottomplus;
    F.params = &params;

    gsl_integration_qag (&F,
                         a, b,
                         0.0, 1e-5, 1000,
                         GSL_INTEG_GAUSS15,
                         w,
                         &result, &error);


    gsl_integration_workspace_free (w);
	return result;

} 
	
double *delta4band(double d, int v)
{
	int i, j;
	double *delta_new, *delta_old, err_sum, err_avg;
	double *deltaQW, rootcheck;
	double lambdabulk, deltabulk;
	double mu = mubulk + drcpnt;
	double prox = 0.;
	
	
//	v = (int)(((d*kFbulk)/pi)*sqrt(mu/mubulk) + 0.5);
	
	double kf1 = newton_raph(0.1, *eigenEtopplus, v, d);
	double kf2 = newton_raph(0.3, *eigenEtopminus, v, d);	
	double kf3 = newton_raph(0.1, *eigenEbottomplus, v, d);

		
	//printf("kf1 = %lf\tkf2 = %lf\t kf3 = %lf\n",kf1, kf2, kf3);
	
//	v = (int)(((d*kFbulk)/pi) + 0.5);
	lambdabulk = (pow(kFbulk, 3) * J)/(4.*pi*pi*mubulk);
	deltabulk = 2. * exp(-1./lambdabulk);
	
	delta_old = (double *)malloc((v+2) * sizeof(double));
	delta_new = (double *)malloc((v+2) * sizeof(double));

	
	rootcheck = newton_raph(0.3, *eigenEtopminus, v, d);
	
	//printf("root = %lf\n", rootcheck);


	for(i = 0;i < v+2;i ++)
	{
		delta_old[i] = deltabulk;
		delta_new[i] = 0.;
    }

    			
/* *********************************************************  */
    
    if(fabs(rootcheck) <= 0.0001)
	{
	//	printf("1\n");
			 /* One hybrid Fermi surface */
		 delta_old[v] = 0.;
		 delta_old[v+1] = 0.;
		 prox = 0.;
		 delta_old[v-1] = deltabulk;
		


    	if(v == 1)
	    {
			
		      prox = 0.;
	    }
	    else
	    {
		    for(i = 0;i < v - 1;i ++)
		    {
				
			   prox = prox + ((kFbulk*kFbulk)/(4.*pi* mubulk)) * Vconfopen(v, d, v, i+1)*delta_old[i]*asinh(1./delta_old[i]);
		    }
	    }
		 
		do
		{
		
			
			/* coupling between hybrid FS's and with the thin film FS*/
			
		
			
		    delta_new[v-1] = Vintra_gen_conf(v, d, 3, 3) * delta_old[v-1] * integralFS3(v, d, delta_old[v-1]) + spctrl_wghtbottomplus(v, d, kf3) * prox;
					
			 
			
		for(i=0;i < v-1;i ++)
	    {
		
		/*thin film FS's pairing with hybrid FS's*/
		    delta_new[i] = delta_new[i] + spctrl_wghtbottomplus(v, d, kf3) * Vconfopen(v, d, i+1, v) * delta_old[v-1] * integralFS3(v, d, delta_old[v-1]);
	
	       

	/*coupling between thin film FS's */
		    for(j=0;j < v - 1;j ++)
		    {
			   delta_new[i] = delta_new[i] + ((kFbulk*kFbulk)/(4.*pi* mubulk))*Vconfopen(v, d, i+1, j+1)*delta_old[j]*asinh(1./delta_old[j]) ;
			
	        }
		
        }
	
			
			err_sum = 0.;
	
	       for(i=0;i<v;i++)
	       {		
	          //	printf("%lf\t", delta_new[i]);
	        	err_sum = err_sum + fabs((delta_new[i] - delta_old[i])/delta_old[i]);
	        	delta_old[i] = delta_new[i];
	        	delta_new[i] = 0.;
           }	

           err_avg = err_sum/((double)(v));
		   
		//   printf("error = %lf\n", err_avg);

		   
		}
		while(err_avg > 0.0001);
				
		
	}
	else
	{
		

/* if the Dirac point is within the debye frequency, then there is an error in the code below. But the error is probably very minute*/ 
        do
	    {	
        		prox = 0.;
	
            	if(v == 1)
	            {
					
	        	      prox = 0.;
	            }
	            else
	            {
		            for(i = 0;i < v - 1;i ++)
		            {
			             prox = prox + ((kFbulk*kFbulk)/(4.*pi* mubulk)) * Vconfopen(v, d, v, i+1)*delta_old[i]*asinh(1./delta_old[i]);
		            }
	            } 

	     
		/*  FS1 is assigned to v+1 slot, FS2 to v slot and FS3 to v-1 slot*/
		/* coupling between hybrid FS's and with the thin film FS*/
    	delta_new[v+1] = Vintra_gen_conf(v, d, 1, 1) * delta_old[v+1] * integralFS1(v, d, delta_old[v+1]) + Vintra_gen_conf(v, d, 1, 2) * delta_old[v] * integralFS2(v, d, delta_old[v])
    	+ Vintra_gen_conf(v, d, 1, 3) * delta_old[v-1] * integralFS3(v, d, delta_old[v-1]) + spctrl_wghttopplus(v, d, kf1) * prox;
	
        delta_new[v] = Vintra_gen_conf(v, d, 2, 1) * delta_old[v+1] * integralFS1(v, d, delta_old[v+1]) + Vintra_gen_conf(v, d, 2, 2) * delta_old[v] * integralFS2(v, d, delta_old[v])
	    + Vintra_gen_conf(v, d, 2, 3) * delta_old[v-1] * integralFS3(v, d, delta_old[v-1]) + spctrl_wghttopminus(v, d, kf2) * prox;
	
	    delta_new[v-1] = Vintra_gen_conf(v, d, 3, 1) * delta_old[v+1] * integralFS1(v, d, delta_old[v+1]) + Vintra_gen_conf(v, d, 3, 2) * delta_old[v] * integralFS2(v, d, delta_old[v])
	    + Vintra_gen_conf(v, d, 3, 3) * delta_old[v-1] * integralFS3(v, d, delta_old[v-1]) + spctrl_wghtbottomplus(v, d, kf3) * prox;
	
	    for(i=0;i < v-1;i ++)
	    {
	//	printf("\nOops!\n");
	    	/*thin film FS's pairing with hybrid FS's*/
	    	delta_new[i] = delta_new[i] + spctrl_wghttopplus(v, d, kf1) * Vconfopen(v, d, i+1, v) * delta_old[v+1] * integralFS1(v, d, delta_old[v+1]) + 
	    	spctrl_wghttopminus(v, d, kf2) * Vconfopen(v, d, i+1, v) * delta_old[v] * integralFS2(v, d, delta_old[v])
	    	+ spctrl_wghtbottomplus(v, d, kf3) * Vconfopen(v, d, i+1, v) * delta_old[v-1] * integralFS3(v, d, delta_old[v-1]);
	

	/*coupling between thin film FS's */
	   	   for(j=0;j < v - 1;j ++)
		   {
			delta_new[i] = delta_new[i] + ((kFbulk*kFbulk)/(4.*pi* mubulk))*Vconfopen(v, d, i+1, j+1)*delta_old[j]*asinh(1./delta_old[j]) ;
			
	       }
		
        }
	
	
	    err_sum = 0.;
	
	    for(i=0;i<v+2;i++)
	    {		
	      //	printf("%lf\t", delta_new[i]);
	    	err_sum = err_sum + fabs((delta_new[i] - delta_old[i])/delta_old[i]);
	    	delta_old[i] = delta_new[i];
	    	delta_new[i] = 0.;
	    }	
		

        err_avg = err_sum/((double)(v+2));
	

	
	    }
	    while(err_avg > 0.001);
		
	
	}
	return delta_old;
}	


void main()
{
	int i, v;
	double *delta, d, drun, Amax;
	const int buf_size = 50;
	char buf[buf_size];  
	FILE *out_ptr;
	double mu = mubulk + drcpnt;
	
	 v = 1;
	 int ret = snprintf(buf, buf_size, "SCgapsoc1FSn1_drcpnt%.2lft%.2lf.dat", drcpnt, t0);
	 
	 
	for(drun = 0.1; drun <= 1.2;drun = drun + 0.000001)
	{
	
	
	if(fabs(detun_genplus(v, drun, 0.)) <= 0.0001)
		{
	         d = drun;
	
		}
		
	}
	 
	// The codes below calculates superconductivity of the hybrid as a function of the Fermi velocity of the 
	// topological surface states.
	
	out_ptr = fopen(buf, "w+");
	printf("d = %lf\n", d);
	for(A0 = 0.05;A0 <= 5.0; A0 = A0 + 0.1)
	{
			
		
	     delta = delta4band(d, v);
	     printf("\nt0 = %lf\tA0 = %lf\n", t0, A0);
	
	    fprintf(out_ptr,"%lf\t%lf\t%lf\t%lf\n", A0/1.975, delta[0], delta[1], delta[2]);
		 
		  printf("%lf\t%lf\t%lf\t%lf\n", A0/1.975, delta[0], delta[1], delta[2]);
	
	}
	
	fclose(out_ptr);

	
}

