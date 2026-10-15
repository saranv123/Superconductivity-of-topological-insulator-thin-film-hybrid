#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<stdlib.h>
#define pi 3.14159


double A0;
double drcpnt = -0.05;
double t0 = 0.2;


double mubulk = 9.47;
double kFbulk = 1.58;
double n0 = 0.132; 
double J = 37.03;
double debye_conv = 0.633;

/*
double detun(int v, double d)
{
	
	double detune;
	
	detune = (mubulk - ((mubulk*pi*pi)/(kFbulk*kFbulk))*(pow(v,2)/pow(d,2)))/((2*t0)/sqrt(d));
	
	return detune;
	
}

*/

double V(int n, double d, int i, int j)
{
	
	double coupl;
	
			if(i == j)
				coupl = (1.5*J)/d;
			else
				coupl = J/d;
	
		
	return coupl;
}
/* let us set open boundary conditions for phonons too*/
double Vconfopen(int n, double d, int i, int j)
{
	int l_max, l;
	double coupl, beta;
	double ldouble, idouble, jdouble, kD, kF;
	double mu = mubulk + drcpnt;
	
	beta = 0.;
	
	idouble = (double)i;
	jdouble = (double)j;
	
	kF = kFbulk * sqrt(mu/mubulk);
	
	kD = pow(2.,1./3.)*kFbulk;
	//kD = kF;
	
	l_max = (int)(((d*kD)/pi) + 0.5);

		
	for(l = 1; l <= l_max; l++)
	{
			ldouble = (double)l;
			beta = beta + pow(((1./pi)*(((ldouble - 0.5)/(pow((ldouble - 0.5), 2.) - pow(idouble - jdouble, 2.))) 
			 - ((ldouble - 0.5)/(pow((ldouble - 0.5),2.) - pow(idouble + jdouble - 1., 2.))))),2.);
			ldouble = 0.;
	
	}
	
	coupl = (2.* beta*J)/d;
	
	return coupl;
}

double Vconfopenphclsd(int n, double d, int i, int j)
{
	int l_max, l;
	double coupl, beta;
	double ldouble, idouble, jdouble, kD, kF;
	double mu = mubulk + drcpnt;
	
	beta = 0.;
	
	idouble = (double)i;
	jdouble = (double)j;
	
	kF = kFbulk * sqrt(mu/mubulk);
	
	//kD = pow(2.,1./3.)*kF;
	kD = kF;
	
	l_max = (int)(((d*kD)/pi) + 0.5);

		
	for(l = 1; l <= l_max; l++)
	{
		if((i + j + l)%2 == 0)
	    {
			ldouble = (double)l;
			beta = beta + pow(((2./pi)*(ldouble/(ldouble*ldouble - pow(idouble + jdouble - 1., 2)))), 2);
			ldouble = 0.;
		}
		else
		{
			ldouble = (double)l;
			beta = beta + pow(((2./pi)*(ldouble/(ldouble*ldouble - pow(idouble - jdouble, 2)))), 2);
			ldouble = 0.;
		}
	}
	
	coupl = (2.* beta*J)/d;
	
	return coupl;
}


double Vconfclsd(int n, double d, int i, int j)
{
	int l_max, l;
	double coupl, beta, kF, kD;
	double ldouble, idouble, jdouble;
	double mu = mubulk + drcpnt;
	
	beta = 0.;
	
	kF = kFbulk * sqrt(mu/mubulk);
	
	kD = kF;
	
	
	
	l_max = (int)((d*kD)/pi);
	
/*	l_max = (int) (((d*kFbulk)/pi)*sqrt(mu/mubulk));*/
  
/*
	l_max = (int) (debye_conv*d);
	*/
	idouble = (double)i;
	jdouble = (double)j;
			

	
	if((i + j)%2 == 0)
	{
		for(l = 1; l <= l_max; l= l+2)
		{
			ldouble = (double)l;
			beta = beta + pow(((2./pi)*((ldouble/ (ldouble*ldouble - pow(idouble - jdouble, 2))) - (ldouble/ (ldouble*ldouble - pow(idouble + jdouble, 2))))),2);
			ldouble = 0.;
		}	  
	}

		else
	{
		for(l = 0; l <= l_max; l= l+2)
		{
			ldouble = (double)l;
			beta = beta + pow(((2./pi)*((ldouble/ (ldouble*ldouble - pow(idouble - jdouble, 2))) - (ldouble/ (ldouble*ldouble - pow(idouble + jdouble, 2))))),2);
			ldouble = 0.;
		}
	}	
	
	/*
	
	for(l = 1; l <= l_max; l++)
	{
		if((i + j + l)%2 != 0)
	    {
			ldouble = (double)l;
			beta = beta + pow(((2./pi)*(ldouble/(ldouble*ldouble - pow(idouble - jdouble, 2)) - ldouble/(ldouble*ldouble - pow(idouble + jdouble, 2)))), 2);
			ldouble = 0.;
		}
		
	} 
	*/
	coupl = (2.* beta*J)/d;
	
	return coupl;
}	



int subband(double d)
{
	int n, m;
	double mu1;
	
	m = 10;
	
	do
	{
        n = m;	    
		
	    mu1 = ((pi*d*2.*mubulk)/(n*kFbulk*kFbulk))*(n0 + ((pi/(6.*d*d*d))*(double)n*((double)n+0.5)*((double)n+1.)));  
	
    	m = (int) ((d/pi)*sqrt(kFbulk*kFbulk*mu1/mubulk));
	   
   	
	}	
	
	while(m != n);
	
	return m;
}	


double *delta_confined(double l)
{
		
    int i, j, v;
	double *delta_old, *delta_new;
	double deltabulk, lambdabulk, err_sum, err_avg;
	double mu = mubulk + drcpnt;
		
	v = (int)(((l*kFbulk)/pi)*sqrt(mu/mubulk) + 0.5);
	err_sum = 0.;
	lambdabulk = (pow(kFbulk, 3) * J)/(4*pi*pi*mubulk);
	deltabulk = 2. * exp(-1./lambdabulk);
	
	delta_old = (double *)malloc((v-1) * sizeof(double));
	delta_new = (double *)malloc((v-1) * sizeof(double));
	
	for(i=0;i< v-1;i++)
	{
		delta_old[i] = deltabulk;
	}
	
	do
	{
		err_sum = 0.;
	for(i=0;i < v - 1;i ++)
	{
		
		for(j=0;j < v - 1;j ++)
		{
			delta_new[i] = delta_new[i] + ((kFbulk*kFbulk)/(4.*pi* mubulk))*Vconfopen(v, l, i+1, j+1)*delta_old[j]*asinh(1./delta_old[j]);
			
	    }
		
    }
	
	for(i=0;i < v - 1;i++)
	{
		err_sum = err_sum + fabs((delta_new[i] - delta_old[i])/delta_old[i]);
		delta_old[i] = delta_new[i];
		delta_new[i] = 0.;
		
	}	
	
	

    err_avg = err_sum/v;
	}
    while(err_avg > 0.001);	
	
	
	
	return delta_old;
}

double *delta_confinedthfmopen(double l)
{
		
    int i, j, v;
	double *delta_old, *delta_new;
	double deltabulk, lambdabulk, err_sum, err_avg;
	double mu = mubulk + drcpnt;
		
	v = (int)(((l*kFbulk)/pi)*sqrt(mu/mubulk) + 0.5);
	err_sum = 0.;
	lambdabulk = (pow(kFbulk, 3) * J)/(4*pi*pi*mubulk);
	deltabulk = 2. * exp(-1./lambdabulk);
	
	delta_old = (double *)malloc(v * sizeof(double));
	delta_new = (double *)malloc(v * sizeof(double));
	
	for(i=0;i<v;i++)
	{
		delta_old[i] = deltabulk;
	}
	
	do
	{
		err_sum = 0.;
	for(i=0;i<v;i ++)
	{
		
		for(j=0;j<v;j ++)
		{
			delta_new[i] = delta_new[i] + ((kFbulk*kFbulk)/(4.*pi* mubulk))*Vconfopen(v, l, i+1, j+1)*delta_old[j]*asinh(1./delta_old[j]);
			
	    }
		
    }
	
	for(i=0;i<v;i++)
	{
		err_sum = err_sum + fabs((delta_new[i] - delta_old[i])/delta_old[i]);
		delta_old[i] = delta_new[i];
		delta_new[i] = 0.;
		
	}	
	
	

    err_avg = err_sum/v;
	}
    while(err_avg > 0.001);	
	
	
	
	return delta_old;
}

double *delta_confinedthfmclsd(double l)
{
		
    int i, j, v;
	double *delta_old, *delta_new;
	double deltabulk, lambdabulk, err_sum, err_avg;
	double mu = mubulk + drcpnt;
		
	v = (int)(((l*kFbulk)/pi)*sqrt(mu/mubulk) + 0.5);
	err_sum = 0.;
	lambdabulk = (pow(kFbulk, 3) * J)/(4*pi*pi*mubulk);
	deltabulk = 2. * exp(-1./lambdabulk);
	
	delta_old = (double *)malloc(v * sizeof(double));
	delta_new = (double *)malloc(v * sizeof(double));
	
	for(i=0;i<v;i++)
	{
		delta_old[i] = deltabulk;
	}
	
	do
	{
		err_sum = 0.;
	for(i=0;i<v;i ++)
	{
		
		for(j=0;j<v;j ++)
		{
			delta_new[i] = delta_new[i] + ((kFbulk*kFbulk)/(4.*pi* mubulk))*Vconfclsd(v, l, i+1, j+1)*delta_old[j]*asinh(1./delta_old[j]);
			
	    }
		
    }
	
	for(i=0;i<v;i++)
	{
		err_sum = err_sum + fabs((delta_new[i] - delta_old[i])/delta_old[i]);
		delta_old[i] = delta_new[i];
		delta_new[i] = 0.;
		
	}	
	
	

    err_avg = err_sum/v;
	}
    while(err_avg > 0.001);	
	
	
	
	return delta_old;
}




double *delta_bulk(double l)
{
		
    int i, j, v;
	double *delta_old, *delta_new;
	double deltabulk, lambdabulk, err_sum, err_avg;
	double mu = mubulk + drcpnt;
		
	v = (int)(((l*kFbulk)/pi)*sqrt(mu/mubulk) + 0.5);
	err_sum = 0.;
	lambdabulk = (pow(kFbulk, 3) * J)/(4*pi*pi*mubulk);
	deltabulk = 2. * exp(-1/lambdabulk);
		
	
	delta_old = (double *)malloc((v) * sizeof(double));
	delta_new = (double *)malloc((v) * sizeof(double));
	
	for(i=0;i<v;i++)
	{
		delta_old[i] = deltabulk;
	}
	
	do
	{
		err_sum = 0.;
	for(i=0;i<v;i ++)
	{
		
		for(j=0;j<v;j ++)
		{
			delta_new[i] = delta_new[i] + ((kFbulk*kFbulk)/(4.*pi* mubulk))*V(v, l, i+1, j+1)*delta_old[j]*asinh(1./delta_old[j]);
			
	    }
		

    }
	
	for(i=0;i < v;i++)
	{
		err_sum = err_sum + fabs((delta_new[i] - delta_old[i])/delta_old[i]);
		delta_old[i] = delta_new[i];
		delta_new[i] = 0.;
	}	

    err_avg = err_sum/v;
	}
    while(err_avg > 0.0001);	
	
	
	
	return delta_old;
}

