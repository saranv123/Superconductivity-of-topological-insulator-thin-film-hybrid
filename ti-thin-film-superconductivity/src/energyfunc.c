#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<stdlib.h>
#include "hybrid_resonant_funcs.c"

double debye = 0.009;

double deriv(double x0, double h, double (*func)(int v, double d, double x), int v, double d)
{
	double ans;								
	
		ans = ((*func)(v, d, x0 + h) - (*func)(v, d, x0))/h;
	
	return(ans);
}


double newton_raph(double m, double (*func)(int v, double d, double x), int v, double d)
{

int i = 0, n = 40;
double h = 1e-6;
double c, s;

c = m; 

while(i <= n)
{

s = c - ((*func)(v, d, c)/(deriv(c, h, func, v, d)));

c = s;

i++;

}

if(fabs((*func)(v, d, s)) <= 0.01)
	return s;
else
	return 0.;

}

double Ethinfilm(int v, double d, double k)
{
	double Eth;
	double mu = mubulk + drcpnt;
	
	Eth = (mubulk* pow(k,2))/(kFbulk*kFbulk) + (mubulk/(kFbulk*kFbulk))*pow((((double)v - 0.5)*pi)/d, 2)  - mu;

	
	return Eth;
}

double ETIplus(double k)
{
	double eti;
	
	eti = A0*sqrt(pow(k,2)) - drcpnt;
	
	return eti;
}

double ETIminus(double k)
{
	double eti;
	
	eti = -A0*sqrt(pow(k,2)) - drcpnt;
	
	return eti;
}
	
double eigenEtopplus(int v, double d, double k)
{
	double E1, E1dimless;
	
	E1 = 0.5*(Ethinfilm(v, d, k) + ETIplus(k)) + sqrt(pow(0.5*(Ethinfilm(v, d, k) - ETIplus(k)), 2) + (pow(t0,2)/d));
	
	E1dimless = E1/debye;
	
	return E1dimless;
}	

double eigenEtopminus(int v, double d, double k)
{
	double E1, E1dimless;
	
	E1 = 0.5*(Ethinfilm(v, d, k) + ETIminus(k)) + sqrt(pow(0.5*(Ethinfilm(v, d, k) - ETIminus(k)), 2) + (pow(t0,2)/d));
	
	E1dimless = E1/debye;
	
	return E1dimless;
}	
	
double eigenEbottomplus(int v, double d, double k)
{
	double E1, E1dimless;
	
	E1 = 0.5*(Ethinfilm(v, d, k) + ETIplus(k)) - sqrt(pow(0.5*(Ethinfilm(v, d, k) - ETIplus(k)), 2) + (pow(t0,2)/d));
	
	E1dimless = E1/debye;
	
	return E1dimless;
}	

double eigenEbottomminus(int v, double d, double k)
{
	double E1, E1dimless;
	
	E1 = 0.5*(Ethinfilm(v, d, k) + ETIminus(k)) - sqrt(pow(0.5*(Ethinfilm(v, d, k) - ETIminus(k)), 2) + (pow(t0,2)/d));
	
	E1dimless = E1/debye;
	
	return E1dimless;
}	

/*effective potential*/

double detun_genplus(int v, double d, double k)
{
	double detun;
	
		detun = (Ethinfilm(v, d, k) - ETIplus(k))/((2.*t0)/sqrt(d));
		
	return detun;
}

double detun_genminus(int v, double d, double k)
{
	double detun;
	
		detun = (Ethinfilm(v, d, k) - ETIminus(k))/((2.*t0)/sqrt(d));
		
	return detun;
}		

double spctrl_wghttopplus(int v, double d, double k)
{
	double spc_wght;
	
	spc_wght = 0.5 * (1. + (detun_genplus(v, d, k)/sqrt(1. + pow(detun_genplus(v, d, k), 2))));
	
	return spc_wght;
}	

double spctrl_wghttopminus(int v, double d, double k)
{
	double spc_wght;
	
	spc_wght = 0.5 * (1. + (detun_genminus(v, d, k)/sqrt(1. + pow(detun_genminus(v, d, k), 2))));
	
	return spc_wght;
}	

double spctrl_wghtbottomplus(int v, double d, double k)
{
	double spc_wght;
	
	spc_wght = 0.5 * (1. - (detun_genplus(v, d, k)/sqrt(1. + pow(detun_genplus(v, d, k), 2))));
	
	return spc_wght;
}

double spctrl_wghtbottomminus(int v, double d, double k)
{
	double spc_wght;
	
	spc_wght = 0.5 * (1. - (detun_genminus(v, d, k)/sqrt(1. + pow(detun_genminus(v, d, k), 2))));
	
	return spc_wght;
	
	
}

double Vintra_gen_conf(int v, double d, int i, int j)
{
    double Vintra = 0.;
	double kf1, kf2, kf3;
	
	kf1 = newton_raph(0.1, *eigenEtopplus, v, d);
	
	kf2 = newton_raph(0.2, *eigenEtopminus, v, d);
	
	kf3 = newton_raph(0.1, *eigenEbottomplus, v, d);
	
	
	if(kf2 == 0.)
	{
		if((i == 3)&&(j == 3))
			Vintra = Vconfopen(v, d, v, v) * spctrl_wghtbottomplus(v, d, kf3) * spctrl_wghtbottomplus(v, d, kf3);
		else
			Vintra = 0.;
	}
	else
	{

	    if((i == 1)&&(j == 1))
	    {				
		   Vintra = Vconfopen(v, d, v, v) * spctrl_wghttopplus(v, d, kf1) * spctrl_wghttopplus(v, d, kf1);
	    }
	    else if(((i == 1)&&(j == 2))||((i == 2)&&(j == 1)))
	    {		
		   Vintra = Vconfopen(v, d, v, v) * spctrl_wghttopplus(v, d, kf1) * spctrl_wghttopminus(v, d, kf2);
	    }
	    else if(((i == 1)&&(j == 3))||((i == 3)&&(j == 1)))
	    {
		   Vintra = Vconfopen(v, d, v, v) * spctrl_wghttopplus(v, d, kf1) * spctrl_wghtbottomplus(v, d, kf3);
	    }
	    else if(((i == 2)&&(j == 3))||((i == 3)&&(j == 2)))
	    {
		   Vintra = Vconfopen(v, d, v, v) * spctrl_wghttopminus(v, d, kf2) * spctrl_wghtbottomplus(v, d, kf3);
	    }
	    if((i == 2)&&(j == 2))
	    {
		   Vintra = Vconfopen(v, d, v, v) * spctrl_wghttopminus(v, d, kf2) * spctrl_wghttopminus(v, d, kf2);
	    }
	    if((i == 3)&&(j == 3))
	    {
		   Vintra = Vconfopen(v, d, v, v) * spctrl_wghtbottomplus(v, d, kf3) * spctrl_wghtbottomplus(v, d, kf3);
	    }
		
	}

	return Vintra;
}

double zplus(int v, double d, double k)
{
	double zplus;
	
	zplus = 0.5*(spctrl_wghttopplus(v, d, k) + spctrl_wghttopminus(v, d, k));
	
	return zplus;
}	

double zminus(int v, double d, double k)
{
	double zminus;
	
	zminus = 0.5*(spctrl_wghtbottomplus(v, d, k) + spctrl_wghtbottomminus(v, d, k));
	
	return zminus;
}	




	

/*integration limit :*/

double uplimit(double v, double d, double (*func)(int v, double d, double x), double start_val)
{
	
	double k, kup;
	
	k = start_val;
	
	while(k < kFbulk)
	{
		if(fabs(func(v, d, k) - 1.) <= 0.01)
		{
			kup = k;
			break;
		}
		k = k + 0.00001;
	}
	
	return kup;
}	
		
double lowlimit(double v, double d, double (*func)(int v, double d, double x), double start_val)
{
	
	double k, kdown = 0.;
	
	k = start_val;
	
	if( k > 0.001)
	{
		while(k > 0.)
	{
	   	if(fabs(func(v, d, k) + 1.) <= 0.01)
		{
			kdown = k;
			break;
		}
		k = k - 0.00001;
	}
	
	}
	else
	{
		while(k < kFbulk)
	    {
	     	if(fabs(func(v, d, k) + 1.) <= 0.01)
	        	{
		    	kdown = k;
		    	break;
	    	    }
		k = k + 0.00001;
	    }
	}


		
	return kdown;
}	
			
