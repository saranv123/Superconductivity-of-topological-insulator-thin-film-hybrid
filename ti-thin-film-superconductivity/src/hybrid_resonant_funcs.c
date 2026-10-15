#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<stdlib.h>
#include "thinfilmSC_funcs.c"



double detun(int v, double d)
{
	
	double detune;
	
	detune = (mubulk - ((mubulk*pi*pi)/(kFbulk*kFbulk))*(pow(v,2)/pow(d,2)))/((2*t0)/sqrt(d));
	
	return detune;
	
}


double drcpnteff(int v, double d)
{
	double epsilon;
	
	epsilon = (drcpnt/(t0/sqrt(d)) + detun(v, d) - sqrt(1. + detun(v, d)*detun(v, d)))*(t0/sqrt(d));
	
	return epsilon;
}	

double effspinorbit(int v, double d)
{
	double spinorbit;
	
	spinorbit = A0*0.5*(1 + (detun(v, d)/sqrt(1. + detun(v, d)*detun(v, d))));
	
	return spinorbit;
}	

double spctrl_wght(int v, double d)
{
	double spctrlwght;
	
	spctrlwght = 0.5*(1. - (detun(v, d)/sqrt(1. + detun(v, d)*detun(v, d))));
	
	return spctrlwght;
}	

double spctrl_wght2(int v, double d)
{
	double spctrlwght;
	
	spctrlwght = 0.5*(1. + (detun(v, d)/sqrt(1. + detun(v, d)*detun(v, d))));
	
	return spctrlwght;
}	


double dnstyofstates(int v, double d)
{
	double N_d;
	
	N_d = fabs(drcpnteff(v, d))/(4.*pi*effspinorbit(v, d)*effspinorbit(v, d));
	
	return N_d;
}	



