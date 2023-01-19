//
//  Metals.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/11/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "Metals.h"

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

@implementation FeSolid

-(id)init {
	if ((self = [super init])) {
		tr =     298.15;
		pr =       1.0;
		k0 =      23.9788;
		k1 =      -0.12875;
		k2 =     662.0648;
		k3 = -111298.07;
		ct =       0.00836712;

		v1 = -0.594e-6;
		v2 =  0.0;
		v3 = 75.8e-06;
		v4 =  0.0;

		h0 = 7788.0;
		s0 =   35.545;
		v0 =    0.709258*0.964596;

	}
	return self;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double h = h0 + k0*(t-tr) + 0.5*ct*(t*t-tr*tr) + 2.0*k1*(sqrt(t)-sqrt(tr)) - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
	double s = s0 + k0*log(t/tr) - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) + ct*(t-tr) - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr)) - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
	return h - t*s + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return h0 + k0*(t-tr) + 0.5*ct*(t*t-tr*tr) + 2.0*k1*(sqrt(t)-sqrt(tr)) - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr))
	+ v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr)) - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
}
-(double)getEntropyFromT:(double)t andP:(double)p {
	return s0 + k0*log(t/tr) - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) + ct*(t-tr) - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr)) - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr))
	-v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
}
-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return k0 + k1/sqrt(t) + ct*t + k2/SQUARE(t) + k3/CUBE(t) -t*v0*2.0*v4*(p-pr);
}
-(double)getDcpDtFromT:(double)t andP:(double)p {
	return - 0.5*k1/pow(t, (double) 1.5) + ct - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t) -v0*2.0*v4*(p-pr);
}
-(double)getVolumeFromT:(double)t andP:(double)p {
	return v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3 + (t-tr)*(t-tr)*v4);
}
-(double)getDvDtFromT:(double)t andP:(double)p {
	return v0*(v3 + 2.0*(t-tr)*v4);
}
-(double)getDvDpFromT:(double)t andP:(double)p {
	return v0*(v1 + 2.0*(p-pr)*v2);
}
-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return v0*2.0*v4;
}
-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return 0.0;
}
-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return v0*2.0*v2;
}

@end

@implementation FeLiquid

-(id)init {
	if ((self = [super init])) {
		tr  =  298.15;
		tr1 = 1200.0;
		pr  =    1.0;

		h0 = 12395.0+32028.0;
		s0 = 80.932;

		v0 =  0.652697;
		v1 = -8.87646e-7;
		v2 =  1.306567e-12;
		v3 =  1.31136e-4;
		v4 =  7.2557e-9;

		cp = 46.024;
	}
	return self;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double h = h0 + cp*(t-tr1);
	double s = s0 + cp*log(t/tr1);
	return h - t*s + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
}
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return h0 + cp*(t-tr1) + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr)) - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
}
-(double)getEntropyFromT:(double)t andP:(double)p {
	return s0 + cp*log(t/tr1) -v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
}
-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return cp -t*v0*2.0*v4*(p-pr);
}
-(double)getDcpDtFromT:(double)t andP:(double)p {
	return -v0*2.0*v4*(p-pr);
}
-(double)getVolumeFromT:(double)t andP:(double)p {
	return v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3 + (t-tr)*(t-tr)*v4);
}
-(double)getDvDtFromT:(double)t andP:(double)p {
	return v0*(v3 + 2.0*(t-tr)*v4);
}
-(double)getDvDpFromT:(double)t andP:(double)p {
	return v0*(v1 + 2.0*(p-pr)*v2);
}
-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return v0*2.0*v4;
}
-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return 0.0;
}
-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return v0*2.0*v2;
}

@end

@implementation NiSolid

-(id)init {
	if ((self = [super init])) {
		tr1   = 298.15;
		tr    = 750.0;
		pr    =   1.0;

		h0    = 3.32124*4184.0;
		s0    = 0.013795*4184.0;

		v0    =  0.659;
		vni1  = -0.536e-6;
		vni3a = 47.900e-6;
		vni3b = 53.900e-6;

		cpa   =     5.071e-3*4184.0;
		cpb   = 2.0*1.157e-6*4184.0;
		cpc   =    -3.136e-2*4184.0;
	}
	return self;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double h = h0 + cpa*(t-tr) + 0.5*cpb*(t*t-tr*tr) + cpc*(1.0/t-1.0/tr);
	double s = s0 + cpa*(log(t)-log(tr)) + cpb*(t-tr) + 0.5*cpc*(1/t/t-1/tr/tr);
	return h - t*s + v0*(1.0+vni3a*(631.0-tr1))*((vni1/2.0)*(p*p-pr*pr) + (1.0-vni1+vni3b*(t-631.0))*(p-pr));
}
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return h0 + cpa*(t-tr) + 0.5*cpb*(t*t-tr*tr) + cpc*(1.0/t-1.0/tr)
	+ v0*(1.0+vni3a*(631.0-tr1))*(0.5*vni1*(p*p-pr*pr) + (1.0-vni1+vni3b*(t-631.0))*(p-pr)) - t*v0*(1.0+vni3a*(631.0-tr1))*vni3b*(p-pr);
}
-(double)getEntropyFromT:(double)t andP:(double)p {
	return s0 + cpa*(log(t)-log(tr)) + cpb*(t-tr) + 0.5*cpc*(1/t/t-1/tr/tr) -v0*(1.0+vni3a*(631.0-tr1))*vni3b*(p-pr);
}
-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return cpa + cpb*t - cpc/(t*t);
}
-(double)getDcpDtFromT:(double)t andP:(double)p {
	return cpb + 2.0*cpc/(t*t*t);
}
-(double)getVolumeFromT:(double)t andP:(double)p {
	return v0*(1.0+vni3a*(631.0-tr1))*(vni1*p + (1.0-vni1+vni3b*(t-631.0)));
}
-(double)getDvDtFromT:(double)t andP:(double)p {
	return v0*(1.0+vni3a*(631.0-tr1))*vni3b;
}
-(double)getDvDpFromT:(double)t andP:(double)p {
	return v0*(1.0+vni3a*(631.0-tr1))*vni1;
}
-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return 0.0;
}
-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return 0.0;
}
-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return 0.0;
}

@end

@implementation NiLiquid

-(id)init {
	if ((self = [super init])) {
		tr  =  298.15;
		tr1 = 1300.0;
		pr  =    1.0;

		h0 = 17479.0+30382.0;
		s0 = 84.680;

		v0 =  0.6286951;
		v1 = -8.87646e-7;
		v2 =  1.306567e-12;
		v3 =  1.2486e-4;
		v4 =  3.2325e-9;

		cp = 38.911;
	}
	return self;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double h = h0 + cp*(t-tr1);
	double s = s0 + cp*log(t/tr1);
	return h - t*s + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
}
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return h0 + cp*(t-tr1) + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr)) - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
}
-(double)getEntropyFromT:(double)t andP:(double)p {
	return s0 + cp*log(t/tr1) - v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
}
-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return cp - t*v0*2.0*v4*(p-pr);
}
-(double)getDcpDtFromT:(double)t andP:(double)p {
	return -v0*2.0*v4*(p-pr);
}
-(double)getVolumeFromT:(double)t andP:(double)p {
	return v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3 + (t-tr)*(t-tr)*v4);
}
-(double)getDvDtFromT:(double)t andP:(double)p {
	return v0*(v3 + 2.0*(t-tr)*v4);
}
-(double)getDvDpFromT:(double)t andP:(double)p {
	return v0*(v1 + 2.0*(p-pr)*v2);
}
-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return v0*2.0*v4;
}
-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return 0.0;
}
-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return v0*2.0*v2;
}

@end
