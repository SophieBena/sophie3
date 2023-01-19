//
//  LiquidpMeltsGenericEM.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/16/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "LiquidpMeltsGenericEM.h"

static const unsigned int returnValueOfG       =  1;
static const unsigned int returnValueOfH       =  2;
static const unsigned int returnValueOfS       =  3;
static const unsigned int returnValueOfCP      =  4;
static const unsigned int returnValueOfDCPDT   =  5;
static const unsigned int returnValueOfV       =  6;
static const unsigned int returnValueOfdVdT    =  7;
static const unsigned int returnValueOfdVdP    =  8;
static const unsigned int returnValueOfd2VdT2  =  9;
static const unsigned int returnValueOfd2VdTdP = 10;
static const unsigned int returnValueOfd2VdP2  = 11;

@implementation LiquidpMeltsGenericEM

-(id)initWithH:(double)hIn
			 S:(double)sIn
			k0:(double)k0In
			k1:(double)k1In
			k2:(double)k2In
			k3:(double)k3In
			l1:(double)l1In
			l2:(double)l2In
			Tt:(double)TtIn
		deltaH:(double)deltaHIn
		  vLiq:(double)vLiqIn
	   dvdtLiq:(double)dvdtLiqIn
	   dvdpLiq:(double)dvdpLiqIn
	d2vdtdpLiq:(double)d2vdtdpLiqIn
	 d2vdp2Liq:(double)d2vdp2LiqIn
	   tFusion:(double)tFusionIn
	   sFusion:(double)sFusionIn
		 cpLiq:(double)cpLiqIn {
	if ((self = [super initWithH:hIn S:sIn k0:k0In k1:k1In k2:k2In k3:k3In l1:l1In l2:l2In Tt:TtIn deltaH:deltaHIn v0:0.0 v1:0.0 v2:0.0 v3:0.0 v4:0.0])) {
		vLiq = vLiqIn;
		dvdtLiq = dvdtLiqIn;
		dvdpLiq = dvdpLiqIn;
		d2vdtdpLiq = d2vdtdpLiqIn;
		d2vdp2Liq = d2vdp2LiqIn;
		tFusion = tFusionIn;
		sFusion = sFusionIn;
		cpLiq = cpLiqIn;

		hLiqAtTfusion = [super getEnthalpyFromT:tFusion andP:1.0];
		sLiqAtTfusion = [super getEntropyFromT:tFusion andP:1.0];
		gLiqAtTfusion = hLiqAtTfusion - tFusionIn*sLiqAtTfusion;
		hLiqAtTfusion += sFusion*tFusion;
		sLiqAtTfusion += sFusion;

		trl = 1673.15;

	}
	return self;
}

-(double)calculateWithT:(double)t andWithP:(double)p forCase:(NSUInteger)returnMode {
	if (p > pr && (dvdpLiq + d2vdtdpLiq*(t-trl)) != 0.0) {
		double d2v0BMdtdp = d2vdtdpLiq;
		double v0BM	      = vLiq + dvdtLiq*(t-trl);
		double dv0BMdt    = dvdtLiq;
		double dv0BMdp    = dvdpLiq + d2v0BMdtdp*(t-trl);
		double K	      = -v0BM/dv0BMdp;
		double dKdt       =      d2v0BMdtdp*v0BM/(dv0BMdp*dv0BMdp)                                       -     dv0BMdt/dv0BMdp;
		double d2Kdt2     = -2.0*d2v0BMdtdp*d2v0BMdtdp*v0BM/(dv0BMdp*dv0BMdp*dv0BMdp)                    + 2.0*d2v0BMdtdp*dv0BMdt/(dv0BMdp*dv0BMdp);
		double d3Kdt3     =  6.0*d2v0BMdtdp*d2v0BMdtdp*d2v0BMdtdp*v0BM/(dv0BMdp*dv0BMdp*dv0BMdp*dv0BMdp) - 6.0*d2v0BMdtdp*d2v0BMdtdp*dv0BMdt/(dv0BMdp*dv0BMdp*dv0BMdp);
		double Kp	      = 5.0;
		int iter          = 0;
		const int iterMax = 1000;
		double v, vLast, fn = 0.0 , dfn;
		double f, dfdt, d2fdt2, d3fdt3, A, B, C, D, E, F, dDdt, dEdt, dFdt, d3vldt3, d2gdt2, d3gdt3;
		double minusIntPdV, d_minusIntPdV_df, d2_minusIntPdV_df2, d3_minusIntPdV_df3;

		/* Newton's method for volume from BM */
		vLast  = v0BM;
		v      = v0BM*0.99;
		while ((fabs(v - vLast) > 100.0*DBL_EPSILON) && (iter < iterMax)) {
			fn = (3.0/2.0)*K*(pow(v0BM/v, (double) 7.0/3.0) - pow(v0BM/v, (double) 5.0/3.0))
			*(1.0-(3.0/4.0)*(4.0-Kp)*(pow(v0BM/v, (double) 2.0/3.0) - 1.0))
			- p;
			dfn = (3.0/2.0)*K*((7.0/3.0)*pow(v0BM/v, (double) 4.0/3.0) - (5.0/3.0)*pow(v0BM/v, (double) 2.0/3.0))*(-v0BM/(v*v))
			*(1.0-(3.0/4.0)*(4.0-Kp)*(pow(v0BM/v, (double) 2.0/3.0) - 1.0))
			+ (3.0/2.0)*K*(pow(v0BM/v, (double) 7.0/3.0) - pow(v0BM/v, (double) 5.0/3.0))
			*(-(3.0/4.0)*(4.0-Kp)*(2.0/3.0)*pow(v0BM/v, (double) -1.0/3.0))*(-v0BM/(v*v));
			vLast = v;
			v += -fn/dfn;
			if (v > v0BM     ) v = v0BM;
			if (v < 0.01*v0BM) v = 0.01*v0BM;
			iter++;
		}
		if (iter >= iterMax) {
			NSLog(@"Convergence error in Birch-Murnaghan Volume routine in LiquidpMeltsGenericEM.");
			NSLog(@"  For %@,  v = %g, dV = %g, f = %g at T = %g and P = %g in %d iterations.", [self phaseName], v, fabs(v - vLast), fn, t, p, iter);
		}

		f = (pow(v0BM/v, (double) 2.0/3.0) - 1.0)/2.0;
		A = 1.0/f + 5.0/(1.0+2.0*f) + (3.0/2.0)*(Kp-4.0)/(1.0+(3.0/2.0)*f*(Kp-4.0));
		B = 1.0/(f*f) + 10.0/((1.0+2.0*f)*(1.0+2.0*f)) + pow((3.0/2.0)*(Kp-4.0)/(1.0+(3.0/2.0)*f*(Kp-4.0)), (double) 2.0);
		C = -2.0/(f*f*f) - 40.0/((1.0+2.0*f)*(1.0+2.0*f)*(1.0+2.0*f)) - 2.0*pow((3.0/2.0)*(Kp-4.0)/(1.0+(3.0/2.0)*f*(Kp-4.0)), (double) 3.0);

		D    = dKdt*dv0BMdt/K;
		dDdt = dv0BMdt*(d2Kdt2/K - dKdt*dKdt/(K*K));
		E    = (dv0BMdt*dKdt + v0BM*d2Kdt2)/K;
		dEdt = 2.0*dv0BMdt*d2Kdt2/K + v0BM*d3Kdt3/K - dv0BMdt*dKdt*dKdt/(K*K) - v0BM*dKdt*d2Kdt2/(K*K);
		F    = - v0BM*dKdt*dKdt*(A - 5.0/(1.0+2.0*f) + B/A)/(A*K*K);
		dFdt = (-dv0BMdt*dKdt*dKdt - 2.0*v0BM*dKdt*d2Kdt2 + 2.0*v0BM*dKdt*dKdt*dKdt/K + v0BM*dKdt*dKdt*dKdt*B/(A*A*K))
		*(A - 5.0/(1.0+2.0*f) + B/A)/(A*K*K)
		- v0BM*dKdt*dKdt*dKdt*(B - 10.0/((1.0+2.0*f)*(1.0+2.0*f)) + C/A - B*B/(A*A))/(A*A*K*K*K);

		minusIntPdV        = (9.0/2.0)*K*(Kp-4.0)*v0BM*f*f*f + (9.0/2.0)*K*v0BM*f*f;
		d_minusIntPdV_df   = (27.0/2.0)*K*(Kp-4.0)*v0BM*f*f + 9.0*K*v0BM*f;
		d2_minusIntPdV_df2 = 27.0*K*(Kp-4.0)*v0BM*f + 9.0*K*v0BM;
		d3_minusIntPdV_df3 = 27.0*K*(Kp-4.0)*v0BM;

		dfdt = -dKdt/(A*K);
		d2fdt2 = -d2Kdt2/(A*K) + dKdt*dKdt/(A*K*K) + dKdt*dKdt*B/(A*A*A*K*K);
		d3fdt3 = -d3Kdt3/(A*K) + 3.0*dKdt*d2Kdt2*(B/(A*A)+1.0)/(A*K*K)
		+ dKdt*dKdt*dKdt*(C/(A*A*A)-3.0*B/(A*A)-2.0-3.0*B*B/(A*A*A*A))/(A*K*K*K);

		d2gdt2  = p*((3.0/(A*pow(1.0+2.0*f, (double) 5.0/2.0)))*(dv0BMdt*dKdt/K + (d2Kdt2*v0BM+dKdt*dv0BMdt)/K - v0BM*dKdt*dKdt*(A-5.0/(1.0+2.0*f)+B/A)/(A*K*K)))
		+ 2.0*(minusIntPdV/(v0BM*K))*dKdt*dv0BMdt + 2.0*(d_minusIntPdV_df/v0BM)*dfdt*dv0BMdt
		+ 2.0*(d_minusIntPdV_df/K)*dfdt*dKdt + d2_minusIntPdV_df2*dfdt*dfdt
		+ (minusIntPdV/K)*d2Kdt2 + d_minusIntPdV_df*d2fdt2;
		d3vldt3 = 3.0*(dDdt + dEdt + dFdt)/(A*pow(1.0+2.0*f, (double) 5.0/2.0))
		- 3.0*(D + E + F)*dKdt*(B/A - 5.0/(1.0+2.0*f))/(K*A*A*pow(1.0+2.0*f, (double) 5.0/2.0));
		d3gdt3  = p*d3vldt3 + 3.0*d2_minusIntPdV_df2*dfdt*dfdt*dKdt/K + 6.0*d_minusIntPdV_df*dv0BMdt*dfdt*dKdt/(K*v0BM)
		+ 3.0*d_minusIntPdV_df*d2fdt2*dKdt/K + d3_minusIntPdV_df3*dfdt*dfdt*dfdt + 3.0*d2_minusIntPdV_df2*dfdt*dfdt*dv0BMdt/v0BM
		+ 3.0*d_minusIntPdV_df*dfdt*d2Kdt2/K + 3.0*d2_minusIntPdV_df2*dfdt*d2fdt2 + 3.0*minusIntPdV*dv0BMdt*d2Kdt2/(v0BM*K)
		+ 3.0*d_minusIntPdV_df*dv0BMdt*d2fdt2/v0BM + minusIntPdV*d3Kdt3/K + d_minusIntPdV_df*d3fdt3;

		switch (returnMode) {
			case 1:
				return p*v - pr*v0BM + minusIntPdV;
				break;
			case 2:
				return p*(v+3.0*v0BM*dKdt/(A*pow(1.0+2.0*f, (double) 5.0/2.0)*K) + dv0BMdt/pow(1.0+2.0*f, (double) 3.0/2.0)) - pr*(v0BM+dv0BMdt)
				       + minusIntPdV + (minusIntPdV/v0BM)*dv0BMdt + (minusIntPdV/K)*dKdt + d_minusIntPdV_df*dfdt;
				break;
			case 3:
				return -(p*(3.0*v0BM*dKdt/(A*pow(1.0+2.0*f, (double) 5.0/2.0)*K) + dv0BMdt/pow(1.0+2.0*f, (double) 3.0/2.0))
						 - pr*dv0BMdt + (minusIntPdV/v0BM)*dv0BMdt + (minusIntPdV/K)*dKdt + d_minusIntPdV_df*dfdt);
				break;
			case 4:
				return -t*d2gdt2;
				break;
			case 5:
				return -t*d3gdt3 - d2gdt2;
				break;
			case 6:
				return v;
				break;
			case 7:
				return 3.0*v0BM*dKdt/(A*pow(1.0+2.0*f, (double) 5.0/2.0)*K) + dv0BMdt/pow(1.0+2.0*f, (double) 3.0/2.0);
				break;
			case 8:
				return -3.0*v0BM/(p*pow(1.0+2.0*f, (double) 5.0/2.0)*A);
				break;
			case 9:
				return (3.0/(A*pow(1.0+2.0*f, (double) 5.0/2.0)))
				       *(dv0BMdt*dKdt/K + (d2Kdt2*v0BM+dKdt*dv0BMdt)/K - v0BM*dKdt*dKdt*(A-5.0/(1.0+2.0*f)+B/A)/(A*K*K));
				break;
			case 10:
				return (3.0*v0BM/(p*p*pow(1.0+2.0*f, (double) 5.0/2.0)*A*A))*(A-B/A+5.0/(1.0+2.0*f));
				break;
			case 11:
				return (3.0/(p*pow(1.0+2.0*f, (double) 5.0/2.0)*A))*(-dv0BMdt + v0BM*dKdt*(B/A-5.0/(1.0+2.0*f))/(K*A));
				break;
			default:
				return 0.0;
				break;
		}

	} else {
		switch (returnMode) {
			case 6:
				return vLiq;
				break;
			default:
				return 0.0;
				break;
		}
	}
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return hLiqAtTfusion + cpLiq*(t-tFusion) - t*(sLiqAtTfusion + cpLiq*log(t/tFusion))
		   + [self calculateWithT:t andWithP:p forCase:returnValueOfG];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return hLiqAtTfusion + cpLiq*(t-tFusion) + [self calculateWithT:t andWithP:p forCase:returnValueOfH];
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	return sLiqAtTfusion + cpLiq*log(t/tFusion) + [self calculateWithT:t andWithP:p forCase:returnValueOfS];
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return cpLiq + [self calculateWithT:t andWithP:p forCase:returnValueOfCP];
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfDCPDT];
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfV];
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfdVdT];
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfdVdP];
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdT2];
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdTdP];
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdP2];
}

@end
