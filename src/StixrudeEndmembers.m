//
//  StixrudeEndmembers.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/9/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "StixrudeEndmembers.h"
#import "DoubleMatrix.h"

@implementation StixrudeEndmembers

+(double)chebvalat:(double)x {
    double c[17] = {
        2.707737068327440945 / 2.0, 0.340068135211091751, -0.12945150184440869e-01, 0.7963755380173816e-03,
        -0.546360009590824e-04, 0.39243019598805e-05, -0.2894032823539e-06, 0.217317613962e-07, -0.16542099950e-08,
        0.1272796189e-09, -0.987963460e-11, 0.7725074e-12, -0.607797e-13, 0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17
    };
    double x2 = 2 * x;
    double c0 = c[17-2];
    double c1 = c[17-1];
    for (NSUInteger i=3; i<18; i++) {
        double tmp = c0;
        c0 = c[17-i] - c1;
        c1 = tmp + c1 * x2;
    }
    return c0 + c1 * x;
}

+(double)debyeIntegral:(double)x {
    double val_infinity = 19.4818182068004875;
    double sqrt_eps = sqrt(DBL_EPSILON);
    double log_eps = log(DBL_EPSILON);
    double xcut = -log_eps;

    if (x <= 0.0) return 0.0;

    if (x < (2.0*sqrt(2.0)*sqrt_eps)) return 1.0 - 3.0*x/8.0 + x*x/20.0;
    else if (x <= 4.0) {
        double t = x*x/8.0 - 1.0;
        double c = [StixrudeEndmembers chebvalat:t];
        return c - 0.375*x;
    } else if (x < -(log(2.0)+log_eps)) {
        NSInteger nexp = (NSInteger)(floor(xcut / x));
        double ex = exp(-x);
        double xk = nexp * x;
        double rk = nexp;
        double sum = 0.0;
        for (NSInteger i=nexp; i>0; i--) {
            double xk_inv = 1.0/xk;
            sum *= ex;
            sum += (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0)/rk;
            rk -= 1.0;
            xk -= x;
        }
        return val_infinity / (x * x * x) - 3.0 * sum * ex;
    } else if (x < xcut) {
        double x3 = x*x*x;
        double sum = 6.0 + 6.0*x + 3.0*x*x + x3;
        return (val_infinity - 3.0*sum*exp(-x))/x3;
    } else return ((val_infinity/x)/x)/x;
}

+(double)debyeThermalEnergyWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double E_th = 3.0*n*r*t*[StixrudeEndmembers debyeIntegral:x];
    return E_th;
}

+(double)debyeHelmholtzWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double F = n*r*t*(3.0*log(1.0 - exp(-x)) - [StixrudeEndmembers debyeIntegral:x]);
    return F;
}

+(double)debyeEntropyWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double S = n*r*(4.0*[StixrudeEndmembers debyeIntegral:x] - 3.0*log(1.0 - exp(-x)));
    return S;
}

+(double)debyeDaDvWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n andDdebyeDv:(double)DthetaTdV {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double dFdV = 3.0*n*r*DthetaTdV*[StixrudeEndmembers debyeIntegral:x]/x;
    return dFdV;
}

+(double)debyeHeatCapacityWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double C_v = 3.0*n*r*(4.0*[StixrudeEndmembers debyeIntegral:x] - 3.0*x/(exp(x) - 1.0));
    return C_v;
}

+(double)debyeD2aDv2WithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n andDdebyeDv:(double)DthetaTdV andD2debyeDv2:(double)D2thetaTdV2 {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double d2FdV2 = -12.0*n*r*DthetaTdV*DthetaTdV*[StixrudeEndmembers debyeIntegral:x]/(t*x*x);
           d2FdV2 += 9.0*n*r*DthetaTdV*DthetaTdV/(t*x*(exp(x)-1.0));
           d2FdV2 += 3.0*n*r*[StixrudeEndmembers debyeIntegral:x]*D2thetaTdV2/x;
    return d2FdV2;
}

+(double)debyeD2aDtDvWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n andDdebyeDv:(double)DthetaTdV {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double r = 8.314472;
    double d2FdTdV = 12.0*n*r*DthetaTdV*[StixrudeEndmembers debyeIntegral:x]/(t*x);
           d2FdTdV += -9.0*n*r*DthetaTdV/(t*(exp(x)-1.0));
    return d2FdTdV;
}

+(double)debyeD3aDv3WithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n andDdebyeDv:(double)DthetaTdV andD2debyeDv2:(double)D2thetaTdV2 andD3debyeDv3:(double)D3thetaTdV3 {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    //double DxDt = -x/t;
    double DxDv = DthetaTdV/t;
    double r = 8.314472;
    double D3_x = [StixrudeEndmembers debyeIntegral:x];
    double dD3_xDx = 3.0*(1.0/(exp(x)-1.0) - D3_x/x);

    //double d2FdV2 = -12.0*n*r*DthetaTdV*DthetaTdV*D3_x/(t*x*x);
    double d3FdV3 = -12.0*n*r*2.0*DthetaTdV*D2thetaTdV2*D3_x/(t*x*x) - 12.0*n*r*DthetaTdV*DthetaTdV*dD3_xDx*DxDv/(t*x*x) + 2.0*12.0*n*r*DthetaTdV*DthetaTdV*D3_x*DxDv/(t*x*x*x);
    //d2FdV2 += 9.0*n*r*DthetaTdV*DthetaTdV/(t*x*(exp(x)-1.0));
    d3FdV3 += 9.0*n*r*2.0*DthetaTdV*D2thetaTdV2/(t*x*(exp(x)-1.0)) - 9.0*n*r*DthetaTdV*DthetaTdV*(exp(x) - 1.0 + x*exp(x))*DxDv/(t*pow(x*(exp(x)-1.0), 2.0));
    //d2FdV2 += 3.0*n*r*[StixrudeEndmembers debyeIntegral:x]*D2thetaTdV2/x;
    d3FdV3 += 3.0*n*r*dD3_xDx*DxDv*D2thetaTdV2/x + 3.0*n*r*D3_x*D3thetaTdV3/x - 3.0*n*r*D3_x*D2thetaTdV2*DxDv/(x*x);

    return d3FdV3;
}

+(double)debyeD3aDv2DtWithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n andDdebyeDv:(double)DthetaTdV andD2debyeDv2:(double)D2thetaTdV2 {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double DxDt = -x/t;
    double r = 8.314472;
    double D3_x = [StixrudeEndmembers debyeIntegral:x];
    double dD3_xDx = 3.0*(1.0/(exp(x)-1.0) - D3_x/x);

    //double d2FdV2 = -12.0*n*r*DthetaTdV*DthetaTdV*D3_x/(t*x*x);
    double d3FdV2dT = -12.0*n*r*DthetaTdV*DthetaTdV*dD3_xDx*DxDt/(t*x*x) + 12.0*n*r*DthetaTdV*DthetaTdV*D3_x/(t*t*x*x) + 2.0*12.0*n*r*DthetaTdV*DthetaTdV*D3_x*DxDt/(t*x*x*x);
    //d2FdV2 += 9.0*n*r*DthetaTdV*DthetaTdV/(t*x*(exp(x)-1.0));
    d3FdV2dT += -9.0*n*r*DthetaTdV*DthetaTdV/(t*t*x*(exp(x)-1.0)) - 9.0*n*r*DthetaTdV*DthetaTdV*DxDt/(t*x*x*(exp(x)-1.0)) - 9.0*n*r*DthetaTdV*DthetaTdV*exp(x)*DxDt/(t*x*pow(exp(x)-1.0, 2.0));
    //d2FdV2 += 3.0*n*r*[StixrudeEndmembers debyeIntegral:x]*D2thetaTdV2/x;
    d3FdV2dT += 3.0*n*r*dD3_xDx*DxDt*D2thetaTdV2/x - 3.0*n*r*D3_x*D2thetaTdV2*DxDt/x/x;
    return d3FdV2dT;
}

+(double)debyeD3aDvDt2WithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n andDdebyeDv:(double)DthetaTdV {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double DxDt = -x/t;
    double r = 8.314472;
    double D3_x = [StixrudeEndmembers debyeIntegral:x];
    double dD3_xDx = 3.0*(1.0/(exp(x)-1.0) - D3_x/x);
    //double d2FdTdV = 12.0*n*r*DthetaTdV*D3_x/(t*x);
    double d3FdVdT2 = 12.0*n*r*DthetaTdV*dD3_xDx*DxDt/(t*x) - 12.0*n*r*DthetaTdV*D3_x/(t*t*x) - 12.0*n*r*DthetaTdV*D3_x*DxDt/(t*x*x);
    //d2FdTdV += -9.0*n*r*DthetaTdV/(t*(exp(x)-1.0));
    d3FdVdT2 += 9.0*n*r*DthetaTdV/(t*t*(exp(x)-1.0)) + 9.0*n*r*DthetaTdV*exp(x)*DxDt/(t*pow(exp(x)-1.0, 2.0));
    return d3FdVdT2;
}

+(double)debyeD3aDt3WithT:(double)t andDebyeT:(double)thetaT andN:(NSInteger)n {
    if (t <= DBL_EPSILON) return 0.0;
    double x = thetaT/t;
    double DxDt = -x/t;
    double r = 8.314472;
    double D3_x = [StixrudeEndmembers debyeIntegral:x];
    double dD3_xDx = 3.0*(1.0/(exp(x)-1.0) - D3_x/x);
    double d3FdT3 = -3.0*n*r*(4.0*dD3_xDx*DxDt - 3.0*DxDt/(exp(x)-1.0) + 3.0*x*exp(x)*DxDt/pow(exp(x)-1.0, 2.0))/t
                  + 3.0*n*r*(4.0*D3_x - 3.0*x/(exp(x) - 1.0))/t/t;
    return d3FdT3;
}

#define FUNCTION 0
#define GRADIENT 1
#define HESSIAN  2
#define TENSOR   3

#define T   0
#define V   1
#define TT  2
#define TV  3
#define VV  4
#define VVV 5
#define VVT 6
#define VTT 7
#define TTT 8

+(double)helmholtz:(int)type returnVariable:(int)variable temperature:(double)t volume:(double)v params:(double *)params {
	double tr, r, a0, n, v0, k00, k0p, theta0, gamma0, q, refS;
	double f, c1, c2, c5, c7, tht, tht0, d0;

	// constants
	tr    = 300.00;
	r     = 8.314472;

	a0	   = params[0]*1000.0;  //  1 J/m
	n	   = params[1];         //  2
	v0	   = params[2]/10.0;    //  3 J/bar-m
	k00	   = params[3]*10000.0; //  4 bar
	k0p	   = params[4];         //  5
	theta0 = params[5];         //  6 K
	gamma0 = params[6];         //  7
	q	   = params[7];         //  8
	refS   = params[8];         // 10

	c1   = 9.0*k00*v0;
	c2   = k0p/2.0 - 2.0;
	c5   = 3.0*gamma0;
	c7   = c5*(-2.0 + 6.0*gamma0 - 3.0*q);

	f    = 0.5*pow(v0/v, 2.0/3.0) - 0.5;
	double arg = 1.0 + c7*f*f + 2.0*c5*f;

	if (arg < 0.0) {
		NSException *stixrudeException = [NSException exceptionWithName:@"Stixrude Internal Error"
																 reason:@"Argument to square root function is negative."
															   userInfo:[NSDictionary dictionaryWithObjectsAndKeys:
																		 [NSNumber numberWithDouble:arg], @"Argument",
																		 [NSNumber numberWithDouble:f], @"f",
																		 [NSNumber numberWithDouble:c7], @"c7",
																		 [NSNumber numberWithDouble:c5], @"c5",
																		 nil]
										  ];
		@throw stixrudeException;
	}

	d0   = theta0*sqrt(arg);
	tht  = theta0/t*d0;
	tht0 = tht*t/tr;

	if (type == FUNCTION) {
        double F_quasiharmonic = [StixrudeEndmembers debyeHelmholtzWithT:t andDebyeT:d0 andN:(NSInteger)n]
                               - [StixrudeEndmembers debyeHelmholtzWithT:tr andDebyeT:d0 andN:(NSInteger)n];
        return a0 + c1*f*f*(0.5 + c2*f) + F_quasiharmonic;
	} else if ((type == GRADIENT) && (variable == T)) {
        return -[StixrudeEndmembers debyeEntropyWithT:t andDebyeT:d0 andN:(NSInteger)n];
	} else if ((type == GRADIENT) && (variable == V)) {
		double dfdv    = 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-v0/v/v);
		double dd0dv   = theta0*0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)/sqrt(arg);
        return 2.0*c1*f*dfdv*(0.5 + c2*f) + c1*f*f*c2*dfdv
               + [StixrudeEndmembers debyeDaDvWithT:t andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv]
               - [StixrudeEndmembers debyeDaDvWithT:tr andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv];
	} else if ((type == HESSIAN) && (variable == TT)) {
        return -[StixrudeEndmembers debyeHeatCapacityWithT:t andDebyeT:d0 andN:(NSInteger)n]/t;
	} else if ((type == HESSIAN) && (variable == VV)) {
		double dfdv      = 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-v0/v/v);
		double d2fdv2    = 0.5*(2.0/3.0)*(1.0/3.0)*pow(v0/v, 2.0/3.0)*(-1.0/v/v) + 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(2.0*v0/v/v/v);
		double dd0dv     = theta0*0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)/sqrt(arg);
		double d2d0dv2   = theta0*0.5*(2.0*c7*dfdv*dfdv + 2.0*c7*f*d2fdv2 + 2.0*c5*d2fdv2)/sqrt(arg)
                         - 0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)*dd0dv/arg;
        return 2.0*c1*dfdv*dfdv*(0.5 + c2*f) + 2.0*c1*f*d2fdv2*(0.5 + c2*f) + 2.0*c1*f*dfdv*c2*dfdv + 2.0*c1*f*dfdv*c2*dfdv + c1*f*f*c2*d2fdv2
               + [StixrudeEndmembers debyeD2aDv2WithT:t andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv andD2debyeDv2:d2d0dv2]
               - [StixrudeEndmembers debyeD2aDv2WithT:tr andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv andD2debyeDv2:d2d0dv2];
	} else if ((type == HESSIAN) && (variable == TV)) {
		double dfdv       = 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-v0/v/v);
        double dd0dv     = theta0*0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)/sqrt(arg);
        return [StixrudeEndmembers debyeD2aDtDvWithT:t andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv];
    } else if ((type == TENSOR) && (variable == TTT)) {
        return [StixrudeEndmembers debyeD3aDt3WithT:t andDebyeT:d0 andN:(NSInteger)n];
    } else if ((type == TENSOR) && (variable == VVV)) {
        double dfdv      = 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-v0/v/v);
        double d2fdv2    = 0.5*(2.0/3.0)*(1.0/3.0)*pow(v0/v, 2.0/3.0)*(-1.0/v/v) + 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(2.0*v0/v/v/v);
        double d3fdv3    = 0.5*(2.0/3.0)*(1.0/3.0)*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-1.0/v/v)*(-v0/v/v) + 0.5*(2.0/3.0)*(1.0/3.0)*pow(v0/v, 2.0/3.0)*(2.0/v/v/v)
                         + 0.5*(2.0/3.0)*(1.0/3.0)*pow(v0/v, 2.0/3.0)*(2.0/v/v/v) + 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-6.0*v0/v/v/v/v);
        double dd0dv     = theta0*0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)/sqrt(arg);
        double d2d0dv2   = theta0*0.5*(2.0*c7*dfdv*dfdv + 2.0*c7*f*d2fdv2 + 2.0*c5*d2fdv2)/sqrt(arg)
                         - 0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)*dd0dv/arg;
        double DargDv   = 2.0*c7*f*dfdv + 2.0*c5*dfdv;
        double D2argDv2 = 2.0*c7*dfdv*dfdv + 2.0*c7*f*d2fdv2 + 2.0*c5*d2fdv2;
        double D3argDv3 = 2.0*c7*2.0*dfdv*d2fdv2 + 2.0*c7*dfdv*d2fdv2 + 2.0*c7*f*d3fdv3 + 2.0*c5*d3fdv3;
        double d3d0dv3   = theta0*0.5*(D3argDv3)/sqrt(arg) - theta0*0.5*(D2argDv2)*DargDv/2.0/pow(arg, 3.0/2.0)
                        - 0.5*(D2argDv2)*dd0dv/arg - 0.5*(DargDv)*d2d0dv2/arg + 0.5*(DargDv)*dd0dv*DargDv/pow(arg, 2.0);
        return 2.0*c1*2.0*dfdv*d2fdv2*(0.5 + c2*f) + 2.0*c1*dfdv*dfdv*c2*dfdv
               + 2.0*c1*dfdv*d2fdv2*(0.5 + c2*f) + 2.0*c1*f*d3fdv3*(0.5 + c2*f) + 2.0*c1*f*d2fdv2*c2*dfdv
               + 2.0*c1*dfdv*dfdv*c2*dfdv + 2.0*c1*f*2.0*dfdv*c2*d2fdv2
               + 2.0*c1*dfdv*dfdv*c2*dfdv + 2.0*c1*f*2.0*dfdv*c2*d2fdv2
               + c1*2.0*f*dfdv*c2*d2fdv2 + c1*f*f*c2*d3fdv3
               + [StixrudeEndmembers debyeD3aDv3WithT:t andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv andD2debyeDv2:d2d0dv2 andD3debyeDv3:d3d0dv3]
               - [StixrudeEndmembers debyeD3aDv3WithT:tr andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv andD2debyeDv2:d2d0dv2 andD3debyeDv3:d3d0dv3];
    } else if ((type == TENSOR) && (variable == VVT)) {
        double dfdv    = 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-v0/v/v);
        double d2fdv2  = 0.5*(2.0/3.0)*(1.0/3.0)*pow(v0/v, 2.0/3.0)*(-1.0/v/v) + 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(2.0*v0/v/v/v);
        double dd0dv   = theta0*0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)/sqrt(arg);
        double d2d0dv2 = theta0*0.5*(2.0*c7*dfdv*dfdv + 2.0*c7*f*d2fdv2 + 2.0*c5*d2fdv2)/sqrt(arg)
                       - 0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)*dd0dv/arg;
        return [StixrudeEndmembers debyeD3aDv2DtWithT:t andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv andD2debyeDv2:d2d0dv2];
    } else if ((type == TENSOR) && (variable == VTT)) {
        double dfdv       = 0.5*(2.0/3.0)*pow(v/v0, 1.0/3.0)*(-v0/v/v);
        double dd0dv     = theta0*0.5*(2.0*c7*f*dfdv + 2.0*c5*dfdv)/sqrt(arg);
        return [StixrudeEndmembers debyeD3aDvDt2WithT:t andDebyeT:d0 andN:(NSInteger)n andDdebyeDv:dd0dv];

    }
	return 0.0;
}

+(NSDictionary *)calculateThermodynamicPropertiesOfEndmemberAtT:(double)t andAtP:(double)p withParameters:(NSArray *)parameters {
	double params[9];

	params[0] = [[parameters objectAtIndex:0] doubleValue]; //  1 kJ/m    a0
	params[1] = [[parameters objectAtIndex:1] doubleValue]; //  2         n
	params[2] = [[parameters objectAtIndex:2] doubleValue]; //  3 cc/m    v0
	params[3] = [[parameters objectAtIndex:3] doubleValue]; //  4 GPa     k00
	params[4] = [[parameters objectAtIndex:4] doubleValue]; //  5         k0p
	params[5] = [[parameters objectAtIndex:5] doubleValue]; //  6 K       theta0
	params[6] = [[parameters objectAtIndex:6] doubleValue]; //  7         gamma0
	params[7] = [[parameters objectAtIndex:7] doubleValue]; //  8         q
	params[8] = [[parameters objectAtIndex:8] doubleValue]; // 10         refS
	double vRef = params[2]/10.0;
	double v    = vRef*0.90;
	double targetP = p;

	p = -[self helmholtz:GRADIENT returnVariable:V temperature:t volume:v params:params];
	double lastGoodVolume = v;
	int iter = 0;
	while ((fabs(p-targetP) > 0.1) && (iter < 100)) {
		double f, df;
		@try {
			p = -[self helmholtz:GRADIENT returnVariable:V temperature:t volume:v params:params];
			f  = p - targetP;
			df = -[self helmholtz:HESSIAN returnVariable:VV temperature:t volume:v params:params];
			lastGoodVolume = v;
			v -= f/df;
			if      (v < 0.10*vRef) v = vRef/4.0;
			else if (v > 10.0*vRef) v = 4.0*vRef;
			iter++;
		}
		@catch (NSException * e) {
			BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
			if (debug) NSLog(@"Caught an exception at t = %g, p = %g, iter = %d  v = %g last good v = %g vRef = %g", t, targetP, iter, v, lastGoodVolume, vRef);
			v = (v+lastGoodVolume)/2.0;
		}
	}
	if (iter == 100) {
		NSLog(@"FAILURE in Stixrude endmembers: Iterations exceeded at t = %g, p = %g, p est = %g", t, targetP, p);
		NSDictionary *resultDictionary = [NSDictionary dictionaryWithObjectsAndKeys:@"FAILURE", @"Status",
										  [NSNumber numberWithDouble:t],   @"Temperature",
										  [NSNumber numberWithDouble:p],   @"Pressure",
										  [NSNumber numberWithDouble:0.0], @"GibbsFreeEnergy",
										  [NSNumber numberWithDouble:0.0], @"HelmholtzFreeEnergy",
										  [NSNumber numberWithDouble:0.0], @"Enthalpy",
										  [NSNumber numberWithDouble:0.0], @"InternalEnergy",
										  [NSNumber numberWithDouble:0.0], @"Entropy",
										  [NSNumber numberWithDouble:0.0], @"Volume",
										  [NSNumber numberWithDouble:0.0], @"HeatCapacityAtConstantVolume",
										  [NSNumber numberWithDouble:0.0], @"HeatCapacityAtConstantPressure",
										  [NSNumber numberWithDouble:0.0], @"ThermalExpansionIsobaric",
										  [NSNumber numberWithDouble:0.0], @"CompressibilityIsothermal",
										  [NSNumber numberWithDouble:0.0], @"dVolume/dT",
										  [NSNumber numberWithDouble:0.0], @"dVolume/dP",
                                          [NSNumber numberWithDouble:0.0], @"d2Volume/dT2",
                                          [NSNumber numberWithDouble:0.0], @"d2Volume/dP2",
                                          [NSNumber numberWithDouble:0.0], @"d2Volume/dTdP",
                                          [NSNumber numberWithDouble:0.0], @"dCp/dT",
										  nil];
		return resultDictionary;
	}

	double a = [self helmholtz:FUNCTION returnVariable:0 temperature:t volume:v params:params];

	double s = -[self helmholtz:GRADIENT returnVariable:T temperature:t volume:v params:params];
	p = -[self helmholtz:GRADIENT returnVariable:V temperature:t volume:v params:params];

	double Cv    = -t*[self helmholtz:HESSIAN returnVariable:TT temperature:t volume:v params:params];
	double beta  = 1.0/v/[self helmholtz:HESSIAN returnVariable:VV temperature:t volume:v params:params];
	double K     = 1.0/beta;
	double alpha = -[self helmholtz:HESSIAN returnVariable:TV temperature:t volume:v params:params]*beta;
	double Cp    = Cv + alpha*alpha*v*t*K;

    double d2adv2   = [self helmholtz:HESSIAN returnVariable:VV temperature:t volume:v params:params];
    double d2advdt  = [self helmholtz:HESSIAN returnVariable:TV temperature:t volume:v params:params];
    double d2adt2   = [self helmholtz:HESSIAN returnVariable:TT temperature:t volume:v params:params];
    double d3adv3   = [self helmholtz:TENSOR returnVariable:VVV temperature:t volume:v params:params];
    double d3adv2dt = [self helmholtz:TENSOR returnVariable:VVT temperature:t volume:v params:params];
    double d3advdt2 = [self helmholtz:TENSOR returnVariable:VTT temperature:t volume:v params:params];
    double d3adt3   = [self helmholtz:TENSOR returnVariable:TTT temperature:t volume:v params:params];

    double dvdp = -1.0/d2adv2;
    double dvdt = -d2advdt/d2adv2;

    double Dd2adv2Dp  = d3adv3*dvdp;
    double Dd2adv2Dt  = d3adv3*dvdt + d3adv2dt;
    double Dd2advdtDt = d3adv2dt*dvdt + d3advdt2;

    double d2vdt2  = - Dd2advdtDt/d2adv2 + (d2advdt/d2adv2/d2adv2)*Dd2adv2Dt;
    double d2vdtdp = Dd2adv2Dt/d2adv2/d2adv2;
    double d2vdp2  = Dd2adv2Dp/d2adv2/d2adv2;
    double dcpdt   = -d2adt2 - t*(d3advdt2*dvdt + d3adt3) - dvdt*dvdt/dvdp - t*2.0*dvdt*d2vdt2/dvdp + t*dvdt*dvdt*d2vdtdp/dvdp/dvdp;

	NSDictionary *resultDictionary = [NSDictionary dictionaryWithObjectsAndKeys:@"OK", @"Status",
									  [NSNumber numberWithDouble:t],         @"Temperature",
									  [NSNumber numberWithDouble:p],         @"Pressure",
									  [NSNumber numberWithDouble:a+p*v],     @"GibbsFreeEnergy",
									  [NSNumber numberWithDouble:a],         @"HelmholtzFreeEnergy",
									  [NSNumber numberWithDouble:a+p*v+t*s], @"Enthalpy",
									  [NSNumber numberWithDouble:a-p*v],     @"InternalEnergy",
									  [NSNumber numberWithDouble:s],         @"Entropy",
									  [NSNumber numberWithDouble:v],         @"Volume",
									  [NSNumber numberWithDouble:Cv],        @"HeatCapacityAtConstantVolume",
									  [NSNumber numberWithDouble:Cp],        @"HeatCapacityAtConstantPressure",
									  [NSNumber numberWithDouble:alpha],     @"ThermalExpansionIsobaric",
									  [NSNumber numberWithDouble:beta],      @"CompressibilityIsothermal",
									  [NSNumber numberWithDouble:alpha*v],   @"dVolume/dT",
									  [NSNumber numberWithDouble:-beta*v],   @"dVolume/dP",
                                      [NSNumber numberWithDouble:d2vdt2],    @"d2Volume/dT2",
                                      [NSNumber numberWithDouble:d2vdp2],    @"d2Volume/dP2",
                                      [NSNumber numberWithDouble:d2vdtdp],   @"d2Volume/dTdP",
                                      [NSNumber numberWithDouble:dcpdt],     @"dCp/dT",
									  nil];
	return resultDictionary;

}

@end
