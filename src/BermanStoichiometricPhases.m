//
//  BermanStoichiometricPhases.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "BermanStoichiometricPhases.h"
#ifndef __APPLE__
#include <dispatch/dispatch.h>
#endif

#define QUARTZ_ADJUSTMENT -1291.0

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))
#define QUINTIC(x) ((x)*(x)*(x)*(x)*(x))

@implementation AegirineBerman

-(id)init {
	if ((self = [super initWithH:-2576800.0
							  S:170.57
							 k0:311.29+303.909-297.499
							 k1:-20.051e2+(-14.1767e2)-(-13.5596e2)
							 k2:-53.503e5+(-43.5654e5)-(-67.0219e5)
							 k3:66.257e7+(35.2523e7)-(75.9082e7)
							 v0:63.997
							 v1:0.0
							 v2:0.0
							 v3:0.0
							 v4:0.0])) {
		[self setPhaseFormula:@"NaFeSi2O6"];
		[self setPhaseName:@"Aegirine"];
	}
	return self;
}

@end

@implementation AenigmatiteBerman

-(id)init {
	if ((self = [super initWithH:-8472805.0
							  S:740.23
							 k0:1092.073
							 k1:-61.5981e2
							 k2:-225.2679e5
							 k3:326.8659e7
							 v0:22.8546
							 v1:0.0
							 v2:0.0
							 v3:0.0
							 v4:0.0])) {
		[self setPhaseFormula:@"Na2Fe5TiSi6O20"];
		[self setPhaseName:@"Aenigmatite"];
	}
	return self;
}

@end

@implementation AkermaniteBerman

-(id)init {
	if ((self = [super initWithH:-3860441.0
							  S:212.0
							 k0:387.06
							 k1:-29.388e2
							 k2:0.0
							 k3:-4.079e7
							 l1:0.0
							 l2:0.0
							 Tt:358.0
						 deltaH:452.0
							 v0:9.252
							 v1:-0.785e-6
							 v2:0.0
							 v3:25.011e-6
							 v4:67.224e-10])) {
		[self setPhaseFormula:@"Ca2MgSi2O7"];
		[self setPhaseName:@"Akermanite"];
	}
	return self;
}

@end

@implementation AndalusiteBerman

-(id)init {
	if ((self = [super initWithH:-2589972.0
							  S:91.434
							 k0:236.48
							 k1:-11.029e2
							 k2:-75.268e5
							 k3:93.644e7
							 v0:5.147
							 v1:-0.770e-6
							 v2:1.923e-12
							 v3:23.443e-6
							 v4:7.189e-10])) {
		[self setPhaseFormula:@"Al2SiO5"];
		[self setPhaseName:@"Andalusite"];
	}
	return self;
}

@end

@implementation ApatiteBerman

-(id)init {
	if ((self = [super initWithH:-6694689.0
							  S:398.74
							 k0:758.81
							 k1:-64.806E2
							 k2:0.0
							 k3:44.794E7
							 v0:16.4025
							 v1:0.0
							 v2:0.0
							 v3:0.0
							 v4:0.0])) {
		[self setPhaseFormula:@"Ca5(PO4)3OH"];
		[self setPhaseName:@"Apatite"];
	}
	return self;
}

@end

@implementation ChromiteBerman

-(id)init {
	if ((self = [super initWithH:-1445490.0
							  S:142.676
							 k0:236.874
							 k1:-16.796E2
							 k2:0.0
							 k3:-16.765E7
							 v0:4.4010
							 v1:0.0
							 v2:0.0
							 v3:0.0
							 v4:0.0])) {
		[self setPhaseFormula:@"FeCr2O4"];
		[self setPhaseName:@"Chromite"];
	}
	return self;
}

@end

@implementation CoesiteBerman

-(id)init {
	if ((self = [super initWithH:-907604.0
							  S:39.424
							 k0: 94.907
							 k1:-7.673e2
							 k2:-5.279e5
							 k3:2.627e7
							 v0:2.064
							 v1:-1.037e-6
							 v2:3.000e-12
							 v3:7.396e-6
							 v4:43.605e-10])) {
		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Coesite"];
	}
	return self;
}

@end


@implementation CorundumBerman

-(id)init {
	if ((self = [super initWithH:-1675700.0
							  S:50.820
							 k0:155.02
							 k1:-8.284E2
							 k2:-38.614E5
							 k3:40.908E7
							 v0:2.558
							 v1:-0.385E-6
							 v2:0.375E-12
							 v3:21.342E-6
							 v4:47.180E-10])) {
		[self setPhaseFormula:@"Al2O3"];
		[self setPhaseName:@"Corundum"];
	}
	return self;
}

@end

@implementation CristobaliteBerman

-(id)init {
	if ((self = [super initWithH:0 S:0 k0:0 k1:0 k2:0 k3:0])) {
		alphaCristobalite = [[BermanProperties alloc] initWithH:-907753.0
															  S:43.394
															 k0:83.51
															 k1:-3.747E2
															 k2:-24.554E5
															 k3:28.007E7
															 v0:2.587
															 v1:-2.515E-6
															 v2:0.0
															 v3:20.824E-6
															 v4:0.0];
        [alphaCristobalite setPhaseFormula:@"SiO2"];
        [alphaCristobalite setPhaseName:@"alpha-Cristobalite"];
		betaCristobalite = [[BermanProperties alloc] initWithH:-906377.0
															 S:46.029
															k0:83.51
															k1:-3.747E2
															k2:-24.554E5
															k3:28.007E7
															v0:2.730
															v1:-1.100e-6
															v2:5.535e-12
															v3:3.189e-6
															v4:0.0];
        [betaCristobalite setPhaseFormula:@"SiO2"];
        [betaCristobalite setPhaseName:@"beta-Cristobalite"];

		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Cristobalite"];
		Tt     = 535.0;     // These are alpha-quartz properties to be applied here
		l1     = -14.216E-2;
		l2     = 44.142E-5;
		Tt     = 535.0;
		deltaH = 0.0;
		dTtdp  = 0.0480;    //  Berman (1988) dt/dp
	}
	return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kCristobaliteBerman_dTtdp = @"dTtdp";
static NSString *kCristobaliteBerman_alphaCristobalite = @"alphaCristobalite";
static NSString *kCristobaliteBerman_betaCristobalite = @"betaCristobalite";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        dTtdp = [aDecoder decodeDoubleForKey:kCristobaliteBerman_dTtdp];
#ifdef __APPLE__
        alphaCristobalite = (BermanProperties *) [aDecoder decodeObjectOfClass:[BermanProperties class] forKey:kCristobaliteBerman_alphaCristobalite];
        betaCristobalite  = (BermanProperties *) [aDecoder decodeObjectOfClass:[BermanProperties class] forKey:kCristobaliteBerman_betaCristobalite];
#else
        alphaCristobalite = (BermanProperties *) [aDecoder decodeObjectForKey:kCristobaliteBerman_alphaCristobalite];
        betaCristobalite  = (BermanProperties *) [aDecoder decodeObjectForKey:kCristobaliteBerman_betaCristobalite];
#endif
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:dTtdp forKey:kCristobaliteBerman_dTtdp];
        [aCoder encodeObject:alphaCristobalite forKey:kCristobaliteBerman_alphaCristobalite];
        [aCoder encodeObject:betaCristobalite forKey:kCristobaliteBerman_betaCristobalite];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

-(BOOL)isAlphaPhaseAtT:(double)t andP:(double)p {
    double cp_t = Tt + dTtdp*(p-1.0);
    if (t > cp_t) return NO;
    return YES;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [self getEnthalpyFromT:t andP:p] - t*[self getEntropyFromT:t andP:p];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getEnthalpyFromT:t andP:p];

	double result = [alphaCristobalite getEnthalpyFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;
	double x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double x4 = l2*l2;

	result += x1*(t-(trl-delt)) + x2*(t*t-SQUARE(trl-delt))/2.0 + x3*(t*t*t-CUBE(trl-delt))/3.0 + x4*(t*t*t*t-QUARTIC(trl-delt))/4.0;
	return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getEntropyFromT:t andP:p];

	double result = [alphaCristobalite getEntropyFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;
	double x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double x4 = l2*l2;

	result += x1*(log(t)-log(trl-delt)) + x2*(t-(trl-delt)) + x3*(t*t-SQUARE(trl-delt))/2.0 + x4*(t*t*t-CUBE(trl-delt))/3.0;
	return result;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getHeatCapacityFromT:t andP:p];

	double result = [alphaCristobalite getHeatCapacityFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;

	result += (t+delt)*SQUARE(l1+l2*(t+delt));
	return result;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getDcpDtFromT:t andP:p];
	double result = [alphaCristobalite getDcpDtFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;

	result += SQUARE(l1+l2*(t+delt)) + (t+delt)*2.0*(l1+l2*(t+delt))*l2;
	return result;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getVolumeFromT:t andP:p];
	double result = [alphaCristobalite getVolumeFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1    = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double x2    = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double Dx2Dp = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double x3    = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp = 3.0*l2*l2*DdeltDp;
	double x4    = l2*l2;
	double Dx4Dp = 0.0;

	result += Dx1Dp*(t-(trl-delt)) + Dx2Dp*(t*t-SQUARE(trl-delt))/2.0 + Dx3Dp*(t*t*t-CUBE(trl-delt))/3.0
	+ Dx4Dp*(t*t*t*t-QUARTIC(trl-delt))/4.0
	- t*( Dx1Dp*(log(t)-log(trl-delt)) + Dx2Dp*(t-(trl-delt))
		 + Dx3Dp*(t*t-SQUARE(trl-delt))/2.0
		 + Dx4Dp*(t*t*t-CUBE(trl-delt))/3.0
		 )
	+ x1*DdeltDp + x2*(trl-delt)*DdeltDp
	+ x3*SQUARE(trl-delt)*DdeltDp
	+ x4*CUBE(trl-delt)*DdeltDp
	- t*(x1*DdeltDp/(trl-delt) + x2*DdeltDp + x3*(trl-delt)*DdeltDp
		 + x4*SQUARE(trl-delt)*DdeltDp
		 );
	return result;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getDvDtFromT:t andP:p];
	double result = [alphaCristobalite getDvDtFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1    = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double x2    = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double Dx2Dp = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double x3    = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp = 3.0*l2*l2*DdeltDp;
	double x4    = l2*l2;
	double Dx4Dp = 0.0;

	result += Dx1Dp + Dx2Dp*t + Dx3Dp*t*t + Dx4Dp*t*t*t
	- ( Dx1Dp*(log(t)-log(trl-delt)) + Dx2Dp*(t-(trl-delt))
	   + Dx3Dp*(t*t-SQUARE(trl-delt))/2.0
	   + Dx4Dp*(t*t*t-CUBE(trl-delt))/3.0
	   ) - t*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
	- (x1*DdeltDp/(trl-delt) + x2*DdeltDp + x3*(trl-delt)*DdeltDp
	   + x4*SQUARE(trl-delt)*DdeltDp
	   );
	return result;
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getDvDpFromT:t andP:p];
	double result = [alphaCristobalite getDvDpFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
	double x2      = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
	double x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp   = 3.0*l2*l2*DdeltDp;
	double D2x3Dp2 = 0.0;
	double x4      = l2*l2;
	double Dx4Dp   = 0.0;
	double D2x4Dp2 = 0.0;

	result += D2x1Dp2*(t-(trl-delt)) + Dx1Dp*DdeltDp
	+ D2x2Dp2*(t*t-SQUARE(trl-delt))/2.0 + Dx2Dp*(trl-delt)*DdeltDp
	+ D2x3Dp2*(t*t*t-CUBE(trl-delt))/3.0
	+ Dx3Dp*SQUARE(trl-delt)*DdeltDp
	+ D2x4Dp2*(t*t*t*t-QUARTIC(trl-delt))/4.0
	+ Dx4Dp*CUBE(trl-delt)*DdeltDp
	- t*(  D2x1Dp2*(log(t)-log(trl-delt)) + Dx1Dp*DdeltDp/(trl-delt)
		 + D2x2Dp2*(t-(trl-delt)) + Dx2Dp*DdeltDp
		 + D2x3Dp2*(t*t-SQUARE(trl-delt))/2.0 + Dx3Dp*(trl-delt)*DdeltDp
		 + D2x4Dp2*(t*t*t-CUBE(trl-delt))/3.0
		 + Dx4Dp*SQUARE(trl-delt)*DdeltDp
		 )
	+ Dx1Dp*DdeltDp
	+ Dx2Dp*(trl-delt)*DdeltDp - x2*DdeltDp*DdeltDp
	+ Dx3Dp*SQUARE(trl-delt)*DdeltDp
	- x3*2.0*(trl-delt)*SQUARE(DdeltDp)
	+ Dx4Dp*CUBE(trl-delt)*DdeltDp
	- x4*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	- t*(  Dx1Dp*DdeltDp/(trl-delt)
		 + x1*SQUARE(DdeltDp)/SQUARE(trl-delt)
		 + Dx2Dp*DdeltDp
		 + Dx3Dp*(trl-delt)*DdeltDp - x3*SQUARE(DdeltDp)
		 + Dx4Dp*SQUARE(trl-delt)*DdeltDp
		 - x4*2.0*(trl-delt)*SQUARE(DdeltDp)
		 );
	return result;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getD2vDt2FromT:t andP:p];
	double result = [alphaCristobalite getD2vDt2FromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double Dx1Dp = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double Dx2Dp = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double Dx3Dp = 3.0*l2*l2*DdeltDp;
	double Dx4Dp = 0.0;

	result += Dx2Dp + 2.0*Dx3Dp*t + 3.0*Dx4Dp*t*t
	- 2.0*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
	- t*(-Dx1Dp/(t*t) + Dx3Dp + 2.0*Dx4Dp*t);
	return result;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getD2vDtDpFromT:t andP:p];
	double result = [alphaCristobalite getD2vDtDpFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
	double Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
	double x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp   = 3.0*l2*l2*DdeltDp;
	double D2x3Dp2 = 0.0;
	double x4      = l2*l2;
	double Dx4Dp   = 0.0;
	double D2x4Dp2 = 0.0;

	result += D2x1Dp2 + D2x2Dp2*t + D2x3Dp2*t*t + D2x4Dp2*t*t*t
	- (  D2x1Dp2*(log(t)-log(trl-delt)) + Dx1Dp*DdeltDp/(trl-delt)
	   + D2x2Dp2*(t-(trl-delt)) + Dx2Dp*DdeltDp
	   + D2x3Dp2*(t*t-SQUARE(trl-delt))/2.0 + Dx3Dp*(trl-delt)*DdeltDp
	   + D2x4Dp2*(t*t*t-CUBE(trl-delt))/3.0
	   + Dx4Dp*SQUARE(trl-delt)*DdeltDp
	   )
	- t*(D2x1Dp2/t + D2x2Dp2 + D2x3Dp2*t + D2x4Dp2*t*t)
	- (  Dx1Dp*DdeltDp/(trl-delt)
	   + x1*SQUARE(DdeltDp)/SQUARE(trl-delt)
	   + Dx2Dp*DdeltDp
	   + Dx3Dp*(trl-delt)*DdeltDp - x3*SQUARE(DdeltDp)
	   + Dx4Dp*SQUARE(trl-delt)*DdeltDp
	   - x4*2.0*(trl-delt)*SQUARE(DdeltDp)
	   );
	return result;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaCristobalite getD2vDp2FromT:t andP:p];
	double result = [alphaCristobalite getD2vDp2FromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;

	double x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
	double D3x1Dp3 = 6.0*l2*l2*DdeltDp*DdeltDp*DdeltDp;

	double Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
	double D3x2Dp3 = 0.0;

	double x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp   = 3.0*l2*l2*DdeltDp;
	double D2x3Dp2 = 0.0;
	double D3x3Dp3 = 0.0;

	double x4      = l2*l2;
	double Dx4Dp   = 0.0;
	double D2x4Dp2 = 0.0;
	double D3x4Dp3 = 0.0;

	result += D3x1Dp3*(t-(trl-delt)) + D2x1Dp2*DdeltDp
	+ D2x1Dp2*DdeltDp
	+ D3x2Dp3*(t*t-SQUARE(trl-delt))/2.0 + D2x2Dp2*(trl-delt)*DdeltDp
	+ D2x2Dp2*(trl-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
	+ D3x3Dp3*(t*t*t-CUBE(trl-delt))/3.0
	+ D2x3Dp2*SQUARE(trl-delt)*DdeltDp
	+ D2x3Dp2*SQUARE(trl-delt)*DdeltDp
	- Dx3Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
	+ D3x4Dp3*(t*t*t*t-QUARTIC(trl-delt))/4.0
	+ D2x4Dp2*CUBE(trl-delt)*DdeltDp
	+ D2x4Dp2*CUBE(trl-delt)*DdeltDp
	- Dx4Dp*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	- t*(  D3x1Dp3*(log(t)-log(trl-delt)) + D2x1Dp2*DdeltDp/(trl-delt)
		 + D2x1Dp2*DdeltDp/(trl-delt)
		 + Dx1Dp*SQUARE(DdeltDp/(trl-delt))
		 + D3x2Dp3*(t-(trl-delt)) + D2x2Dp2*DdeltDp
		 + D2x2Dp2*DdeltDp
		 + D3x3Dp3*(t*t-SQUARE(trl-delt))/2.0
		 + D2x3Dp2*(trl-delt)*DdeltDp
		 + D2x3Dp2*(trl-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
		 + D3x4Dp3*(t*t*t-CUBE(trl-delt))/3.0
		 + D2x4Dp2*SQUARE(trl-delt)*DdeltDp
		 + D2x4Dp2*SQUARE(trl-delt)*DdeltDp
		 - Dx4Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
		 )
	+ D2x1Dp2*DdeltDp
	+ D2x2Dp2*(trl-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
	- Dx2Dp*DdeltDp*DdeltDp
	+ D2x3Dp2*SQUARE(trl-delt)*DdeltDp
	- Dx3Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
	- Dx3Dp*2.0*(trl-delt)*SQUARE(DdeltDp) + x3*2.0*CUBE(DdeltDp)
	+ D2x4Dp2*CUBE(trl-delt)*DdeltDp
	- Dx4Dp*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	- Dx4Dp*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	+ x4*6.0*(trl-delt)*CUBE(DdeltDp)
	- t*(  D2x1Dp2*DdeltDp/(trl-delt)
		 + Dx1Dp*SQUARE(DdeltDp/(trl-delt))
		 + Dx1Dp*SQUARE(DdeltDp/(trl-delt))
		 + x1*2.0*CUBE(DdeltDp/(trl-delt))
		 + D2x2Dp2*DdeltDp
		 + D2x3Dp2*(trl-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
		 - Dx3Dp*SQUARE(DdeltDp)
		 + D2x4Dp2*SQUARE(trl-delt)*DdeltDp
		 - Dx4Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
		 - Dx4Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
		 + x4*2.0*CUBE(DdeltDp)
		 );
	return result;
}

@end

@implementation FayaliteBerman

-(id)init {
	if ((self = [super initWithH:-1479360.0
							  S:150.930
							 k0:248.93
							 k1:-19.239e2
							 k2:0.0
							 k3:-13.910e7
							 v0:4.630
							 v1:-0.730E-6
							 v2:0.0
							 v3:26.546E-6
							 v4:79.482E-10])) {
		[self setPhaseFormula:@"Fe2SiO4"];
		[self setPhaseName:@"Fayalite"];
	}
	return self;
}

@end

@implementation ForsteriteBerman

-(id)init {
	if ((self = [super initWithH:-2174420.0
							  S:94.010
							 k0:238.64
							 k1:-20.013E2
							 k2:0.0
							 k3:-11.624E7
							 v0:4.366
							 v1:-0.791E-6
							 v2:1.351E-12
							 v3:29.464E-6
							 v4:88.633E-10])) {
		[self setPhaseFormula:@"Mg2SiO4"];
		[self setPhaseName:@"Forsterite"];
	}
	return self;
}

@end


@implementation GehleniteBerman

-(id)init {
	if ((self = [super initWithH:-3988158.0
							  S:198.600
							 k0:373.09
							 k1:-22.768e2
							 k2:-47.785e5
							 k3:47.791e7
							 v0:9.033
							 v1:-0.996e-6
							 v2:2.488e-12
							 v3:24.926e-6
							 v4:5.644e-10])) {
		[self setPhaseFormula:@"Ca2Al2SiO7"];
		[self setPhaseName:@"Gehlenite"];
		d0 = -221.74;
		d1 = 0.0;
		d2 = 172.91e5;
		d3 = 36.950e-2;
		d4 = -146.90e-6;
		d5 = 0.0;
	}
	return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kGehleniteBerman_d0 = @"d0";
static NSString *kGehleniteBerman_d1 = @"d1";
static NSString *kGehleniteBerman_d2 = @"d2";
static NSString *kGehleniteBerman_d3 = @"d3";
static NSString *kGehleniteBerman_d4 = @"d4";
static NSString *kGehleniteBerman_d5 = @"d5";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        d0 = [aDecoder decodeDoubleForKey:kGehleniteBerman_d0];
        d1 = [aDecoder decodeDoubleForKey:kGehleniteBerman_d1];
        d2 = [aDecoder decodeDoubleForKey:kGehleniteBerman_d2];
        d3 = [aDecoder decodeDoubleForKey:kGehleniteBerman_d3];
        d4 = [aDecoder decodeDoubleForKey:kGehleniteBerman_d4];
        d5 = [aDecoder decodeDoubleForKey:kGehleniteBerman_d5];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:d0 forKey:kGehleniteBerman_d0];
        [aCoder encodeDouble:d1 forKey:kGehleniteBerman_d1];
        [aCoder encodeDouble:d2 forKey:kGehleniteBerman_d2];
        [aCoder encodeDouble:d3 forKey:kGehleniteBerman_d3];
        [aCoder encodeDouble:d4 forKey:kGehleniteBerman_d4];
        [aCoder encodeDouble:d5 forKey:kGehleniteBerman_d5];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

/*
-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double td = (t < 1600.0) ? t : 1600.0;
	double dhdis = d0*(td-698.0) + 2.0*d1*(sqrt(td)-sqrt(698.0)) - d2*(1.0/td-1.0/698.0) + d3*(td*td-698.0*698.0)/2.0 + d4*(td*td*td-698.0*698.0*698.0)/3.0;
	double dsdis = d0*(log(td)-log(698.0)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(698.0)) - d2*(1.0/(td*td)-1.0/(698.0*698.0))/2.0 + d3*(td-698.0) + d4*((td*td)-(698.0*698.0))/2.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	return [super getGibbsFreeEnergyFromT:t andP:p] + dhdis - t*dsdis + dvdis*(p-pr);
}
*/
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	double td = (t < 1600.0) ? t : 1600.0;
	double dhdis = d0*(td-698.0) + 2.0*d1*(sqrt(td)-sqrt(698.0)) - d2*(1.0/td-1.0/698.0) + d3*(td*td-698.0*698.0)/2.0 + d4*(td*td*td-698.0*698.0*698.0)/3.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	dhdis += (t < 1600.0 && d5 != 0.0) ? -(p-pr)*t*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
	return [super getEnthalpyFromT:t andP:p] + dhdis + dvdis*(p-pr);
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	double td = (t < 1600.0) ? t : 1600.0;
	double dsdis = d0*(log(td)-log(698.0)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(698.0)) - d2*(1.0/(td*td)-1.0/(698.0*698.0))/2.0 + d3*(td-698.0) + d4*((td*td)-(698.0*698.0))/2.0;
	dsdis += (t < 1600.0 && d5 != 0.0) ? -(p-pr)*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
	return [super getEntropyFromT:t andP:p] + dsdis;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	double cpdis = (t < 1600.0) ? d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t : 0.0;
	cpdis += (t < 1600.0 && d5 != 0.0) ? -(p-pr)*t*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
	return [super getHeatCapacityFromT:t andP:p] + cpdis;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	double dcpdisdt = (t < 1600.0) ? -d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t : 0.0;
	dcpdisdt += (t < 1600.0 && d5 != 0.0) ? -(p-pr)*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5
	-(p-pr)*t*(1.5*0.5*d1/pow(t,(double) 2.5) + 6.0*d2/(t*t*t*t) + 2.0*d4)/d5 : 0.0;
	return [super getDcpDtFromT:t andP:p] + dcpdisdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	double td = (t < 1600.0) ? t : 1600.0;
	double dhdis = d0*(td-698.0) + 2.0*d1*(sqrt(td)-sqrt(698.0)) - d2*(1.0/td-1.0/698.0) + d3*(td*td-698.0*698.0)/2.0 + d4*(td*td*td-698.0*698.0*698.0)/3.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	return [super getVolumeFromT:t andP:p] + dvdis;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	double dvdisdt = (t < 1600.0 && d5 != 0.0) ? (d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*(t*t))/d5 : 0.0;
	return [super getDvDtFromT:t andP:p] + dvdisdt;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	double d2vdisdt2 = (t < 1600.0 && d5 != 0.0) ? (-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
	return [super getD2vDt2FromT:t andP:p] + d2vdisdt2;
}

@end

@implementation HematiteBerman

-(id)init {
	if ((self = [super initWithH:-825627.0
							  S:87.437
							 k0:146.86
							 k1:0.0
							 k2:-55.768e5
							 k3:52.563e7
							 l1:-7.403e-2
							 l2:27.921e-5
							 Tt:955.0
						 deltaH:1287.0
							 v0:3.027
							 v1:-0.479e-6
							 v2:0.304e-12
							 v3:38.310e-6
							 v4:1.650e-10])) {
		[self setPhaseFormula:@"Fe2O3"];
		[self setPhaseName:@"Hematite"];
	}
	return self;
}

@end

@implementation IlmeniteBerman

-(id)init {
	if ((self = [super initWithH:-1231947.0
							  S:108.628
							 k0:150.0
							 k1:-4.416e2
							 k2:-33.237e5
							 k3:34.815e7
							 v0:3.170
							 v1:-0.584e-6
							 v2:1.230e-12
							 v3:27.248e-6
							 v4:29.968e-10])) {
		[self setPhaseFormula:@"FeTiO3"];
		[self setPhaseName:@"Ilmenite"];
	}
	return self;
}

@end

@implementation KalsiliteBerman

-(id)init {
	if ((self = [super initWithH:-2109563.55-(-1.35*1000.0*4.184)/4.0-3572.76/4.0
							  S:133.9653
							 k0:186.0
							 k1:0.0
							 k2:-131.067E5
							 k3:213.893E7
							 l1:-7.096454E-2
							 l2:21.682E-5
							 Tt:800.15
						 deltaH:1154.0
							 v0:6.043478-(-0.00001*1000.0*4.184)/4.0
							 v1:-2.0500e-6
							 v2:5.2000e-12
							 v3:31.802e-6
							 v4:213.0e-10])) {
		[self setPhaseFormula:@"KAlSiO4"];
		[self setPhaseName:@"Kalsilite"];
	}
	return self;
}

@end

@implementation KyaniteBerman

-(id)init {
	if ((self = [super initWithH:-2594220.0
							  S:82.430
							 k0:262.68
							 k1:-20.014e2
							 k2:-19.997e5
							 k3:-6.318e7
							 v0:4.412
							 v1:-0.646e-6
							 v2:0.0
							 v3:23.973e-6
							 v4:0.0])) {
		[self setPhaseFormula:@"Al2SiO5"];
		[self setPhaseName:@"Kyanite"];
	}
	return self;
}

@end

@implementation LimeBerman

-(id)init {
	if ((self = [super initWithH:-635090.
							  S:37.750
							 k0:58.79
							 k1:-1.339e2
							 k2:-11.471e5
							 k3:10.298e7
							 v0:1.676
							 v1:-1.022e-6
							 v2:2.565e-12
							 v3:34.610e-6
							 v4:67.406e-10])) {
		[self setPhaseFormula:@"CaO"];
		[self setPhaseName:@"Lime"];
	}
	return self;
}

@end

@implementation MagnetiteBerman

-(id)init {
	if ((self = [super initWithH:-1117403.0
							  S:146.114
							 k0:207.93
							 k1:0.0
							 k2:-72.433e5
							 k3:66.436e7
							 l1:-19.502e-2
							 l2:61.037e-5
							 Tt:848.0
						 deltaH:1565.0
							 v0:4.452
							 v1:-0.582e-6
							 v2:1.751e-12
							 v3:30.291e-6
							 v4:138.470e-10])) {
		[self setPhaseFormula:@"Fe3O4"];
		[self setPhaseName:@"Magnetite"];
	}
	return self;
}

@end

@implementation MuscoviteBerman

-(id)init {
	if ((self = [super initWithH:-5976740.0
							  S:293.157
							 k0:651.49
							 k1:-38.732e2
							 k2:-185.232e5
							 k3:274.247e7
							 v0:14.087
							 v1:-1.717e-6
							 v2:4.295e-12
							 v3:33.527e-6
							 v4:0.0])) {
		[self setPhaseFormula:@"KAl2Si3AlO10(OH)2"];
		[self setPhaseName:@"Muscovite"];
	}
	return self;
}

@end

@implementation NephelineBerman

-(id)init {
	if ((self = [super initWithH:-2093004.0
							  S:124.641
							 k0:205.24
							 k1:-7.599E2
							 k2:-108.383E5
							 k3:208.182E7
							 l1:-50.249E-2
							 l2:165.95E-5
							 Tt:467.0
						 deltaH:241.0
							 v0:5.4131
							 v1:-2.0500e-6
							 v2:5.2000e-12
							 v3:31.802e-6
							 v4:213.0e-10])) {
		[self setPhaseFormula:@"NaAlSiO4"];
		[self setPhaseName:@"Nepheline"];
	}
	return self;
}

@end


@implementation PericlaseBerman

-(id)init {
	if ((self = [super initWithH:-601500.0
							  S:26.951
							 k0:61.10965
							 k1:-2.96199e2
							 k2:-6.21154e5
							 k3:0.5844612e7
							 v0:1.125
							 v1:-0.6223329e-6
							 v2:1.5114e-12
							 v3:37.4774182e-6
							 v4:3.556e-10])) {
		[self setPhaseFormula:@"MgO"];
		[self setPhaseName:@"Periclase"];
	}
	return self;
}

@end

@implementation PerovskiteBerman

-(id)init {
	if ((self = [super initWithH:-1660630.0
							  S:93.64
							 k0:150.49
							 k1:-6.213E2
							 k2:0.0
							 k3:-43.010E7
							 l1:0.0
							 l2:0.0
							 Tt:1530.0
						 deltaH:550.0*4.184
							 v0:3.3626
							 v1:0.0
							 v2:0.0
							 v3:0.0
							 v4:0.0])) {
		[self setPhaseFormula:@"CaTiO3"];
		[self setPhaseName:@"Perovskite"];
	}
	return self;
}

@end

@implementation PhlogopiteBerman

-(id)init {
	if ((self = [super initWithH:-6207342.0
							  S:334.158
							 k0:610.38
							 k1:-20.838e2
							 k2:-215.330e5
							 k3:284.104e7
							 v0:14.977
							 v1:-1.697e-6
							 v2:0.0
							 v3:34.447e-6
							 v4:0.0])) {
		[self setPhaseFormula:@"KMg3AlSi3O10(OH)2"];
		[self setPhaseName:@"Phlogopite"];
	}
	return self;
}

@end

@interface QuartzBermanSettings : NSObject
@property (atomic, assign, getter = isQuartzCorrectionUsed) BOOL quartzCorrectionUsed;
@end

@implementation QuartzBermanSettings
- (instancetype)init {
    self = [super init];
    if (self) {
        _quartzCorrectionUsed = YES;
    }
    return self;
}
@end

@implementation QuartzBerman

+(QuartzBermanSettings *)quartzSettings {
    static QuartzBermanSettings *quartzSettings = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        quartzSettings = [[QuartzBermanSettings alloc] init];
    });
    return quartzSettings;
}

+(void)enableQuartzCorrectionUsed {
    @synchronized (self.quartzSettings) {
        if (!self.quartzSettings.isQuartzCorrectionUsed) {
            self.quartzSettings.quartzCorrectionUsed = YES;
        }
    }
}

+(void)disableQuartzCorrectionUsed {
    @synchronized (self.quartzSettings) {
        if (self.quartzSettings.isQuartzCorrectionUsed) {
            self.quartzSettings.quartzCorrectionUsed = NO;
        }
    }
}

-(id)init {
	if ((self = [super initWithH:0 S:0 k0:0 k1:0 k2:0 k3:0])) {
		alphaQuartz = [[BermanProperties alloc] initWithH:-910700.0
														S:41.460
													   k0:80.01
													   k1:-2.403E2
													   k2:-35.467E5
													   k3:49.157E7
													   v0:2.269
													   v1:-2.434E-6
													   v2:10.137E-12
													   v3:23.895E-6
													   v4:0.0];
        [alphaQuartz setPhaseName:@"alpha-Quartz"];
        [alphaQuartz setPhaseFormula:@"SiO2"];
		betaQuartz = [[BermanProperties alloc] initWithH:-908627.0
													   S:44.207
													  k0:80.01
													  k1:-2.403E2
													  k2:-35.467E5
													  k3:49.157E7
													  v0:2.370
													  v1:-1.238e-6
													  v2:7.087e-13
													  v3:0.0
													  v4:0.0];
        [betaQuartz setPhaseName:@"beta-Quartz"];
        [betaQuartz setPhaseFormula:@"SiO2"];

		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Quartz"];
		Tt     = 848.0;     // These are alpha-quartz properties to be applied here
		l1     = -9.187E-2;
		l2     = 24.607E-5;
		Tt     = 848.0;
		deltaH = 0.0;
		dTtdp  = 0.0237;    //  Berman (1988) dt/dp
		trl    = 373.0;

	}
	return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kQuartzBerman_dTtdp = @"dTtdp";
static NSString *kQuartzBerman_alphaQuartz = @"alphaQuartz";
static NSString *kQuartzBerman_betaQuartz = @"betaQuartz";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        dTtdp = [aDecoder decodeDoubleForKey:kQuartzBerman_dTtdp];
#ifdef __APPLE__
        alphaQuartz = [aDecoder decodeObjectOfClass:[BermanProperties class] forKey:kQuartzBerman_alphaQuartz];
        betaQuartz = [aDecoder decodeObjectOfClass:[BermanProperties class] forKey:kQuartzBerman_betaQuartz];
#else
        alphaQuartz = [aDecoder decodeObjectForKey:kQuartzBerman_alphaQuartz];
        betaQuartz = [aDecoder decodeObjectForKey:kQuartzBerman_betaQuartz];
#endif
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:dTtdp forKey:kQuartzBerman_dTtdp];
        [aCoder encodeObject:alphaQuartz forKey:kQuartzBerman_alphaQuartz];
        [aCoder encodeObject:betaQuartz forKey:kQuartzBerman_betaQuartz];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

-(BOOL)isAlphaPhaseAtT:(double)t andP:(double)p {
    double cp_t = Tt + dTtdp*(p-1.0);
    if (t > cp_t) return NO;
    return YES;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [self getEnthalpyFromT:t andP:p] - t*[self getEntropyFromT:t andP:p];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
    if (t > cp_t) {
        double result = [betaQuartz getEnthalpyFromT:t andP:p];
        if (QuartzBerman.quartzSettings.isQuartzCorrectionUsed) result += (QUARTZ_ADJUSTMENT);
        return  result;
    }

	double result = [alphaQuartz getEnthalpyFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;
	double x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double x4 = l2*l2;

	result += x1*(t-(trl-delt)) + x2*(t*t-SQUARE(trl-delt))/2.0 + x3*(t*t*t-CUBE(trl-delt))/3.0 + x4*(t*t*t*t-QUARTIC(trl-delt))/4.0;
    if (QuartzBerman.quartzSettings.isQuartzCorrectionUsed) result += (QUARTZ_ADJUSTMENT);
	return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getEntropyFromT:t andP:p];

	double result = [alphaQuartz getEntropyFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;
	double x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double x4 = l2*l2;

	result += x1*(log(t)-log(trl-delt)) + x2*(t-(trl-delt)) + x3*(t*t-SQUARE(trl-delt))/2.0 + x4*(t*t*t-CUBE(trl-delt))/3.0;
	return result;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getHeatCapacityFromT:t andP:p];

	double result = [alphaQuartz getHeatCapacityFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;

	result += (t+delt)*SQUARE(l1+l2*(t+delt));
	return result;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getDcpDtFromT:t andP:p];
	double result = [alphaQuartz getDcpDtFromT:t andP:p]; // without lambda transition
	double delt = Tt - cp_t;

	result += SQUARE(l1+l2*(t+delt)) + (t+delt)*2.0*(l1+l2*(t+delt))*l2;
	return result;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getVolumeFromT:t andP:p];
	double result = [alphaQuartz getVolumeFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1    = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double x2    = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double Dx2Dp = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double x3    = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp = 3.0*l2*l2*DdeltDp;
	double x4    = l2*l2;
	double Dx4Dp = 0.0;

	result += Dx1Dp*(t-(trl-delt)) + Dx2Dp*(t*t-SQUARE(trl-delt))/2.0 + Dx3Dp*(t*t*t-CUBE(trl-delt))/3.0
	+ Dx4Dp*(t*t*t*t-QUARTIC(trl-delt))/4.0
	- t*( Dx1Dp*(log(t)-log(trl-delt)) + Dx2Dp*(t-(trl-delt))
		 + Dx3Dp*(t*t-SQUARE(trl-delt))/2.0
		 + Dx4Dp*(t*t*t-CUBE(trl-delt))/3.0
		 )
	+ x1*DdeltDp + x2*(trl-delt)*DdeltDp
	+ x3*SQUARE(trl-delt)*DdeltDp
	+ x4*CUBE(trl-delt)*DdeltDp
	- t*(x1*DdeltDp/(trl-delt) + x2*DdeltDp + x3*(trl-delt)*DdeltDp
		 + x4*SQUARE(trl-delt)*DdeltDp
		 );
	return result;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getDvDtFromT:t andP:p];
	double result = [alphaQuartz getDvDtFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1    = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double x2    = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double Dx2Dp = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double x3    = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp = 3.0*l2*l2*DdeltDp;
	double x4    = l2*l2;
	double Dx4Dp = 0.0;

	result += Dx1Dp + Dx2Dp*t + Dx3Dp*t*t + Dx4Dp*t*t*t
	- ( Dx1Dp*(log(t)-log(trl-delt)) + Dx2Dp*(t-(trl-delt))
	   + Dx3Dp*(t*t-SQUARE(trl-delt))/2.0
	   + Dx4Dp*(t*t*t-CUBE(trl-delt))/3.0
	   ) - t*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
	- (x1*DdeltDp/(trl-delt) + x2*DdeltDp + x3*(trl-delt)*DdeltDp
	   + x4*SQUARE(trl-delt)*DdeltDp
	   );
	return result;
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getDvDpFromT:t andP:p];
	double result = [alphaQuartz getDvDpFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
	double x2      = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
	double Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
	double x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp   = 3.0*l2*l2*DdeltDp;
	double D2x3Dp2 = 0.0;
	double x4      = l2*l2;
	double Dx4Dp   = 0.0;
	double D2x4Dp2 = 0.0;

	result += D2x1Dp2*(t-(trl-delt)) + Dx1Dp*DdeltDp
	+ D2x2Dp2*(t*t-SQUARE(trl-delt))/2.0 + Dx2Dp*(trl-delt)*DdeltDp
	+ D2x3Dp2*(t*t*t-CUBE(trl-delt))/3.0
	+ Dx3Dp*SQUARE(trl-delt)*DdeltDp
	+ D2x4Dp2*(t*t*t*t-QUARTIC(trl-delt))/4.0
	+ Dx4Dp*CUBE(trl-delt)*DdeltDp
	- t*(  D2x1Dp2*(log(t)-log(trl-delt)) + Dx1Dp*DdeltDp/(trl-delt)
		 + D2x2Dp2*(t-(trl-delt)) + Dx2Dp*DdeltDp
		 + D2x3Dp2*(t*t-SQUARE(trl-delt))/2.0 + Dx3Dp*(trl-delt)*DdeltDp
		 + D2x4Dp2*(t*t*t-CUBE(trl-delt))/3.0
		 + Dx4Dp*SQUARE(trl-delt)*DdeltDp
		 )
	+ Dx1Dp*DdeltDp
	+ Dx2Dp*(trl-delt)*DdeltDp - x2*DdeltDp*DdeltDp
	+ Dx3Dp*SQUARE(trl-delt)*DdeltDp
	- x3*2.0*(trl-delt)*SQUARE(DdeltDp)
	+ Dx4Dp*CUBE(trl-delt)*DdeltDp
	- x4*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	- t*(  Dx1Dp*DdeltDp/(trl-delt)
		 + x1*SQUARE(DdeltDp)/SQUARE(trl-delt)
		 + Dx2Dp*DdeltDp
		 + Dx3Dp*(trl-delt)*DdeltDp - x3*SQUARE(DdeltDp)
		 + Dx4Dp*SQUARE(trl-delt)*DdeltDp
		 - x4*2.0*(trl-delt)*SQUARE(DdeltDp)
		 );
	return result;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getD2vDt2FromT:t andP:p];
	double result = [alphaQuartz getD2vDt2FromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double Dx1Dp = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double Dx2Dp = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double Dx3Dp = 3.0*l2*l2*DdeltDp;
	double Dx4Dp = 0.0;

	result += Dx2Dp + 2.0*Dx3Dp*t + 3.0*Dx4Dp*t*t
	- 2.0*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
	- t*(-Dx1Dp/(t*t) + Dx3Dp + 2.0*Dx4Dp*t);
	return result;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getD2vDtDpFromT:t andP:p];
	double result = [alphaQuartz getD2vDtDpFromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;
	double x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
	double Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
	double x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp   = 3.0*l2*l2*DdeltDp;
	double D2x3Dp2 = 0.0;
	double x4      = l2*l2;
	double Dx4Dp   = 0.0;
	double D2x4Dp2 = 0.0;

	result += D2x1Dp2 + D2x2Dp2*t + D2x3Dp2*t*t + D2x4Dp2*t*t*t
	- (  D2x1Dp2*(log(t)-log(trl-delt)) + Dx1Dp*DdeltDp/(trl-delt)
	   + D2x2Dp2*(t-(trl-delt)) + Dx2Dp*DdeltDp
	   + D2x3Dp2*(t*t-SQUARE(trl-delt))/2.0 + Dx3Dp*(trl-delt)*DdeltDp
	   + D2x4Dp2*(t*t*t-CUBE(trl-delt))/3.0
	   + Dx4Dp*SQUARE(trl-delt)*DdeltDp
	   )
	- t*(D2x1Dp2/t + D2x2Dp2 + D2x3Dp2*t + D2x4Dp2*t*t)
	- (  Dx1Dp*DdeltDp/(trl-delt)
	   + x1*SQUARE(DdeltDp)/SQUARE(trl-delt)
	   + Dx2Dp*DdeltDp
	   + Dx3Dp*(trl-delt)*DdeltDp - x3*SQUARE(DdeltDp)
	   + Dx4Dp*SQUARE(trl-delt)*DdeltDp
	   - x4*2.0*(trl-delt)*SQUARE(DdeltDp)
	   );
	return result;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	double cp_t = Tt + dTtdp*(p-1.0);
	if (t > cp_t) return [betaQuartz getD2vDp2FromT:t andP:p];
	double result = [alphaQuartz getD2vDp2FromT:t andP:p]; // without lambda transition

	double delt = Tt - cp_t;
	double DdeltDp = - dTtdp;

	double x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
	double Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp + 3.0*l2*l2*delt*delt*DdeltDp;
	double D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
	double D3x1Dp3 = 6.0*l2*l2*DdeltDp*DdeltDp*DdeltDp;

	double Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
	double D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
	double D3x2Dp3 = 0.0;

	double x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
	double Dx3Dp   = 3.0*l2*l2*DdeltDp;
	double D2x3Dp2 = 0.0;
	double D3x3Dp3 = 0.0;

	double x4      = l2*l2;
	double Dx4Dp   = 0.0;
	double D2x4Dp2 = 0.0;
	double D3x4Dp3 = 0.0;

	result += D3x1Dp3*(t-(trl-delt)) + D2x1Dp2*DdeltDp
	+ D2x1Dp2*DdeltDp
	+ D3x2Dp3*(t*t-SQUARE(trl-delt))/2.0 + D2x2Dp2*(trl-delt)*DdeltDp
	+ D2x2Dp2*(trl-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
	+ D3x3Dp3*(t*t*t-CUBE(trl-delt))/3.0
	+ D2x3Dp2*SQUARE(trl-delt)*DdeltDp
	+ D2x3Dp2*SQUARE(trl-delt)*DdeltDp
	- Dx3Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
	+ D3x4Dp3*(t*t*t*t-QUARTIC(trl-delt))/4.0
	+ D2x4Dp2*CUBE(trl-delt)*DdeltDp
	+ D2x4Dp2*CUBE(trl-delt)*DdeltDp
	- Dx4Dp*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	- t*(  D3x1Dp3*(log(t)-log(trl-delt)) + D2x1Dp2*DdeltDp/(trl-delt)
		 + D2x1Dp2*DdeltDp/(trl-delt)
		 + Dx1Dp*SQUARE(DdeltDp/(trl-delt))
		 + D3x2Dp3*(t-(trl-delt)) + D2x2Dp2*DdeltDp
		 + D2x2Dp2*DdeltDp
		 + D3x3Dp3*(t*t-SQUARE(trl-delt))/2.0
		 + D2x3Dp2*(trl-delt)*DdeltDp
		 + D2x3Dp2*(trl-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
		 + D3x4Dp3*(t*t*t-CUBE(trl-delt))/3.0
		 + D2x4Dp2*SQUARE(trl-delt)*DdeltDp
		 + D2x4Dp2*SQUARE(trl-delt)*DdeltDp
		 - Dx4Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
		 )
	+ D2x1Dp2*DdeltDp
	+ D2x2Dp2*(trl-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
	- Dx2Dp*DdeltDp*DdeltDp
	+ D2x3Dp2*SQUARE(trl-delt)*DdeltDp
	- Dx3Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
	- Dx3Dp*2.0*(trl-delt)*SQUARE(DdeltDp) + x3*2.0*CUBE(DdeltDp)
	+ D2x4Dp2*CUBE(trl-delt)*DdeltDp
	- Dx4Dp*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	- Dx4Dp*3.0*SQUARE(trl-delt)*SQUARE(DdeltDp)
	+ x4*6.0*(trl-delt)*CUBE(DdeltDp)
	- t*(  D2x1Dp2*DdeltDp/(trl-delt)
		 + Dx1Dp*SQUARE(DdeltDp/(trl-delt))
		 + Dx1Dp*SQUARE(DdeltDp/(trl-delt))
		 + x1*2.0*CUBE(DdeltDp/(trl-delt))
		 + D2x2Dp2*DdeltDp
		 + D2x3Dp2*(trl-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
		 - Dx3Dp*SQUARE(DdeltDp)
		 + D2x4Dp2*SQUARE(trl-delt)*DdeltDp
		 - Dx4Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
		 - Dx4Dp*2.0*(trl-delt)*SQUARE(DdeltDp)
		 + x4*2.0*CUBE(DdeltDp)
		 );
	return result;
}

@end

@implementation RutileBerman

-(id)init {
	if ((self = [super initWithH:-944750.0
							  S:50.460
							 k0:77.84
							 k1:0.0
							 k2:-33.678E5
							 k3:40.294E7
							 v0:1.882
							 v1:-0.454E-6
							 v2:0.584E-12
							 v3:25.716E-6
							 v4:15.409E-10])) {
		[self setPhaseFormula:@"TiO2"];
		[self setPhaseName:@"Rutile"];
	}
	return self;
}

@end

@implementation SanidineBerman

-(id)init {
	if ((self = [super initWithH:-3970791.0+3400.0  // 3400 is the rhyolite sanidine adjustment
							  S:214.145
							 k0:381.37
							 k1:-19.410E2
							 k2:-120.373E5
							 k3:183.643E7
							 v0:10.869
							 v1:-1.805E-6
							 v2:5.112E-12
							 v3:15.145E-6
							 v4:54.850E-10])) {
		[self setPhaseFormula:@"KAlSi3O8"];
		[self setPhaseName:@"sanidine"];
		d0 = 282.98;
		d1 = -4.83e3;
		d2 = 36.21e5;
		d3 = -15.733e-2;
		d4 = 34.770e-6;
		d5 = 41.063e4;
	}
	return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kSanidineBerman_d0 = @"d0";
static NSString *kSanidineBerman_d1 = @"d1";
static NSString *kSanidineBerman_d2 = @"d2";
static NSString *kSanidineBerman_d3 = @"d3";
static NSString *kSanidineBerman_d4 = @"d4";
static NSString *kSanidineBerman_d5 = @"d5";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        d0 = [aDecoder decodeDoubleForKey:kSanidineBerman_d0];
        d1 = [aDecoder decodeDoubleForKey:kSanidineBerman_d1];
        d2 = [aDecoder decodeDoubleForKey:kSanidineBerman_d2];
        d3 = [aDecoder decodeDoubleForKey:kSanidineBerman_d3];
        d4 = [aDecoder decodeDoubleForKey:kSanidineBerman_d4];
        d5 = [aDecoder decodeDoubleForKey:kSanidineBerman_d5];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:d0 forKey:kSanidineBerman_d0];
        [aCoder encodeDouble:d1 forKey:kSanidineBerman_d1];
        [aCoder encodeDouble:d2 forKey:kSanidineBerman_d2];
        [aCoder encodeDouble:d3 forKey:kSanidineBerman_d3];
        [aCoder encodeDouble:d4 forKey:kSanidineBerman_d4];
        [aCoder encodeDouble:d5 forKey:kSanidineBerman_d5];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

/*
-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double td = (t < 1436.0) ? t : 1436.0;
	double dhdis = d0*(td-tr) + 2.0*d1*(sqrt(td)-sqrt(tr)) - d2*(1.0/td-1.0/tr) + d3*(td*td-tr*tr)/2.0 + d4*(td*td*td-tr*tr*tr)/3.0;
	double dsdis = d0*(log(td)-log(tr)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(tr)) - d2*(1.0/(td*td)-1.0/(tr*tr))/2.0 + d3*(td-tr) + d4*((td*td)-(tr*tr))/2.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	return [super getGibbsFreeEnergyFromT:t andP:p] + dhdis - t*dsdis + dvdis*(p-pr);
}
*/
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	double td = (t < 1436.0) ? t : 1436.0;
	double dhdis = d0*(td-tr) + 2.0*d1*(sqrt(td)-sqrt(tr)) - d2*(1.0/td-1.0/tr) + d3*(td*td-tr*tr)/2.0 + d4*(td*td*td-tr*tr*tr)/3.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	dhdis += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
	return [super getEnthalpyFromT:t andP:p] + dhdis + dvdis*(p-pr);
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	double td = (t < 1436.0) ? t : 1436.0;
	double dsdis = d0*(log(td)-log(tr)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(tr)) - d2*(1.0/(td*td)-1.0/(tr*tr))/2.0 + d3*(td-tr) + d4*((td*td)-(tr*tr))/2.0;
	dsdis += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
	return [super getEntropyFromT:t andP:p] + dsdis;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	double cpdis = (t < 1436.0) ? d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t : 0.0;
	cpdis += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
	return [super getHeatCapacityFromT:t andP:p] + cpdis;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	double dcpdisdt = (t < 1436.0) ? -d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t : 0.0;
	dcpdisdt += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5
	-(p-pr)*t*(1.5*0.5*d1/pow(t,(double) 2.5) + 6.0*d2/(t*t*t*t) + 2.0*d4)/d5 : 0.0;
	return [super getDcpDtFromT:t andP:p] + dcpdisdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	double td = (t < 1436.0) ? t : 1436.0;
	double dhdis = d0*(td-tr) + 2.0*d1*(sqrt(td)-sqrt(tr)) - d2*(1.0/td-1.0/tr) + d3*(td*td-tr*tr)/2.0 + d4*(td*td*td-tr*tr*tr)/3.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	return [super getVolumeFromT:t andP:p] + dvdis;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	double dvdisdt = (t < 1436.0 && d5 != 0.0) ? (d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*(t*t))/d5 : 0.0;
	return [super getDvDtFromT:t andP:p] + dvdisdt;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	double d2vdisdt2 = (t < 1436.0 && d5 != 0.0) ? (-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
	return [super getD2vDt2FromT:t andP:p] + d2vdisdt2;
}

@end

@implementation SillimaniteBerman

-(id)init {
	if ((self = [super initWithH:-2586091.0
							  S:95.930
							 k0:256.73
							 k1:-18.872E2
							 k2:-29.774E5
							 k3:25.096E7
							 v0:4.983
							 v1:-0.753E-6
							 v2:0.0E-12
							 v3:13.431E-6
							 v4:0.0E-10])) {
		[self setPhaseFormula:@"Al2SiO5"];
		[self setPhaseName:@"Sillimanite"];
	}
	return self;
}

@end

@implementation SpheneBerman

-(id)init {
	if ((self = [super initWithH:-2596652.0
							  S:129.290
							 k0:234.62
							 k1:-10.403E2
							 k2:-51.183E5
							 k3:59.146E7
							 v0:5.565
							 v1:-0.590E-6
							 v2:0.0
							 v3:25.200E-6
							 v4:0.0])) {
		[self setPhaseFormula:@"CaTiSiO5"];
		[self setPhaseName:@"Sphene"];
	}
	return self;
}

@end

@implementation TridymiteBerman

-(id)init {
	if ((self = [super initWithH:0 S:0 k0:0 k1:0 k2:0 k3:0])) {
		alphaTridymite = [[BermanProperties alloc] initWithH:-907750.0
														   S:43.770
														  k0:75.37
														  k1:0.0
														  k2:-59.581E5
														  k3:95.825E7
														  l1:42.670E-2
														  l2:-144.575E-5
														  Tt:383.0
													  deltaH:130.0
														  v0:2.675
														  v1:-2.508E-6
														  v2:0.0
														  v3:19.339E-6
														  v4:0.0];
        [alphaTridymite setPhaseFormula:@"SiO2"];
        [alphaTridymite setPhaseName:@"alpha-Tridymite"];
		betaTridymite = [[BermanProperties alloc] initWithH:-907045.0
														  S:45.524
														 k0:75.37
														 k1:0.0
														 k2:-59.581E5
														 k3:95.825E7
														 v0:2.737
														 v1:-0.740e-6
														 v2:3.735e-12
														 v3:4.829e-6
														 v4:0.0];
        [betaTridymite setPhaseFormula:@"SiO2"];
        [betaTridymite setPhaseName:@"beta-Tridymite"];

		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Tridymite"];
		Tt     = 383.0;     // These are alpha-quartz properties to be applied here
	}
	return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kTridymiteBerman_alphaTridymite = @"alphaTridymite";
static NSString *kTridymiteBerman_betaTridymite = @"betaTridymite";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
#ifdef __APPLE__
        alphaTridymite = [aDecoder decodeObjectOfClass:[BermanProperties class] forKey:kTridymiteBerman_alphaTridymite];
        betaTridymite = [aDecoder decodeObjectOfClass:[BermanProperties class] forKey:kTridymiteBerman_betaTridymite];
#else
        alphaTridymite = [aDecoder decodeObjectForKey:kTridymiteBerman_alphaTridymite];
        betaTridymite = [aDecoder decodeObjectForKey:kTridymiteBerman_betaTridymite];
#endif
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeObject:alphaTridymite forKey:kTridymiteBerman_alphaTridymite];
        [aCoder encodeObject:betaTridymite forKey:kTridymiteBerman_betaTridymite];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

-(BOOL)isAlphaPhaseAtT:(double)t andP:(double)p {
    if (t > Tt) return NO;
    return YES;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [self getEnthalpyFromT:t andP:p] - t*[self getEntropyFromT:t andP:p];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getEnthalpyFromT:t andP:p];
	else        return [alphaTridymite getEnthalpyFromT:t andP:p];
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getEntropyFromT:t andP:p];
	else        return [alphaTridymite getEntropyFromT:t andP:p];
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getHeatCapacityFromT:t andP:p];
	else        return [alphaTridymite getHeatCapacityFromT:t andP:p];
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getDcpDtFromT:t andP:p];
	else        return [alphaTridymite getDcpDtFromT:t andP:p];
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getVolumeFromT:t andP:p];
	else        return [alphaTridymite getVolumeFromT:t andP:p];
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getDvDtFromT:t andP:p];
	else        return [alphaTridymite getDvDtFromT:t andP:p];
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getDvDpFromT:t andP:p];
	else        return [alphaTridymite getDvDpFromT:t andP:p];
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getD2vDt2FromT:t andP:p];
	else        return [alphaTridymite getD2vDt2FromT:t andP:p];
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getD2vDtDpFromT:t andP:p];
	else        return [alphaTridymite getD2vDtDpFromT:t andP:p];
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	if (t > Tt) return [betaTridymite getD2vDp2FromT:t andP:p];
	else        return [alphaTridymite getD2vDp2FromT:t andP:p];
}

@end

@implementation WhitlockiteBerman

-(id)init {
	if ((self = [super initWithH:-4097169.0
							  S:235.978
							 k0:402.997
							 k1:-28.0835E2
							 k2:0.0
							 k3:-32.6230E7
							 l1:2.5427E-2
							 l2:19.255E-5
							 Tt:1373.0
						 deltaH:14059.0
							 v0:9.7620
							 v1:0.0
							 v2:0.0
							 v3:0.0
							 v4:0.0])) {
		[self setPhaseFormula:@"Ca3(PO4)2"];
		[self setPhaseName:@"Whitlockite"];
	}
	return self;
}

@end

// Additional Berman Phases

@implementation AlbiteBerman

#define NS 2

-(id)init {
    if ((self = [super initWithH:-3921618.201
                               S:224.412
                              k0:393.63574
                              k1:-2415.498
                              k2:-7892826.0
                              k3:1070636032.0
                              v0:10.083
                              v1:-1.94469e-06
                              v2:4.8611e-13
                              v3:2.63072e-05
                              v4:3.2407e-09])) {
        [self setPhaseFormula:@"NaAlSi3O8"];
        [self setPhaseName:@"Albite"];
        tOld = -9999.0;
        pOld = -9999.0;
        for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
        gDis    = 0.0;
        hDis    = 0.0;
        sDis    = 0.0;
        cpDis   = 0.0;
        dcpdt   = 0.0;
        vDis    = 0.0;
        dvdt    = 0.0;
        dvdp    = 0.0;
        d2vdt2  = 0.0;
        d2vdtdp = 0.0;
        d2vdp2  = 0.0;
    }
    return self;
}

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

#define MAX_ITER 200

static const double a0   =     5.479;
static const double b	 =  6854.0;
static const double aod0 =    41.620;
static const double bod  = -9301.0;
static const double cod  = 43600.0;
static const double d0   =    -2.171;
static const double d1   =    -3.043;
static const double d2   =    -0.001569;
static const double d3   =     0.000002109;
static const double tc   =  1251.0;
static const double tod  =   824.1;

/* q = sOrd[0], qod = sOrd[1] */

#define S   -(0.5*a0*sOrd[0]*sOrd[0] + 0.5*aod0*sOrd[1]*sOrd[1] \
+ (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[0]*sOrd[1])
#define H   - 0.5*a0*tc*sOrd[0]*sOrd[0] + 0.25 *b*sOrd[0]*sOrd[0]*sOrd[0]*sOrd[0] \
- 0.5*aod0*tod*sOrd[1]*sOrd[1] + 0.25*bod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (cod/6.0)*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (d0 - d2*t*t - 2.0*d3*t*t*t)*sOrd[0]*sOrd[1]
#define V   (H)/335282.925
#define G   (H) - t*(S) + (p-1.0)*(V)

/*----------------------------------------------------------------------------*/

#define DGDS0 (-a0*tc*sOrd[0] + b*sOrd[0]*sOrd[0]*sOrd[0] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[1])*(1.0+(p-1.0)/335282.925) \
+ t*(a0*sOrd[0] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[1])
#define DGDS1 (-aod0*tod*sOrd[1] + bod*sOrd[1]*sOrd[1]*sOrd[1] \
+ cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[0])*(1.0+(p-1.0)/335282.925) \
+ t*(aod0*sOrd[1] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[0])
#define DGDT  -(S) + (-2.0*d2*t-6.0*d3*t*t)*sOrd[0]*sOrd[1]*(p-1.0)/335282.925
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDS0S0 (-a0*tc + 3.0*b*sOrd[0]*sOrd[0])*(1.0+(p-1.0)/335282.925) + t*a0
#define D2GDS0S1 (d0-d2*t*t-2.0*d3*t*t*t)*(1.0+(p-1.0)/335282.925) \
+ t*(d1+2.0*d2*t+3.0*d3*t*t)
#define D2GDS0DT -2.0*(d2*t+3.0*d3*t*t)*sOrd[1]*(p-1.0)/335282.925 \
+ a0*sOrd[0] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[1]
#define D2GDS0DP (-a0*tc*sOrd[0] + b*sOrd[0]*sOrd[0]*sOrd[0] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[1])/335282.925

#define D2GDS1S1 (-aod0*tod + 3.0*bod*sOrd[1]*sOrd[1] + 5.0*cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]) \
*(1.0+(p-1.0)/335282.925) + t*aod0
#define D2GDS1DT -2.0*(d2*t+3.0*d3*t*t)*sOrd[0]*(p-1.0)/335282.925 \
+ aod0*sOrd[1] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[0]
#define D2GDS1DP (-aod0*tod*sOrd[1] + bod*sOrd[1]*sOrd[1]*sOrd[1] \
+ cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[0])/335282.925

#define D2GDT2   -2.0*(d2+6.0*d3*t)*sOrd[0]*sOrd[1]*(p-1.0)/335282.925 \
+ (2.0*d2+6.0*d3*t)*sOrd[0]*sOrd[1]
#define D2GDTDP  -2.0*(d2*t+3.0*d3*t*t)*sOrd[0]*sOrd[1]/335282.925
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDS0S0S0 6.0*b*sOrd[0]*(1.0+(p-1.0)/335282.925)
#define D3GDS0S0S1 0.0
#define D3GDS0S0DT a0
#define D3GDS0S0DP (-a0*tc + 3.0*b*sOrd[0]*sOrd[0])/335282.925
#define D3GDS0S1S1 0.0
#define D3GDS0S1DT (-2.0*d2*t-6.0*d3*t*t)*(1.0+(p-1.0)/335282.925) \
+ d1 + 2.0*d2*t + 3.0*d3*t*t + t*(2.0*d2+6.0*d3*t)
#define D3GDS0S1DP (d0-d2*t*t-2.0*d3*t*t*t)/335282.925
#define D3GDS0DT2  -2.0*(d2+6.0*d3*t)*sOrd[1]*(p-1.0)/335282.925 \
+ (2.0*d2+6.0*d3*t)*sOrd[1]
#define D3GDS0DTDP -2.0*(d2*t+3.0*d3*t*t)*sOrd[1]/335282.925
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 (6.0*bod*sOrd[1] + 20.0*cod*sOrd[1]*sOrd[1]*sOrd[1]) \
*(1.0+(p-1.0)/335282.925)
#define D3GDS1S1DT aod0
#define D3GDS1S1DP (-aod0*tod + 3.0*bod*sOrd[1]*sOrd[1] \
+ 5.0*cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1])/335282.925
#define D3GDS1DT2  -2.0*(d2+6.0*d3*t)*sOrd[0]*(p-1.0)/335282.925 \
+ (2.0*d2+6.0*d3*t)*sOrd[0]
#define D3GDS1DTDP -2.0*(d2*t+3.0*d3*t*t)*sOrd[0]/335282.925
#define D3GDS1DP2  0.0

#define D3GDT3     -12.0*d3*sOrd[0]*sOrd[1]*(p-1.0)/335282.925 + 6.0*d3*sOrd[0]*sOrd[1]
#define D3GDT2DP   -2.0*(d2+6.0*d3*t)*sOrd[0]*sOrd[1]/335282.925
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define fillD2GDS2 \
d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT \
d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP;

#define fillD3GDS3 \
d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2;

/**
 Matrix inversion routine converted from Numerical recipies to have zero array indexing
 */

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

-(void)gaussj:(double [NS][NS])a {
    int indxc[NS], indxr[NS], ipiv[NS];
    int i, icol = -1, irow = -1, j, k, l,ll;
    double big, dum, pivinv, temp;

    for (j=0; j<NS; j++) ipiv[j]=0;
    for (i=0; i<NS; i++) {
        big=0.0;
        for (j=0; j<NS; j++)
            if (ipiv[j] != 1)
                for (k=0; k<NS; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) return;
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=0; l<NS; l++) SWAP(a[irow][l],a[icol][l])
                }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) return;
        pivinv = 1.0/a[icol][icol];
        a[icol][icol] = 1.0;
        for (l=0; l<NS; l++) a[icol][l] *= pivinv;
        for (ll=0; ll<NS; ll++)
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l=0; l<NS; l++) a[ll][l] -= a[icol][l]*dum;
            }
    }
    for (l=(NS-1); l>=0; l--) {
        if (indxr[l] != indxc[l])
            for (k=0; k<NS; k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
}

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives
 */

-(void)order:(double)t
           p:(double)p
           s:(double [NS])sOrd        // s[NS]                BINARY MASK: 0000000001
          dt:(double [NS])dt          // ds[NS]/dt            BINARY MASK: 0000000100
          dp:(double [NS])dp          // ds[NS]/dp            BINARY MASK: 0000001000
         dt2:(double [NS])dt2         // d2s[NS]/dt2          BINARY MASK: 0010000000
         dtp:(double [NS])dtp         // d2s[NS]/dtp          BINARY MASK: 0100000000
         dp2:(double [NS])dp2         // d2s[NS]/dp2          BINARY MASK: 1000000000
{
    int i, j, k, l, iter=0;
    double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
    d3gdsdt2[NS], temp[NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS];

    /* look-up or compute the current ordering state */
    if ( (t != tOld) || (p != pOld) ) {
        double dgds[NS], sNew[NS];

        for (i=0; i<NS; i++) sOld[i] = 2.0;

        sNew[0] = 0.6;
        sNew[1] = 0.9;

        while (   ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) || (fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON))
               && (iter < MAX_ITER)) {
            double sOrd[NS], deltaS[NS];

            for (i=0; i<NS; i++) sOrd[i] = sNew[i];

            dgds[0] = DGDS0;
            dgds[1] = DGDS1;

            invd2gds2[0][0] = D2GDS0S0;
            invd2gds2[0][1] = D2GDS0S1;
            invd2gds2[1][0] = invd2gds2[0][1];
            invd2gds2[1][1] = D2GDS1S1;

            for (i=0; i<NS; i++) sOld[i] = sOrd[i];

            [self gaussj:invd2gds2];

            for (i=0; i<NS; i++) {
                for(j=0; j<NS; j++) sOrd[i] += - invd2gds2[i][j]*dgds[j];
                deltaS[i] = sOrd[i] - sOld[i];
            }

            for (i=0; i<NS; i++) sNew[i] = sOrd[i];
            iter++;
        }
        tOld = t;
        pOld = p;
    }

    /* s */
    for (i=0; i<NS; i++) sOrd[i] = sOld[i];

    /* dsdt */
    fillD2GDSDT

    for (i=0; i<NS; i++) {
        dt[i] = 0.0;
        for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }

    /* dsdp */
    fillD2GDSDP

    for (i=0; i<NS; i++) {
        dp[i] = 0.0;
        for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }

    /* d2sdt2 */
    fillD2GDSDT
    fillD3GDS3
    fillD3GDS2DT
    fillD3GDSDT2

    for (i=0; i<NS; i++) {
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdt2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dt[j][k]*dt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dt[k]*dt[l];
            }
        }
        dt2[i] = 0.0;
        for (j=0; j<NS; j++) dt2[i] += - invd2gds2[i][j]*temp[j];
    }

    /* d2sdtdp */
    fillD2GDSDT
    fillD2GDSDP
    fillD3GDS3
    fillD3GDS2DT
    fillD3GDS2DP
    fillD3GDSDTDP

    for (i=0; i<NS; i++) {
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<NS; k++) {
                temp[j] += d3gds2dt[j][k]*dp[k] + d3gds2dp[j][k]*dt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dt[k]*dp[l];
            }
        }
        dtp[i] = 0.0;
        for (j=0; j<NS; j++) dtp[i] += - invd2gds2[i][j]*temp[j];
    }

    /* d2sdp2 */
    fillD2GDSDP
    fillD3GDS3
    fillD3GDS2DP
    fillD3GDSDP2

    for (i=0; i<NS; i++) {
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdp2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dp[j][k]*dp[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dp[k]*dp[l];
            }
        }
        dp2[i] = 0.0;
        for (j=0; j<NS; j++) dp2[i] += - invd2gds2[i][j]*temp[j];
    }
}

-(void)albiteSalje:(double)p
                 t:(double)t
{
    double sOrd[NS], dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS], d2sdp2[NS];

    [self order:t p:p s:sOrd dt:dsdt dp:dsdp dt2:d2sdt2 dtp:d2sdtdp dp2:d2sdp2];
    int skip = (fabs(sOrd[0]) < sqrt(DBL_EPSILON)) && (fabs(sOrd[1]) < sqrt(DBL_EPSILON));

    if (skip) {
        gDis    = 0.0;
        hDis    = 0.0;
        sDis    = 0.0;
        cpDis   = 0.0;
        dcpdt   = 0.0;
        vDis    = 0.0;
        dvdt    = 0.0;
        dvdp    = 0.0;
        d2vdt2  = 0.0;
        d2vdtdp = 0.0;
        d2vdp2  = 0.0;
        return;
    } else {
        int i, j, k;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d2gdt2, d2gdtdp, d2gdp2;
        double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gds2dp[NS][NS],
        d3gdsdt2[NS], d3gdsdtdp[NS], d3gdsdp2[NS], d3gdt3, d3gdt2dp, d3gdtdp2, d3gdp3;
        double temp;

        /* ---------- */
        gDis = (G);
        /* ---------- */
        hDis = gDis - t*(DGDT);
        /* ---------- */
        sDis = -(DGDT);
        /* ---------- */
        fillD2GDS2
        fillD2GDSDT
        d2gdt2  = D2GDT2;

        cpDis = d2gdt2;
        for (i=0; i<NS; i++) {
            cpDis += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) cpDis += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        temp = cpDis;
        cpDis *= -t;
        /* ---------- */
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2
        d3gdt3 = D3GDT3;

        dcpdt = d3gdt3;
        for (i=0; i<NS; i++) {
            dcpdt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
            for (j=0; j<NS; j++) {
                dcpdt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j] + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) dcpdt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
            }
        }
        dcpdt = -t*dcpdt - temp;
        /* ---------- */
        vDis = (V);
        /* ---------- */
        fillD2GDSDP
        d2gdtdp = D2GDTDP;

        dvdt = d2gdtdp;
        for (i=0; i<NS; i++) {
            dvdt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NS; j++) dvdt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }
        /* ---------- */
        d2gdp2 = D2GDP2;

        dvdp = d2gdp2;
        for (i=0; i<NS; i++) {
            dvdp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NS; j++) dvdp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }
        /* ---------- */
        fillD3GDS2DP
        fillD3GDSDTDP
        d3gdt2dp = D3GDT2DP;

        d2vdt2 = d3gdt2dp;
        for (i=0; i<NS; i++) {
            d2vdt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i] + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
            for (j=0; j<NS; j++) {
                d2vdt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) d2vdt2
                    += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
            }
        }
        /* ---------- */
        fillD3GDSDP2
        d3gdtdp2 = D3GDTDP2;

        d2vdtdp = d3gdtdp2;
        for (i=0; i<NS; i++) {
            d2vdtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i] + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
            for (j=0; j<NS; j++) {
                d2vdtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j] + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
                + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j] + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) d2vdtdp
                    += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
            }
        }
        /* ---------- */
        d3gdp3 = D3GDP3;

        d2vdp2 = d3gdp3;
        for (i=0; i<NS; i++) {
            d2vdp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
            for (j=0; j<NS; j++) {
                d2vdp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j] + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) d2vdp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
            }
        }
    }

}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getEnthalpyFromT:t andP:p] + hDis;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getEntropyFromT:t andP:p] + sDis;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getHeatCapacityFromT:t andP:p] + cpDis;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getDcpDtFromT:t andP:p] + dcpdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getVolumeFromT:t andP:p] + vDis;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getDvDtFromT:t andP:p] + dvdt;
}

-(double)getDvDpFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getDvDpFromT:t andP:p] + dvdp;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getD2vDt2FromT:t andP:p] + d2vdt2;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getD2vDtDpFromT:t andP:p] + d2vdtdp;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
    if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
    return [super getD2vDp2FromT:t andP:p] + d2vdp2;
}

#undef NS
#undef S
#undef H
#undef V
#undef G
#undef DGDS0
#undef DGDS1
#undef DGDT
#undef DGDP
#undef D2GDS0S0
#undef D2GDS0S1
#undef D2GDS0DT
#undef D2GDS0DP
#undef D2GDS1S1
#undef D2GDS1DT
#undef D2GDS1DP
#undef D2GDT2
#undef D2GDTDP
#undef D2GDP2
#undef D3GDS0S0S0
#undef D3GDS0S0S1
#undef D3GDS0S0DT
#undef D3GDS0S0DP
#undef D3GDS0S1S1
#undef D3GDS0S1DT
#undef D3GDS0S1DP
#undef D3GDS0DT2
#undef D3GDS0DTDP
#undef D3GDS0DP2
#undef D3GDS1S1S1
#undef D3GDS1S1DT
#undef D3GDS1S1DP
#undef D3GDS1DT2
#undef D3GDS1DTDP
#undef D3GDS1DP2
#undef D3GDT3
#undef D3GDT2DP
#undef D3GDTDP2
#undef D3GDP3
#undef fillD2GDS2
#undef fillD2GDSDT
#undef fillD2GDSDP
#undef fillD3GDS3
#undef fillD3GDS2DT
#undef fillD3GDS2DP
#undef fillD3GDSDT2
#undef fillD3GDSDTDP
#undef fillD3GDSDP2

@end

@implementation High_AlbiteBerman

-(id)init {
    if ((self = [super initWithH:-3921618.2
                               S:224.412
                              k0:393.63574
                              k1:-2415.498
                              k2:-7892826.0
                              k3:1070636032.0
                              v0:10.083
                              v1:-1.94469e-06
                              v2:4.8611e-13
                              v3:2.63072e-05
                              v4:3.2407e-09])) {
        [self setPhaseFormula:@"NaAlSi3O8"];
        [self setPhaseName:@"High_Albite"];
    }
    return self;
}

@end

@implementation Low_AlbiteBerman

-(id)init {
    if ((self = [super initWithH:-3935100.1
                               S:207.443
                              k0:393.63574
                              k1:-2415.498
                              k2:-7892826.0
                              k3:1070636032.0
                              v0:10.043
                              v1:-1.94469e-06
                              v2:4.8611e-13
                              v3:2.63072e-05
                              v4:3.2407e-09])) {
        [self setPhaseFormula:@"NaAlSi3O8"];
        [self setPhaseName:@"Low_Albite"];
    }
    return self;
}

@end

@implementation AlmandineBerman

-(id)init {
    if ((self = [super initWithH:-5265502.33
                               S:339.927
                              k0:573.96191
                              k1:-1483.127
                              k2:-29291968.0
                              k3:5022076928.0
                              v0:11.511
                              v1:-5.57783e-07
                              v2:3.211e-14
                              v3:1.86125e-05
                              v4:7.4539e-09])) {
        [self setPhaseFormula:@"Si3Fe3Al2O12"];
        [self setPhaseName:@"Almandine"];
    }
    return self;
}

@end

@implementation AnorthiteBerman

-(id)init {
    if ((self = [super initWithH:-4228730.47
                               S:200.186
                              k0:439.36938
                              k1:-3734.149
                              k2:0.0
                              k3:-317023232.0
                              v0:10.075
                              v1:-1.27243e-06
                              v2:3.1762e-13
                              v3:1.09181e-05
                              v4:4.1985e-09])) {
        [self setPhaseFormula:@"Al2CaSi2O8"];
        [self setPhaseName:@"Anorthite"];
    }
    return self;
}

@end

@implementation AnthophylliteBerman

-(id)init {
    if ((self = [super initWithH:-12069031.71
                               S:535.196
                              k0:1219.31348
                              k1:-5766.535
                              k2:-34766144.0
                              k3:4400898048.0
                              v0:26.56
                              v1:-1.25901e-06
                              v2:0.0
                              v3:2.70597e-05
                              v4:3.1325e-09])) {
        [self setPhaseFormula:@"Mg7Si8O24H2"];
        [self setPhaseName:@"Anthophyllite"];
    }
    return self;
}

@end

@implementation AntigoriteBerman

-(id)init {
    if ((self = [super initWithH:-71364156.29
                               S:3602.996
                              k0:7394.51172
                              k1:0.0
                              k2:-548363008.0
                              k3:87284121600.0
                              v0:174.246
                              v1:-1.97778e-06
                              v2:4.9442e-13
                              v3:2.49647e-05
                              v4:3.9444e-09])) {
        [self setPhaseFormula:@"Mg48Si34O99H62O48"];
        [self setPhaseName:@"Antigorite"];
    }
    return self;
}

@end

@implementation BruciteBerman

-(id)init {
    if ((self = [super initWithH:-925936.6
                               S:63.064
                              k0:136.83899
                              k1:-537.13
                              k2:-4361949.0
                              k3:552694528.0
                              v0:2.468
                              v1:-2.02285e-06
                              v2:6.7261e-13
                              v3:3.28538e-05
                              v4:1.094e-09])) {
        [self setPhaseFormula:@"MgO2H2"];
        [self setPhaseName:@"Brucite"];
    }
    return self;
}

@end

@implementation Ca_Al_PyroxeneBerman

-(id)init {
    if ((self = [super initWithH:-3298766.68
                               S:140.751
                              k0:310.69775
                              k1:-1671.627
                              k2:-7455263.0
                              k3:948781568.0
                              v0:6.356
                              v1:-8.70044e-07
                              v2:2.1712e-13
                              v3:2.22505e-05
                              v4:5.2863e-09])) {
        [self setPhaseFormula:@"CaAl2SiO6"];
        [self setPhaseName:@"Ca-Al_Pyroxene"];
    }
    return self;
}

@end

@implementation CalciteBerman

-(id)init {
    if ((self = [super initWithH:-1206819.09
                               S:91.725
                              k0:178.18748
                              k1:-1657.697
                              k2:-482722.0
                              k3:166604928.0
                              v0:3.69
                              v1:-1.40019e-06
                              v2:0.0
                              v3:8.90682e-06
                              v4:2.27402e-08])) {
        [self setPhaseFormula:@"CaCO3"];
        [self setPhaseName:@"Calcite"];
    }
    return self;
}

@end

@implementation ChrysotileBerman

-(id)init {
    if ((self = [super initWithH:-4363356.18
                               S:220.135
                              k0:610.02197
                              k1:-5581.152
                              k2:-1857300.0
                              k3:195472000.0
                              v0:10.72
                              v1:-1.81028e-06
                              v2:4.5243e-13
                              v3:2.71513e-05
                              v4:6.7351e-09])) {
        [self setPhaseFormula:@"Mg3Si2O9H4"];
        [self setPhaseName:@"Chrysotile"];
    }
    return self;
}

@end

@implementation ClinochloreBerman

-(id)init {
    if ((self = [super initWithH:-8909589.84
                               S:435.154
                              k0:1214.28442
                              k1:-11217.129
                              k2:0.0
                              k3:-1256253184.0
                              v0:21.147
                              v1:-1.81947e-06
                              v2:0.0
                              v3:2.64515e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"Mg5Al2Si3O18H8"];
        [self setPhaseName:@"Clinochlore"];
    }
    return self;
}

@end

@implementation CordieriteBerman

-(id)init {
    if ((self = [super initWithH:-9158726.56
                               S:417.97
                              k0:954.3865
                              k1:-7962.274
                              k2:-2317258.0
                              k3:-370214090.0
                              v0:23.311
                              v1:-1.15825e-06
                              v2:0.0
                              v3:3.00287e-06
                              v4:1.8017e-09])) {
        [self setPhaseFormula:@"Mg2Al4Si5O18"];
        [self setPhaseName:@"Cordierite"];
    }
    return self;
}

@end

@implementation DiasporeBerman

-(id)init {
    if ((self = [super initWithH:-999377.89
                               S:35.308
                              k0:143.23953
                              k1:-1540.409
                              k2:-323098.0
                              k3:64632256.0
                              v0:1.776
                              v1:-5.99292e-07
                              v2:0.0
                              v3:2.9718e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"AlO2H"];
        [self setPhaseName:@"Diaspore"];
    }
    return self;
}

@end

@implementation DiopsideBerman

-(id)init {
    if ((self = [super initWithH:-3200583.5
                               S:142.5
                              k0:305.41333
                              k1:-1604.931
                              k2:-7165973.0
                              k3:921837568.0
                              v0:6.62
                              v1:-8.72477e-07
                              v2:1.7069e-13
                              v3:2.77952e-05
                              v4:8.3082e-09])) {
        [self setPhaseFormula:@"MgCaSi2O6"];
        [self setPhaseName:@"Diopside"];
    }
    return self;
}

@end

@implementation DolomiteBerman

-(id)init {
    if ((self = [super initWithH:-2325248.43
                               S:154.89
                              k0:328.48267
                              k1:-2554.425
                              k2:-4688477.0
                              k3:790382336.0
                              v0:6.432
                              v1:-1.07019e-06
                              v2:0.0
                              v3:1.42302e-05
                              v4:3.63778e-08])) {
        [self setPhaseFormula:@"MgCaC2O6"];
        [self setPhaseName:@"Dolomite"];
        d0 = -9.42;
        d1 = 0.0;
        d2 = 3.85e5;
        d3 = 1.732e-2;
        d4 = 5.0200e-6;
        d5 = 0.0;
    }
    return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kDolomiteBerman_d0 = @"d0";
static NSString *kDolomiteBerman_d1 = @"d1";
static NSString *kDolomiteBerman_d2 = @"d2";
static NSString *kDolomiteBerman_d3 = @"d3";
static NSString *kDolomiteBerman_d4 = @"d4";
static NSString *kDolomiteBerman_d5 = @"d5";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        d0 = [aDecoder decodeDoubleForKey:kDolomiteBerman_d0];
        d1 = [aDecoder decodeDoubleForKey:kDolomiteBerman_d1];
        d2 = [aDecoder decodeDoubleForKey:kDolomiteBerman_d2];
        d3 = [aDecoder decodeDoubleForKey:kDolomiteBerman_d3];
        d4 = [aDecoder decodeDoubleForKey:kDolomiteBerman_d4];
        d5 = [aDecoder decodeDoubleForKey:kDolomiteBerman_d5];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:d0 forKey:kDolomiteBerman_d0];
        [aCoder encodeDouble:d1 forKey:kDolomiteBerman_d1];
        [aCoder encodeDouble:d2 forKey:kDolomiteBerman_d2];
        [aCoder encodeDouble:d3 forKey:kDolomiteBerman_d3];
        [aCoder encodeDouble:d4 forKey:kDolomiteBerman_d4];
        [aCoder encodeDouble:d5 forKey:kDolomiteBerman_d5];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

/*
 -(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double td = (t < 1423.0) ? t : 1423.0;
	double dhdis = d0*(td-298.0) + 2.0*d1*(sqrt(td)-sqrt(298.0)) - d2*(1.0/td-1.0/298.0) + d3*(td*td-298.0*298.0)/2.0 + d4*(td*td*td-298.0*298.0*298.0)/3.0;
	double dsdis = d0*(log(td)-log(298.0)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(298.0)) - d2*(1.0/(td*td)-1.0/(298.0*298.0))/2.0 + d3*(td-298.0) + d4*((td*td)-(298.0*298.0))/2.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	return [super getGibbsFreeEnergyFromT:t andP:p] + dhdis - t*dsdis + dvdis*(p-pr);
 }
 */
-(double)getEnthalpyFromT:(double)t andP:(double)p {
    double td = (t < 1423.0) ? t : 1423.0;
    double dhdis = d0*(td-298.0) + 2.0*d1*(sqrt(td)-sqrt(298.0)) - d2*(1.0/td-1.0/298.0) + d3*(td*td-298.0*298.0)/2.0 + d4*(td*td*td-298.0*298.0*298.0)/3.0;
    double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
    dhdis += (t < 1423.0 && d5 != 0.0) ? -(p-pr)*t*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
    return [super getEnthalpyFromT:t andP:p] + dhdis + dvdis*(p-pr);
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    double td = (t < 1423.0) ? t : 1423.0;
    double dsdis = d0*(log(td)-log(298.0)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(298.0)) - d2*(1.0/(td*td)-1.0/(298.0*298.0))/2.0 + d3*(td-298.0) + d4*((td*td)-(298.0*298.0))/2.0;
    dsdis += (t < 1423.0 && d5 != 0.0) ? -(p-pr)*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
    return [super getEntropyFromT:t andP:p] + dsdis;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
    double cpdis = (t < 1423.0) ? d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t : 0.0;
    cpdis += (t < 1423.0 && d5 != 0.0) ? -(p-pr)*t*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
    return [super getHeatCapacityFromT:t andP:p] + cpdis;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
    double dcpdisdt = (t < 1423.0) ? -d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t : 0.0;
    dcpdisdt += (t < 1423.0 && d5 != 0.0) ? -(p-pr)*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5
    -(p-pr)*t*(1.5*0.5*d1/pow(t,(double) 2.5) + 6.0*d2/(t*t*t*t) + 2.0*d4)/d5 : 0.0;
    return [super getDcpDtFromT:t andP:p] + dcpdisdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
    double td = (t < 1423.0) ? t : 1423.0;
    double dhdis = d0*(td-298.0) + 2.0*d1*(sqrt(td)-sqrt(298.0)) - d2*(1.0/td-1.0/298.0) + d3*(td*td-298.0*298.0)/2.0 + d4*(td*td*td-298.0*298.0*298.0)/3.0;
    double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
    return [super getVolumeFromT:t andP:p] + dvdis;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
    double dvdisdt = (t < 1423.0 && d5 != 0.0) ? (d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*(t*t))/d5 : 0.0;
    return [super getDvDtFromT:t andP:p] + dvdisdt;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
    double d2vdisdt2 = (t < 1423.0 && d5 != 0.0) ? (-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
    return [super getD2vDt2FromT:t andP:p] + d2vdisdt2;
}

@end

@implementation ClinoenstatiteBerman

-(id)init {
    if ((self = [super initWithH:-1545926.25
                               S:66.325
                              k0:139.95824
                              k1:-497.034
                              k2:-4400237.0
                              k3:535708928.0
                              v0:3.131
                              v1:-7.496e-07
                              v2:4.48e-14
                              v3:2.1915e-05
                              v4:7.492e-09])) {
        [self setPhaseFormula:@"MgSiO3"];
        [self setPhaseName:@"Clinoenstatite"];
    }
    return self;
}

@end

@implementation OrthoenstatiteBerman

-(id)init {
    if ((self = [super initWithH:-1545552.25
                               S:66.17
                              k0:166.5795
                              k1:-1200.588
                              k2:-2270560.0
                              k3:279150336.0
                              v0:3.133
                              v1:-7.49346e-07
                              v2:4.467e-14
                              v3:2.46558e-05
                              v4:7.467e-09])) {
        [self setPhaseFormula:@"MgSiO3"];
        [self setPhaseName:@"Orthoenstatite"];
    }
    return self;
}

@end

@implementation ProtoenstatiteBerman

-(id)init {
    if ((self = [super initWithH:-1543958.64
                               S:67.438
                              k0:166.5795
                              k1:-1200.588
                              k2:-2270560.0
                              k3:279150336.0
                              v0:3.242
                              v1:-7.496e-07
                              v2:4.48e-14
                              v3:1.68317e-05
                              v4:1.1665e-08])) {
        [self setPhaseFormula:@"MgSiO3"];
        [self setPhaseName:@"Protoenstatite"];
    }
    return self;
}

@end

@implementation FerrosiliteBerman

-(id)init {
    if ((self = [super initWithH:-1194374.61
                               S:95.882
                              k0:169.06058
                              k1:-1192.971
                              k2:-2097138.0
                              k3:292532960.0
                              v0:3.296
                              v1:-9.89988e-07
                              v2:0.0
                              v3:3.18076e-05
                              v4:7.585e-09])) {
        [self setPhaseFormula:@"SiFeO3"];
        [self setPhaseName:@"Ferrosilite"];
    }
    return self;
}

@end

@implementation GrossularBerman

-(id)init {
    if ((self = [super initWithH:-6632859.38
                               S:255.15
                              k0:573.43042
                              k1:-2039.405
                              k2:-18887168.0
                              k3:2319311872.0
                              v0:12.538
                              v1:-6.53914e-07
                              v2:1.635e-13
                              v3:1.89942e-05
                              v4:7.9756e-09])) {
        [self setPhaseFormula:@"Ca3Al2Si3O12"];
        [self setPhaseName:@"Grossular"];
    }
    return self;
}

@end

@implementation JadeiteBerman

-(id)init {
    if ((self = [super initWithH:-3025118.24
                               S:133.574
                              k0:311.29297
                              k1:-2005.101
                              k2:-5350264.0
                              k3:662566912.0
                              v0:6.034
                              v1:-8.59802e-07
                              v2:2.1488e-13
                              v3:2.31184e-05
                              v4:2.5785e-09])) {
        [self setPhaseFormula:@"NaAlSi2O6"];
        [self setPhaseName:@"Jadeite"];
    }
    return self;
}

@end

@implementation KaoliniteBerman

-(id)init {
    if ((self = [super initWithH:-4120327.22
                               S:203.7
                              k0:523.23291
                              k1:-4426.668
                              k2:-2244259.0
                              k3:92308768.0
                              v0:9.952
                              v1:-1.2e-06
                              v2:0.0
                              v3:3.2e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"Al2Si2O9H4"];
        [self setPhaseName:@"Kaolinite"];
    }
    return self;
}

@end

@implementation LawsoniteBerman

-(id)init {
    if ((self = [super initWithH:-4865665.95
                               S:229.176
                              k0:728.6731
                              k1:-8248.055
                              k2:0.0
                              k3:850556416.0
                              v0:10.144
                              v1:-7.68809e-07
                              v2:1.9223e-13
                              v3:2.62834e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"CaAl2Si2O10H4"];
        [self setPhaseName:@"Lawsonite"];
    }
    return self;
}

@end

@implementation MagnesiteBerman

-(id)init {
    if ((self = [super initWithH:-1114505.11
                               S:65.09
                              k0:162.29663
                              k1:-1109.267
                              k2:-4882646.0
                              k3:874661888.0
                              v0:2.803
                              v1:-8.90086e-07
                              v2:1.533e-13
                              v3:1.84362e-05
                              v4:4.15968e-08])) {
        [self setPhaseFormula:@"MgCO3"];
        [self setPhaseName:@"Magnesite"];
    }
    return self;
}

@end

@implementation MargariteBerman

-(id)init {
    if ((self = [super initWithH:-6236602.96
                               S:265.084
                              k0:699.79565
                              k1:-5587.105
                              k2:-6807750.0
                              k3:734321664.0
                              v0:12.958
                              v1:-1.1549e-06
                              v2:2.8862e-13
                              v3:2.10185e-05
                              v4:1.24556e-08])) {
        [self setPhaseFormula:@"CaAl4Si2O12H2"];
        [self setPhaseName:@"Margarite"];
    }
    return self;
}

@end

@implementation MeioniteBerman

-(id)init {
    if ((self = [super initWithH:-13849722.94
                               S:730.0
                              k0:1511.34497
                              k1:-13243.328
                              k2:0.0
                              k3:-751613952.0
                              v0:34.036
                              v1:-1.11027e-06
                              v2:0.0
                              v3:9.34e-06
                              v4:0.0])) {
        [self setPhaseFormula:@"Ca4Al6Si6O27C"];
        [self setPhaseName:@"Meionite"];
    }
    return self;
}

@end

@implementation MerwiniteBerman

-(id)init {
    if ((self = [super initWithH:-4537497.44
                               S:251.777
                              k0:453.61963
                              k1:-3250.015
                              k2:0.0
                              k3:-344227072.0
                              v0:9.847
                              v1:-5.50543e-07
                              v2:1.3811e-13
                              v3:2.93756e-05
                              v4:8.7235e-09])) {
        [self setPhaseFormula:@"Ca3MgSi2O8"];
        [self setPhaseName:@"Merwinite"];
    }
    return self;
}

@end

@implementation MonticelliteBerman

-(id)init {
    if ((self = [super initWithH:-2250027.48
                               S:108.3
                              k0:226.34225
                              k1:-1542.741
                              k2:-1179709.0
                              k3:-23285888.0
                              v0:5.148
                              v1:-9.03655e-07
                              v2:2e-13
                              v3:2.78632e-05
                              v4:7.6339e-09])) {
        [self setPhaseFormula:@"CaMgSiO4"];
        [self setPhaseName:@"Monticellite"];
    }
    return self;
}

@end

@implementation ParagoniteBerman

-(id)init {
    if ((self = [super initWithH:-5944208.14
                               S:277.699
                              k0:577.56934
                              k1:-1472.801
                              k2:-32214368.0
                              k3:5050077184.0
                              v0:13.216
                              v1:-1.97323e-06
                              v2:4.9335e-13
                              v3:3.94241e-05
                              v4:5.9701e-09])) {
        [self setPhaseFormula:@"NaAl3Si3O12H2"];
        [self setPhaseName:@"Paragonite"];
    }
    return self;
}

@end

@implementation Potassium_FeldsparBerman

-(id)init {
    if ((self = [super initWithH:-3970790.781
                               S:214.145
                              k0:381.37231
                              k1:-1941.045
                              k2:-12037252.0
                              k3:1836425472.0
                              v0:10.869
                              v1:-1.80451e-06
                              v2:5.112e-13
                              v3:1.51451e-05
                              v4:5.485e-09])) {
        [self setPhaseFormula:@"KAlSi3O8"];
        [self setPhaseName:@"Potassium_Feldspar"];
        d0 = 282.98;
        d1 = -4.83e3;
        d2 = 36.21e5;
        d3 = -15.733e-2;
        d4 = 34.770e-6;
        d5 = 41.063e4;
    }
    return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kPotassiumFeldsparBerman_d0 = @"d0";
static NSString *kPotassiumFeldsparBerman_d1 = @"d1";
static NSString *kPotassiumFeldsparBerman_d2 = @"d2";
static NSString *kPotassiumFeldsparBerman_d3 = @"d3";
static NSString *kPotassiumFeldsparBerman_d4 = @"d4";
static NSString *kPotassiumFeldsparBerman_d5 = @"d5";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        d0 = [aDecoder decodeDoubleForKey:kPotassiumFeldsparBerman_d0];
        d1 = [aDecoder decodeDoubleForKey:kPotassiumFeldsparBerman_d1];
        d2 = [aDecoder decodeDoubleForKey:kPotassiumFeldsparBerman_d2];
        d3 = [aDecoder decodeDoubleForKey:kPotassiumFeldsparBerman_d3];
        d4 = [aDecoder decodeDoubleForKey:kPotassiumFeldsparBerman_d4];
        d5 = [aDecoder decodeDoubleForKey:kPotassiumFeldsparBerman_d5];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:d0 forKey:kPotassiumFeldsparBerman_d0];
        [aCoder encodeDouble:d1 forKey:kPotassiumFeldsparBerman_d1];
        [aCoder encodeDouble:d2 forKey:kPotassiumFeldsparBerman_d2];
        [aCoder encodeDouble:d3 forKey:kPotassiumFeldsparBerman_d3];
        [aCoder encodeDouble:d4 forKey:kPotassiumFeldsparBerman_d4];
        [aCoder encodeDouble:d5 forKey:kPotassiumFeldsparBerman_d5];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", self.className];
}

#pragma mark -
#pragma mark instance methods

/*
 -(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	double td = (t < 1436.0) ? t : 1436.0;
	double dhdis = d0*(td-298.0) + 2.0*d1*(sqrt(td)-sqrt(298.0)) - d2*(1.0/td-1.0/298.0) + d3*(td*td-298.0*298.0)/2.0 + d4*(td*td*td-298.0*298.0*298.0)/3.0;
	double dsdis = d0*(log(td)-log(298.0)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(298.0)) - d2*(1.0/(td*td)-1.0/(298.0*298.0))/2.0 + d3*(td-298.0) + d4*((td*td)-(298.0*298.0))/2.0;
	double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
	return [super getGibbsFreeEnergyFromT:t andP:p] + dhdis - t*dsdis + dvdis*(p-pr);
 }
 */
-(double)getEnthalpyFromT:(double)t andP:(double)p {
    double td = (t < 1436.0) ? t : 1436.0;
    double dhdis = d0*(td-298.0) + 2.0*d1*(sqrt(td)-sqrt(298.0)) - d2*(1.0/td-1.0/298.0) + d3*(td*td-298.0*298.0)/2.0 + d4*(td*td*td-298.0*298.0*298.0)/3.0;
    double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
    dhdis += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
    return [super getEnthalpyFromT:t andP:p] + dhdis + dvdis*(p-pr);
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    double td = (t < 1436.0) ? t : 1436.0;
    double dsdis = d0*(log(td)-log(298.0)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(298.0)) - d2*(1.0/(td*td)-1.0/(298.0*298.0))/2.0 + d3*(td-298.0) + d4*((td*td)-(298.0*298.0))/2.0;
    dsdis += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*(d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t)/d5 : 0.0;
    return [super getEntropyFromT:t andP:p] + dsdis;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
    double cpdis = (t < 1436.0) ? d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*t*t : 0.0;
    cpdis += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
    return [super getHeatCapacityFromT:t andP:p] + cpdis;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
    double dcpdisdt = (t < 1436.0) ? -d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t : 0.0;
    dcpdisdt += (t < 1436.0 && d5 != 0.0) ? -(p-pr)*(-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5
    -(p-pr)*t*(1.5*0.5*d1/pow(t,(double) 2.5) + 6.0*d2/(t*t*t*t) + 2.0*d4)/d5 : 0.0;
    return [super getDcpDtFromT:t andP:p] + dcpdisdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
    double td = (t < 1436.0) ? t : 1436.0;
    double dhdis = d0*(td-298.0) + 2.0*d1*(sqrt(td)-sqrt(298.0)) - d2*(1.0/td-1.0/298.0) + d3*(td*td-298.0*298.0)/2.0 + d4*(td*td*td-298.0*298.0*298.0)/3.0;
    double dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;
    return [super getVolumeFromT:t andP:p] + dvdis;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
    double dvdisdt = (t < 1436.0 && d5 != 0.0) ? (d0 + d1/sqrt(t) + d2/(t*t) + d3*t + d4*(t*t))/d5 : 0.0;
    return [super getDvDtFromT:t andP:p] + dvdisdt;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
    double d2vdisdt2 = (t < 1436.0 && d5 != 0.0) ? (-d1*0.5/pow(t, (double) 1.5) - 2.0*d2/(t*t*t) + d3 + 2.0*d4*t)/d5 : 0.0;
    return [super getD2vDt2FromT:t andP:p] + d2vdisdt2;
}

@end

@implementation MicroclineBerman

-(id)init {
    if ((self = [super initWithH:-3970790.78
                               S:214.145
                              k0:381.37231
                              k1:-1941.045
                              k2:-12037252.0
                              k3:1836425472.0
                              v0:10.869
                              v1:-1.80451e-06
                              v2:5.112e-13
                              v3:1.51451e-05
                              v4:5.485e-09])) {
        [self setPhaseFormula:@"KAlSi3O8"];
        [self setPhaseName:@"Microcline"];
    }
    return self;
}

@end

@implementation PrehniteBerman

-(id)init {
    if ((self = [super initWithH:-6198606.33
                               S:288.634
                              k0:716.05176
                              k1:-6404.609
                              k2:-2182456.0
                              k3:268502784.0
                              v0:14.016
                              v1:-1.42694e-06
                              v2:0.0
                              v3:1.46841e-06
                              v4:1.12835e-07])) {
        [self setPhaseFormula:@"Ca2Al2Si3O12H2"];
        [self setPhaseName:@"Prehnite"];
    }
    return self;
}

@end

@implementation PyropeBerman

-(id)init {
    if ((self = [super initWithH:-6286547.62
                               S:266.359
                              k0:640.71997
                              k1:-4542.07
                              k2:-4701900.0
                              k3:0.0
                              v0:11.316
                              v1:-5.76209e-07
                              v2:4.42e-14
                              v3:2.25187e-05
                              v4:3.7044e-09])) {
        [self setPhaseFormula:@"Mg3Al2Si3O12"];
        [self setPhaseName:@"Pyrope"];
    }
    return self;
}

@end

@implementation PyrophylliteBerman

-(id)init {
    if ((self = [super initWithH:-5640780.66
                               S:239.4
                              k0:665.93066
                              k1:-5897.422
                              k2:-4979913.0
                              k3:661808896.0
                              v0:12.76
                              v1:-1.35449e-06
                              v2:0.0
                              v3:1.26374e-05
                              v4:3.81661e-08])) {
        [self setPhaseFormula:@"Al2Si4O12H2"];
        [self setPhaseName:@"Pyrophyllite"];
    }
    return self;
}

@end

@implementation Mg_Al_SpinelBerman

-(id)init {
    if ((self = [super initWithH:-2300312.74
                               S:84.535
                              k0:235.89993
                              k1:-1766.578
                              k2:-1710415.0
                              k3:40616928.0
                              v0:3.977
                              v1:-4.89417e-07
                              v2:0.0
                              v3:2.16914e-05
                              v4:5.0528e-09])) {
        [self setPhaseFormula:@"MgAl2O4"];
        [self setPhaseName:@"Mg_Al_Spinel"];
    }
    return self;
}

@end

@implementation TalcBerman

-(id)init {
    if ((self = [super initWithH:-5897386.72
                               S:261.24
                              k0:664.10522
                              k1:-5187.172
                              k2:-2147218.0
                              k3:-327371776.0
                              v0:13.61
                              v1:-1.6989e-06
                              v2:5.665e-13
                              v3:2.94468e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"Mg3Si4O12H2"];
        [self setPhaseName:@"Talc"];
    }
    return self;
}

@end

@implementation TremoliteBerman

-(id)init {
    if ((self = [super initWithH:-12305578.13
                               S:551.15
                              k0:1229.36007
                              k1:-6401.9
                              k2:-32089890.0
                              k3:4208807760.0
                              v0:27.268
                              v1:-1.39177e-06
                              v2:3.4809e-13
                              v3:2.43739e-05
                              v4:9.8338e-09])) {
        [self setPhaseFormula:@"Ca2Mg5Si8O24H2"];
        [self setPhaseName:@"Tremolite"];
    }
    return self;
}

@end

@implementation WollastoniteBerman

-(id)init {
    if ((self = [super initWithH:-1631500.0
                               S:81.81
                              k0:149.07266
                              k1:-690.295
                              k2:-3659348.0
                              k3:484349440.0
                              v0:3.983
                              v1:-1.24492e-06
                              v2:3.1132e-13
                              v3:2.818e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"CaSiO3"];
        [self setPhaseName:@"Wollastonite"];
    }
    return self;
}

@end

@implementation PseudowollastoniteBerman

-(id)init {
    if ((self = [super initWithH:-1627427.11
                               S:85.279
                              k0:141.15611
                              k1:-417.232
                              k2:-5857595.0
                              k3:940734976.0
                              v0:4.016
                              v1:-1.24492e-06
                              v2:3.1132e-13
                              v3:2.818e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"CaSiO3"];
        [self setPhaseName:@"Pseudowollastonite"];
    }
    return self;
}

@end

@implementation ZoisiteBerman

-(id)init {
    if ((self = [super initWithH:-6889488.28
                               S:297.576
                              k0:749.17041
                              k1:-6509.281
                              k2:-2380525.0
                              k3:124858368.0
                              v0:13.588
                              v1:-5.1516e-07
                              v2:1.2879e-13
                              v3:3.46696e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"Ca2Al3Si3O13H"];
        [self setPhaseName:@"Zoisite"];
    }
    return self;
}

@end

@implementation ClinozoisiteBerman

-(id)init {
    if ((self = [super initWithH:-6894967.74
                               S:287.076
                              k0:749.17037
                              k1:-6509.283
                              k2:-2380525.0
                              k3:124858370.0
                              v0:13.673
                              v1:-5.1516e-07
                              v2:1.2879e-13
                              v3:3.46696e-05
                              v4:0.0])) {
        [self setPhaseFormula:@"Ca2Al3Si3O13H"];
        [self setPhaseName:@"Clinozoisite"];
    }
    return self;
}

@end

@implementation Oxygen_GasBerman

-(id)init {
    if ((self = [super initWithH:205.033
                               S:0.0
                              k0:23.10248
                              k1:804.888
                              k2:1762835.0
                              k3:0.0
                              v0:0.0
                              v1:0.0
                              v2:0.0
                              v3:0.0
                              v4:0.0])) {
        [self setPhaseFormula:@"O2"];
        [self setPhaseName:@"Oxygen_Gas"];
    }
    return self;
}

@end

@implementation Sulfur_GasBerman

-(id)init {
    if ((self = [super initWithH:228.07
                               S:0.0
                              k0:36.4845
                              k1:0.0
                              k2:-376560.0
                              k3:0.0
                              v0:0.0
                              v1:0.0
                              v2:0.0
                              v3:0.0
                              v4:0.0])) {
        [self setPhaseFormula:@"S2"];
        [self setPhaseName:@"Sulfur_Gas"];
    }
    return self;
}

@end

@implementation Hydrogen_GasBerman

-(id)init {
    if ((self = [super initWithH:130.574
                               S:0.0
                              k0:27.2797
                              k1:0.0
                              k2:-50208.0
                              k3:0.0
                              v0:0.0
                              v1:0.0
                              v2:0.0
                              v3:0.0
                              v4:0.0])) {
        [self setPhaseFormula:@"H2"];
        [self setPhaseName:@"Hydrogen_Gas"];
    }
    return self;
}

@end
