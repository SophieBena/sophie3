//
//  StixrudeProperties.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/10/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "StixrudeProperties.h"
#import "StixrudeEndmembers.h"
#import "DoubleVector.h"

@implementation StixrudeProperties

#pragma mark -
#pragma mark Initialization methods

-(id)initWithParameters:(NSArray *)paramsIn andLandauTerms:(NSArray *)landauIn {
	if ((self = [super init])) {
		parameters = [NSArray arrayWithArray:paramsIn];

		currentT = -9999.0;
		currentP = -9999.0;
		currentResults = nil;
		TC0 = 0.0;
		VD  = 0.0;
		SD  = 1.0;

		if (landauIn) {
			TC0 = [[landauIn objectAtIndex:0] doubleValue];
			VD  = [[landauIn objectAtIndex:1] doubleValue];
			SD  = [[landauIn objectAtIndex:2] doubleValue];
		}
	}
	return self;
}

-(id)initWithParameters:(NSArray *)paramsIn {
	return [self initWithParameters:paramsIn andLandauTerms:nil];
}

#pragma mark -
#pragma mark Overloaded PhaseBase method

-(void)setPhaseFormula:(NSString *)formulaString {
	[super setPhaseFormula:formulaString];
	double *composition = [[self formulaAsElementArray] pointerToDouble];
	numberOfFeAtoms = composition[26];
}

#pragma mark -
#pragma mark Support methods

-(void)doCalculationWithT:(double)t andWithP:(double)p {
	if ((t != currentT) || (p != currentP)) {
		currentResults = [StixrudeEndmembers calculateThermodynamicPropertiesOfEndmemberAtT:t andAtP:p withParameters:parameters];
		currentT = t;
		currentP = p;
		/*
		if ([[self phaseName] isEqualToString:@"Forsterite"]) {
			NSLog(@"%@, t = %g, p = %g", [self phaseName], t, p);
			NSLog(@"... a     = %20.13g", [[currentResults objectForKey:@"HelmholtzFreeEnergy"] doubleValue]);
			NSLog(@"... g     = %20.13g", [[currentResults objectForKey:@"GibbsFreeEnergy"] doubleValue]);
			NSLog(@"... h     = %20.13g", [[currentResults objectForKey:@"Enthalpy"] doubleValue]);
			NSLog(@"... s     = %20.13g", [[currentResults objectForKey:@"Entropy"] doubleValue]);
			NSLog(@"... Cp    = %20.13g", [[currentResults objectForKey:@"HeatCapacityAtConstantPressure"] doubleValue]);
			NSLog(@"... Cv    = %20.13g", [[currentResults objectForKey:@"HeatCapacityAtConstantVolume"] doubleValue]);
			NSLog(@"... v     = %20.13g", [[currentResults objectForKey:@"Volume"] doubleValue]);
			NSLog(@"... dvdt  = %20.13g", [[currentResults objectForKey:@"dVolume/dT"] doubleValue]);
			NSLog(@"... dvdp  = %20.13g", [[currentResults objectForKey:@"dVolume/dP"] doubleValue]);
			NSLog(@"... alpha = %20.13g", [[currentResults objectForKey:@"ThermalExpansionIsobaric"] doubleValue]);
			NSLog(@"... beta  = %20.13g", [[currentResults objectForKey:@"CompressibilityIsothermal"] doubleValue]);
		}
		*/
	}
}

#pragma mark -
#pragma mark Stoichiometric phase protocol methods

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double g = [[currentResults objectForKey:@"GibbsFreeEnergy"] doubleValue];

	g += - 8.3143*t*numberOfFeAtoms*log(5.0);

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return g;

	double Q2 = sqrt(1.0 - t/TC);
	return g + SD*( (t-TC)*Q2 + TC0*Q2*Q2*Q2/3.0);
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double h = [[currentResults objectForKey:@"Enthalpy"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return h;

	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	return h + SD*( (t-TC)*Q2 + TC0*Q2*Q2*Q2/3.0) - t*SD*(Q2 + (t-TC)*dQ2dt + TC0*Q2*Q2*dQ2dt);
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double s = [[currentResults objectForKey:@"Entropy"] doubleValue];

	s += 8.3143*numberOfFeAtoms*log(5.0);

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return s;

	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	return s - SD*(Q2 + (t-TC)*dQ2dt + TC0*Q2*Q2*dQ2dt);
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double cp =  [[currentResults objectForKey:@"HeatCapacityAtConstantPressure"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return cp;

	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	double d2Q2dt2 = (1.0/TC)*dQ2dt/Q2/Q2/2.0;
	return cp - t*SD*(dQ2dt + dQ2dt + (t-TC)*d2Q2dt2 + 2.0*TC0*Q2*dQ2dt*dQ2dt + TC0*Q2*Q2*d2Q2dt2);
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	double dcpdt = [[currentResults objectForKey:@"dCp/dT"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return dcpdt;

	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	double d2Q2dt2 = (1.0/TC)*dQ2dt/Q2/Q2/2.0;
	double d3Q2dt3 = - 2.0*(1.0/TC)*dQ2dt*dQ2dt/Q2/Q2/Q2/2.0 + (1.0/TC)*d2Q2dt2/Q2/Q2/2.0;
	return dcpdt -   SD*(2.0*dQ2dt + (t-TC)*d2Q2dt2 + 2.0*TC0*Q2*dQ2dt*dQ2dt + TC0*Q2*Q2*d2Q2dt2)
	             - t*SD*(3.0*d2Q2dt2 + (t-TC)*d3Q2dt3 + 2.0*TC0*dQ2dt*dQ2dt*dQ2dt + 6.0*TC0*Q2*dQ2dt*d2Q2dt2 + TC0*Q2*Q2*d3Q2dt3);
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double v = [[currentResults objectForKey:@"Volume"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return v;

	double dTCdp = VD/SD;
	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dp = (t/TC/TC)*dTCdp/Q2/2.0;
	return v + SD*( (t-TC)*dQ2dp - dTCdp*Q2 + TC0*Q2*Q2*dQ2dp);
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double dvdt = [[currentResults objectForKey:@"dVolume/dT"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return dvdt;

	double dTCdp = VD/SD;
	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	double dQ2dp = (t/TC/TC)*dTCdp/Q2/2.0;
	double d2Q2dtdp = (1.0/TC/TC)*dTCdp/Q2/2.0 - (t/TC/TC)*dTCdp*dQ2dt/Q2/Q2/2.0;
	return dvdt + SD*( dQ2dp + (t-TC)*d2Q2dtdp - dTCdp*dQ2dt + 2.0*TC0*Q2*dQ2dt*dQ2dp + TC0*Q2*Q2*d2Q2dtdp);

}

-(double)getDvDpFromT:(double)t andP:(double)p {
	[self doCalculationWithT:t andWithP:p];
	double dvdp = [[currentResults objectForKey:@"dVolume/dP"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return dvdp;

	double dTCdp = VD/SD;
	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dp = (t/TC/TC)*dTCdp/Q2/2.0;
	double d2Q2dp2 = -2.0*(t/TC/TC/TC)*dTCdp*dTCdp/Q2/2.0 - (t/TC/TC)*dTCdp*dQ2dp/Q2/Q2/2.0;
	return dvdp + SD*( -dTCdp*dQ2dp + (t-TC)*d2Q2dp2 - dTCdp*dQ2dp + 2.0*TC0*Q2*dQ2dp*dQ2dp + TC0*Q2*Q2*d2Q2dp2);

}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	double d2vdt2 = [[currentResults objectForKey:@"d2Volume/dT2"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return d2vdt2;

	double dTCdp = VD/SD;
	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	double dQ2dp = (t/TC/TC)*dTCdp/Q2/2.0;
	double d2Q2dt2 = (1.0/TC)*dQ2dt/Q2/Q2/2.0;
	double d2Q2dtdp = (1.0/TC/TC)*dTCdp/Q2/2.0 - (t/TC/TC)*dTCdp*dQ2dt/Q2/Q2/2.0;
	double d3Q2dt2dp = - 2.0*(1.0/TC/TC)*dTCdp*dQ2dt/Q2/Q2/2.0  - (t/TC/TC)*dTCdp*d2Q2dt2/Q2/Q2/2.0 + 2.0*(t/TC/TC)*dTCdp*dQ2dt*dQ2dt/Q2/Q2/Q2/2.0;
	return d2vdt2 + SD*( 2.0*d2Q2dtdp + (t-TC)*d3Q2dt2dp - dTCdp*d2Q2dt2 + 2.0*TC0*dQ2dt*dQ2dt*dQ2dp + 2.0*TC0*Q2*d2Q2dt2*dQ2dp + 2.0*TC0*Q2*dQ2dt*d2Q2dtdp
						+ 2.0*TC0*Q2*dQ2dt*d2Q2dtdp + TC0*Q2*Q2*d3Q2dt2dp);

}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	double d2vdtdp = [[currentResults objectForKey:@"d2Volume/dTdP"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return d2vdtdp;

	double dTCdp = VD/SD;
	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dt = -(1.0/TC)/Q2/2.0;
	double dQ2dp = (t/TC/TC)*dTCdp/Q2/2.0;
	double d2Q2dtdp = (1.0/TC/TC)*dTCdp/Q2/2.0 - (t/TC/TC)*dTCdp*dQ2dt/Q2/Q2/2.0;
	double d2Q2dp2 = -2.0*(t/TC/TC/TC)*dTCdp*dTCdp/Q2/2.0 - (t/TC/TC)*dTCdp*dQ2dp/Q2/Q2/2.0;
	double d3Q2dtdp2 = -2.0*(1.0/TC/TC/TC)*dTCdp*dTCdp/Q2/2.0 + 2.0*(t/TC/TC/TC)*dTCdp*dTCdp*dQ2dt/Q2/Q2/2.0
	- (1.0/TC/TC)*dTCdp*dQ2dp/Q2/Q2/2.0 - (t/TC/TC)*dTCdp*d2Q2dtdp/Q2/Q2/2.0 + 2.0*(t/TC/TC)*dTCdp*dQ2dp*dQ2dt/Q2/Q2/Q2/2.0;
	return d2vdtdp + SD*( -dTCdp*d2Q2dtdp + d2Q2dp2 + (t-TC)*d3Q2dtdp2 - dTCdp*d2Q2dtdp + 2.0*TC0*dQ2dt*dQ2dp*dQ2dp + 4.0*TC0*Q2*dQ2dp*d2Q2dtdp
						 + 2.0*TC0*Q2*dQ2dt*d2Q2dp2 + TC0*Q2*Q2*d3Q2dtdp2);

}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	double d2vdp2 = [[currentResults objectForKey:@"d2Volume/dP2"] doubleValue];

	double TC = TC0 + VD*p/SD;
	if (t >= TC) return d2vdp2;

	double dTCdp = VD/SD;
	double Q2 = sqrt(1.0 - t/TC);
	double dQ2dp = (t/TC/TC)*dTCdp/Q2/2.0;
	double d2Q2dp2 = -2.0*(t/TC/TC/TC)*dTCdp*dTCdp/Q2/2.0 - (t/TC/TC)*dTCdp*dQ2dp/Q2/Q2/2.0;
	double d3Q2dp3 = 6.0*(t/TC/TC/TC/TC)*dTCdp*dTCdp*dTCdp/Q2/2.0 + 2.0*(t/TC/TC/TC)*dTCdp*dTCdp*dQ2dp/Q2/Q2/2.0
	+ 2.0*(t/TC/TC/TC)*dTCdp*dTCdp*dQ2dp/Q2/Q2/2.0 - (t/TC/TC)*dTCdp*d2Q2dp2/Q2/Q2/2.0 + 2.0*(t/TC/TC)*dTCdp*dQ2dp*dQ2dp/Q2/Q2/Q2/2.0;
	return d2vdp2 + SD*( -dTCdp*d2Q2dp2 - dTCdp*d2Q2dp2 + (t-TC)*d3Q2dp3 - dTCdp*d2Q2dp2 + 2.0*TC0*dQ2dp*dQ2dp*dQ2dp + 4.0*TC0*Q2*dQ2dp*d2Q2dp2
						+ 2.0*TC0*Q2*dQ2dp*d2Q2dp2 + TC0*Q2*Q2*d3Q2dp3);

}

@end
