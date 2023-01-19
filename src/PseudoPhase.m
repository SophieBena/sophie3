//
//  PseudoPhase.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 8/20/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "PseudoPhase.h"

@implementation PseudoPhase

- (id)init
{
    self = [super init];
    if (self) {
        g       = 0.0;
        h       = 0.0;
        s       = 0.0;
        cp      = 0.0;
        dcpdt   = 0.0;
        v       = 0.0;
        dvdt    = 0.0;
        dvdp    = 0.0;
        d2vdt2  = 0.0;
        d2vdtdp = 0.0;
        d2vdp2  = 0.0;
        gLast   = 1.23456789;
    }

    return self;
}

@synthesize g, h, s, cp, dcpdt, v, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2;
@synthesize gLast;

#pragma mark -
#pragma mark Stoichiometric phase protocol methods

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return g;
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return h;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	return s;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return cp;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	return dcpdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	return v;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	return dvdt;
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	return dvdp;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return d2vdt2;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return d2vdtdp;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return d2vdp2;
}

@end
