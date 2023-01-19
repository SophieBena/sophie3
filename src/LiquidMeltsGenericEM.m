//
//  LiquidMeltsGenericEM.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/18/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "LiquidMeltsGenericEM.h"


@implementation LiquidMeltsGenericEM

@synthesize vLiq, dvdtLiq, dvdpLiq, d2vdtdpLiq, d2vdp2Liq, tFusion, sFusion, cpLiq;

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

#pragma mark -
#pragma mark NSSecureCoding protocol methods

static NSString *kvLiq = @"vLiq";
static NSString *kdvdtLiq = @"dvdtLiq";
static NSString *kdvdpLiq = @"dvdpLiq";
static NSString *kd2vdtdpLiq = @"d2vdtdpLiq";
static NSString *kd2vdp2Liq = @"d2vdp2Liq";
static NSString *ktFusion = @"tFusion";
static NSString *ksFusion = @"sFusion";
static NSString *kcpLiq = @"cpLiq";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        vLiq = [aDecoder decodeDoubleForKey:kvLiq];
        dvdtLiq = [aDecoder decodeDoubleForKey:kdvdtLiq];
        dvdpLiq = [aDecoder decodeDoubleForKey:kdvdpLiq];
        d2vdtdpLiq = [aDecoder decodeDoubleForKey:kd2vdtdpLiq];
        d2vdp2Liq = [aDecoder decodeDoubleForKey:kd2vdp2Liq];
        tFusion = [aDecoder decodeDoubleForKey:ktFusion];
        sFusion = [aDecoder decodeDoubleForKey:ksFusion];
        cpLiq = [aDecoder decodeDoubleForKey:kcpLiq];

        hLiqAtTfusion = [super getEnthalpyFromT:tFusion andP:1.0];
        sLiqAtTfusion = [super getEntropyFromT:tFusion andP:1.0];
        gLiqAtTfusion = hLiqAtTfusion - tFusion*sLiqAtTfusion;
        hLiqAtTfusion += sFusion*tFusion;
        sLiqAtTfusion += sFusion;

        trl = 1673.15;
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:vLiq forKey:kvLiq];
        [aCoder encodeDouble:dvdtLiq forKey:kdvdtLiq];
        [aCoder encodeDouble:dvdpLiq forKey:kdvdpLiq];
        [aCoder encodeDouble:d2vdtdpLiq forKey:kd2vdtdpLiq];
        [aCoder encodeDouble:d2vdp2Liq forKey:kd2vdp2Liq];
        [aCoder encodeDouble:tFusion forKey:ktFusion];
        [aCoder encodeDouble:sFusion forKey:ksFusion];
        [aCoder encodeDouble:cpLiq forKey:kcpLiq];
    } // else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", [self className]];
}

#pragma mark -
#pragma mark Override superclass methods

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return hLiqAtTfusion + cpLiq*(t-tFusion) - t*(sLiqAtTfusion + cpLiq*log(t/tFusion))
	+ (vLiq + dvdtLiq*(t-trl))*(p-pr) + 0.5*(dvdpLiq + (t-trl)*d2vdtdpLiq)*(p*p-pr*pr)
	- (dvdpLiq + (t-trl)*d2vdtdpLiq)*pr*(p-pr)
	+ d2vdp2Liq*( (p*p*p-pr*pr*pr)/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 );
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return hLiqAtTfusion + cpLiq*(t-tFusion)
	+ (vLiq + dvdtLiq*(t-trl))*(p-pr)
	+ 0.5*(dvdpLiq + (t-trl)*d2vdtdpLiq)*(p*p-pr*pr)
	- (dvdpLiq + (t-trl)*d2vdtdpLiq)*pr*(p-pr)
	+ d2vdp2Liq*( (p*p*p-pr*pr*pr)/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 )
	- t*(dvdtLiq*(p-pr) + 0.5*d2vdtdpLiq*(p-pr)*(p-pr));
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	return sLiqAtTfusion + cpLiq*log(t/tFusion) - (dvdtLiq*(p-pr) + 0.5*d2vdtdpLiq*(p-pr)*(p-pr));
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return cpLiq;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	return 0.0;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	return vLiq + dvdtLiq*(t-trl) + (dvdpLiq + d2vdtdpLiq*(t-trl))*(p-pr) + d2vdp2Liq*(0.5*p*p - pr*(p-pr));
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	return dvdtLiq + d2vdtdpLiq*(p-pr);
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	return dvdpLiq + d2vdtdpLiq*(t-trl) + d2vdp2Liq*(p-pr);
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return 0.0;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return d2vdtdpLiq;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return d2vdp2Liq;
}

@end
