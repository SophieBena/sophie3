//
//  LiquidMeltsSiO2.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/18/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "LiquidMeltsSiO2.h"


@implementation LiquidMeltsSiO2

@synthesize tr, pr, trl, h0_sio2, s0_sio2, al_sio2, bl_sio2, cl_sio2, dl_sio2, tg_sio2, cp_sio2, vLiq, dvdtLiq, dvdpLiq, d2vdtdpLiq, d2vdp2Liq;

-(id)init {
	if ((self = [super init])) {
		tr = 298.15;
		pr = 1.0;
		trl = 1673.0;
		h0_sio2 = -901554.0;
		s0_sio2 = 48.475;
		al_sio2 = 127.200;
		bl_sio2 = -10.777e-3;
		cl_sio2 = 4.3127e5;
		dl_sio2 = -1463.8;
		tg_sio2 = 1480.0;
		cp_sio2 = 81.373;
		vLiq = 2.690;
		dvdtLiq = 0.0;
		dvdpLiq = -1.89e-5;
		d2vdtdpLiq = 1.3e-8;
		d2vdp2Liq = 3.6e-10;
        _gibbsFreeEnergyReferenceStateUsed = NO;
        self.phaseName = @"SiO2";
        self.phaseFormula = @"SiO2";
	}
	return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol methods

static NSString *ktr = @"tr";
static NSString *kpr = @"pr";
static NSString *ktrl = @"trl";
static NSString *kh0_sio2 = @"h0_sio2";
static NSString *ks0_sio2 = @"s0_sio2";
static NSString *kal_sio2 = @"al_sio2";
static NSString *kbl_sio2 = @"bl_sio2";
static NSString *kcl_sio2 = @"cl_sio2";
static NSString *kdl_sio2 = @"dl_sio2";
static NSString *ktg_sio2 = @"tg_sio2";
static NSString *kcp_sio2 = @"cp_sio2";
static NSString *kvLiq = @"vLiq";
static NSString *kdvdtLiq = @"dvdtLiq";
static NSString *kdvdpLiq = @"dvdpLiq";
static NSString *kd2vdtdpLiq = @"d2vdtdpLiq";
static NSString *kd2vdp2Liq = @"d2vdp2Liq";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        tr  = [aDecoder decodeDoubleForKey:ktr];
        pr = [aDecoder decodeDoubleForKey:kpr];
        trl = [aDecoder decodeDoubleForKey:ktrl];
        h0_sio2 = [aDecoder decodeDoubleForKey:kh0_sio2];
        s0_sio2 = [aDecoder decodeDoubleForKey:ks0_sio2];
        al_sio2 = [aDecoder decodeDoubleForKey:kal_sio2];
        bl_sio2 = [aDecoder decodeDoubleForKey:kbl_sio2];
        cl_sio2 = [aDecoder decodeDoubleForKey:kcl_sio2];
        dl_sio2 = [aDecoder decodeDoubleForKey:kdl_sio2];
        tg_sio2 = [aDecoder decodeDoubleForKey:ktg_sio2];
        cp_sio2 = [aDecoder decodeDoubleForKey:kcp_sio2];
        vLiq = [aDecoder decodeDoubleForKey:kvLiq];
        dvdtLiq = [aDecoder decodeDoubleForKey:kdvdtLiq];
        dvdpLiq = [aDecoder decodeDoubleForKey:kdvdpLiq];
        d2vdtdpLiq = [aDecoder decodeDoubleForKey:kd2vdtdpLiq];
        d2vdp2Liq = [aDecoder decodeDoubleForKey:kd2vdp2Liq];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:tr forKey:ktr];
        [aCoder encodeDouble:pr forKey:kpr];
        [aCoder encodeDouble:trl forKey:ktrl];
        [aCoder encodeDouble:h0_sio2 forKey:kh0_sio2];
        [aCoder encodeDouble:s0_sio2 forKey:ks0_sio2];
        [aCoder encodeDouble:al_sio2 forKey:kal_sio2];
        [aCoder encodeDouble:bl_sio2 forKey:kbl_sio2];
        [aCoder encodeDouble:cl_sio2 forKey:kcl_sio2];
        [aCoder encodeDouble:dl_sio2 forKey:kdl_sio2];
        [aCoder encodeDouble:tg_sio2 forKey:ktg_sio2];
        [aCoder encodeDouble:cp_sio2 forKey:kcp_sio2];
        [aCoder encodeDouble:vLiq forKey:kvLiq];
        [aCoder encodeDouble:dvdtLiq forKey:kdvdtLiq];
        [aCoder encodeDouble:dvdpLiq forKey:kdvdpLiq];
        [aCoder encodeDouble:d2vdtdpLiq forKey:kd2vdtdpLiq];
        [aCoder encodeDouble:d2vdp2Liq forKey:kd2vdp2Liq];
    } // else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", [self className]];
}

#pragma mark -
#pragma mark Override superclass methods

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [self getEnthalpyFromT:t andP:p] - t*[self getEntropyFromT:t andP:p];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	double hl = 0.0;
	if (t >= tg_sio2) {
		hl = h0_sio2 + al_sio2*(tg_sio2 - tr) + bl_sio2*((tg_sio2*tg_sio2)-tr*tr)/2.0 - cl_sio2*(1.0/tg_sio2-1.0/tr)
		   + 2.0 *dl_sio2*(sqrt(tg_sio2)-sqrt(tr)) + (t - tg_sio2)*cp_sio2;
	} else {
		hl = h0_sio2 + al_sio2*(t-tr) + bl_sio2*(t*t-tr*tr)/2.0 - cl_sio2*(1.0/t - 1.0/tr) + 2.0*dl_sio2*(sqrt(t)-sqrt(tr));
	}
    double result = hl + (vLiq + dvdtLiq*(t-trl))*(p-pr) + 0.5*(dvdpLiq + (t-trl)*d2vdtdpLiq)*(p*p-pr*pr)
                       - (dvdpLiq + (t-trl)*d2vdtdpLiq)*pr*(p-pr)
                       + d2vdp2Liq*( (p*p*p-pr*pr*pr)/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 )
                       - t*(dvdtLiq*(p-pr) + 0.5*d2vdtdpLiq*(p-pr)*(p-pr));
    if (self.isGibbsFreeEnergyReferenceStateUsed) result += 298.15*self.entropyFromRobieEtAl1979;
    return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	double sl = 0.0;
	if (t >= tg_sio2) {
		sl = s0_sio2 + al_sio2*log(tg_sio2/tr) + bl_sio2*(tg_sio2 - tr) - (cl_sio2/2.0)*(1.0/(tg_sio2*tg_sio2) - 1.0/(tr*tr))
		   - 2.0*dl_sio2*(1.0/sqrt(tg_sio2)-1.0/sqrt(tr)) + cp_sio2*log(t/tg_sio2);
	} else {
		sl = s0_sio2 + al_sio2*log(t/tr) + bl_sio2*(t - tr) - (cl_sio2/2.0)*(1.0/(t*t) - 1.0/(tr*tr))
		   - 2.0*dl_sio2*(1.0/sqrt(t)-1.0/sqrt(tr));
	}

	return sl - (dvdtLiq*(p-pr) + 0.5*d2vdtdpLiq*(p-pr)*(p-pr));
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	if (t >= tg_sio2) {
		return cp_sio2;
	} else {
		return al_sio2 + bl_sio2*t + cl_sio2/(t*t) + dl_sio2/sqrt(t);
	}
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	if (t >= tg_sio2) {
		return 0.0;
	} else {
		return bl_sio2 - 2.0*cl_sio2/(t*t*t) - 0.5*dl_sio2/pow(t, (double) 1.5);
	}
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
