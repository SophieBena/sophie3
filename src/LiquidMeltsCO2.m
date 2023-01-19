//
//  LiquidMeltsCO2.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/26/12.
//  Copyright 2012 OFM Research Inc.. All rights reserved.
//

#import "LiquidMeltsCO2.h"

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

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

@implementation LiquidMeltsCO2

-(id)init {
    if ((self = [super init])) {
        _enthalpyCorrectionForCO2 = -6.3093193811701e-01*1000.0;
        _entropyCorrectionForCO2  = -1.0939331414050e+02;
        _volumeForCO2  =  4.0157994267547e+00;
        _dvdtForCO2    =  1.2131890000000e+00/1000.0;
        _dvdpForCO2    = -4.2673870000000e-01/10000.0;
        _d2vdtdpForCO2 = 0.0;
        _d2vdp2ForCO2  = 0.0;
        duanCO2 = [[DuanCO2 alloc] init];
        _gibbsFreeEnergyReferenceStateUsed = NO;
        self.phaseName = @"CO2";
        self.phaseFormula = @"CO2";
    }
    return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol methods

static NSString *kenthalpyCorrectionForCO2 = @"enthalpyCorrectionForCO2";
static NSString *kentropyCorrectionForCO2 = @"entropyCorrectionForCO2";
static NSString *kvolumeForCO2 = @"volumeForCO2";
static NSString *kdvdtForCO2 = @"dvdtForCO2";
static NSString *kdvdpForCO2 = @"dvdpForCO2";
static NSString *kd2vdtdpForCO2 = @"d2vdtdpForCO2";
static NSString *kd2vdp2ForCO2 = @"d2vdp2ForCO2";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        _enthalpyCorrectionForCO2 = [aDecoder decodeDoubleForKey:kenthalpyCorrectionForCO2];
        _entropyCorrectionForCO2 = [aDecoder decodeDoubleForKey:kentropyCorrectionForCO2];
        _volumeForCO2 = [aDecoder decodeDoubleForKey:kvolumeForCO2];
        _dvdtForCO2 = [aDecoder decodeDoubleForKey:kdvdtForCO2];
        _dvdpForCO2 = [aDecoder decodeDoubleForKey:kdvdpForCO2];
        _d2vdtdpForCO2 = [aDecoder decodeDoubleForKey:kd2vdtdpForCO2];
        _d2vdp2ForCO2 = [aDecoder decodeDoubleForKey:kd2vdp2ForCO2];

        duanCO2 = [[DuanCO2 alloc] init];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:_enthalpyCorrectionForCO2 forKey:kenthalpyCorrectionForCO2];
        [aCoder encodeDouble:_entropyCorrectionForCO2 forKey:kentropyCorrectionForCO2];
        [aCoder encodeDouble:_volumeForCO2 forKey:kvolumeForCO2];
        [aCoder encodeDouble:_dvdtForCO2 forKey:kdvdtForCO2];
        [aCoder encodeDouble:_dvdpForCO2 forKey:kdvdpForCO2];
        [aCoder encodeDouble:_d2vdtdpForCO2 forKey:kd2vdtdpForCO2];
        [aCoder encodeDouble:_d2vdp2ForCO2 forKey:kd2vdp2ForCO2];
    } // else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", [self className]];
}

#pragma mark -
#pragma mark Stoichiometric Phase protocol methods

-(double)calculateWithT:(double)t andWithP:(double)p forCase:(NSUInteger)returnMode {
	double vCO2       = [self volumeForCO2];
    double dvCO2dt    = [self dvdtForCO2];
    double dvCO2dp    = [self dvdpForCO2];
    double d2vCO2dtdp = [self d2vdtdpForCO2];
    double d2vCO2dp2  = [self d2vdp2ForCO2];
	double result = 0.0;

	double tRef = 1673.15;
    double pRef = 1.0;

    double vl       = vCO2 + dvCO2dt*(t-tRef) + dvCO2dp*(p-pRef) + d2vCO2dtdp*(t-tRef)*(p-pRef) + d2vCO2dp2*(p-pRef)*(p-pRef);
    double dvldt    = dvCO2dt + d2vCO2dtdp*(p-pRef);
    double dvldp    = dvCO2dp + 2.0*d2vCO2dp2*(p-pRef);
    double d2vldt2  = 0.0;
    double d2vldp2  = 2.0*d2vCO2dp2;
    double d2vldtdp = d2vCO2dtdp;

    double duanGatOneBar     = [duanCO2 getGibbsFreeEnergyFromT:t andP:1.0];
    double duanHatOneBar     = [duanCO2 getEnthalpyFromT:t andP:1.0];
    double duanSatOneBar     = [duanCO2 getEntropyFromT:t andP:1.0];
    double duanCpatOneBar    = [duanCO2 getHeatCapacityFromT:t andP:1.0];
    double duandCpDtatOneBar = [duanCO2 getDcpDtFromT:t andP:1.0];

    double gCO2     = duanGatOneBar + [self enthalpyCorrectionForCO2] - t*[self entropyCorrectionForCO2];
    double hCO2     = duanHatOneBar + [self enthalpyCorrectionForCO2];
    double sCO2     = duanSatOneBar + [self entropyCorrectionForCO2];
    double cpCO2    = duanCpatOneBar;
    double dcpdtCO2 = duandCpDtatOneBar;

	switch (returnMode) {
		case 1:
			result = gCO2 + vCO2*(p-pRef) + dvCO2dt*(t-tRef)*(p-pRef) + dvCO2dp*(p-pRef)*(p-pRef)/2.0 + d2vCO2dtdp*(t-tRef)*(p-pRef)*(p-pRef)/2.0 + d2vCO2dp2*(p-pRef)*(p-pRef)*(p-pRef)/3.0;
			break;
		case 2:
			result = hCO2 + vCO2*(p-pRef) + dvCO2dt*(t-tRef)*(p-pRef) + dvCO2dp*(p-pRef)*(p-pRef)/2.0 + d2vCO2dtdp*(t-tRef)*(p-pRef)*(p-pRef)/2.0 + d2vCO2dp2*(p-pRef)*(p-pRef)*(p-pRef)/3.0
            - t*dvCO2dt*(p-pRef) - t*d2vCO2dtdp*(p-pRef)*(p-pRef)/2.0;
			break;
		case 3:
			result = sCO2 - dvCO2dt*(p-pRef) - d2vCO2dtdp*(p-pRef)*(p-pRef)/2.0;
			break;
		case 4:
			result = cpCO2;
			break;
		case 5:
			result = dcpdtCO2;
			break;
        case 6:
            result = vl;
            break;
        case 7:
            result = dvldt;
            break;
        case 8:
            result = dvldp;
            break;
        case 9:
            result = d2vldt2;
            break;
        case 10:
            result = d2vldtdp;
            break;
        case 11:
            result = d2vldp2;
            break;
		default:
			break;
	}
	return result;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    double result = [self calculateWithT:t andWithP:p forCase:returnValueOfG];
    if (self.isGibbsFreeEnergyReferenceStateUsed) result += 298.15*self.entropyFromRobieEtAl1979;
	return result;
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
    double result = [self calculateWithT:t andWithP:p forCase:returnValueOfH];
    if (self.isGibbsFreeEnergyReferenceStateUsed) result += 298.15*self.entropyFromRobieEtAl1979;
	return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfS];
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfCP];
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfDCPDT];
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	return  [self calculateWithT:t andWithP:p forCase:returnValueOfV];
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfdVdT];
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	return  [self calculateWithT:t andWithP:p forCase:returnValueOfdVdP];
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdT2];
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdTdP];
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdP2];;
}

@end
