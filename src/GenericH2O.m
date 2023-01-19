//
//  GenericH2O.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 3/1/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "GenericH2O.h"
#import "FluidDuan.h"
#import "DEWH2O.h"
#import "HoltenJPCRD2014.h"
#import "Wagner2002.h"
#import "ZhangDuan2009.h"

typedef enum {
    NotApplicable,
    DuanAndZhang2006,
    ZhangAndDuan2005,
    HoltenEtAl2014,
    WagnerEtAl2002
} RegionType;

@interface GenericH2O() {
    DuanH2O *duanH2O;
    PhaseBase<StoichiometricPhaseProtocol> *dewH2O;
    HoltenJPCRD2014 *holtenJPCRD2014;
    Wagner2002 *wagner2002;
    RegionType forceRegion;
    double HliqReferenceWagner2002, SliqReferenceWagner2002, CpliqReferenceWagner2002;
    double widthOfWagnerToDZ2006SmoothInterval, centerOfWagnerToDZ2006SmoothInterval;
    double widthOfWagnerToZD2005SmoothInterval, centerOfWagnerToZD2005SmoothInterval;
    double widthOfDZ2006ToZD2005SmoothInterval, centerOfDZ2006ToZD2005SmoothInterval;
    double lowTofZD2005AdjustmentRegion, highTofZD2005AdjustmentRegion, pOfZD2005AdjustmentRegion, offsetsForZD2005AdjustmentRegion[11];
}

@end

@implementation GenericH2O

static const double HliqReference = -285830.0;      // J/mol, 298.15 K, 1 bar  NIST, CODATA
//static const double HgasReference = -241826.0;    // J/mol, 298.15 K, 1 bar  NIST, CODATA
static const double SliqReference =  69.950;        // J/K-mol, 298.15, 1 bar  NIST, CODATA
//static const double SgasReference = 188.835;      // J/K-mol, 298.15, 1 bar  NIST, CODATA
//static const double CpLiqReference = 75.3749;     // J/K-mol, 298.15, 1 bar  NIST, CODATA
static const double GliqReference = -56687.0*4.184; // J/mol, 298.15 K, 1 bar Helgeson and Kirkham, 1974

static const unsigned int valueOfG       =  0;
static const unsigned int valueOfH       =  1;
static const unsigned int valueOfS       =  2;
static const unsigned int valueOfCP      =  3;
static const unsigned int valueOfDCPDT   =  4;
static const unsigned int valueOfV       =  5;
static const unsigned int valueOfdVdT    =  6;
static const unsigned int valueOfdVdP    =  7;
static const unsigned int valueOfd2VdT2  =  8;
static const unsigned int valueOfd2VdTdP =  9;
static const unsigned int valueOfd2VdP2  = 10;

#pragma mark -
#pragma mark Private methods

- (instancetype)init {
    self = [super init];
    if (self) {
        self.phaseName    = @"Water";
        self.phaseFormula = @"H2O";
        duanH2O = [[DuanH2O alloc] init];
        dewH2O = [[DEWH2O alloc] init];
        holtenJPCRD2014 = [[HoltenJPCRD2014 alloc] init];
        wagner2002 = [[Wagner2002 alloc] init];
        forceRegion = NotApplicable;
        widthOfWagnerToDZ2006SmoothInterval  =  50.0;  // K
        centerOfWagnerToDZ2006SmoothInterval = 673.15; // K
        lowTofZD2005AdjustmentRegion         = 298.15; // K
        highTofZD2005AdjustmentRegion        = 398.15; // K
        pOfZD2005AdjustmentRegion            = -1.0;
        widthOfWagnerToZD2005SmoothInterval  =  100.0; // bars
        centerOfWagnerToZD2005SmoothInterval = 1000.0; // bars
        widthOfDZ2006ToZD2005SmoothInterval  =  100.0; // bars
        centerOfDZ2006ToZD2005SmoothInterval = 1000.0; // bars

        HliqReferenceWagner2002  = [wagner2002 getEnthalpyFromT:298.15 andP:1.0];
        SliqReferenceWagner2002  = [wagner2002 getEntropyFromT:298.15 andP:1.0];
        CpliqReferenceWagner2002 = [wagner2002 getHeatCapacityFromT:298.15 andP:1.0];
    }
    return self;
}

- (RegionType)classifyRegionForT:(double)t andForP:(double)p {
    RegionType region = NotApplicable;

    if      ((t > 673.15) && (p < 1000.0)) region = DuanAndZhang2006;
    else if ((t > 298.15) && (p > 1000.0)) region = ZhangAndDuan2005;
    else if  (t < 298.15)                  region = HoltenEtAl2014;
    else if  (t < 673.15)                  region = WagnerEtAl2002;
    else                                   region = ZhangAndDuan2005;

    return region;
}

#pragma mark -
#pragma mark Public methods

static NSString *region0 = @"No EOS is applicable at the specified temperature and pressure";
static NSString *region1 = @"MELTS H2O-CO2 from Duan and Zhang (2006)";
static NSString *region2 = @"DEW H2O from Zhang and Duan (2005)";
static NSString *region3 = @"Supercooled H2O from Holten et al. (2014)";
static NSString *region4 = @"Steam Properties from Wagner et al. (2002)";

- (void)forceModeChoiceAutomatic  { forceRegion = NotApplicable;    }
- (void)forceModeChoiceTo:(NSString *)region {
    if      ([region isEqualToString:region0]) forceRegion = NotApplicable;
    else if ([region isEqualToString:region1]) forceRegion = DuanAndZhang2006;
    else if ([region isEqualToString:region2]) forceRegion = ZhangAndDuan2005;
    else if ([region isEqualToString:region3]) forceRegion = HoltenEtAl2014;
    else if ([region isEqualToString:region4]) forceRegion = WagnerEtAl2002;
}

- (void)setModeToDuanAndZhang2006 { forceRegion = DuanAndZhang2006; }
- (void)setModeToZhangAndDuan2005 { forceRegion = ZhangAndDuan2005; }
- (void)setModeToHoltenEtAl2014   { forceRegion = HoltenEtAl2014;   }
- (void)setModeToWagnerEtAl2002   { forceRegion = WagnerEtAl2002;   }

- (void)useZhangAndDuan2005forDEW { dewH2O = [[DEWH2O alloc] init]; }
- (void)useDuanAndZhang2006forDEW { dewH2O = [[DuanH2O alloc] init]; }
- (void)useZhangAndDuan2009forDEW { dewH2O = [[ZhangDuan2009 alloc] init]; }

- (double)tTransHoltenWagner { return 298.15; }
- (double)tTransWagnerDZ2006 { return 673.15; }
- (double)tTransHoltenZD2005 { return 298.15; }
- (double)pTransWagnerZD2005 { return 1000.0; }
- (double)pTransDZ2006ZD2005 { return 1000.0; }

- (NSString *)eosForRegionAtT:(double)t andP:(double)p {
    RegionType region = (forceRegion == NotApplicable) ? [self classifyRegionForT:t andForP:p] : forceRegion;
    NSString *result = region0;
    switch (region) {
        case DuanAndZhang2006:
            result = region1;
            break;
        case ZhangAndDuan2005:
            result = region2;
            break;
        case HoltenEtAl2014:
            result = region3;
            break;
        case WagnerEtAl2002:
            result = region4;
            break;
        default:
            break;
    }
    return result;
}

- (double)lowerTemperatureLimitAtPinBars:(double)p {
    double TK = 0.0;
    if (p >= 2000.0) TK = [holtenJPCRD2014 homogeneousIceNucleationTemperatureForPressureInBars:p];
    else {
        double TKold = 200.0;
        double pOld = [holtenJPCRD2014 homogeneousIceNucleationPressureForTemperatureInK:TKold];
        TK = 220.0;
        double pGuess = [holtenJPCRD2014 homogeneousIceNucleationPressureForTemperatureInK:TK];
        NSUInteger iter = 0;
        while ((fabs(pGuess-p) > 0.1) && (iter < 200)) {
            double deriv = (pGuess - pOld)/(TK-TKold);
            TKold = TK;
            pOld = pGuess;
            TK = TK + (p-pGuess)/deriv;
            if (TK < 175.0 ) TK = 175.0;
            if (TK > 235.15) TK = 235.15;
            pGuess = [holtenJPCRD2014 homogeneousIceNucleationPressureForTemperatureInK:TK];
            iter++;
        }
    }
    return TK;
}

#pragma mark -
#pragma mark Generic smoothing method

- (double)getProperty:(int)valueOfProperty forClassInstance:(PhaseBase<StoichiometricPhaseProtocol>*)instance atT:(double)t andP:(double)p {
    double result = 0.0;
    switch (valueOfProperty) {
        case valueOfG:
        {
            result = [instance getGibbsFreeEnergyFromT:t andP:p];
            if ([instance isMemberOfClass:[holtenJPCRD2014 class]] || [instance isMemberOfClass:[wagner2002 class]]) {
                result +=     HliqReference - HliqReferenceWagner2002;
                result += -t*(SliqReference - SliqReferenceWagner2002);
            }
            result += GliqReference - (HliqReference - 298.15*SliqReference);
            break;
        }
        case valueOfH: {
            result = [instance getEnthalpyFromT:t andP:p];
            if ([instance isMemberOfClass:[holtenJPCRD2014 class]] || [instance isMemberOfClass:[wagner2002 class]]) {
                result += HliqReference - HliqReferenceWagner2002;
            }
            break;
        }
        case valueOfS: {
            result = [instance getEntropyFromT:t andP:p];
            if ([instance isMemberOfClass:[holtenJPCRD2014 class]] || [instance isMemberOfClass:[wagner2002 class]]) {
                result += SliqReference - SliqReferenceWagner2002;
            }
            break;
        }
        case valueOfCP:
            result = [instance getHeatCapacityFromT:t andP:p];
            break;
        case valueOfDCPDT:
            result = [instance getDcpDtFromT:t andP:p];
            break;
        case valueOfV:
            result = [instance getVolumeFromT:t andP:p];
            break;
        case valueOfdVdT:
            result = [instance getDvDtFromT:t andP:p];
            break;
        case valueOfdVdP:
            result = [instance getDvDpFromT:t andP:p];
            break;
        case valueOfd2VdT2:
            result = [instance getD2vDt2FromT:t andP:p];
            break;
        case valueOfd2VdTdP:
            result = [instance getD2vDtDpFromT:t andP:p];
            break;
        case valueOfd2VdP2:
            result = [instance getD2vDp2FromT:t andP:p];
            break;
        default:
            NSAssert(NO, @"Inconsistent method invocation: GenericH2O:getProperty:forClassInstance:atT:andP:");
            break;
    }
    return result;
}

- (double)smoothedPropertyAtT:(double)t andP:(double)p forMethodThatRetrieves:(int)valueOfProperty {
    double result = 0.0;
    if (forceRegion != NotApplicable) {
        switch (forceRegion) {
            case DuanAndZhang2006: {
                result = [self getProperty:valueOfProperty forClassInstance:duanH2O atT:t andP:p];
                break;
            }
            case ZhangAndDuan2005: {
                result = [self getProperty:valueOfProperty forClassInstance:dewH2O atT:t andP:p];
                break;
            }
            case HoltenEtAl2014: {
                result = [self getProperty:valueOfProperty forClassInstance:holtenJPCRD2014 atT:t andP:p];
                break;
            }
            case WagnerEtAl2002: {
                result += [self getProperty:valueOfProperty forClassInstance:wagner2002 atT:t andP:p];
                break;
            }
            default:
                break;
        }
        return  result;
    }
    RegionType region = [self classifyRegionForT:t andForP:p];
    switch (region) {
        case DuanAndZhang2006: {
            double weightT = 1.0/2.0 - tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
            double weightP = 1.0/2.0 + tanh((p-centerOfDZ2006ToZD2005SmoothInterval)/widthOfDZ2006ToZD2005SmoothInterval)/2.0;
            if (weightP <= 0.001) weightP = 0.0;
            if (weightT <= 0.001) weightT = 0.0;
            if (weightP  > 0.001) result += weightP*[self getProperty:valueOfProperty forClassInstance:dewH2O atT:t andP:p];
            if (weightT  > 0.001) result += weightT*(1.0-weightP)*[self getProperty:valueOfProperty forClassInstance:wagner2002 atT:t andP:p];
            result += (1.0-weightP-weightT*(1.0-weightP))*[self getProperty:valueOfProperty forClassInstance:duanH2O atT:t andP:p];
            break;
        }
        case ZhangAndDuan2005: {
            if (t < centerOfWagnerToDZ2006SmoothInterval) {
                double weightT = 1.0/2.0 + tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
                double weightP = 1.0/2.0 - tanh((p-centerOfWagnerToZD2005SmoothInterval)/widthOfWagnerToZD2005SmoothInterval)/2.0;
                if (weightP <= 0.001) weightP = 0.0;
                if (weightT <= 0.001) weightT = 0.0;
                if (weightP  > 0.001) result += weightP*(1.0-weightT)*[self getProperty:valueOfProperty forClassInstance:wagner2002 atT:t andP:p];
                if (weightT  > 0.001) result += weightP*weightT*[self getProperty:valueOfProperty forClassInstance:duanH2O atT:t andP:p];
                if (t < highTofZD2005AdjustmentRegion) {
                    if (p != pOfZD2005AdjustmentRegion) {
                        double holtenAtLowT = [self getProperty:valueOfProperty forClassInstance:holtenJPCRD2014 atT:lowTofZD2005AdjustmentRegion andP:p];
                        double ZD2005AtLowT = [self getProperty:valueOfProperty forClassInstance:dewH2O atT:lowTofZD2005AdjustmentRegion andP:p];
                        offsetsForZD2005AdjustmentRegion[valueOfProperty] = holtenAtLowT - ZD2005AtLowT;
                        pOfZD2005AdjustmentRegion = p;
                    }
                    double adjust = -offsetsForZD2005AdjustmentRegion[valueOfProperty]*(t-lowTofZD2005AdjustmentRegion)
                                     /(highTofZD2005AdjustmentRegion-lowTofZD2005AdjustmentRegion)
                                  + offsetsForZD2005AdjustmentRegion[valueOfProperty];
                    result += (1.0-weightP)*adjust;
                }
                result += (1.0-weightP)*[self getProperty:valueOfProperty forClassInstance:dewH2O atT:t andP:p];
            } else {
                double weightT = 1.0/2.0 - tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
                double weightP = 1.0/2.0 - tanh((p-centerOfDZ2006ToZD2005SmoothInterval)/widthOfDZ2006ToZD2005SmoothInterval)/2.0;
                if (weightP <= 0.001) weightP = 0.0;
                if (weightT <= 0.001) weightT = 0.0;
                if (weightP  > 0.001) result += weightP*(1.0-weightT)*[self getProperty:valueOfProperty forClassInstance:duanH2O atT:t andP:p];
                if (weightT  > 0.001) result += weightP*weightT*[self getProperty:valueOfProperty forClassInstance:wagner2002 atT:t andP:p];
                result += (1.0-weightP)*[self getProperty:valueOfProperty forClassInstance:dewH2O atT:t andP:p];
            }
            break;
        }
        case HoltenEtAl2014: {
            result = [self getProperty:valueOfProperty forClassInstance:holtenJPCRD2014 atT:t andP:p];
            break;
        }
        case WagnerEtAl2002: {
            double weightT = 1.0/2.0 + tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
            double weightP = 1.0/2.0 + tanh((p-centerOfWagnerToZD2005SmoothInterval)/widthOfWagnerToZD2005SmoothInterval)/2.0;
            if (weightP <= 0.001) weightP = 0.0;
            if (weightT <= 0.001) weightT = 0.0;
            if (weightP  > 0.001) {
                if (t < highTofZD2005AdjustmentRegion) {
                    if (p != pOfZD2005AdjustmentRegion) {
                        double holtenAtLowT = [self getProperty:valueOfProperty forClassInstance:holtenJPCRD2014 atT:lowTofZD2005AdjustmentRegion andP:p];
                        double ZD2005AtLowT = [self getProperty:valueOfProperty forClassInstance:dewH2O atT:lowTofZD2005AdjustmentRegion andP:p];
                        offsetsForZD2005AdjustmentRegion[valueOfProperty] = holtenAtLowT - ZD2005AtLowT;
                        pOfZD2005AdjustmentRegion = p;
                    }
                    double adjust = -offsetsForZD2005AdjustmentRegion[valueOfProperty]*(t-lowTofZD2005AdjustmentRegion)
                                     /(highTofZD2005AdjustmentRegion-lowTofZD2005AdjustmentRegion)
                                  + offsetsForZD2005AdjustmentRegion[valueOfProperty];
                    result += weightP*adjust;
                }
                result += weightP*[self getProperty:valueOfProperty forClassInstance:dewH2O atT:t andP:p];
            }
            if (weightT  > 0.001) result += weightT*(1.0-weightP)*[self getProperty:valueOfProperty forClassInstance:duanH2O atT:t andP:p];
            result += (1.0-weightP-weightT*(1.0-weightP))*[self getProperty:valueOfProperty forClassInstance:wagner2002 atT:t andP:p];
            break;
        }
        default:
            break;
    }
    return  result;
}

#pragma mark -
#pragma mark Stoichiometric Proptocol methods

- (double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfG];
}

- (double)getEnthalpyFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfH];
}

- (double)getEntropyFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfS];
}

- (double)getHeatCapacityFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfCP];
}

- (double)getDcpDtFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfDCPDT];
}

- (double)getVolumeFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfV];
}

- (double)getDvDtFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfdVdT];
}

- (double)getDvDpFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfdVdP];
}

- (double)getD2vDt2FromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfd2VdT2];
}

- (double)getD2vDtDpFromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfd2VdTdP];
}

- (double)getD2vDp2FromT:(double)t andP:(double)p {
    return [self smoothedPropertyAtT:t andP:p forMethodThatRetrieves:valueOfd2VdP2];
}

@end
