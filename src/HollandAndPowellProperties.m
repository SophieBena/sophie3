//
//  HollandAndPowellProperties.m
//  PhaseMELTSobjC
//
//  Created by Mark Ghiorso on 6/6/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import "HollandAndPowellProperties.h"
#ifndef __APPLE__
#include <dispatch/dispatch.h>
#endif

@interface HollandAndPowellPropertiesSettings : NSObject
@property (atomic, assign, getter = isGibbsFreeEnergyReferenceStateUsed) BOOL gibbsFreeEnergyReferenceStateUsed;
@end

@implementation HollandAndPowellPropertiesSettings
- (instancetype)init {
    self = [super init];
    if (self) {
        _gibbsFreeEnergyReferenceStateUsed = NO;
    }
    return self;
}
@end

@implementation HollandAndPowellProperties

+(HollandAndPowellPropertiesSettings *)hollandAndPowellPropertySettings {
    static HollandAndPowellPropertiesSettings *bermanPropertySettings = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        bermanPropertySettings = [[HollandAndPowellPropertiesSettings alloc] init];
    });
    return bermanPropertySettings;
}

+(void)enableGibbsFreeEnergyReferenceStateUsed {
    @synchronized (self.hollandAndPowellPropertySettings) {
        if (!self.hollandAndPowellPropertySettings.isGibbsFreeEnergyReferenceStateUsed) {
            self.hollandAndPowellPropertySettings.gibbsFreeEnergyReferenceStateUsed = YES;
        }
    }
}

+(void)disableGibbsFreeEnergyReferenceStateUsed {
    @synchronized (self.hollandAndPowellPropertySettings) {
        if (self.hollandAndPowellPropertySettings.isGibbsFreeEnergyReferenceStateUsed) {
            self.hollandAndPowellPropertySettings.gibbsFreeEnergyReferenceStateUsed = NO;
        }
    }
}

-(id)initWithH:(double)hIn
             S:(double)sIn
             a:(double)aIn
             b:(double)bIn
             c:(double)cIn
             d:(double)dIn
           Tc0:(double)Tc0In
          Smax:(double)SmaxIn
          Vmax:(double)VmaxIn
            v0:(double)v0In
            a0:(double)a0In
             K:(double)KIn {
    if ((self = [super init])) {
        h    = hIn;
        s    = sIn;
        a    = aIn;
        b    = bIn;
        c    = cIn;
        d    = dIn;
        Tc0  = Tc0In;
        Smax = SmaxIn;
        Vmax = VmaxIn;
        v0   = v0In;
        a0   = a0In;
        K    = KIn;

        tr = 298.15;
        pr = 1.0;

        _useHandPexcessVolumes = NO;

        if (Tc0 > 0.0) {
            Q2Tr = sqrt(1.0-tr/Tc0);
            spTr = Smax*Q2Tr;
            hpTr = Smax*Tc0*(Q2Tr - Q2Tr*Q2Tr*Q2Tr/3.0);
            vpTr = Vmax*Q2Tr;
            h  += hpTr;
            s  += spTr;
            v0 += vpTr;
        } else {
            Q2Tr = 0.0;
            spTr = 0.0;
            hpTr = 0.0;
        }

    }
    return self;
}

-(void)setTr:(double)trIn {
    tr = trIn;
}

-(void)setPr:(double)prIn {
    pr = prIn;
}

// Methods inherited from the stoichiometric phase protocol.

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    double htPr = h + a*(t-tr) + b*(t*t-tr*tr)/2.0 - c*(1.0/t-1.0/tr) + 2.0*d*(sqrt(t)-sqrt(tr));
    double stPr = s + a*log(t/tr) + b*(t-tr) - c*(1.0/t/t-1.0/tr/tr)/2.0 - 2.0*d*(1.0/sqrt(t)-1.0/sqrt(tr));
    double result = htPr - t*stPr;
    double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
    if(Tc > 0.0) {
        if (t >= Tc) result += 2.0*Smax*Tc/3.0 - t*Smax;
        else {
            double Q2 = sqrt(1.0 - t/Tc);
            result += 2.0*Smax*Tc*(Q2*Q2*Q2/6.0 - Q2/2.0 + 1.0/3.0) - t*Smax*(1.0-Q2);
        }
    }
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double intVdp = vT*kT*(pow(1.0+4.0*p/kT, 3.0/4.0) - pow(1.0+4.0/kT, 3.0/4.0))/3.0;
    result += intVdp;
    if (HollandAndPowellProperties.hollandAndPowellPropertySettings.isGibbsFreeEnergyReferenceStateUsed) result += 298.15*self.entropyFromRobieEtAl1979;
    return result;
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
    double result = h + a*(t-tr) + b*(t*t-tr*tr)/2.0 - c*(1.0/t-1.0/tr) + 2.0*d*(sqrt(t)-sqrt(tr));
    double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
    if(Tc > 0.0) {
        if (t >= Tc) result += 2.0*Smax*Tc/3.0;
        else {
            double Q2 = sqrt(1.0 - t/Tc);
            result += 2.0*Smax*Tc*(Q2*Q2*Q2/6.0 - Q2/2.0 + 1.0/3.0);
        }
    }
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double intVdp = vT*kT*(pow(1.0+4.0*p/kT, 3.0/4.0) - pow(1.0+4.0/kT, 3.0/4.0))/3.0;
    double DintVdpDt = intVdp*(DvTdT/vT + DkTdT/kT) - vT*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0) + vT*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0);
    result += intVdp - t*DintVdpDt;
    if (HollandAndPowellProperties.hollandAndPowellPropertySettings.isGibbsFreeEnergyReferenceStateUsed) result += 298.15*self.entropyFromRobieEtAl1979;
    return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    double result = s + a*log(t/tr) + b*(t-tr) - c*(1.0/t/t-1.0/tr/tr)/2.0 - 2.0*d*(1.0/sqrt(t)-1.0/sqrt(tr));
    double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
    if(Tc > 0.0) {
        if(t >= Tc) result += Smax;
        else {
            double Q2 = sqrt(1.0 - t/Tc);
            result += Smax*(1.0-Q2);
        }
    }
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double intVdp = vT*kT*(pow(1.0+4.0*p/kT, 3.0/4.0) - pow(1.0+4.0/kT, 3.0/4.0))/3.0;
    double DintVdpDt = intVdp*(DvTdT/vT + DkTdT/kT) - vT*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0) + vT*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0);
    result += -DintVdpDt;
    return result;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
    double result = a + b*t + c/t/t + d/sqrt(t);
    double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
    if (t < Tc) result += t*Smax/2.0/sqrt(Tc)/sqrt(Tc-t);
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double D2vTdT2 = v0*5.0*a0/pow(t,3.0/2.0);
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double intVdp = vT*kT*(pow(1.0+4.0*p/kT, 3.0/4.0) - pow(1.0+4.0/kT, 3.0/4.0))/3.0;
    double DintVdpDt = intVdp*(DvTdT/vT + DkTdT/kT) - vT*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0) + vT*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0);
    double D2intVdpDt2 = DintVdpDt*(DvTdT/vT + DkTdT/kT)
                       + intVdp*(-DvTdT*DvTdT/vT/vT + D2vTdT2/vT - DkTdT*DkTdT/kT/kT)
                       - DvTdT*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       + vT*p*DkTdT*DkTdT/kT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       - vT*p*p*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0*p/kT, 5.0/4.0)
                       + DvTdT*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       - vT*DkTdT*DkTdT/kT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       + vT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0/kT, 5.0/4.0);
    result += -t*D2intVdpDt2;
    return result;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
    double result = b - 2.0*c/t/t/t - d/2.0/t/sqrt(t);
    double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
    if (t < Tc) result += Smax/2.0/sqrt(Tc)/sqrt(Tc-t) + t*Smax/2.0/sqrt(Tc)/2.0/pow(Tc-t, 3.0/2.0);
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double D2vTdT2 = v0*5.0*a0/pow(t,3.0/2.0);
    double D3vTdT3 = -3.0*v0*5.0*a0/2.0/pow(t,5.0/2.0);
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double intVdp = vT*kT*(pow(1.0+4.0*p/kT, 3.0/4.0) - pow(1.0+4.0/kT, 3.0/4.0))/3.0;
    double DintVdpDt = intVdp*(DvTdT/vT + DkTdT/kT) - vT*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0) + vT*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0);
    double D2intVdpDt2 = DintVdpDt*(DvTdT/vT + DkTdT/kT)
                       + intVdp*(-DvTdT*DvTdT/vT/vT + D2vTdT2/vT - DkTdT*DkTdT/kT/kT)
                       - DvTdT*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       + vT*p*DkTdT*DkTdT/kT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       - vT*p*p*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0*p/kT, 5.0/4.0)
                       + DvTdT*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       - vT*DkTdT*DkTdT/kT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       + vT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0/kT, 5.0/4.0);
    double D3intVdpDt3 = D2intVdpDt2*(DvTdT/vT + DkTdT/kT)
                       + DintVdpDt*(D2vTdT2/vT - DvTdT*DvTdT/vT/vT - DkTdT*DkTdT/kT/kT)
                       + DintVdpDt*(-DvTdT*DvTdT/vT/vT + D2vTdT2/vT - DkTdT*DkTdT/kT/kT)
                       + intVdp*(-2.0*DvTdT*D2vTdT2/vT/vT + 2.0*DvTdT*DvTdT*DvTdT/vT/vT/vT + D3vTdT3/vT - D2vTdT2*DvTdT/vT/vT + 2.0*DkTdT*DkTdT*DkTdT/kT/kT/kT)
                       - D2vTdT2*p*DkTdT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       + DvTdT*p*DkTdT*DkTdT/kT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       - p*p*DvTdT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0*p/kT, 5.0/4.0)
                       + DvTdT*p*DkTdT*DkTdT/kT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       - 2.0*vT*p*DkTdT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0*p/kT, 1.0/4.0)
                       + p*p*vT*DkTdT*DkTdT*DkTdT/kT/kT/kT/kT/pow(1.0+4.0*p/kT, 5.0/4.0)
                       - DvTdT*p*p*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0*p/kT, 5.0/4.0)
                       + 3.0*vT*p*p*DkTdT*DkTdT*DkTdT/kT/kT/kT/kT/pow(1.0+4.0*p/kT, 5.0/4.0)
                       - 5.0*vT*p*p*p*DkTdT*DkTdT*DkTdT/kT/kT/kT/kT/kT/pow(1.0+4.0*p/kT, 9.0/4.0)
                       + D2vTdT2*DkTdT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       - DvTdT*DkTdT*DkTdT/kT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       + DvTdT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0/kT, 5.0/4.0)
                       - DvTdT*DkTdT*DkTdT/kT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       + 2.0*vT*DkTdT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       - vT*DkTdT*DkTdT*DkTdT/kT/kT/kT/kT/pow(1.0+4.0/kT, 1.0/4.0)
                       + DvTdT*DkTdT*DkTdT/kT/kT/kT/pow(1.0+4.0/kT, 5.0/4.0)
                       - 3.0*vT*DkTdT*DkTdT*DkTdT/kT/kT/kT/kT/pow(1.0+4.0/kT, 5.0/4.0)
                       + 5.0*vT*DkTdT*DkTdT*DkTdT/kT/kT/kT/kT/kT/pow(1.0+4.0/kT, 9.0/4.0);
    result += -D2intVdpDt2  -t*D3intVdpDt3;
    return result;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double vTP = vT*pow(1.0-4.0*p/(kT+4.0*p), 1.0/4.0);
    if (Tc0 > 0.0) {
        double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
        if (t < Tc) {
            double Q2 = sqrt(1.0 - t/Tc);
            if (self.useHandPexcessVolumes) vTP += Vmax*(1.0-Q2);
            else vTP += 2.0*Vmax*(Q2*Q2*Q2/6.0 - Q2/2.0 + 1.0/3.0);
        } else vTP += 2.0*Vmax/3.0;
    }
    return vTP;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double vTP = vT*pow(1.0-4.0*p/(kT+4.0*p), 1.0/4.0);
    double DvTPdT = vTP*(DvTdT/vT + p*DkTdT/kT/(kT+4.0*p));
    if (Tc0 > 0.0) {
        double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
        if (t < Tc) {
            double Q2 = sqrt(1.0 - t/Tc);
            double DQ2dT = -1.0/2.0/sqrt(1.0 - t/Tc)/Tc;
            if (self.useHandPexcessVolumes) DvTPdT += -Vmax*DQ2dT;
            else DvTPdT += 2.0*Vmax*(Q2*Q2*DQ2dT/2.0 - DQ2dT/2.0);
        }
    }
    return DvTPdT;
}

-(double)getDvDpFromT:(double)t andP:(double)p {
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double vTP = vT*pow(1.0-4.0*p/(kT+4.0*p), 1.0/4.0);
    double DvTPdP = -vTP/(kT+4.0*p);
    if (Tc0 > 0.0) {
        double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
        if (t < Tc) {
            double Q2 = sqrt(1.0 - t/Tc);
            double DQ2dP = t*Vmax/Smax/2.0/sqrt(1.0 - t/Tc)/Tc/Tc;
            if (self.useHandPexcessVolumes) DvTPdP += -Vmax*DQ2dP;
            else DvTPdP += 2.0*Vmax*(Q2*Q2*DQ2dP/2.0 - DQ2dP/2.0);
        }
    }
    return DvTPdP;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double D2vTdT2 = v0*5.0*a0/pow(t,3.0/2.0);
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double vTP = vT*pow(1.0-4.0*p/(kT+4.0*p), 1.0/4.0);
    double DvTPdT = vTP*(DvTdT/vT + p*DkTdT/kT/(kT+4.0*p));
    double D2vDt2 = DvTPdT*DvTPdT/vTP + vTP*(-DvTdT*DvTdT/vT + D2vTdT2)/vT -vTP*2.0*p*(kT+2.0*p)*DkTdT*DkTdT/pow((kT+4.0*p)*kT, 2.0);
    if (Tc0 > 0.0) {
        double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
        if (t < Tc) {
            double Q2 = sqrt(1.0 - t/Tc);
            double DQ2dT = -1.0/2.0/sqrt(1.0 - t/Tc)/Tc;
            double D2Q2dT2 = -1.0/4.0/pow(1.0 - t/Tc, 3.0/2.0)/Tc/Tc;
            if (self.useHandPexcessVolumes) D2vDt2 += -Vmax*D2Q2dT2;
            else D2vDt2 += 2.0*Vmax*(Q2*DQ2dT*DQ2dT + Q2*Q2*D2Q2dT2/2.0 - D2Q2dT2/2.0);
        }
    }
    return D2vDt2;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double DvTdT = v0*(a0 - 10.0*a0/sqrt(t));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double DkTdT = -K*1.5e-4;
    double vTP = vT*pow(1.0-4.0*p/(kT+4.0*p), 1.0/4.0);
    double DvTPdP = -vTP/(kT+4.0*p);
    double DvTPdT = vTP*(DvTdT/vT + p*DkTdT/kT/(kT+4.0*p));
    double D2vDtDp = DvTPdP*(DvTPdT/vTP - DkTdT/(kT+4.0*p));
    if (Tc0 > 0.0) {
        double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
        if (t < Tc) {
            double Q2 = sqrt(1.0 - t/Tc);
            double DQ2dP = t*Vmax/Smax/2.0/sqrt(1.0 - t/Tc)/Tc/Tc;
            double DQ2dT = -1.0/2.0/sqrt(1.0 - t/Tc)/Tc;
            double D2Q2dTdP = Vmax/Smax/2.0/sqrt(1.0 - t/Tc)/Tc/Tc + t*Vmax/Smax/4.0/pow(1.0 - t/Tc, 3.0/2.0)/Tc/Tc/Tc;
            if (self.useHandPexcessVolumes) D2vDtDp += -Vmax*D2Q2dTdP;
            else D2vDtDp += 2.0*Vmax*(Q2*DQ2dP*DQ2dT + Q2*Q2*D2Q2dTdP/2.0 - D2Q2dTdP/2.0);
        }
    }
    return D2vDtDp;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
    double vT = v0*(1.0 + a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double kT = K*(1.0 - 1.5e-4*(t-tr));
    double vTP = vT*pow(1.0-4.0*p/(kT+4.0*p), 1.0/4.0);
    double DvTPdP = -vTP/(kT+4.0*p);
    double D2vDp2 = DvTPdP*(DvTPdP/vTP - 4.0/(kT+4.0*p));
    if (Tc0 > 0.0) {
        double Tc = (Tc0 > 0.0) ? Tc0 + Vmax*(p-pr)/Smax : 0.0;
        if (t < Tc) {
            double Q2 = sqrt(1.0 - t/Tc);
            double DQ2dP = t*Vmax/Smax/2.0/sqrt(1.0 - t/Tc)/Tc/Tc;
            double D2Q2dP2 = -t*t*Vmax*Vmax/Smax/Smax/4.0/pow(1.0 - t/Tc, 3.0/2.0)/Tc/Tc/Tc/Tc
                           - 2.0*t*Vmax*Vmax/Smax/Smax/2.0/sqrt(1.0 - t/Tc)/Tc/Tc/Tc;
            if (self.useHandPexcessVolumes) D2vDp2 += -Vmax*D2Q2dP2;
            else D2vDp2 += 2.0*Vmax*(Q2*DQ2dP*DQ2dP + Q2*Q2*D2Q2dP2/2.0  - D2Q2dP2/2.0);
        }
    }
    return D2vDp2;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kh    = @"h";
static NSString *ks    = @"s";
static NSString *ka    = @"a";
static NSString *kb    = @"b";
static NSString *kc    = @"c";
static NSString *kd    = @"d";
static NSString *kTc0  = @"Tc0";
static NSString *kSmax = @"Smax";
static NSString *kVmax = @"Vmax";
static NSString *kv0   = @"v0";
static NSString *ka0   = @"a0";
static NSString *kK    = @"K";
static NSString *ktr   = @"tr";
static NSString *kpr   = @"pr";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        h    = [aDecoder decodeDoubleForKey:kh];
        s    = [aDecoder decodeDoubleForKey:ks];
        a    = [aDecoder decodeDoubleForKey:ka];
        b    = [aDecoder decodeDoubleForKey:kb];
        c    = [aDecoder decodeDoubleForKey:kc];
        d    = [aDecoder decodeDoubleForKey:kd];
        Tc0  = [aDecoder decodeDoubleForKey:kTc0];
        Smax = [aDecoder decodeDoubleForKey:kSmax];
        Vmax = [aDecoder decodeDoubleForKey:kVmax];
        v0   = [aDecoder decodeDoubleForKey:kv0];
        a0   = [aDecoder decodeDoubleForKey:ka0];
        K    = [aDecoder decodeDoubleForKey:kK];
        tr   = [aDecoder decodeDoubleForKey:ktr];
        pr   = [aDecoder decodeDoubleForKey:kpr];
        if (Tc0 > 0.0) {
            Q2Tr = sqrt(1.0-tr/Tc0);
            spTr = Smax*Q2Tr;
            hpTr = Smax*Tc0*(Q2Tr - Q2Tr*Q2Tr*Q2Tr/3.0);
            vpTr = Vmax*Q2Tr;
            h  += hpTr;
            s  += spTr;
            v0 += vpTr;
        } else {
            Q2Tr = 0.0;
            spTr = 0.0;
            hpTr = 0.0;
        }
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:h    forKey:kh];
        [aCoder encodeDouble:s    forKey:ks];
        [aCoder encodeDouble:a    forKey:ka];
        [aCoder encodeDouble:b    forKey:kb];
        [aCoder encodeDouble:c    forKey:kc];
        [aCoder encodeDouble:d    forKey:kd];
        [aCoder encodeDouble:Tc0  forKey:kTc0];
        [aCoder encodeDouble:Smax forKey:kSmax];
        [aCoder encodeDouble:Vmax forKey:kVmax];
        [aCoder encodeDouble:v0   forKey:kv0];
        [aCoder encodeDouble:a0   forKey:ka0];
        [aCoder encodeDouble:K    forKey:kK];
        [aCoder encodeDouble:tr   forKey:ktr];
        [aCoder encodeDouble:pr   forKey:kpr];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class BermanProperties only supports NSKeyedArchiver coders."];
}

@synthesize h, s, a, b, c, d, Tc0, Smax, Vmax, v0, a0, K;
@end
