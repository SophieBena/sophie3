//
//  Wagner2002.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 3/12/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "Wagner2002.h"
#ifdef BUILD_WITH_INLINE_FREESTEAM
#import "steam.h"
#import "region4.h"
#import "backwards.h"
#import "b23.h"
#import "derivs.h"
#import "zeroin.h"
#import "region3.h"
#import "solver2.h"
#import "steam_ph.h"
#import "steam_ps.h"
#import "steam_Ts.h"
#import "steam_pT.h"
#import "steam_pv.h"
#import "steam_Tx.h"
#import "region1.h"
#import "viscosity.h"
#import "thcond.h"
#import "surftens.h"
#else
#import <FreeSteam21/FreeSteam21.h>
#endif

@interface Wagner2002() {
    double EPS, EPS2;
}

@end

@implementation Wagner2002

static const double MolesPerKg = 55.508435;

- (instancetype)init {
    self = [super init];
    if (self) {
        self.phaseName    = @"Wagner Water";
        self.phaseFormula = @"H2O";
        EPS = sqrt(DBL_EPSILON);
        EPS2 = sqrt(sqrt(DBL_EPSILON));

    }
    return self;
}

#pragma mark -
#pragma mark Stoichiometric Proptocol methods

- (double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t);
    double result = (freesteam_h(S) - t*freesteam_s(S))/MolesPerKg;
    return result;
}

- (double)getEnthalpyFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t);
    double result = freesteam_h(S)/MolesPerKg;
    return result;
}

- (double)getEntropyFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t);
    double result = freesteam_s(S)/MolesPerKg;
    return result;
}

- (double)getHeatCapacityFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t);
    double result = freesteam_cp(S)/MolesPerKg;
    return result;
}

- (double)getDcpDtFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t*(1.0+EPS));
    double result  = freesteam_cp(S);
    S = freesteam_set_pT(p*1.0e5, t*(1.0-EPS));
    result -= freesteam_cp(S);
    result *= 1.0/t/2.0/EPS/MolesPerKg; // e5
    return  result;
}

- (double)getVolumeFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t);
    double result = 1.0e5*freesteam_v(S)/MolesPerKg;
    return  result;
}

- (double)getDvDtFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t*(1.0+EPS));
    double result  = freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5, t*(1.0-EPS));
    result -= freesteam_v(S);
    result *= 1.0e5/t/2.0/EPS/MolesPerKg;
    return  result;
}

- (double)getDvDpFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5*(1.0+EPS), t);
    double result  = freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-EPS), t);
    result -= freesteam_v(S);
    result *= 1.0e5/p/2.0/EPS/MolesPerKg;
    return  result;
}

- (double)getD2vDt2FromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5, t*(1.0+2.0*EPS2));
    double result  = -(1.0/12.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5, t*(1.0+EPS2));
    result +=  (4.0/3.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5, t);
    result += -(5.0/2.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5, t*(1.0-EPS2));
    result +=  (4.0/3.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5, t*(1.0-2.0*EPS2));
    result += -(1.0/12.0)*freesteam_v(S);
    result *= 1.0e5/t/t/EPS2/EPS2/MolesPerKg;
    return  result;
}

- (double)getD2vDtDpFromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5*(1.0+2.0*EPS2), t*(1.0+2.0*EPS2));
    double result  = -44.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0+2.0*EPS2), t*(1.0-2.0*EPS2));
    result +=  44.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-2.0*EPS2), t*(1.0+2.0*EPS2));
    result +=  44.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-2.0*EPS2), t*(1.0-2.0*EPS2));
    result += -44.0*freesteam_v(S);

    S = freesteam_set_pT(p*1.0e5*(1.0+EPS2), t*(1.0+EPS2));
    result +=  74.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0+EPS2), t*(1.0-EPS2));
    result += -74.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-EPS2), t*(1.0+EPS2));
    result += -74.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-EPS2), t*(1.0-EPS2));
    result +=  74.0*freesteam_v(S);

    S = freesteam_set_pT(p*1.0e5*(1.0+2.0*EPS2), t*(1.0+EPS2));
    result +=  63.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0+EPS2), t*(1.0+2.0*EPS2));
    result +=  63.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-2.0*EPS2), t*(1.0-EPS2));
    result +=  63.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-EPS2), t*(1.0-2.0*EPS2));
    result +=  63.0*freesteam_v(S);

    S = freesteam_set_pT(p*1.0e5*(1.0+EPS2), t*(1.0-2.0*EPS2));
    result += -63.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0+2.0*EPS2), t*(1.0-EPS2));
    result += -63.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-2.0*EPS2), t*(1.0+EPS2));
    result += -63.0*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-EPS2), t*(1.0+2.0*EPS2));
    result += -63.0*freesteam_v(S);

    result *= 1.0e5/600.0/p/EPS2/t/EPS2/MolesPerKg;
    return  result;
}

- (double)getD2vDp2FromT:(double)t andP:(double)p {
    SteamState S = freesteam_set_pT(p*1.0e5*(1.0+2.0*EPS2), t);
    double result  = -(1.0/12.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0+EPS2), t);
    result +=  (4.0/3.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5, t);
    result += -(5.0/2.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-EPS2), t);
    result +=  (4.0/3.0)*freesteam_v(S);
    S = freesteam_set_pT(p*1.0e5*(1.0-2.0*EPS2), t);
    result += -(1.0/12.0)*freesteam_v(S);
    result *= 1.0e5/p/p/EPS2/EPS2/MolesPerKg;
    return  result;
}

@end
