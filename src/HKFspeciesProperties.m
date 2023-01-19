//
//  HKFspeciesProperties.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 11/2/15.
//  Copyright © 2015 Mark Ghiorso. All rights reserved.
//

#import "HKFspeciesProperties.h"
#import "DEWDielectricConstant.h"

@interface HKFspeciesProperties () {
    DEWDielectricConstant *dielectricConstantFunctions;
    double tLast, pLast;
}
@end

#define CapPSI 2600.0   // bars, Tanger and Helgeson (1988)
#define THETA   228.0   // K, Helgeson and Kirkham (1976)

#define VERY_SMALL 0.01 // K or bars
#define R   8.3144598   // J/K/mol
#define MW 18.01528     // g/mol

#define PR 1.0          // bar
#define TR 298.15       // K

#define epsilonTRPR 78.47
#define YatTRPR     (-0.0000579865)

@implementation HKFspeciesProperties

#pragma mark -
#pragma mark Class initializers

- (instancetype)init {
    self = [super init];
    if (self) {
        dielectricConstantFunctions = [[DEWDielectricConstant alloc] init];

        _deltaGibbsFreeEnergyOfFormationInTheReferenceState = 0.0;
        _deltaEnthalpyOfFormationInTheReferenceState = 0.0;
        _entropyInTheReferenceState = 0.0;
        _volumeInTheReferenceState = 0.0;
        _heatCapacityInTheReferenceState = 0.0;
        _a1HKF = 0.0;
        _a2HKF = 0.0;
        _a3HKF = 0.0;
        _a4HKF = 0.0;
        _c1HKF = 0.0;
        _c2HKF = 0.0;
        _omegaHKF = 0.0;
        _charge = 0.0;

        tLast = -1.0;
        pLast = -1.0;
    }
    return self;
}

- (instancetype)initWithGibbsFreeEnergyOfFormation:(double)gIn
                        andWithEnthalpyOfFormation:(double)hIn
                           andWithReferenceEntropy:(double)sIn
                            andWithReferenceVolume:(double)vIn
                      andWithReferenceHeatCapacity:(double)cpIn
                             andWithHKFa1Parameter:(double)a1In
                             andWithHKFa2Parameter:(double)a2In
                             andWithHKFa3Parameter:(double)a3In
                             andWithHKFa4Parameter:(double)a4In
                             andWithHKFc1Parameter:(double)c1In
                             andWithHKFc2Parameter:(double)c2In
                             andWithHKFomegaParameter:(double)omegaIn
                                     andWithCharge:(double)chargeIn {
    self = [super init];
    if (self) {
        dielectricConstantFunctions = [[DEWDielectricConstant alloc] init];

        _deltaGibbsFreeEnergyOfFormationInTheReferenceState = gIn;
        _deltaEnthalpyOfFormationInTheReferenceState = hIn;
        _entropyInTheReferenceState = sIn;
        _volumeInTheReferenceState = vIn;
        _heatCapacityInTheReferenceState = cpIn;
        _a1HKF = a1In;
        _a2HKF = a2In;
        _a3HKF = a3In;
        _a4HKF = a4In;
        _c1HKF = c1In;
        _c2HKF = c2In;
        _omegaHKF = omegaIn;
        _charge = chargeIn;

        tLast = -1.0;
        pLast = -1.0;
    }
    return self;
}

#pragma mark -
#pragma mark Generic private functions

- (void)updateDensityPropertiesWithT:(double)t andP:(double)p {
    [dielectricConstantFunctions loadDensityPropertiesWithT:t andWithP:p];
    tLast = t;
    pLast = p;
}

#pragma mark -
#pragma mark Born term solvent functions

// Shock et al., 1992, J. Chem. Soc. Faraday Trans. 88(6) 803-826
#define agP   (-2.037662)
#define agPP  ( 5.747000e-3)
#define agPPP (-6.557892e-6)

#define bgP   (6.107361)
#define bgPP  (-1.074337e-2)
#define bgPPP (1.268348e-5)

#define ag1   ( 3.66666e-16)
#define ag2   (-1.504956e-10)
#define ag3   ( 5.01799e-14)

- (double)gSolventAtT:(double)t andP:(double)p {
    if ((fabs(t-tLast) < VERY_SMALL) || (fabs(p-pLast) < VERY_SMALL)) [self updateDensityPropertiesWithT:t andP:p];
    double rho = dielectricConstantFunctions.rhoLast;
    double tc = t - 273.15;

    if (rho >= 1.0) return 0.0;
    else {
        double aG = agP + agPP*tc + agPPP*tc*tc;
        double bG = bgP + bgPP*tc + bgPPP*tc*tc;

        double f = 0.0;
        if ((p <= 1000) && (tc >= 155.0) && (tc <= 355.0)) {
            f += (pow((tc-155.0)/300.0, 4.8) + ag1*pow((tc-155.0)/300.0, 16.0))
               * (ag2*pow(1000.0-p, 3.0) + ag3*pow(1000.0-p, 4.0));
        }

        return aG * pow(1.0-rho, bG) - f;
    }
}

- (double)DgSolventDtAtT:(double)t andP:(double)p {
    if ((fabs(t-tLast) < VERY_SMALL) || (fabs(p-pLast) < VERY_SMALL)) [self updateDensityPropertiesWithT:t andP:p];
    double rho    = dielectricConstantFunctions.rhoLast;
    double DrhoDt = dielectricConstantFunctions.DrhoDtLast;
    double tc = t - 273.15;

    if (rho >= 1.0) return 0.0;
    else {
        double  aG   = agP + agPP*tc +     agPPP*tc*tc;
        double daGdt =       agPP    + 2.0*agPPP*tc;
        double  bG   = bgP + bgPP*tc +     bgPPP*tc*tc;
        double dbGdt =       bgPP    + 2.0*bgPPP*tc;

        double g = aG * pow(1.0-rho, bG);
        double dgdt = daGdt*pow(1.0-rho, bG) + g*(dbGdt*log(1.0-rho) - bG*DrhoDt/(1.0-rho));

        double dfdt = 0.0;
        if ((p <= 1000) && (tc >= 155.0) && (tc <= 355.0)) {
            // Note that Sverjensky et al. (2014) has an entirely independent function for this derivative
            // which I have not adopted, prefering internal consistency to accuracy
            dfdt += ((4.8/300.0)*pow((tc-155.0)/300.0, 3.8) + ag1*(16.0/300.0)*pow((tc-155.0)/300.0, 15.0))
                  * (ag2*pow(1000.0-p, 3.0) + ag3*pow(1000.0-p, 4.0));
        }

        return dgdt - dfdt;
    }
}

- (double)DgSolventDpAtT:(double)t andP:(double)p {
    if ((fabs(t-tLast) < VERY_SMALL) || (fabs(p-pLast) < VERY_SMALL)) [self updateDensityPropertiesWithT:t andP:p];
    double rho    = dielectricConstantFunctions.rhoLast;
    double DrhoDp = dielectricConstantFunctions.DrhoDpLast;
    double tc = t - 273.15;

    if (rho >= 1.0) return 0.0;
    else {
        double aG = agP + agPP*tc + agPPP*tc*tc;
        double bG = bgP + bgPP*tc + bgPPP*tc*tc;

       double dgdp = -aG*bG*DrhoDp*pow(1.0-rho, bG-1.0);

        double dfdp = 0.0;
        if ((p <= 1000) && (tc >= 155.0) && (tc <= 355.0)) {
            dfdp += (pow((tc-155.0)/300.0, 4.8) + ag1*pow((tc-155.0)/300.0, 16.0))
                  * (-3.0*ag2*pow(1000.0-p, 2.0) - 4.0*ag3*pow(1000.0-p, 3.0));
        }

        return dgdp - dfdp;
    }
}

- (double)D2gSolventDt2AtT:(double)t andP:(double)p {
    if ((fabs(t-tLast) < VERY_SMALL) || (fabs(p-pLast) < VERY_SMALL)) [self updateDensityPropertiesWithT:t andP:p];
    double   rho    = dielectricConstantFunctions.rhoLast;
    double  DrhoDt  = dielectricConstantFunctions.DrhoDtLast;
    double D2rhoDt2 = dielectricConstantFunctions.D2rhoDt2Last;
    double tc = t - 273.15;

    if (rho >= 1.0) return 0.0;
    else {
        double   aG    = agP + agPP*tc +     agPPP*tc*tc;
        double  daGdt  =       agPP    + 2.0*agPPP*tc;
        double d2aGdt2 =                 2.0*agPPP;
        double   bG    = bgP + bgPP*tc +     bgPPP*tc*tc;
        double  dbGdt  =       bgPP    + 2.0*bgPPP*tc;
        double d2bGdt2 =                 2.0*bgPPP;

        double g = aG * pow(1.0-rho, bG);
        double d2gdt2 = d2aGdt2*pow(1.0-rho, bG)
                      + 2.0*daGdt*pow(1.0-rho, bG)*(dbGdt*log(1.0-rho) - bG*DrhoDt/(1.0-rho))
                      + g*pow(dbGdt*log(1.0-rho) - bG*DrhoDt/(1.0-rho), 2.0)
                      + g*(d2bGdt2*log(1.0-rho) - 2.0*dbGdt*DrhoDt/(1.0-rho) -bG*pow(DrhoDt/(1.0-rho), 2.0) - bG*D2rhoDt2/(1.0-rho));

        double d2fdt2 = 0.0;
        if ((p <= 1000) && (tc >= 155.0) && (tc <= 355.0)) {
            d2fdt2 += ((3.8/300.0)*(4.8/300.0)*pow((tc-155.0)/300.0, 2.8) + ag1*(15.0/300.0)*(16.0/300.0)*pow((tc-155.0)/300.0, 14.0))
                    * (ag2*pow(1000.0-p, 3.0) + ag3*pow(1000.0-p, 4.0));
        }

        return d2gdt2 - d2fdt2;
    }
}

- (double)D2gSolventDtDpAtT:(double)t andP:(double)p {
    if ((fabs(t-tLast) < VERY_SMALL) || (fabs(p-pLast) < VERY_SMALL)) [self updateDensityPropertiesWithT:t andP:p];
    double   rho     = dielectricConstantFunctions.rhoLast;
    double  DrhoDt   = dielectricConstantFunctions.DrhoDtLast;
    double  DrhoDp   = dielectricConstantFunctions.DrhoDpLast;
    double D2rhoDtDp = dielectricConstantFunctions.D2rhoDtDpLast;
    double tc = t - 273.15;

    if (rho >= 1.0) return 0.0;
    else {
        double  aG   = agP + agPP*tc +     agPPP*tc*tc;
        double daGdt =       agPP    + 2.0*agPPP*tc;
        double  bG   = bgP + bgPP*tc +     bgPPP*tc*tc;
        double dbGdt =       bgPP    + 2.0*bgPPP*tc;

        double g = aG * pow(1.0-rho, bG);
        double d2gdtdp = -bG*daGdt*DrhoDp*pow(1.0-rho, bG-1.0)
                       + g*(-dbGdt*DrhoDp/(1.0-rho) - bG*DrhoDt*DrhoDp/pow(1.0-rho, 2.0) - bG*D2rhoDtDp/(1.0-rho))
                       - aG*bG*DrhoDp*pow(1.0-rho, bG-1.0)*(dbGdt*log(1.0-rho) - bG*DrhoDt/(1.0-rho));

        double d2fdtdp = 0.0;
        if ((p <= 1000) && (tc >= 155.0) && (tc <= 355.0)) {
            d2fdtdp += ((4.8/300.0)*pow((tc-155.0)/300.0, 3.8) + ag1*(16.0/300.0)*pow((tc-155.0)/300.0, 15.0))
                     * (-3.0*ag2*pow(1000.0-p, 2.0) - 4.0*ag3*pow(1000.0-p, 3.0));
        }

        return d2gdtdp - d2fdtdp;
    }
}

- (double)D2gSolventDp2AtT:(double)t andP:(double)p {
    if ((fabs(t-tLast) < VERY_SMALL) || (fabs(p-pLast) < VERY_SMALL)) [self updateDensityPropertiesWithT:t andP:p];
    double   rho    = dielectricConstantFunctions.rhoLast;
    double  DrhoDp  = dielectricConstantFunctions.DrhoDpLast;
    double D2rhoDp2 = dielectricConstantFunctions.D2rhoDp2Last;
    double tc = t - 273.15;

    if (rho >= 1.0) return 0.0;
    else {
        double aG = agP + agPP*tc + agPPP*tc*tc;
        double bG = bgP + bgPP*tc + bgPPP*tc*tc;

        double d2gdp2 = aG*bG*(bG-1.0)*pow(1.0-rho, bG-2.0)*DrhoDp*DrhoDp
                      - aG*bG*pow(1.0-rho, bG-1.0)*D2rhoDp2;

        double d2fdp2 = 0.0;
        if ((p <= 1000) && (tc >= 155.0) && (tc <= 355.0)) {
            d2fdp2 += (pow((tc-155.0)/300.0, 4.8) + ag1*pow((tc-155.0)/300.0, 16.0))
                    * (6.0*ag2*(1000.0-p) + 12.0*ag3*pow(1000.0-p, 2.0));
        }

        return d2gdp2 - d2fdp2;
    }
}

#pragma mark -
#pragma mark Omega and its derivatives

#define CALtoJOULES (4.184)

// Shock et al., 1992, J. Chem. Soc. Faraday Trans. 88(6) 803-826
#define ETA         (1.66027e5) // Angstroms-cal/mole
#define RErefHplus  (3.082)     // Angstroms
#define CALtoJOULES (4.184)

- (double)omegaAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    if ((z == 0.0) || (p > 6000.0)) return omegaRef;

    NSAssert((omegaRef/(ETA*CALtoJOULES) + z/RErefHplus) != 0.0,
             @"Omega reference value prohibits calculation of effective electrostatic radius");
    double reRef = z*z/(omegaRef/(ETA*CALtoJOULES) + z/RErefHplus);

    double g   = [self gSolventAtT:t andP:p];

    double omega = (ETA*CALtoJOULES)*z*(z/(reRef + fabs(z)*g) - 1.0/(RErefHplus + g));
    return omega;
}

- (double)DomegaDtAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    if ((z == 0.0) || (p > 6000.0)) return 0.0;

    NSAssert((omegaRef/(ETA*CALtoJOULES) + z/RErefHplus) != 0.0,
             @"Omega reference value prohibits calculation of effective electrostatic radius");
    double reRef = z*z/(omegaRef/(ETA*CALtoJOULES) + z/RErefHplus);

    double g    = [self gSolventAtT:t andP:p];
    double dgdt = [self DgSolventDtAtT:t andP:p];

    double DomegaDt = (ETA*CALtoJOULES)*z*(-z*fabs(z)*dgdt/pow(reRef + fabs(z)*g, 2.0) + dgdt/pow(RErefHplus + g, 2.0));
    return DomegaDt;
}

- (double)DomegaDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    if ((z == 0.0) || (p > 6000.0)) return 0.0;

    NSAssert((omegaRef/(ETA*CALtoJOULES) + z/RErefHplus) != 0.0,
             @"Omega reference value prohibits calculation of effective electrostatic radius");
    double reRef = z*z/(omegaRef/(ETA*CALtoJOULES) + z/RErefHplus);

    double g    = [self gSolventAtT:t andP:p];
    double dgdp = [self DgSolventDpAtT:t andP:p];

    double DomegaDp = (ETA*CALtoJOULES)*z*(-z*fabs(z)*dgdp/pow(reRef + fabs(z)*g, 2.0) + dgdp/pow(RErefHplus + g, 2.0));
    return DomegaDp;
}

- (double)D2omegaDt2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    if ((z == 0.0) || (p > 6000.0)) return 0.0;

    NSAssert((omegaRef/(ETA*CALtoJOULES) + z/RErefHplus) != 0.0,
             @"Omega reference value prohibits calculation of effective electrostatic radius");
    double reRef = z*z/(omegaRef/(ETA*CALtoJOULES) + z/RErefHplus);

    double g      = [self gSolventAtT:t andP:p];
    double dgdt   = [self DgSolventDtAtT:t andP:p];
    double d2gdt2 = [self D2gSolventDt2AtT:t andP:p];

    double D2omegaDt2 = (ETA*CALtoJOULES)*z*(
                               -z*fabs(z)*d2gdt2/pow(reRef + fabs(z)*g, 2.0) + 2.0*z*z*z*dgdt*dgdt/pow(reRef + fabs(z)*g, 3.0)
                               + d2gdt2/pow(RErefHplus + g, 2.0) - 2.0*dgdt*dgdt/pow(RErefHplus + g, 3.0)
                               );
    return D2omegaDt2;
}

- (double)D2omegaDtDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    if ((z == 0.0) || (p > 6000.0)) return 0.0;

    NSAssert((omegaRef/(ETA*CALtoJOULES) + z/RErefHplus) != 0.0,
             @"Omega reference value prohibits calculation of effective electrostatic radius");
    double reRef = z*z/(omegaRef/(ETA*CALtoJOULES) + z/RErefHplus);

    double g       = [self gSolventAtT:t andP:p];
    double dgdt    = [self DgSolventDtAtT:t andP:p];
    double dgdp    = [self DgSolventDpAtT:t andP:p];
    double d2gdtdp = [self D2gSolventDtDpAtT:t andP:p];

    double D2omegaDtDp = (ETA*CALtoJOULES)*z*(
                                -z*fabs(z)*d2gdtdp/pow(reRef + fabs(z)*g, 2.0) + 2.0*z*z*z*dgdt*dgdp/pow(reRef + fabs(z)*g, 3.0)
                                + d2gdtdp/pow(RErefHplus + g, 2.0) - 2.0*dgdt*dgdp/pow(RErefHplus + g, 3.0)
                                );
    return D2omegaDtDp;
}

- (double)D2omegaDp2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    if ((z == 0.0) || (p > 6000.0)) return 0.0;

    NSAssert((omegaRef/(ETA*CALtoJOULES) + z/RErefHplus) != 0.0,
             @"Omega reference value prohibits calculation of effective electrostatic radius");
    double reRef = z*z/(omegaRef/(ETA*CALtoJOULES) + z/RErefHplus);

    double g      = [self gSolventAtT:t andP:p];
    double dgdp   = [self DgSolventDpAtT:t andP:p];
    double d2gdp2 = [self D2gSolventDp2AtT:t andP:p];

    double D2omegaDp2 = (ETA*CALtoJOULES)*z*(
                               -z*fabs(z)*d2gdp2/pow(reRef + fabs(z)*g, 2.0) + 2.0*z*z*z*dgdp*dgdp/pow(reRef + fabs(z)*g, 3.0)
                               + d2gdp2/pow(RErefHplus + g, 2.0) - 2.0*dgdp*dgdp/pow(RErefHplus + g, 3.0)
                               );
    return D2omegaDp2;
}

#pragma mark -
#pragma mark Some higher order numerical derivatives of the omega function

#define EPS 1.0e-6

- (double)D3omegaDt2DpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    [dielectricConstantFunctions loadDensityPropertiesWithT:t*(1.0+EPS) andWithP:p]; // must be forced to recompute rho and derivatives
    double D2omegaDtDpPlus = [self D2omegaDtDpAtT:t*(1.0+EPS) andP:p withOmegaRef:omegaRef andCharge:z];
    [dielectricConstantFunctions loadDensityPropertiesWithT:t*(1.0-EPS) andWithP:p]; // must be forced to recompute rho and derivatives
    double D2omegaDtDpMinus = [self D2omegaDtDpAtT:t*(1.0-EPS) andP:p withOmegaRef:omegaRef andCharge:z];
    return (D2omegaDtDpPlus-D2omegaDtDpMinus)/2.0/(t*EPS);
}

- (double)D3omegaDtDp2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    [dielectricConstantFunctions loadDensityPropertiesWithT:t andWithP:p*(1.0+EPS)]; // must be forced to recompute rho and derivatives
    double D2omegaDtDpPlus = [self D2omegaDtDpAtT:t andP:p*(1.0+EPS) withOmegaRef:omegaRef andCharge:z];
    [dielectricConstantFunctions loadDensityPropertiesWithT:t andWithP:p*(1.0-EPS)]; // must be forced to recompute rho and derivatives
    double D2omegaDtDpMinus = [self D2omegaDtDpAtT:t andP:p*(1.0-EPS) withOmegaRef:omegaRef andCharge:z];
    return (D2omegaDtDpPlus-D2omegaDtDpMinus)/2.0/(p*EPS);
}

- (double)D3omegaDt3AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    [dielectricConstantFunctions loadDensityPropertiesWithT:t*(1.0+EPS) andWithP:p]; // must be forced to recompute rho and derivatives
    double D2omegaDt2Plus = [self D2omegaDt2AtT:t*(1.0+EPS) andP:p withOmegaRef:omegaRef andCharge:z];
    [dielectricConstantFunctions loadDensityPropertiesWithT:t*(1.0-EPS) andWithP:p]; // must be forced to recompute rho and derivatives
    double D2omegaDt2Minus = [self D2omegaDt2AtT:t*(1.0-EPS) andP:p withOmegaRef:omegaRef andCharge:z];
    return (D2omegaDt2Plus-D2omegaDt2Minus)/2.0/(t*EPS);
}

- (double)D3omegaDp3AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    [dielectricConstantFunctions loadDensityPropertiesWithT:t andWithP:p*(1.0+EPS)]; // must be forced to recompute rho and derivatives
    double D2omegaDp2Plus = [self D2omegaDp2AtT:t andP:p*(1.0+EPS) withOmegaRef:omegaRef andCharge:z];
    [dielectricConstantFunctions loadDensityPropertiesWithT:t andWithP:p*(1.0-EPS)]; // must be forced to recompute rho and derivatives
    double D2omegaDp2Minus = [self D2omegaDp2AtT:t andP:p*(1.0-EPS) withOmegaRef:omegaRef andCharge:z];
    return (D2omegaDp2Plus-D2omegaDp2Minus)/2.0/(p*EPS);
}

#pragma mark -
#pragma mark StoichiometricPhaseProtocol protocol

- (double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    double a1 = self.a1HKF;
    double a2 = self.a2HKF;
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double c1 = self.c1HKF;
    double c2 = self.c2HKF;
    double gTrPr = self.deltaGibbsFreeEnergyOfFormationInTheReferenceState;
    double sTrPr = self.entropyInTheReferenceState;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double g = gTrPr - sTrPr*(t-TR) - c1*(t*log(t/TR) - t + TR)
             - c2*(((1.0/(t-THETA) - 1.0/(TR-THETA))*(THETA-t)/THETA - t*log(TR*(t-THETA)/t/(TR-THETA))/THETA/THETA))
             + a1*(p-PR) + a2*log((CapPSI+p)/(CapPSI+PR))
             + (1.0/(t-THETA))*(a3*(p-PR) + a4*log((CapPSI+p)/(CapPSI+PR)));

    double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
    if (z == 0) g += omegaRef*(YatTRPR*(t-TR) + 1.0/epsilon - 1.0/epsilonTRPR);
    else {
        double omega   = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        g += omega*(1.0/epsilon - 1.0) - omegaRef*(1.0/epsilonTRPR - 1.0) + omegaRef*YatTRPR*(t-TR);
    }
    return g;
}
- (double)getEnthalpyFromT:(double)t andP:(double)p {
    double a1 = self.a1HKF;
    double a2 = self.a2HKF;
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double c1 = self.c1HKF;
    double c2 = self.c2HKF;
    double hTrPr = self.deltaEnthalpyOfFormationInTheReferenceState;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double h = hTrPr + c1*(t - TR) - c2*(1.0/(t-THETA) - 1.0/(TR-THETA))
             + a1*(p-PR) + a2*log((CapPSI+p)/(CapPSI+PR))
             + ((2.0*t-THETA)/pow(t-THETA, 2.0))*(a3*(p-PR) + a4*log((CapPSI+p)/(CapPSI+PR)));

    double Y = [dielectricConstantFunctions YfromT:t andP:p];
    double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
    if (z == 0) h += omegaRef*(Y*t - YatTRPR*TR + 1.0/epsilon - 1.0/epsilonTRPR);
    else {
        double  omega   = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDt = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        h += omega*(1.0/epsilon - 1.0) + omega*t*Y - t*(1.0/epsilon - 1.0)*DomegaDt
           - omegaRef*(1.0/epsilonTRPR - 1.0) - omegaRef*YatTRPR*TR;
    }
    return h;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double c1 = self.c1HKF;
    double c2 = self.c2HKF;
    double sTrPr = self.entropyInTheReferenceState;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double s = sTrPr + c1*log(t/TR) - (c2/THETA)*(1.0/(t-THETA) - 1.0/(TR-THETA) + (1.0/THETA)*log(TR*(t-THETA)/t/(TR-THETA)))
             + pow(1.0/(t-THETA), 2.0)*(a3*(p-PR) + a4*log((CapPSI+p)/(CapPSI+PR)));

    double Y = [dielectricConstantFunctions YfromT:t andP:p];
    if (z == 0) s += omegaRef*(Y - YatTRPR);
    else {
        double epsilon  = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega   = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDt = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        s += omega*Y - (1.0/epsilon - 1.0)*DomegaDt - omegaRef*YatTRPR;
    }
    return s;
}

- (double)getHeatCapacityFromT:(double)t andP:(double)p {
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double c1 = self.c1HKF;
    double c2 = self.c2HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double cp = c1 + c2/pow(t-THETA, 2.0) - (2.0*t/pow(t-THETA, 3.0))*(a3*(p-PR) + a4*log((CapPSI+p)/(CapPSI+PR)));

    double X = [dielectricConstantFunctions XfromT:t andP:p];
    if (z == 0.0) cp += omegaRef*t*X;
    else {
        double Y = [dielectricConstantFunctions YfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double   omega    = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double  DomegaDt  = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDt2 = [self D2omegaDt2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        cp += omega*t*X + 2.0*t*Y*DomegaDt - t*(1.0/epsilon - 1.0)*D2omegaDt2;

    }
    return cp;
}

- (double)getDcpDtFromT:(double)t andP:(double)p {
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double c2 = self.c2HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double DcpDt = -2.0*c2/pow(t-THETA, 3.0) + (6.0*t/pow(t-THETA, 4.0))*(a3*(p-PR) + a4*log((CapPSI+p)/(CapPSI+PR)));

    double X    = [dielectricConstantFunctions XfromT:t andP:p];
    double dXdT = [dielectricConstantFunctions dXdTfromT:t andP:p];
    if (z == 0.0) DcpDt += omegaRef*X + omegaRef*t*dXdT;
    else {
        double Y = [dielectricConstantFunctions YfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double   omega    = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double  DomegaDt  = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDt2 = [self D2omegaDt2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D3omegaDt3 = [self D3omegaDt3AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        DcpDt += omega*X + omega*t*dXdT + 2.0*Y*DomegaDt + 2.0*t*X*DomegaDt + 2.0*t*Y*D2omegaDt2
            - (1.0/epsilon - 1.0)*D2omegaDt2 + t*(Y/epsilon)*D2omegaDt2 - t*(1.0/epsilon - 1.0)*D3omegaDt3;

    }
    return DcpDt;
}

- (double)getVolumeFromT:(double)t andP:(double)p {
    double a1 = self.a1HKF;
    double a2 = self.a2HKF;
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double v = a1 + a2/(CapPSI+p) + (a3 + a4/(CapPSI+p))/(t-THETA);

    double Q = [dielectricConstantFunctions QfromT:t andP:p];
    if (z == 0.0) v -= omegaRef*Q;
    else {
        double epsilon  = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega   = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDp = [self DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        v += -omega*Q + (1.0/epsilon - 1.0)*DomegaDp;
    }
    return v;
}

- (double)getDvDtFromT:(double)t andP:(double)p {
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double dvdt = - (a3 + a4/(CapPSI+p))/pow(t-THETA, 2.0);

    double U = [dielectricConstantFunctions UfromT:t andP:p];
    if (z == 0.0) dvdt -= omegaRef*U;
    else {
        double Q = [dielectricConstantFunctions QfromT:t andP:p];
        double Y = [dielectricConstantFunctions YfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega      = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDp    = [self DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDt    = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDtDp = [self D2omegaDtDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        dvdt += -omega*U -Q*DomegaDt -Y*DomegaDp + (1.0/epsilon - 1.0)*D2omegaDtDp;
    }
    return dvdt;
}

- (double)getDvDpFromT:(double)t andP:(double)p {
    double a2 = self.a2HKF;
    double a4 = self.a4HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double dvdp = - (a2 + a4/(t-THETA))/pow(CapPSI+p, 2.0);

    double N = [dielectricConstantFunctions NfromT:t andP:p];
    if (z == 0.0) dvdp -= omegaRef*N;
    else {
        double Q = [dielectricConstantFunctions QfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega     = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDp   = [self DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDp2 = [self D2omegaDp2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        dvdp += -omega*N -2.0*Q*DomegaDp + (1.0/epsilon - 1.0)*D2omegaDp2;
    }
    return dvdp;
}

- (double)getD2vDt2FromT:(double)t andP:(double)p {
    double a3 = self.a3HKF;
    double a4 = self.a4HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double d2vdt2 =  2.0*(a3 + a4/(CapPSI+p))/pow(t-THETA, 3.0);

    double dUdT = [dielectricConstantFunctions dUdTfromT:t andP:p];
    if (z == 0.0) d2vdt2 -= omegaRef*dUdT;
    else {
        double U = [dielectricConstantFunctions UfromT:t andP:p];
        double Q = [dielectricConstantFunctions QfromT:t andP:p];
        double Y = [dielectricConstantFunctions YfromT:t andP:p];
        double X = [dielectricConstantFunctions XfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega       = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDp     = [self DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDt     = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDtDp  = [self D2omegaDtDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDt2   = [self D2omegaDt2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D3omegaDt2Dp = [self D3omegaDt2DpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        d2vdt2 += -DomegaDt*U - omega*dUdT - U*DomegaDt - Q*D2omegaDt2 - X*DomegaDp - Y*D2omegaDtDp
                - (Y/epsilon)*D2omegaDtDp + (1.0/epsilon - 1.0)*D3omegaDt2Dp;
    }
    return d2vdt2;
}

- (double)getD2vDtDpFromT:(double)t andP:(double)p {
    double a4 = self.a4HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double d2vdtdp = (a4/pow(CapPSI+p, 2.0))/pow(t-THETA, 2.0);

    double dUdP = [dielectricConstantFunctions dUdPfromT:t andP:p];
    if (z == 0.0) d2vdtdp -= omegaRef*dUdP;
    else {
        double Q = [dielectricConstantFunctions QfromT:t andP:p];
        double N = [dielectricConstantFunctions NfromT:t andP:p];
        double U = [dielectricConstantFunctions UfromT:t andP:p];
        double Y = [dielectricConstantFunctions YfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega       = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDp     = [self DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDt     = [self DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDp2   = [self D2omegaDp2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDtDp  = [self D2omegaDtDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D3omegaDtDp2 = [self D3omegaDtDp2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        d2vdtdp += -DomegaDp*U - omega*dUdP  - N*DomegaDt  - Q*D2omegaDtDp - U*DomegaDp - Y*D2omegaDp2
                 - (Q/epsilon)*D2omegaDtDp + (1.0/epsilon - 1.0)*D3omegaDtDp2;
    }
    return d2vdtdp;
    return 0.0;
}

- (double)getD2vDp2FromT:(double)t andP:(double)p {
    double a2 = self.a2HKF;
    double a4 = self.a4HKF;
    double omegaRef = self.omegaHKF;
    double z = self.charge;

    double d2vdp2 =  2.0*(a2 + a4/(t-THETA))/pow(CapPSI+p, 3.0);

    double dNdP = [dielectricConstantFunctions dNdPfromT:t andP:p];
    if (z == 0.0) d2vdp2 -= omegaRef*dNdP;
    else {
        double N = [dielectricConstantFunctions NfromT:t andP:p];
        double Q = [dielectricConstantFunctions QfromT:t andP:p];
        double epsilon = [dielectricConstantFunctions epsilonFromT:t andP:p];
        double  omega     = [self omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double DomegaDp   = [self DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D2omegaDp2 = [self D2omegaDp2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        double D3omegaDp3 = [self D3omegaDp3AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
        d2vdp2 += -DomegaDp*N - omega*dNdP - 2.0*N*DomegaDp - 2.0*Q*D2omegaDp2
                - (Q/epsilon)*D2omegaDp2 + (1.0/epsilon - 1.0)*D3omegaDp3;
    }
    return d2vdp2;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kdeltaGibbsFreeEnergyOfFormationInTheReferenceState = @"deltaGibbsFreeEnergyOfFormationInTheReferenceState";
static NSString *kdeltaEnthalpyOfFormationInTheReferenceState        = @"deltaEnthalpyOfFormationInTheReferenceState";
static NSString *kentropyInTheReferenceState                         = @"entropyInTheReferenceState";
static NSString *kvolumeInTheReferenceState                          = @"volumeInTheReferenceState";
static NSString *kheatCapacityInTheReferenceState                    = @"heatCapacityInTheReferenceState";

static NSString *ka1HKF    = @"a1HKF";
static NSString *ka2HKF    = @"a2HKF";
static NSString *ka3HKF    = @"a3HKF";
static NSString *ka4HKF    = @"a4HKF";
static NSString *kc1HKF    = @"c1HKF";
static NSString *kc2HKF    = @"c2HKF";
static NSString *komegaHKF = @"omegaHKF";
static NSString *kcharge   = @"charge";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super init])) {
        dielectricConstantFunctions = [[DEWDielectricConstant alloc] init];

        _deltaGibbsFreeEnergyOfFormationInTheReferenceState = (double) [aDecoder decodeDoubleForKey:kdeltaGibbsFreeEnergyOfFormationInTheReferenceState];
        _deltaEnthalpyOfFormationInTheReferenceState = (double) [aDecoder decodeDoubleForKey:kdeltaEnthalpyOfFormationInTheReferenceState];
        _entropyInTheReferenceState = (double) [aDecoder decodeDoubleForKey:kentropyInTheReferenceState];
        _volumeInTheReferenceState = (double) [aDecoder decodeDoubleForKey:kvolumeInTheReferenceState];
        _heatCapacityInTheReferenceState = (double) [aDecoder decodeDoubleForKey:kheatCapacityInTheReferenceState];
        _a1HKF = (double) [aDecoder decodeDoubleForKey:ka1HKF];
        _a2HKF = (double) [aDecoder decodeDoubleForKey:ka2HKF];
        _a3HKF = (double) [aDecoder decodeDoubleForKey:ka3HKF];
        _a4HKF = (double) [aDecoder decodeDoubleForKey:ka4HKF];
        _c1HKF = (double) [aDecoder decodeDoubleForKey:kc1HKF];
        _c2HKF = (double) [aDecoder decodeDoubleForKey:kc2HKF];
        _omegaHKF = (double) [aDecoder decodeDoubleForKey:komegaHKF];
        _charge = (double) [aDecoder decodeDoubleForKey:kcharge];

        tLast = -1.0;
        pLast = -1.0;
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:_deltaGibbsFreeEnergyOfFormationInTheReferenceState forKey:kdeltaGibbsFreeEnergyOfFormationInTheReferenceState];
        [aCoder encodeDouble:_deltaEnthalpyOfFormationInTheReferenceState forKey:kdeltaEnthalpyOfFormationInTheReferenceState];
        [aCoder encodeDouble:_entropyInTheReferenceState forKey:kentropyInTheReferenceState];
        [aCoder encodeDouble:_volumeInTheReferenceState forKey:kvolumeInTheReferenceState];
        [aCoder encodeDouble:_heatCapacityInTheReferenceState forKey:kheatCapacityInTheReferenceState];
        [aCoder encodeDouble:_a1HKF forKey:ka1HKF];
        [aCoder encodeDouble:_a2HKF forKey:ka2HKF];
        [aCoder encodeDouble:_a3HKF forKey:ka3HKF];
        [aCoder encodeDouble:_a4HKF forKey:ka4HKF];
        [aCoder encodeDouble:_c1HKF forKey:kc1HKF];
        [aCoder encodeDouble:_c2HKF forKey:kc2HKF];
        [aCoder encodeDouble:_omegaHKF forKey:komegaHKF];
        [aCoder encodeDouble:_charge forKey:kcharge];

    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class HKFspeciesProperties only supports NSKeyedArchiver coders."];
}

+ (BOOL)supportsSecureCoding {
    return YES;
}

@end
