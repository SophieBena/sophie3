//
//  HollandAndPowellProperties.h
//  PhaseMELTSobjC
//
//  Created by Mark Ghiorso on 6/6/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"

@interface HollandAndPowellProperties : PhaseBase <NSSecureCoding,StoichiometricPhaseProtocol> {
@protected
    double h, s, a, b, c, d, Tc0, Smax, Vmax, v0, a0, K;
    double tr, pr;
    double Q2Tr, hpTr, spTr, vpTr;
}

/**
 These class functions set the behavior of all instances of the class, including those previously instantiated.
 Enables or Disables the convention of using the Gibbs free energy of formation from the elements as the
 reference state. The default is to use the enthalpy of formation from the elements as the refence state.

 Default, Disabled: G(T,P) = deltaH(Tr,Pr) + integral(Tr->T) Cp -      T S(Tr,Pr) - T integral(Tr->T) Cp/T + integral(Pr->P) V
 Enabled: G(T,P) = deltaG(Tr,Pr) + integral(Tr->T) Cp - (T-Tr) S(Tr,Pr) - T integral(Tr->T) Cp/T + integral(Pr->P) V

 where deltaG(Tr,Pr) = deltaH(Tr,Pr) - Tr S(Tr,Pr) + Tr Selements(Tr,Pr)

 or, equivalently:

 Default, Disabled: G(T,P) = deltaH(Tr,Pr) + integral(Tr->T) Cp - T S(Tr,Pr) - T integral(Tr->T) Cp/T + integral(Pr->P) V
 Enabled: G(T,P) = deltaH(Tr,Pr) + integral(Tr->T) Cp - T S(Tr,Pr) - T integral(Tr->T) Cp/T + integral(Pr->P) V
 + Tr Selements(Tr,Pr)

 so, the only difference is the constant + Tr Selements(Tr,Pr)
 */
+(void)enableGibbsFreeEnergyReferenceStateUsed;
+(void)disableGibbsFreeEnergyReferenceStateUsed;

@property (readwrite, nonatomic) double h;
@property (readwrite, nonatomic) double s;
@property (readwrite, nonatomic) double a;
@property (readwrite, nonatomic) double b;
@property (readwrite, nonatomic) double c;
@property (readwrite, nonatomic) double d;
@property (readwrite, nonatomic) double Tc0;
@property (readwrite, nonatomic) double Smax;
@property (readwrite, nonatomic) double Vmax;
@property (readwrite, nonatomic) double v0;
@property (readwrite, nonatomic) double a0;
@property (readwrite, nonatomic) double K;

@property (readwrite) BOOL useHandPexcessVolumes;

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
             K:(double)KIn;

/**
 Set the reference temperature for the Cp integration (default 298.15 K)
 */
-(void)setTr:(double)trIn;
/**
 Set the reference pressure (default 1.0 bar)
 */
-(void)setPr:(double)prIn;

@end
