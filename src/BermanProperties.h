//
//  BermanProperties.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"


/**
 Computes thermodynamic properties using a Berman-type database model.
 */

@interface BermanProperties : PhaseBase <NSSecureCoding,StoichiometricPhaseProtocol> {
	@protected
    double h, s, k0, k1, k2, k3, l1, l2, Tt, deltaH, v0, v1, v2, v3, v4;
	double tr, pr, trl;
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
@property (readwrite, nonatomic) double k0;
@property (readwrite, nonatomic) double k1;
@property (readwrite, nonatomic) double k2;
@property (readwrite, nonatomic) double k3;
@property (readwrite, nonatomic) double l1;
@property (readwrite, nonatomic) double l2;
@property (readwrite, nonatomic) double Tt;
@property (readwrite, nonatomic) double deltaH;
@property (readwrite, nonatomic) double v0;
@property (readwrite, nonatomic) double v1;
@property (readwrite, nonatomic) double v2;
@property (readwrite, nonatomic) double v3;
@property (readwrite, nonatomic) double v4;

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
			v0:(double)v0In
			v1:(double)v1In
			v2:(double)v2In
			v3:(double)v3In
			v4:(double)v4In;

-(id)initWithH:(double)hIn
			 S:(double)sIn
			k0:(double)k0In
			k1:(double)k1In
			k2:(double)k2In
			k3:(double)k3In
			v0:(double)v0In
			v1:(double)v1In
			v2:(double)v2In
			v3:(double)v3In
			v4:(double)v4In;

-(id)initWithH:(double)hIn
			 S:(double)sIn
			k0:(double)k0In
			k1:(double)k1In
			k2:(double)k2In
			k3:(double)k3In;

/**
 Set the reference temperature for the Cp integration (default 298.15 K)
 */
-(void)setTr:(double)trIn;
/**
 Set the reference pressure (default 1.0 bar)
 */
-(void)setPr:(double)prIn;
/**
 Set the reference temperature for the lambda Cp correction (default 298.15 K)
 */
-(void)setTrl:(double)trlIn;

@end
