//
//  StoichiometricPhaseProtocol.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@protocol StoichiometricPhaseProtocol

/**
 T (K), P (bars) => Gibbs free energy (J)
 */
-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => enthalpy (J)
 */
-(double)getEnthalpyFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => entropy (J/K)
 */
-(double)getEntropyFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => isobaric heat capacity (J/K)
 */
-(double)getHeatCapacityFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => d(isobaric heat capacity)/dT (J/K^2)
 */
-(double)getDcpDtFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => volume (J/bar)
 */
-(double)getVolumeFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => d(volume)/dT (J/bar-K)
 */
-(double)getDvDtFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => d(volume)/dP (J/bar^2)
 */
-(double)getDvDpFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => d^2(volume)/dT^2 (J/bar-K^2)
 */
-(double)getD2vDt2FromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => d^2(volume)/dTdP (J/bar^2-K)
 */
-(double)getD2vDtDpFromT:(double)t andP:(double)p;
/**
 T (K), P (bars) => d^2(volume)/dP^2 (J/bar^3)
 */
-(double)getD2vDp2FromT:(double)t andP:(double)p;

@optional

/**
 T (K), P (bars) => chemical potential (J)
 */
-(double)getChemicalPotentialFromT:(double)t andP:(double)p;
/**
 => formulae of the phase
 */
-(NSString *)getFormulaFromInternalVariables;

@end
