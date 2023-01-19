//
//  SolutionPhaseProtocol.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/14/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@class DoubleVector;
@class DoubleMatrix;
@class DoubleTensor;

@protocol SolutionPhaseProtocol

/**
 By default all thermodynamic quantities are calculated to include endmember properties (yesForMixing = NO).
 Calling this function with YES as an argument converts output to mixing quantities.
 */
-(void)setResultsToMixingQuantities:(BOOL)yesForMixing;

/**
 Retrieves the number of endmember components in the system
 */
-(NSUInteger)numberOfSolutionComponents;

/**
 Retrieves superclass instance of PhaseBase object for component at specified index
 */
-(id)componentAtIndex:(NSUInteger)index;

/**
 Moles of endmember components => validity of input values
 */
-(BOOL)testPermissibleValuesOfComponents:(double *)m;

/**
 Moles of elements (standard order) => Moles of end member components of the phase
 */
-(DoubleVector *)convertElementsToMoles:(double *)e;
/**
 Moles of elements (standard order) => Total moles of end member components of the phase
 */
-(double)convertElementsToTotalMoles:(double *)e;
/**
 Moles of elements (standard order) => Total mass of the phase (g)
 */
-(double)convertElementsToTotalMass:(double *)e;


/**
 Moles of endmember components => Mole fractions of endmember components
 */
-(DoubleVector *)convertMolesToMoleFractions:(double *)m;
/**
 Moles of endmember components => Moles of elements (standard order)
 */
-(DoubleVector *)convertMolesToElements:(double *)m;
/**
 Moles of endmember components => Molar sum
 */
-(double)totalMolesFromMolesOfComponents:(double *)m;

/**
 Moles of components, T (K), P (bars) => activities of endmember components
 */
-(DoubleVector *)getActivityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => chemical potentials of endmember components (J)
 */
-(DoubleVector *)getChemicalPotentialFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(activities of endmember components)/d(Moles of components)
 */
-(DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;

/**
 Moles of components, T (K), P (bars) => Gibbs free energy (J)
 */
-(double)getGibbsFreeEnergyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(Gibbs free energy)/d(Moles of components) (J)
 */
-(DoubleVector *)getDgDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d^2(Gibbs free energy)/d(Moles of components)^2 (J)
 */
-(DoubleMatrix *)getD2gDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d^3(Gibbs free energy)/d(Moles of components)^3 (J)
 */
-(DoubleTensor *)getD3gDm3FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;

/**
 Moles of components, T (K), P (bars) => enthalpy (J)
 */
-(double)getEnthalpyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;

/**
 Moles of components, T (K), P (bars) => entropy (J/K)
 */
-(double)getEntropyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(entropy)/d(Moles of components) (J/K)
 */
-(DoubleVector *)getDsDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d^2(entropy)/d(Moles of components)^2 (J/K)
 */
-(DoubleMatrix *)getD2sDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => isobaric heat capacity (J/K)
 */
-(double)getHeatCapacityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(isobaric heat capacity)/dT (J/K^2)
 */
-(double)getDcpDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(isobaric heat capacity)/d(Moles of components) (J/K)
 */
-(DoubleVector *)getDCpDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => volume (J/bar)
 */
-(double)getVolumeFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(volume)/d(Moles of components) (J/bar)
 */
-(DoubleVector *)getDvDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d^2(volume)/d(Moles of components)^2 (J/bar)
 */
-(DoubleMatrix *)getD2vDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(volume)/dT (J/bar-K)
 */
-(double)getDvDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d(volume)/dP (J/bar^2)
 */
-(double)getDvDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d2(volume)/dT^2 (J/bar-K^2)
 */
-(double)getD2vDt2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d2(volume)/dTdP (J/bar^2-K)
 */
-(double)getD2vDtDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d2(volume)/dP^2 (J/bar^3)
 */
-(double)getD2vDp2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d2(volume)/d(Moles of components)dT (J/bar-K)
 */
-(DoubleVector *)getD2vDmDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => d2(volume)/d(Moles of components)dP (J/bar^2)
 */
-(DoubleVector *)getD2vDmDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Moles of components, T (K), P (bars) => formulae as an NSString object
 */
-(NSString *)getFormulaFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Retrieves the number of species (dependent endmembers with positive mole fractions) in the soluton
 */
-(NSUInteger)numberOfSolutionSpecies;
/**
 Retrieves the name of the solution species at the specified index
 */
-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index;
/**
 Moles of solution species => moles of endmember components
 */
-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies;
/**
 Retrieves an elemental stoichiometry vector for the species at the specified index
 */
-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index;
/**
 Moles of components, T (K), P (bars) => chemical potentials of solution species (J)
 */
-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
/**
 Chemical potentials of endmembers, T (K), P (bars) => chemical affinity and composition of the phase, etc.

 @param chemicalPotentials
   chemical potentials of endmember components in the solution.  A zero entry indicates a component is absent.

 @param t
   temperature in K

 @param p
   pressure in bars

 @return
   NSArray output structure:
   - (0):      NSNumber object wrapping a double - chemical affinity
   - (1...NA): NSnumber object wrapping a double - X[0] - X[NA-1], mole fraction compositional variables
   - (NA+1):   NSNumber object wapping a BOOL - convergence flag
   - (NA+2):   NSNumber object wrapping an NSUInteger - iteration count
   - (NA+3):   NSNumber object wrapping a double - number of atoms in formula unit to scale affinity
   - (NA+4):   NSNumber object wrapping a double - approximate error in calculated affinity
 */
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p;

@optional

-(void)incrementInstanceCountOfPhase;
-(void)decrementInstanceCountOfPhase;
-(NSString *)nameOfPhaseWithComposition:(double *)refMoles;
-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p;

@end
