//
//  StixrudeSolutionPhase.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/10/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "StixrudeSolutionPhase.h"
#import "StixrudeProperties.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

#define NA 6
#define NR 5
#define NS 16

@implementation StixrudeSolutionPhase

@synthesize numberOfIterationsAllowedInSaturationMethod;
@synthesize minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod;

#pragma mark -
#pragma mark instance methods

-(id)initWithPhaseName:(NSString *)name
		   withSpecies:(NSArray *)speciesProperties
	withSpeciesWeights:(NSArray *)speciesWeights
	withSpeciesXfactor:(NSArray *)speciesXfactor
			withWarray:(NSArray *)wIn
			 withSites:(NSArray *)sites
			WithNAtoms:(NSUInteger)nAtomIn
  withElementConvertor:(void (^)(double *, double *))convertor
 andWithFormulaDisplay:(NSString *(^)(double *))display
andWithSaturationStateAdjuster:(void (^)(double *, double *))adjuster {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(StixrudeSolutionPhase) with %@ ...", name);

		[self setPhaseName:name];
		computeMixingQuantities = NO;

		endmembers = [NSArray arrayWithArray:speciesProperties];
		na = [endmembers count];
		nr = na - 1;
		if (na > NA) {
			[NSException raise:@"StixrudeSolutionPhase initialization ERROR" format:@"Too many endmember species"];
			return nil;
		}
		for (NSUInteger i=0; i<na; i++) dArray[i] = [[speciesWeights objectAtIndex:i] doubleValue];
		for (NSUInteger i=0; i<na; i++) adjustmentCoefficient[i] = [[speciesXfactor objectAtIndex:i] doubleValue];

		if ([wIn count] < na*(na-1)/2) {
			[NSException raise:@"StixrudeSolutionPhase initialization ERROR" format:@"Too few interaction parameters"];
			return nil;
		}

		for (NSUInteger i=0; i<NA; i++) for (NSUInteger j=0; j<NA; j++) wArray[i][j] = 0.0;
		for (NSUInteger i=1, k=0; i<na; i++) {
			for (NSUInteger j=0; j<i; j++) {
				wArray[i][j] = [[wIn objectAtIndex:k++] doubleValue];
				wArray[j][i] = wArray[i][j];
			}
		}

		ns = [sites count];
		if (ns > NS) {
			[NSException raise:@"StixrudeSolutionPhase initialization ERROR" format:@"Too many sites specified"];
			return nil;
		}
		for (NSUInteger i=0; i<ns; i++) {
			siteMultiplicity[i] = [[[sites objectAtIndex:i] objectAtIndex:0] doubleValue];
			NSArray *coeff = [[sites objectAtIndex:i] objectAtIndex:1];
			if ([coeff count] != na) return nil;
			for (NSUInteger j=0; j<na; j++) {
				siteSpeciesCoefficients[i][j] = [[coeff objectAtIndex:j] doubleValue];
			}
		}

		nAtoms = (double) nAtomIn;
		convertMolesOfElementsToMolesOfSpecies = convertor;
		convertMolesOfSpeciesToDisplayFormula = display;
        correctActivityCoefficients = adjuster;

        numberOfIterationsAllowedInSaturationMethod = 75;
        minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod = 0.1;
	}
	return self;
}

#define R 8.3143

-(NSUInteger)numberOfSolutionSpecies {
	return na;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	return [[endmembers objectAtIndex:index] phaseName];
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
    DoubleVector *mComponentsWrapper = [[DoubleVector alloc] initWithSize:na];
	double *mComponents = [mComponentsWrapper pointerToDouble];
	for (NSUInteger i=0; i<na; i++) mComponents[i] = mSpecies[i];
	return mComponentsWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	return [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
	return [[endmembers objectAtIndex:index] formulaAsElementArray];
}

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
	NSMutableArray *results = [NSMutableArray arrayWithCapacity:na+1];
	double deltaMu[NA], xNz[NA], x[NA], gamma[NA], xLast[NA], affinity = 0.0, reducedAdjustmentCoefficient[NA], gammaLast[NA];
	NSUInteger i, j, nz = 0, index[NA];

	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");

	/* Test code
	double mTest[NA], mRef[NA];
	for (i=0; i<na; i++) mRef[i] = (double) (i+1);
	double gRef = [self getGibbsFreeEnergyFromMolesOfComponents:mRef andT:t andP:p];
	double *dgdmRef = [self getDgDmFromMolesOfComponents:mRef andT:t andP:p];
	double **d2gdm2Ref = [self getD2gDm2FromMolesOfComponents:mRef andT:t andP:p];

	for (i=0; i<na; i++) {
		for (j=0; j<na; j++) mTest[j] = mRef[j];
		mTest[i] = (1.0+sqrt(DBL_EPSILON))*mTest[i];
		double gTest = [self getGibbsFreeEnergyFromMolesOfComponents:mTest andT:t andP:p];
		double dgdm = (gTest-gRef)/(sqrt(DBL_EPSILON)*mTest[i]);
		NSLog(@"...Derivative test for i = %d, dgdm est = %20.13e, dgdm = %20.13e, delta = %20.13e, rel error = %20.13e",
			  i, dgdm, dgdmRef[i], dgdm-dgdmRef[i], (dgdm-dgdmRef[i])/dgdmRef[i]);
	}

	for (i=0; i<na; i++) {
		for (j=0; j<na; j++) {
			for (int k=0; k<na; k++) mTest[k] = mRef[k];
			mTest[j] = (1.0+sqrt(DBL_EPSILON))*mTest[j];
			double *dgdmTest = [self getDgDmFromMolesOfComponents:mTest andT:t andP:p];
			double d2gdm2 = (dgdmTest[i]-dgdmRef[i])/(sqrt(DBL_EPSILON)*mTest[j]);
			NSLog(@"...2nd derivative test for i = %d, j = %d, d2gdm2 est = %20.13e, d2gdm2 = %20.13e, delta = %20.13e, rel error = %20.13e",
				  i, j, d2gdm2, d2gdm2Ref[i][j], d2gdm2-d2gdm2Ref[i][j], (d2gdm2-d2gdm2Ref[i][j])/d2gdm2Ref[i][j]);
		}
	}
	*/

	// Compute solid -> liquid delta mus and deflate composition space
	for (i=0; i<na; i++) {
		if (chemicalPotentials[i] != 0.0) {
			deltaMu[nz] = chemicalPotentials[i] - [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
			index[nz] = i;
			gamma[nz] = 1.0;
			reducedAdjustmentCoefficient[nz] = adjustmentCoefficient[i];
			nz++;
		}
		x[i] = 0.0;
		xLast[i] = 0.0;
	}

	// There are no non-zero chemical potentials, so the phase can never form.  Return an affinity of zero and zero the mole fractions
	if (nz == 0) {
		[results addObject:[NSNumber numberWithDouble:affinity]];
		for (i=0; i<na; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];
		[results addObject:[NSNumber numberWithBool:YES]];          // convergence flag
		[results addObject:[NSNumber numberWithUnsignedInteger:0]]; // iteration count
		[results addObject:[NSNumber numberWithDouble:nAtoms]];     // number of atoms used to scale affinity
		[results addObject:[NSNumber numberWithDouble:0.0]];        // likely error in affinity

		if (debugV) NSLog(@"... Terminated. Trivial case.");
		return [NSArray arrayWithArray:results];
	}

	NSUInteger count = 0;
	BOOL converged = NO;
	double affinityLast;
	do {
		affinityLast = affinity;
		// Solve for mole ractions in the deflated composition space
		if (nz == 1) {
			xNz[0] = 1.0;
			affinity = -(deltaMu[0]-R*t*log(gamma[0]));
		} else {
			double sum = 1.0;
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] = exp(((deltaMu[i]-R*t*reducedAdjustmentCoefficient[i]*log(gamma[i]))
							  -(deltaMu[nz-2]-R*t*reducedAdjustmentCoefficient[nz-2]*log(gamma[nz-2])))/(reducedAdjustmentCoefficient[i]*R*t));
				sum += xNz[i];
			}
			xNz[nz-2] = exp(((deltaMu[nz-2]-R*t*reducedAdjustmentCoefficient[nz-2]*log(gamma[nz-2]))
							 -(deltaMu[nz-1]-R*t*reducedAdjustmentCoefficient[nz-1]*log(gamma[nz-1])))/(reducedAdjustmentCoefficient[nz-2]*R*t));

			xNz[nz-2] /= 1.0 + xNz[nz-2]*sum;
			xNz[nz-1] = 1.0 - xNz[nz-2];
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] *= xNz[nz-2];
				xNz[nz-1] -= xNz[i];
			}
			for (i=0; i<nz; i++) if (xNz[i] <= DBL_EPSILON) xNz[i] = DBL_EPSILON;

			/* compute the chemical affinity (choice of mu[] is arbitrary) */
			affinity = -(deltaMu[0]-R*t*reducedAdjustmentCoefficient[0]*log(gamma[0])) + R*t*reducedAdjustmentCoefficient[0]*log(xNz[0]);
		}

		// Reinflate the solution
		for (i=0; i<nz; i++) x[index[i]] = xNz[i];

		// Determine activity coefficients
        DoubleVector *aWrapper = [self getActivityFromMolesOfComponents:x andT:t andP:p];
		double *a = [aWrapper pointerToDouble];

		for (i=0, j=0; i<na; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];
		if (count > 25) { // pure empiricism
			for (i=0; i<nz; i++) gamma[i] = (gamma[i]+gammaLast[i])/2.0;
		}
		for (i=0; i<nz; i++) gammaLast[i] = gamma[i];
		correctActivityCoefficients(gamma, x);
		if (debugV) {
			for (i=0, j=0; i<na; i++) if (x[i] != 0.0) NSLog(@"... a of %@ is %20.13e with X %20.13e and g %20.13e",
                      [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0], a[i], x[i], gamma[j++]);
            NSLog(@"... a of %@ is %@ with X %@ and g %@ Aff %20.13e delta Aff %20.13e", @"sum                 ",
                  @"                    ", @"                    ", @"                    ", affinity, affinity-affinityLast);
		}
		converged = (fabs(affinity-affinityLast) < minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod);
		count++;

	} while (count < numberOfIterationsAllowedInSaturationMethod && !converged);

	if (debugS) {
		NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
			  converged ? @"YES" : @"NO", [self phaseName], count, fabs(affinity-affinityLast), nAtoms);
		for (i=0, j=0; i<na; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
														 [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
														 gamma[j++], x[i]);
	}
	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");


	[results addObject:[NSNumber numberWithDouble:affinity]];                    // affinity in J
	for (i=0; i<na; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];   // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                     // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];              // iteration count
	[results addObject:[NSNumber numberWithDouble:nAtoms]];                      // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]]; // likely error in affinity

	return [NSArray arrayWithArray:results];
}

// Begin objective-C <SolutionPhaseProtocol> functions

// --> SolutionPhaseProtocol public function
-(void)setResultsToMixingQuantities:(BOOL)yesForMixing {
	computeMixingQuantities = yesForMixing;
}

// --> SolutionPhaseProtocol public function
-(NSUInteger)numberOfSolutionComponents {
	return na;
}

// --> SolutionPhaseProtocol public function
-(id)componentAtIndex:(NSUInteger)index {
	return [endmembers objectAtIndex:index];
}

// Test routine for permissible component numbers

// --> SolutionPhaseProtocol public function
-(BOOL)testPermissibleValuesOfInternalVariables:(double	*)r {
	BOOL result = YES;
	for (NSUInteger i=0; i<nr; i++) result &= (r[i] >= 0.0);
	return result;
}

// --> SolutionPhaseProtocol public function
-(BOOL)testPermissibleValuesOfComponents:(double *)m {
	BOOL result = YES;
	for (NSUInteger i=0; i<na; i++) result &= (m[i] >= 0.0);
	return result;
}

// Component conversion routines

// --> SolutionPhaseProtocol public function
-(DoubleVector *)convertElementsToMoles:(double *)e {
    DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:na];
	double *m = [mWrapper pointerToDouble];
	convertMolesOfElementsToMolesOfSpecies(e, m);
	return mWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)convertElementsToTotalMoles:(double *)e {
	DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:na];
	double *m = [mWrapper pointerToDouble];
	convertMolesOfElementsToMolesOfSpecies(e, m);
	double result = 0.0;
	for (NSUInteger i=0; i<na; i++) result += m[i];
	return result;
}

// --> SolutionPhaseProtocol public function
-(double)convertElementsToTotalMass:(double *)e {
	DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:na];
	double *m = [mWrapper pointerToDouble];
	convertMolesOfElementsToMolesOfSpecies(e, m);
	double result = 0.0;
	for (NSUInteger i=0; i<na; i++) result += m[i]*[[endmembers objectAtIndex:i] mw];
	return result;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)convertMolesToMoleFractions:(double *)m {
    DoubleVector *xWrapper = [[DoubleVector alloc] initWithSize:na];
	double *x = [xWrapper pointerToDouble];
	double mTotal = [self totalMolesFromMolesOfComponents:m];
	for (NSUInteger i=0; i<na; i++) x[i] = m[i]/mTotal;
	return xWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)convertMolesToElements:(double *)m {
    DoubleVector *eWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
	double *e = [eWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) {
        DoubleVector *componentToElementsWrapper = [(StixrudeProperties *)[endmembers objectAtIndex:i] formulaAsElementArray];
		double *componentToElements = [componentToElementsWrapper pointerToDouble];
		for (NSUInteger j=1; j<107; j++) e[j] += m[i]*componentToElements[j];
	}
	return eWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)totalMolesFromMolesOfComponents:(double *)m {
	double result = 0.0;
	for (int i=0; i<na; i++) result += m[i];
	return result;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getActivityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *activityWrapper = [self getDgDmFromMolesOfComponents:m andT:t andP:p];
    double *activity = [activityWrapper pointerToDouble];
	if (!computeMixingQuantities) for (NSUInteger i=0; i<na; i++) activity[i] -= [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
	for (NSUInteger i=0; i<na; i++) activity[i] = (m[i] != 0.0) ? exp(activity[i]/(R*t)) : 0.0;
	return activityWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getChemicalPotentialFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	return [self getDgDmFromMolesOfComponents:m andT:t andP:p];
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *dadmWrapper = [[DoubleMatrix alloc] initWithRowSize:na andWithColumnSize:nr andInitialValue:0.0];
	double **dadm = [dadmWrapper pointerToPointerToDouble];
    DoubleVector *muWrapper = [self getDgDmFromMolesOfComponents:m andT:t andP:p];
	double *mu = [muWrapper pointerToDouble];
    DoubleMatrix *dmudmWrapper = [self getD2gDm2FromMolesOfComponents:m andT:t andP:p];
	double **dmudm = [dmudmWrapper pointerToPointerToDouble];
	if (!computeMixingQuantities) for (NSUInteger i=0; i<na; i++) mu[i] -= [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
	for (NSUInteger i=0; i<na; i++) if (m[i] != 0.0) {
		for (NSUInteger j=0; j<na; j++) if (m[j] != 0.0) dadm[i][j] = dmudm[i][j]*exp(mu[i]/(R*t))/(R*t);
	}
	return dadmWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getGibbsFreeEnergyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];;
	if (mTotal == 0.0) return 0.0;

	double vlWeight = 0.0;
	for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

	double g = 0.0;
	for (NSUInteger i=0; i<na; i++) {
		double phiI = (m[i]/mTotal)*dArray[i]/vlWeight;
		for (NSUInteger j=i+1; j<na; j++) {
			double phiJ = (m[j]/mTotal)*dArray[j]/vlWeight;
			g += phiI*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j]);
		}
	}

	for (NSUInteger i=0; i<ns; i++) {
		double xOnSite = 0.0;
		for (NSUInteger j=0; j<na; j++) xOnSite += siteSpeciesCoefficients[i][j]*m[j]/mTotal;
		g += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*xOnSite*log(xOnSite) : 0.0;
	}

	g *= mTotal;

	if (computeMixingQuantities) return g;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		g += m[i]*[component getGibbsFreeEnergyFromT:t andP:p];
	}
	return g;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDgDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleVector *dgdmWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
	double *dgdm = [dgdmWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dgdm[i] = 0.0;

	if (mTotal == 0.0) return dgdmWrapper;

	double vlWeight = 0.0;
	for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

	for (NSUInteger k=0; k<na; k++) {
		double g = 0.0;
		double dvlWeightdm = dArray[k]/mTotal - vlWeight/mTotal;
		for (NSUInteger i=0; i<na; i++) {
			double phiI = (m[i]/mTotal)*dArray[i]/vlWeight;
			double dphiIdm = -(m[i]/mTotal/mTotal)*dArray[i]/vlWeight - (m[i]/mTotal)*dArray[i]*dvlWeightdm/vlWeight/vlWeight;
			if (i == k) dphiIdm += dArray[k]/vlWeight/mTotal;

			for (NSUInteger j=i+1; j<na; j++) {
				double phiJ = (m[j]/mTotal)*dArray[j]/vlWeight;
				double dphiJdm = -(m[j]/mTotal/mTotal)*dArray[j]/vlWeight - (m[j]/mTotal)*dArray[j]*dvlWeightdm/vlWeight/vlWeight;
				if (j == k) dphiJdm += dArray[k]/vlWeight/mTotal;

				g += phiI*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j]);
				dgdm[k] += dphiIdm*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
				         + phiI*dphiJdm*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
				         + phiI*phiJ*2.0*dvlWeightdm*wArray[i][j]/(dArray[i]+dArray[j]);
			}
		}

		for (NSUInteger i=0; i<ns; i++) {
			double xOnSite = 0.0;
			for (NSUInteger j=0; j<na; j++) xOnSite += siteSpeciesCoefficients[i][j]*m[j]/mTotal;
			double dxOnSitedm = siteSpeciesCoefficients[i][k]/mTotal - xOnSite/mTotal;
			g += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*xOnSite*log(xOnSite) : 0.0;
			dgdm[k] += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*(dxOnSitedm*log(xOnSite) + dxOnSitedm) : 0.0;
		}

		dgdm[k] = g + mTotal*dgdm[k];
	}

	if (computeMixingQuantities) return dgdmWrapper;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dgdm[i] += [component getGibbsFreeEnergyFromT:t andP:p];
	}

	return dgdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getD2gDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleMatrix *d2gdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:na andWithColumnSize:na andInitialValue:0.0];
	double **d2gdm2 = [d2gdm2Wrapper pointerToPointerToDouble];
    for (NSUInteger i=0; i<na; i++) for (NSUInteger j=0; j<na; j++) d2gdm2[i][j] = 0.0;

	if (mTotal == 0.0) return d2gdm2Wrapper;

	double vlWeight = 0.0;
	for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

	for (NSUInteger k=0; k<na; k++) {
		double dvlWeightdmk = dArray[k]/mTotal - vlWeight/mTotal;

		for (NSUInteger l=k; l<na; l++) {
			double dvlWeightdml = dArray[l]/mTotal - vlWeight/mTotal;
			double d2vlWeightdm2 = -dArray[k]/mTotal/mTotal - dArray[l]/mTotal/mTotal + 2.0*vlWeight/mTotal/mTotal;

			double g = 0.0;
			double dgdmk = 0.0;
			double dgdml = 0.0;
			for (NSUInteger i=0; i<na; i++) {
				double phiI = (m[i]/mTotal)*dArray[i]/vlWeight;
				double dphiIdmk = -(m[i]/mTotal/mTotal)*dArray[i]/vlWeight - (m[i]/mTotal)*dArray[i]*dvlWeightdmk/vlWeight/vlWeight;
				if (i == k) dphiIdmk += dArray[k]/vlWeight/mTotal;

				double dphiIdml = -(m[i]/mTotal/mTotal)*dArray[i]/vlWeight - (m[i]/mTotal)*dArray[i]*dvlWeightdml/vlWeight/vlWeight;
				if (i == l) dphiIdml += dArray[l]/vlWeight/mTotal;

				double d2phiIdm2 = 2.0*(m[i]/mTotal/mTotal/mTotal)*dArray[i]/vlWeight
				                 + (m[i]/mTotal/mTotal)*dArray[i]*dvlWeightdml/vlWeight/vlWeight
				                 + (m[i]/mTotal/mTotal)*dArray[i]*dvlWeightdmk/vlWeight/vlWeight
				                 - (m[i]/mTotal)*dArray[i]*d2vlWeightdm2/vlWeight/vlWeight
				                 + 2.0*(m[i]/mTotal)*dArray[i]*dvlWeightdmk*dvlWeightdml/vlWeight/vlWeight/vlWeight;
				if (i == k) d2phiIdm2 += - dArray[k]*dvlWeightdml/vlWeight/vlWeight/mTotal - dArray[k]/vlWeight/mTotal/mTotal;
				if (i == l) d2phiIdm2 += - dArray[l]*dvlWeightdmk/vlWeight/vlWeight/mTotal - dArray[l]/vlWeight/mTotal/mTotal;

				for (NSUInteger j=i+1; j<na; j++) {
					double phiJ = (m[j]/mTotal)*dArray[j]/vlWeight;
					double dphiJdmk = -(m[j]/mTotal/mTotal)*dArray[j]/vlWeight - (m[j]/mTotal)*dArray[j]*dvlWeightdmk/vlWeight/vlWeight;
					if (j == k) dphiJdmk += dArray[k]/vlWeight/mTotal;

					double dphiJdml = -(m[j]/mTotal/mTotal)*dArray[j]/vlWeight - (m[j]/mTotal)*dArray[j]*dvlWeightdml/vlWeight/vlWeight;
					if (j == l) dphiJdml += dArray[l]/vlWeight/mTotal;


					double d2phiJdm2 = 2.0*(m[j]/mTotal/mTotal/mTotal)*dArray[j]/vlWeight
					                 + (m[j]/mTotal/mTotal)*dArray[j]*dvlWeightdml/vlWeight/vlWeight
					                 + (m[j]/mTotal/mTotal)*dArray[j]*dvlWeightdmk/vlWeight/vlWeight
					                 - (m[j]/mTotal)*dArray[j]*d2vlWeightdm2/vlWeight/vlWeight
									 + 2.0*(m[j]/mTotal)*dArray[j]*dvlWeightdmk*dvlWeightdml/vlWeight/vlWeight/vlWeight;
					if (j == k) d2phiJdm2 += - dArray[k]*dvlWeightdml/vlWeight/vlWeight/mTotal - dArray[k]/vlWeight/mTotal/mTotal;
					if (j == l) d2phiJdm2 += - dArray[l]*dvlWeightdmk/vlWeight/vlWeight/mTotal - dArray[l]/vlWeight/mTotal/mTotal;

					g += phiI*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j]);
					dgdmk += dphiIdmk*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
						   + phiI*dphiJdmk*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
						   + phiI*phiJ*2.0*dvlWeightdmk*wArray[i][j]/(dArray[i]+dArray[j]);
					dgdml += dphiIdml*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
						   + phiI*dphiJdml*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
						   + phiI*phiJ*2.0*dvlWeightdml*wArray[i][j]/(dArray[i]+dArray[j]);

					d2gdm2[k][l] += d2phiIdm2*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
					              + dphiIdmk*dphiJdml*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
					              + dphiIdmk*phiJ*2.0*dvlWeightdml*wArray[i][j]/(dArray[i]+dArray[j])

					              + dphiIdml*dphiJdmk*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
					              + phiI*d2phiJdm2*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j])
					              + phiI*dphiJdmk*2.0*dvlWeightdml*wArray[i][j]/(dArray[i]+dArray[j])

					              + dphiIdml*phiJ*2.0*dvlWeightdmk*wArray[i][j]/(dArray[i]+dArray[j])
					              + phiI*dphiJdml*2.0*dvlWeightdmk*wArray[i][j]/(dArray[i]+dArray[j])
					              + phiI*phiJ*2.0*d2vlWeightdm2*wArray[i][j]/(dArray[i]+dArray[j]);
				}
			}

			for (NSUInteger i=0; i<ns; i++) {
				double xOnSite = 0.0;
				for (NSUInteger j=0; j<na; j++) xOnSite += siteSpeciesCoefficients[i][j]*m[j]/mTotal;
				double dxOnSitedmk = siteSpeciesCoefficients[i][k]/mTotal - xOnSite/mTotal;
				double dxOnSitedml = siteSpeciesCoefficients[i][l]/mTotal - xOnSite/mTotal;
				double d2xOnSitedm2 = - siteSpeciesCoefficients[i][k]/mTotal/mTotal - siteSpeciesCoefficients[i][l]/mTotal/mTotal + 2.0*xOnSite/mTotal/mTotal;

				g += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*xOnSite*log(xOnSite) : 0.0;
				dgdmk += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*(dxOnSitedmk*log(xOnSite) + dxOnSitedmk) : 0.0;
				dgdml += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*(dxOnSitedml*log(xOnSite) + dxOnSitedml) : 0.0;

				d2gdm2[k][l] += (xOnSite != 0.0) ? R*t*siteMultiplicity[i]*(d2xOnSitedm2*log(xOnSite) + dxOnSitedmk*dxOnSitedml/xOnSite + d2xOnSitedm2) : 0.0;
			}

			d2gdm2[k][l] = dgdml + dgdmk + mTotal*d2gdm2[k][l];
			d2gdm2[l][k] = d2gdm2[k][l];
		}
	}

	return d2gdm2Wrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleTensor *)getD3gDm3FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleTensor *d3gdm3Wrapper = [[DoubleTensor alloc] initWithFirstSize:na andWithSecondSize:na andWithThirdSize:na andInitialValue:0.0];
	[NSException raise:@"StixrudeSolutionPhase ERROR" format:@"Method getD3gDm3FromMolesOfComponents NYI"];
	return d3gdm3Wrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getEnthalpyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];;
	if (mTotal == 0.0) return 0.0;

	double vlWeight = 0.0;
	for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

	double h = 0.0;
	for (NSUInteger i=0; i<na; i++) {
		double phiI = (m[i]/mTotal)*dArray[i]/vlWeight;
		for (NSUInteger j=i+1; j<na; j++) {
			double phiJ = (m[j]/mTotal)*dArray[j]/vlWeight;
			h += phiI*phiJ*2.0*vlWeight*wArray[i][j]/(dArray[i]+dArray[j]);
		}
	}

	h *= mTotal;

	if (computeMixingQuantities) return h;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		h += m[i]*[component getEnthalpyFromT:t andP:p];
	}
	return h;
}

// --> SolutionPhaseProtocol public function
-(double)getEntropyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];;
	if (mTotal == 0.0) return 0.0;

	double s = 0.0;
	for (NSUInteger i=0; i<ns; i++) {
		double xOnSite = 0.0;
		for (NSUInteger j=0; j<na; j++) xOnSite += siteSpeciesCoefficients[i][j]*m[j]/mTotal;
		s += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*xOnSite*log(xOnSite) : 0.0;
	}
	s *= mTotal;

	if (computeMixingQuantities) return s;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		s += m[i]*[component getEntropyFromT:t andP:p];
	}
	return s;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDsDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleVector *dsdmWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
	double *dsdm = [dsdmWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dsdm[i] = 0.0;

	if (mTotal == 0.0) return dsdmWrapper;

	double s = 0.0;
	for (NSUInteger k=0; k<na; k++) {

		for (NSUInteger i=0; i<ns; i++) {
			double xOnSite = 0.0;
			for (NSUInteger j=0; j<na; j++) xOnSite += siteSpeciesCoefficients[i][j]*m[j]/mTotal;
			double dxOnSitedm = siteSpeciesCoefficients[i][k]/mTotal - xOnSite/mTotal;
			s += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*xOnSite*log(xOnSite) : 0.0;
			dsdm[k] += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*(dxOnSitedm*log(xOnSite) + dxOnSitedm) : 0.0;
		}

		dsdm[k] = s + mTotal*dsdm[k];
	}

	if (computeMixingQuantities) return dsdmWrapper;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dsdm[i] += [component getEntropyFromT:t andP:p];
	}
	return dsdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getD2sDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleMatrix *d2sdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:na andWithColumnSize:na andInitialValue:0.0];
	double **d2sdm2 = [d2sdm2Wrapper pointerToPointerToDouble];
    for (NSUInteger i=0; i<na; i++) for (NSUInteger j=0; j<na; j++) d2sdm2[i][j] = 0.0;

	if (mTotal == 0.0) return d2sdm2Wrapper;

	double vlWeight = 0.0;
	for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

	double s = 0.0;
	for (NSUInteger k=0; k<na; k++) {
		for (NSUInteger l=k; l<na; l++) {
			double dsdmk = 0.0;
			double dsdml = 0.0;

			for (NSUInteger i=0; i<ns; i++) {
				double xOnSite = 0.0;
				for (NSUInteger j=0; j<na; j++) xOnSite += siteSpeciesCoefficients[i][j]*m[j]/mTotal;
				double dxOnSitedmk = siteSpeciesCoefficients[i][k]/mTotal - xOnSite/mTotal;
				double dxOnSitedml = siteSpeciesCoefficients[i][l]/mTotal - xOnSite/mTotal;
				double d2xOnSitedm2 = - siteSpeciesCoefficients[i][k]/mTotal/mTotal - siteSpeciesCoefficients[i][l]/mTotal/mTotal + 2.0*xOnSite/mTotal/mTotal;

				s += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*xOnSite*log(xOnSite) : 0.0;
				dsdmk += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*(dxOnSitedmk*log(xOnSite) + dxOnSitedmk) : 0.0;
				dsdml += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*(dxOnSitedml*log(xOnSite) + dxOnSitedml) : 0.0;

				d2sdm2[k][l] += (xOnSite != 0.0) ? -R*siteMultiplicity[i]*(d2xOnSitedm2*log(xOnSite) + dxOnSitedmk*dxOnSitedml/xOnSite + d2xOnSitedm2) : 0.0;
			}

			d2sdm2[k][l] = dsdml + dsdmk + mTotal*d2sdm2[k][l];
			d2sdm2[l][k] = d2sdm2[k][l];
		}
	}

	return d2sdm2Wrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getHeatCapacityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double cp = 0.0;
	if (computeMixingQuantities) return cp;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		cp += m[i]*[component getHeatCapacityFromT:t andP:p];
	}
	return cp;
}

// --> SolutionPhaseProtocol public function
-(double)getDcpDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double dcpdt = 0.0;
	if (computeMixingQuantities) return dcpdt;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dcpdt += m[i]*[component getDcpDtFromT:t andP:p];
	}
	return dcpdt;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDCpDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dcpdmWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
	double *dcpdm = [dcpdmWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dcpdm[i] = 0.0;

	if (computeMixingQuantities) return dcpdmWrapper;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dcpdm[i] += [component getHeatCapacityFromT:t andP:p];
	}
	return dcpdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getVolumeFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double v = 0.0;
	if (computeMixingQuantities) return v;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		v += m[i]*[component getVolumeFromT:t andP:p];
	}
	return v;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDvDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dvdmWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
	double *dvdm = [dvdmWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dvdm[i] = 0.0;

	if (computeMixingQuantities) return dvdmWrapper;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdm[i] += [component getVolumeFromT:t andP:p];
	}
	return dvdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getD2vDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *d2vdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:na andWithColumnSize:na andInitialValue:0.0];
	return d2vdm2Wrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getDvDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double dvdt = 0.0;
	if (computeMixingQuantities) return dvdt;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdt += m[i]*[component getDvDtFromT:t andP:p];
	}
	return dvdt;
}

// --> SolutionPhaseProtocol public function
-(double)getDvDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double dvdp = 0.0;
	if (computeMixingQuantities) return dvdp;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdp += m[i]*[component getDvDpFromT:t andP:p];
	}
	return dvdp;
}

// --> SolutionPhaseProtocol public function
-(double)getD2vDt2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double d2vdt2 = 0.0;
	if (computeMixingQuantities) return d2vdt2;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		d2vdt2 += m[i]*[component getD2vDt2FromT:t andP:p];
	}
	return d2vdt2;
}

// --> SolutionPhaseProtocol public function
-(double)getD2vDtDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double d2vdtdp = 0.0;
	if (computeMixingQuantities) return d2vdtdp;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		d2vdtdp += m[i]*[component getD2vDtDpFromT:t andP:p];
	}
	return d2vdtdp;
}

// --> SolutionPhaseProtocol public function
-(double)getD2vDp2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double d2vdp2 = 0.0;
	if (computeMixingQuantities) return d2vdp2;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		d2vdp2 += m[i]*[component getD2vDp2FromT:t andP:p];
	}
	return d2vdp2;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getD2vDmDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dvdmdtWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
	double *dvdmdt = [dvdmdtWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dvdmdt[i] = 0.0;

	if (computeMixingQuantities) return dvdmdtWrapper;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdmdt[i] += [component getDvDtFromT:t andP:p];
	}
	return dvdmdtWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getD2vDmDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dvdmdpWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
	double *dvdmdp = [dvdmdpWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dvdmdp[i] = 0.0;

	if (computeMixingQuantities) return dvdmdpWrapper;

	for (NSUInteger i=0; i<na; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdmdp[i] += [component getDvDpFromT:t andP:p];
	}
	return dvdmdpWrapper;
}

// --> SolutionPhaseProtocol public function
-(NSString *)getFormulaFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	return convertMolesOfSpeciesToDisplayFormula(m);
}

@end
