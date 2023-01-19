//
//  StixrudeSolutions.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/17/11.
//  Garnet changed to a three-site formulation by Bob Myhill on 07/16/18.
//  Opx, Cpx and Garnet expanded to reciprocal solutions by Mark Ghiorso on 06/12/20
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "StixrudeSolutions.h"
#import "StixrudeStoichiometricPhases.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

typedef void (^Convertor)(double *, double *);
typedef	NSString *(^Display)(double *);
typedef void (^Adjuster)(double *, double *);

#define NA 6
#define ND 6
#define NS 16

typedef void (^TypeSiteFunction)(double m[NA], double sites[NS]);
typedef void (^TypeDsiteFunction)(double m[NA], double Dsites[NS][NA]);
typedef void (^TypeD2siteFunction)(double m[NA], double D2sites[NS][NA][NA]);

@implementation StixrudeReciprocalSolutionPhase

-(id)initWithPhaseName:(NSString *)name
           withSpecies:(NSArray *)speciesProperties
    withSpeciesWeights:(NSArray *)speciesWeights
    withSpeciesXfactor:(NSArray *)speciesXfactor
            withWarray:(NSArray *)wIn
             withSites:(NSArray *)oldsites
            WithNAtoms:(NSUInteger)nAtomIn
  withElementConvertor:(void (^)(double *, double *))convertor
 andWithFormulaDisplay:(NSString *(^)(double *))display
andWithSaturationStateAdjuster:(void (^)(double *, double *))adjuster
andWithDependentSpecies:(NSArray *)dependentSpeciesProperties
  andWithSiteFunctions:(void (^)(double m[NA], double site[NS]))siteFunction
andWithSiteFirstDerivative:(void (^)(double m[NA], double Dsite[NS][NA]))DsiteFunction
andWithSiteSecondDerivative:(void (^)(double m[NA], double D2site[NS][NA][NA]))D2siteFunction {
    if ((self = [super initWithPhaseName:name
                             withSpecies:speciesProperties
                      withSpeciesWeights:speciesWeights
                      withSpeciesXfactor:speciesXfactor
                              withWarray:wIn withSites:oldsites
                              WithNAtoms:nAtomIn
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
        dependentEndmembers = [NSArray arrayWithArray:dependentSpeciesProperties];
        nd = [dependentEndmembers count];
        siteMoleFraction = siteFunction;
        DsiteMoleFraction = DsiteFunction;
        D2siteMoleFraction = D2siteFunction;
    }
    return self;
}

// --> Override the SolutionPhaseProtocol function
-(NSUInteger)numberOfSolutionSpecies {
    return na + nd;
}

// --> Override the SolutionPhaseProtocol function
-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
    if (index < na) {
        return [[endmembers objectAtIndex:index] phaseName];
    } else {
        return [[dependentEndmembers objectAtIndex:index-na] objectAtIndex:0];
    }
}

// --> Override the SolutionPhaseProtocol function
-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
    DoubleVector *mComponentsWrapper = [[DoubleVector alloc] initWithSize:na];
    double *mComponents = [mComponentsWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) mComponents[i] = mSpecies[i];
    // Add dependent species mappings to components
    for (NSUInteger i=0; i<nd; i++) {
        NSArray *coeffArray = [dependentEndmembers objectAtIndex:i];
        for (NSUInteger j=0; j<na; j++) {
            double coeff = [[coeffArray objectAtIndex:j+1] doubleValue];
            if (coeff != 0.0) mComponents[j] += coeff*mSpecies[na+i];
        }
    }
    return mComponentsWrapper;
}

// --> Override the SolutionPhaseProtocol function
-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *compMuWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
    double *compMu = [compMuWrapper pointerToDouble];
    DoubleVector *speciesMuWrapper = [[DoubleVector alloc] initWithSize:na+nd];
    double *speciesMu = [speciesMuWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) speciesMu[i] = compMu[i];
    // Add dependent species mappings to components
    for (NSUInteger i=0; i<nd; i++) {
        NSArray *coeffArray = [dependentEndmembers objectAtIndex:i];
        speciesMu[na+i] = 0.0;
        for (NSUInteger j=0; j<na; j++) {
            double coeff = [[coeffArray objectAtIndex:j+1] doubleValue];
            if (coeff != 0.0) speciesMu[na+i] += coeff*compMu[j];
        }
    }
    return  speciesMuWrapper;
}

// --> Override the SolutionPhaseProtocol function
-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
    if (index < na) return [[endmembers objectAtIndex:index] formulaAsElementArray];
    DoubleVector *elementWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
    double *element = [elementWrapper pointerToDouble];
    // Add dependent species mappings to components
    NSArray *coeffArray = [dependentEndmembers objectAtIndex:index-na];
    for (NSUInteger i=0; i<na; i++) {
        double coeff = [[coeffArray objectAtIndex:i+1] doubleValue];
        if (coeff != 0.0) {
            DoubleVector *compElementWrapper = [[endmembers objectAtIndex:i] formulaAsElementArray];
            double *compElement = [compElementWrapper pointerToDouble];
            for (NSUInteger j=0; j<107; j++) element[j] += coeff*compElement[j];
        }
    }
    return elementWrapper;
}

#define R 8.3143

// --> Override the SolutionPhaseProtocol function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
    NSMutableArray *results = [NSMutableArray arrayWithCapacity:na+1];
    double mu0[NA+ND], deltaMu[NA+ND], xNz[NA+ND], x[NA+ND], gamma[NA+ND], xLast[NA+ND], affinity = 0.0,
        reducedAdjustmentCoefficient[NA+ND], gammaLast[NA+ND], aSpecies[ND];
    BOOL addme[ND];
    NSUInteger i, j, nz = 0, index[NA+ND];

    BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
    BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
    if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");

    // Compute solid -> liquid delta mus and deflate composition space
    for (i=0; i<na; i++) {
        if (chemicalPotentials[i] != 0.0) {
            mu0[i] = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
            deltaMu[nz] = chemicalPotentials[i] - mu0[i];
            index[nz] = i;
            gamma[nz] = 1.0;
            reducedAdjustmentCoefficient[nz] = adjustmentCoefficient[i];
            nz++;
        }
        x[i] = 0.0;
        xLast[i] = 0.0;
    }

    // Dependent species
    for (i=0; i<nd; i++) {
        NSArray *coeffArray = [dependentEndmembers objectAtIndex:i];
        addme[i] = true;
        for (j=0; j<na; j++) {
            double coeff = [[coeffArray objectAtIndex:1+j] doubleValue];
            if ((coeff != 0.0) && (chemicalPotentials[j] == 0.0)) addme[i] = false;
        }
        if (addme[i]) {
            mu0[na+i] = [[coeffArray objectAtIndex:na+1] doubleValue];
            deltaMu[nz] = -mu0[na+i] ;
            for (NSUInteger j=0; j<na; j++) {
                double coeff = [[coeffArray objectAtIndex:1+j] doubleValue];
                mu0[na+i]   += coeff*mu0[j];
                deltaMu[nz] += coeff*(chemicalPotentials[j] - mu0[j]);
            }
            index[nz] = na+i;
            gamma[nz] = 1.0;
            reducedAdjustmentCoefficient[nz] = 1.0;
            nz++;
        } else mu0[na+i] = 0.0;
        x[na+i] = 0.0;
        xLast[na+i] = 0.0;
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

        // Reduce species -> components
        DoubleVector *xReducedWrapper = [self convertMolesOfSpeciesToMolesOfComponents:x];
        double *xReduced = [xReducedWrapper pointerToDouble];
        if (![self testPermissibleValuesOfComponents:xReduced]) NSLog(@"Composition estimate is infeasible.");
        // Determine activity coefficients
        DoubleVector *aWrapper = [self getActivityFromMolesOfComponents:xReduced andT:t andP:p];
        double *a = [aWrapper pointerToDouble];
        // Compute activity coefficients
        for (i=0, j=0; i<na; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];
        DoubleVector *muWrapper = [self getChemicalPotentialFromMolesOfComponents:xReduced andT:t andP:p];
        double *mu = [muWrapper pointerToDouble];
        for (i=0; i<nd; i++) if (addme[i]) {
            double temp = 0.0;
            for (NSUInteger k=0; k<na; k++) {
                double coeff = [[[dependentEndmembers objectAtIndex:i] objectAtIndex:k+1] doubleValue];
                temp += (coeff != 0.0) ? coeff*mu[k] : 0.0;
            }
            aSpecies[i] = exp((temp-mu0[na+i])/(R*t));
            gamma[j++] = aSpecies[i]/x[na+i];
        }

        if (count > 25) { // pure empiricism
            for (i=0; i<nz; i++) gamma[i] = (gamma[i]+gammaLast[i])/2.0;
        }
        for (i=0; i<nz; i++) gammaLast[i] = gamma[i];
        correctActivityCoefficients(gamma, x);
        if (debugV) {
            for (i=0, j=0; i<na; i++) if (x[i] != 0.0)
                NSLog(@"... a of %@ is %20.13e with X %20.13e and g %20.13e",
                      [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                      a[i], x[i], gamma[j++]);
            for (i=0; i<nd; i++) if (x[na+i] != 0.0)
                NSLog(@"... a of %@ is %20.13e with X %20.13e and g %20.13e",
                      [[[dependentEndmembers objectAtIndex:i] objectAtIndex:0]  stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                      aSpecies[i], x[i+na], gamma[j++]);
            NSLog(@"... a of %@ is %@ with X %@ and g %@ Aff %20.13e delta Aff %20.13e", @"sum                 ",
                  @"                    ", @"                    ", @"                    ", affinity, affinity-affinityLast);
        }
        converged = (fabs(affinity-affinityLast) < minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod);
        count++;

    } while (count < numberOfIterationsAllowedInSaturationMethod && !converged);

    if (debugS) {
        NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
              converged ? @"YES" : @"NO", [self phaseName], count, fabs(affinity-affinityLast), nAtoms);
        for (i=0, j=0; i<na; i++) if (x[i] != 0.0)
            NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
                  [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                  gamma[j++], x[i]);
        for (i=0; i<nd; i++) if (x[na+i] != 0.0)
            NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
                  [[[dependentEndmembers objectAtIndex:i] objectAtIndex:0] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                  gamma[j++], x[na+i]);
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

// --> Override the SolutionPhaseProtocol function
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

    siteMoleFraction(m, sites);
    for (NSUInteger i=0; i<ns; i++) {
        double xOnSite = sites[i];
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

// --> Override the SolutionPhaseProtocol function
-(DoubleVector *)getDgDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleVector *dgdmWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
    double *dgdm = [dgdmWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dgdm[i] = 0.0;

    if (mTotal == 0.0) return dgdmWrapper;

    double vlWeight = 0.0;
    for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

    siteMoleFraction(m, sites);
    DsiteMoleFraction(m, Dsites);
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
            double xOnSite = sites[i];
            double dxOnSitedm = Dsites[i][k];
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

// --> Override the SolutionPhaseProtocol function
-(DoubleMatrix *)getD2gDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleMatrix *d2gdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:na andWithColumnSize:na andInitialValue:0.0];
    double **d2gdm2 = [d2gdm2Wrapper pointerToPointerToDouble];
    for (NSUInteger i=0; i<na; i++) for (NSUInteger j=0; j<na; j++) d2gdm2[i][j] = 0.0;

    if (mTotal == 0.0) return d2gdm2Wrapper;

    double vlWeight = 0.0;
    for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

    siteMoleFraction(m, sites);
    DsiteMoleFraction(m, Dsites);
    D2siteMoleFraction(m, D2sites);
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
                double xOnSite = sites[i];
                double dxOnSitedmk = Dsites[i][k];
                double dxOnSitedml = Dsites[i][l];
                double d2xOnSitedm2 = D2sites[i][k][l];

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

// --> Override the SolutionPhaseProtocol function
-(double)getEntropyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double mTotal = [self totalMolesFromMolesOfComponents:m];;
    if (mTotal == 0.0) return 0.0;

    double s = 0.0;
    siteMoleFraction(m, sites);
    for (NSUInteger i=0; i<ns; i++) {
        double xOnSite = sites[i];
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

// --> Override the SolutionPhaseProtocol function
-(DoubleVector *)getDsDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleVector *dsdmWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
    double *dsdm = [dsdmWrapper pointerToDouble];
    for (NSUInteger i=0; i<na; i++) dsdm[i] = 0.0;

    if (mTotal == 0.0) return dsdmWrapper;

    double s = 0.0;
    siteMoleFraction(m, sites);
    DsiteMoleFraction(m, Dsites);
    for (NSUInteger k=0; k<na; k++) {

        for (NSUInteger i=0; i<ns; i++) {
            double xOnSite = sites[i];
            double dxOnSitedm = Dsites[i][k];
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

// --> Override the SolutionPhaseProtocol function
-(DoubleMatrix *)getD2sDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double mTotal = [self totalMolesFromMolesOfComponents:m];
    DoubleMatrix *d2sdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:na andWithColumnSize:na andInitialValue:0.0];
    double **d2sdm2 = [d2sdm2Wrapper pointerToPointerToDouble];
    for (NSUInteger i=0; i<na; i++) for (NSUInteger j=0; j<na; j++) d2sdm2[i][j] = 0.0;

    if (mTotal == 0.0) return d2sdm2Wrapper;

    double vlWeight = 0.0;
    for (NSUInteger i=0; i<na; i++) vlWeight += dArray[i]*m[i]/mTotal;

    double s = 0.0;
    siteMoleFraction(m, sites);
    DsiteMoleFraction(m, Dsites);
    D2siteMoleFraction(m, D2sites);
    for (NSUInteger k=0; k<na; k++) {
        for (NSUInteger l=k; l<na; l++) {
            double dsdmk = 0.0;
            double dsdml = 0.0;

            for (NSUInteger i=0; i<ns; i++) {
                double xOnSite = sites[i];
                double dxOnSitedmk = Dsites[i][k];
                double dxOnSitedml = Dsites[i][l];
                double d2xOnSitedm2 = D2sites[i][k][l];

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

@end

// -------------------------------------------------------------------------------------------------------------

@implementation FeldsparStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[20];
		m[1] = e[11];
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mCa = m[0]/mTotal;
		double mNa = m[1]/mTotal;
		double mAl = mCa*2.0 + mNa;
		double mSi = mCa*2.0 + mNa*3.0;
		return (NSString *) [NSString stringWithFormat:@"Na%4.2fCa%4.2fAl%4.2fSi%4.2fO8", mNa, mCa, mAl, mSi];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Feldspar"
                             withSpecies:[NSArray arrayWithObjects:[[AnorthiteStixrude alloc] init], [[AlbiteStixrude alloc] init], nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:26000.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Ca on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil], nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Na on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil], nil],
                                          nil]
                              WithNAtoms:13
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
				])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation OlivineStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12]/2.0;
		m[1] = e[26]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = 2.0*m[0]/mTotal;
		double mFe = 2.0*m[1]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fSiO4", mMg, mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Olivine"
                             withSpecies:[NSArray arrayWithObjects:[[ForsteriteStixrude alloc] init], [[FayaliteStixrude alloc] init], nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:7600.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Mg on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil], nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Fe on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil], nil],
                                          nil]
                              WithNAtoms:7
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation WadsleyiteStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12]/2.0;
		m[1] = e[26]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = 2.0*m[0]/mTotal;
		double mFe = 2.0*m[1]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fSiO4", mMg, mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Wadsleyite"
                             withSpecies:[NSArray arrayWithObjects:[[MgWadsleyiteStixrude alloc] init], [[FeWadsleyiteStixrude alloc] init], nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:16500.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Mg on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil], nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Fe on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil], nil],
                                          nil]
                              WithNAtoms:7
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
				])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation RingwooditeStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12]/2.0;
		m[1] = e[26]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = 2.0*m[0]/mTotal;
		double mFe = 2.0*m[1]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fSiO4", mMg, mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Ringwoodite"
                             withSpecies:[NSArray arrayWithObjects:[[MgRingwooditeStixrude alloc] init], [[FeRingwooditeStixrude alloc] init], nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:9100.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Mg on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil], nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Fe on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil], nil],
                                          nil]
                              WithNAtoms:7
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation PerovskiteStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12];
		m[1] = e[26];
		m[2] = e[13]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg =     m[0]/mTotal;
		double mFe =     m[1]/mTotal;
		double mAl = 2.0*m[2]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fAl%4.2fSi%4.2fO3", mMg, mFe, mAl, mMg+mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<3; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Perovskite"
                             withSpecies:[NSArray arrayWithObjects:[[MgPerovskiteStixrude alloc] init], // MgSiO3
                                          [[FePerovskiteStixrude alloc] init],                          // FeSiO3
                                          [[AlPerovskiteStixrude alloc] init],                          // AlAlO3
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:0.39],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // MgSiO3 - FeSiO3
                                          [NSNumber numberWithDouble:116000.0],                     // MgSiO3 - AlAlO3
                                          [NSNumber numberWithDouble:0.0],                          // FeSiO3 - AlAlO3
                                          nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Si on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          nil]
                              WithNAtoms:5
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation PostPerovskiteStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12];
		m[1] = e[26];
		m[2] = e[13]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg =     m[0]/mTotal;
		double mFe =     m[1]/mTotal;
		double mAl = 2.0*m[2]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fAl%4.2fSi%4.2fO3", mMg, mFe, mAl, mMg+mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<3; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Post-Perovskite"
                             withSpecies:[NSArray arrayWithObjects:[[MgPostPerovskiteStixrude alloc] init], // MgSiO3
                                          [[FePostPerovskiteStixrude alloc] init],                          // FeSiO3
                                          [[AlPostPerovskiteStixrude alloc] init],                          // AlAlO3
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // MgSiO3 - FeSiO3
                                          [NSNumber numberWithDouble:60000.0],                      // MgSiO3 - AlAlO3
                                          [NSNumber numberWithDouble:0.0],                          // FeSiO3 - AlAlO3
                                          nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Si on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          nil]
                              WithNAtoms:5
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
				])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation OrthopyroxeneStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = (e[12]-e[20]-e[13]/2.0)/2.0;
		m[1] = e[26]/2.0;
		m[2] = e[13]/2.0;
		m[3] = e[20];
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2] + m[3];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = 2.0*m[0]/mTotal + m[2]/mTotal + m[3]/mTotal;
		double mFe = 2.0*m[1]/mTotal;
		double mAl = 2.0*m[2]/mTotal;
		double mCa =     m[3]/mTotal;
		double mSi = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + m[2]/mTotal + 2.0*m[3]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Ca%4.2fMg%4.2fFe%4.2fAl%4.2fSi%4.2fO6", mCa, mMg, mFe, mAl, mSi];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<7; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

    TypeSiteFunction siteFunction = ^(double m[NA], double x[NS]) {
        double mTotal = m[0] + m[1] + m[2] + m[3]; // order En, Fs, MaTs, Di
        double xEn = m[0]/mTotal;
        double xFs = m[1]/mTotal;
        double xMaTs = m[2]/mTotal;
        double xDi = m[3]/mTotal;
        double denom = 2.0*xEn + 2.0*xFs + xMaTs + xDi;
        // sites XCaM2, XMgM2, XFeM2, XMgM1, XFeM1, XAlM1
        x[0] = xDi;
        x[1] = 1.0 - xDi - 2.0*xFs*(xEn + xFs + xMaTs)/denom;
        x[2] = 2.0*xFs*(xEn + xFs + xMaTs)/denom;
        x[3] = 1.0 - xMaTs - 2.0*xFs*(xEn + xFs + xDi)/denom;
        x[4] = 2.0*xFs*(xEn + xFs + xDi)/denom;
        x[5] = xMaTs;
    };

    TypeDsiteFunction DsiteFunction = ^(double m[NA], double dx[NS][NA]) {
        double mTotal = m[0] + m[1] + m[2] + m[3];
        if (mTotal == 0.0) mTotal = 1.0;
        double xEn   = m[0]/mTotal;
        double xFs   = m[1]/mTotal;
        double xMaTs = m[2]/mTotal;
        double xDi   = m[3]/mTotal;
        double dxEn[4], dxFs[4], dxMaTs[4], dxDi[4];
        for (NSUInteger i=0; i<4;i++) {
            dxEn[i]   = (i == 0) ? (1.0-m[0])/mTotal : -m[0]/mTotal;
            dxFs[i]   = (i == 1) ? (1.0-m[1])/mTotal : -m[1]/mTotal;
            dxMaTs[i] = (i == 2) ? (1.0-m[2])/mTotal : -m[2]/mTotal;
            dxDi[i]   = (i == 3) ? (1.0-m[3])/mTotal : -m[3]/mTotal;
        }
        double denom = 2.0*xEn + 2.0*xFs + xMaTs + xDi;
        double dinvDenom[4];
        for (NSUInteger i=0; i<4; i++) {
            dinvDenom[i] = -(2.0*dxEn[i] + 2.0*dxFs[i] + dxMaTs[i] + dxDi[i])/denom/denom;
            // x[0] = xDi;
            dx[0][i] = dxDi[i];
            // 1.0 - xDi - 2.0*xFs*(xEn + xFs + xMaTs)*invDenom
            dx[1][i] = - dxDi[i] - 2.0*dxFs[i]*(xEn + xFs + xMaTs)/denom - 2.0*xFs*(dxEn[i] + dxFs[i] + dxMaTs[i])/denom
                     - 2.0*xFs*(xEn + xFs + xMaTs)*dinvDenom[i];
            // 2.0*xFs*(xEn + xFs + xMaTs)/denom;
            dx[2][i] = 2.0*dxFs[i]*(xEn + xFs + xMaTs)/denom + 2.0*xFs*(dxEn[i] + dxFs[i] + dxMaTs[i])/denom
                     + 2.0*xFs*(xEn + xFs + xMaTs)*dinvDenom[i];
            // 1.0 - xMaTs - 2.0*xFs*(xEn + xFs + xDi)/denom
            dx[3][i] = - dxMaTs[i] - 2.0*dxFs[i]*(xEn + xFs + xDi)/denom - 2.0*xFs*(dxEn[i] + dxFs[i] + dxDi[i])/denom
                     - 2.0*xFs*(xEn + xFs + xDi)*dinvDenom[i];
            // 2.0*xFs*(xEn + xFs + xDi)/denom
            dx[4][i] = 2.0*dxFs[i]*(xEn + xFs + xDi)/denom + 2.0*xFs*(dxEn[i] + dxFs[i] + dxDi[i])/denom
                     + 2.0*xFs*(xEn + xFs + xDi)*dinvDenom[i];
            // xMaTs
            dx[5][i] = dxMaTs[i];
        }
    };

    TypeD2siteFunction D2siteFunction = ^(double m[NA], double d2x[NS][NA][NA]) {
        double mTotal = m[0] + m[1] + m[2] + m[3];
        if (mTotal == 0.0) mTotal = 1.0;
        double xEn   = m[0]/mTotal;
        double xFs   = m[1]/mTotal;
        double xMaTs = m[2]/mTotal;
        double xDi   = m[3]/mTotal;
        double dxEn[4], dxFs[4], dxMaTs[4], dxDi[4], d2xEn[4][4], d2xFs[4][4], d2xMaTs[4][4], d2xDi[4][4];
        for (NSUInteger i=0; i<4; i++) {
            dxEn[i]   = (i == 0) ? (1.0-m[0])/mTotal : -m[0]/mTotal;
            dxFs[i]   = (i == 1) ? (1.0-m[1])/mTotal : -m[1]/mTotal;
            dxMaTs[i] = (i == 2) ? (1.0-m[2])/mTotal : -m[2]/mTotal;
            dxDi[i]   = (i == 3) ? (1.0-m[3])/mTotal : -m[3]/mTotal;
        }
        for (NSUInteger j=0; j<4; j++) {
            d2xEn[0][j]   = (j == 0) ? -1.0/mTotal - (1.0-m[0])/mTotal/mTotal : -(1.0-m[0])/mTotal/mTotal;
            d2xFs[0][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xMaTs[0][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xDi[0][j]   = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xEn[1][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xFs[1][j]   = (j == 1) ? -1.0/mTotal - (1.0-m[1])/mTotal/mTotal : -(1.0-m[1])/mTotal/mTotal;
            d2xMaTs[1][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xDi[1][j]   = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xEn[2][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xFs[2][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xMaTs[2][j] = (j == 2) ? -1.0/mTotal - (1.0-m[2])/mTotal/mTotal : -(1.0-m[2])/mTotal/mTotal;
            d2xDi[2][j]   = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xEn[3][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xFs[3][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xMaTs[3][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xDi[3][j]   = (j == 3) ? -1.0/mTotal - (1.0-m[3])/mTotal/mTotal : -(1.0-m[3])/mTotal/mTotal;
        }
        double denom = 2.0*xEn + 2.0*xFs + xMaTs + xDi;
        double dinvDenom[4], d2invDenom[4][4];
        for (NSUInteger i=0; i<4; i++) {
            dinvDenom[i] = -(2.0*dxEn[i] + 2.0*dxFs[i] + dxMaTs[i] + dxDi[i])/denom/denom;
            for (NSUInteger j=0; j<4; j++) {
                d2invDenom[i][j] = -(2.0*d2xEn[i][j] + 2.0*d2xFs[i][j] + d2xMaTs[i][j] + d2xDi[i][j])/denom/denom
                + 2.0*(2.0*dxEn[i] + 2.0*dxFs[i] + dxMaTs[i] + dxDi[i])*(2.0*dxEn[j] + 2.0*dxFs[j] + dxMaTs[j] + dxDi[j])/denom/denom/denom;
                // x[0] = xDi;
                d2x[0][i][j] = d2xDi[i][j];
                // 1.0 - xDi - 2.0*xFs*(xEn + xFs + xMaTs)*invDenom
                d2x[1][i][j] = - d2xDi[i][j]
                    - 2.0*d2xFs[i][j]*(xEn + xFs + xMaTs)/denom - 2.0*dxFs[i]*(dxEn[j] + dxFs[j] + dxMaTs[j])/denom
                    - 2.0*dxFs[i]*(xEn + xFs + xMaTs)*dinvDenom[j]
                    - 2.0*dxFs[j]*(dxEn[i] + dxFs[i] + dxMaTs[i])/denom - 2.0*xFs*(d2xEn[i][j] + d2xFs[i][j] + d2xMaTs[i][j])/denom
                    - 2.0*xFs*(dxEn[i] + dxFs[i] + dxMaTs[i])*dinvDenom[j]
                    - 2.0*dxFs[j]*(xEn + xFs + xMaTs)*dinvDenom[i] - 2.0*xFs*(dxEn[j] + dxFs[j] + dxMaTs[j])*dinvDenom[i]
                    - 2.0*xFs*(xEn + xFs + xMaTs)*d2invDenom[i][j];
                // 2.0*xFs*(xEn + xFs + xMaTs)/denom;
                d2x[2][i][j] = 2.0*d2xFs[i][j]*(xEn + xFs + xMaTs)/denom + 2.0*dxFs[i]*(dxEn[j] + dxFs[j] + dxMaTs[j])/denom
                    + 2.0*dxFs[i]*(xEn + xFs + xMaTs)*dinvDenom[j]
                    + 2.0*dxFs[j]*(dxEn[i] + dxFs[i] + dxMaTs[i])/denom + 2.0*xFs*(d2xEn[i][j] + d2xFs[i][j] + d2xMaTs[i][j])/denom
                    + 2.0*xFs*(dxEn[i] + dxFs[i] + dxMaTs[i])*dinvDenom[j]
                    + 2.0*dxFs[j]*(xEn + xFs + xMaTs)*dinvDenom[i] + 2.0*xFs*(dxEn[j] + dxFs[j] + dxMaTs[j])*dinvDenom[i]
                    + 2.0*xFs*(xEn + xFs + xMaTs)*d2invDenom[i][j];
                // 1.0 - xMaTs - 2.0*xFs*(xEn + xFs + xDi)/denom
                d2x[3][i][j] = - d2xMaTs[i][j]
                    - 2.0*d2xFs[i][j]*(xEn + xFs + xDi)/denom - 2.0*dxFs[i]*(dxEn[j] + dxFs[j] + dxDi[j])/denom
                    - 2.0*dxFs[i]*(xEn + xFs + xDi)*dinvDenom[j]
                    - 2.0*dxFs[j]*(dxEn[i] + dxFs[i] + dxDi[i])/denom - 2.0*xFs*(d2xEn[i][j] + d2xFs[i][j] + d2xDi[i][j])/denom
                    - 2.0*xFs*(dxEn[i] + dxFs[i] + dxDi[i])*dinvDenom[j]
                    - 2.0*dxFs[j]*(xEn + xFs + xDi)*dinvDenom[i] - 2.0*xFs*(dxEn[j] + dxFs[j] + dxDi[j])*dinvDenom[i]
                    - 2.0*xFs*(xEn + xFs + xDi)*d2invDenom[i][j];
                // 2.0*xFs*(xEn + xFs + xDi)/denom
                d2x[4][i][j] = 2.0*d2xFs[i][j]*(xEn + xFs + xDi)/denom + 2.0*dxFs[i]*(dxEn[j] + dxFs[j] + dxDi[j])/denom
                    + 2.0*dxFs[i]*(xEn + xFs + xDi)*dinvDenom[j]
                    + 2.0*dxFs[j]*(dxEn[i] + dxFs[i] + dxDi[i])/denom + 2.0*xFs*(d2xEn[i][j] + d2xFs[i][j] + d2xDi[i][j])/denom
                    + 2.0*xFs*(dxEn[i] + dxFs[i] + dxDi[i])*dinvDenom[j]
                    + 2.0*dxFs[j]*(xEn + xFs + xDi)*dinvDenom[i] + 2.0*xFs*(dxEn[j] + dxFs[j] + dxDi[j])*dinvDenom[i]
                    + 2.0*xFs*(xEn + xFs + xDi)*d2invDenom[i][j];
                // xMaTs
                d2x[5][i][j] = d2xMaTs[i][j];
            }
        }
    };

    if ((self = [super initWithPhaseName:@"Orthopyroxene"
                             withSpecies:[NSArray arrayWithObjects:[[EnstatiteStixrude alloc] init], // Mg2Si2O6
                                          [[FerrosiliteStixrude alloc] init],                        // Fe2Si2O6
                                          [[MgTschermaksStixrude alloc] init],                       // MgAlSiAlO6
                                          [[OrthoDiopsideStixrude alloc] init],                      // CaMgSi2O6
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // Fe2Si2O6   - Mg2Si2O6
                                          [NSNumber numberWithDouble:0.0],                          // MgAlSiAlO6 - Mg2Si2O6
                                          [NSNumber numberWithDouble:0.0],                          // MgAlSiAlO6 - Fe2Si2O6
                                          [NSNumber numberWithDouble:32100.0],                      // CaMgSi2O6  - Mg2Si2O6
                                          [NSNumber numberWithDouble:0.0],                          // CaMgSi2O6  - Fe2Si2O6
                                          [NSNumber numberWithDouble:48000.0],                      // CaMgSi2O6  - MgAlSiAlO6
                                          nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Ca on M2 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on M2 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on M2 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on M1 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on M1 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on M1 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          nil]
                              WithNAtoms:10
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 andWithDependentSpecies:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:@"FeAlSiAlO6",
                                           [NSNumber numberWithDouble:-0.5],
                                           [NSNumber numberWithDouble:0.5],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:0.0], // Excess enthalpy
                                           nil],
                                          [NSArray arrayWithObjects:@"CaAlSiAlO6",
                                           [NSNumber numberWithDouble:-1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:48000.0-32100.0], // Excess enthalpy
                                           nil],
                                          [NSArray arrayWithObjects:@"CaFeSi2O6",
                                           [NSNumber numberWithDouble:-0.5],
                                           [NSNumber numberWithDouble:0.5],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:-32100.0/2.0], // Excess enthalpy
                                           nil],
                                          nil]
                   andWithSiteFunctions:siteFunction
                 andWithSiteFirstDerivative:DsiteFunction
                 andWithSiteSecondDerivative:D2siteFunction
				])) {
		// specific no initializers
	}
	return self;
}

// --> Override the SolutionPhaseProtocol function
-(BOOL)testPermissibleValuesOfComponents:(double *)m {
    // [0] Mg2Si2O6 [1] Fe2Si2O6 [2] MgAl2SiO6 [3] CaMgSi2O6
    double mTotal = m[0] + m[1] + m[2] + m[3];
    if (mTotal  < 0.0) return NO;
    if (mTotal == 0.0) mTotal = 1.0;
    double mMg = 2.0*m[0]/mTotal + m[2]/mTotal + m[3]/mTotal;
    double mFe = 2.0*m[1]/mTotal;
    double mAl = 2.0*m[2]/mTotal;
    double mCa =     m[3]/mTotal;
    double mSi = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + m[2]/mTotal + 2.0*m[3]/mTotal;

	BOOL result = YES;
	result &= ((mMg >= 0.0) && (mMg <= 2.0));
	result &= ((mFe >= 0.0) && (mFe <= 2.0));
	result &= ((mAl >= 0.0) && (mAl <= 2.0));
	result &= ((mCa >= 0.0) && (mCa <= 1.0));
	result &= ((mSi >= 1.0) && (mSi <= 2.0));
	return result;
}

@end

@implementation ClinopyroxeneStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[20] - e[26] - (e[13]-e[11])/2.0;
		m[1] = e[26];
		m[2] = (e[12] - (e[20] - e[26] - (e[13]-e[11])/2.0))/2.0;
		m[3] = (e[13]-e[11])/2.0;
		m[4] = e[11];
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2] + m[3] + m[4];
		if (mTotal == 0.0) mTotal = 1.0;
		double mNa = m[4]/mTotal;
		double mMg = m[0]/mTotal + 2.0*m[2]/mTotal;
		double mFe = m[1]/mTotal;
		double mAl = 2.0*m[3]/mTotal + m[4]/mTotal;
		double mCa = m[0]/mTotal + m[1]/mTotal + m[3]/mTotal;
		double mSi = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + 2.0*m[2]/mTotal + m[3]/mTotal + 2.0*m[4]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Na%4.2fCa%4.2fMg%4.2fFe%4.2fAl%4.2fSi%4.2fO6", mNa, mCa, mMg, mFe, mAl, mSi];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<8; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

    TypeSiteFunction siteFunction = ^(double m[NA], double x[NS]) {
        double mTotal = m[0] + m[1] + m[2] + m[3] + m[4]; // order Di, Hd, En, CaTs, Jd
        if (mTotal == 0.0) mTotal = 1.0;
        double xDi   = m[0]/mTotal;
        double xHd   = m[1]/mTotal;
        double xEn   = m[2]/mTotal;
        double xCaTs = m[3]/mTotal;
        double xJd   = m[4]/mTotal;
        double denom = xDi + xHd + 2.0*xEn;
        // sites XNaM2, XCaM2, XMgM2, XFeM2, XMgM1, XFeM1, XAlM1
        x[0] = xJd;
        x[1] = xDi + xHd + xCaTs;
        x[2] = xEn*(1.0-xHd/denom);
        x[3] = xEn*xHd/denom;
        x[4] = (xDi + xHd + xEn)*(1.0-xHd/denom);
        x[5] = (xDi + xHd + xEn)*xHd/denom;
        x[6] = xCaTs + xJd;
    };

    TypeDsiteFunction DsiteFunction = ^(double m[NA], double dx[NS][NA]) {
        double mTotal = m[0] + m[1] + m[2] + m[3] + m[4];
        if (mTotal == 0.0) mTotal = 1.0;
        double xDi   = m[0]/mTotal;
        double xHd   = m[1]/mTotal;
        double xEn   = m[2]/mTotal;
        //double xCaTs = m[3]/mTotal;
        //double xJd   = m[4]/mTotal;
        double dxDi[5], dxHd[5], dxEn[5], dxCaTs[5], dxJd[5];
        for (NSUInteger i=0; i<5;i++) {
            dxDi[i]   = (i == 0) ? (1.0-m[0])/mTotal : -m[0]/mTotal;
            dxHd[i]   = (i == 1) ? (1.0-m[1])/mTotal : -m[1]/mTotal;
            dxEn[i]   = (i == 2) ? (1.0-m[2])/mTotal : -m[2]/mTotal;
            dxCaTs[i] = (i == 3) ? (1.0-m[3])/mTotal : -m[3]/mTotal;
            dxJd[i]   = (i == 4) ? (1.0-m[4])/mTotal : -m[4]/mTotal;
        }
        double denom = xDi + xHd + 2.0*xEn;
        double dinvDenom[5];
        for (NSUInteger i=0; i<5; i++) {
            dinvDenom[i] = -(dxDi[i] + dxHd[i] + 2.0*dxEn[i])/denom/denom;
            // xJd
            dx[0][i] = dxJd[i];
            // xDi + xHd + xCaTs
            dx[1][i] = dxDi[i] + dxHd[i] + dxCaTs[i];
            // xEn*(1.0-xHd/denom)
            dx[2][i] = dxEn[i]*(1.0-xHd/denom) - xEn*dxHd[i]/denom - xEn*xHd*dinvDenom[i];
            // xEn*xHd/denom
            dx[3][i] = dxEn[i]*xHd/denom + xEn*dxHd[i]/denom + xEn*xHd*dinvDenom[i];
            // (xDi + xHd + xEn)*(1.0-xHd/denom)
            dx[4][i] = (dxDi[i] + dxHd[i] + dxEn[i])*(1.0-xHd/denom) - (xDi + xHd + xEn)*dxHd[i]/denom - (xDi + xHd + xEn)*xHd*dinvDenom[i];
            // (xDi + xHd + xEn)*xHd/denom
            dx[5][i] = (dxDi[i] + dxHd[i] + dxEn[i])*xHd/denom + (xDi + xHd + xEn)*dxHd[i]/denom + (xDi + xHd + xEn)*xHd*dinvDenom[i];
            // xCaTs + xJd
            dx[6][i] = dxCaTs[i] + dxJd[i];
        }
    };

    TypeD2siteFunction D2siteFunction = ^(double m[NA], double d2x[NS][NA][NA]) {
        double mTotal = m[0] + m[1] + m[2] + m[3] + m[4];
        if (mTotal == 0.0) mTotal = 1.0;
        double xDi   = m[0]/mTotal;
        double xHd   = m[1]/mTotal;
        double xEn   = m[2]/mTotal;
        // double xCaTs = m[3]/mTotal;
        // double xJd   = m[4]/mTotal;
        double dxDi[5], dxHd[5], dxEn[5], dxCaTs[5], dxJd[5], d2xDi[5][5], d2xHd[5][5], d2xEn[5][5], d2xCaTs[5][5], d2xJd[5][5];
        for (NSUInteger i=0; i<5; i++) {
            dxDi[i]   = (i == 0) ? (1.0-m[0])/mTotal : -m[0]/mTotal;
            dxHd[i]   = (i == 1) ? (1.0-m[1])/mTotal : -m[1]/mTotal;
            dxEn[i]   = (i == 2) ? (1.0-m[2])/mTotal : -m[2]/mTotal;
            dxCaTs[i] = (i == 3) ? (1.0-m[3])/mTotal : -m[3]/mTotal;
            dxJd[i]   = (i == 4) ? (1.0-m[4])/mTotal : -m[4]/mTotal;
        }
        for (NSUInteger j=0; j<5; j++) {
            d2xDi[0][j]   = (j == 0) ? -1.0/mTotal - (1.0-m[0])/mTotal/mTotal : -(1.0-m[0])/mTotal/mTotal;
            d2xHd[0][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xEn[0][j]   = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xCaTs[0][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xJd[0][j]   = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xDi[1][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xHd[1][j]   = (j == 1) ? -1.0/mTotal - (1.0-m[1])/mTotal/mTotal : -(1.0-m[1])/mTotal/mTotal;
            d2xEn[1][j]   = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xCaTs[1][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xJd[1][j]   = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xDi[2][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xHd[2][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xEn[2][j]   = (j == 2) ? -1.0/mTotal - (1.0-m[2])/mTotal/mTotal : -(1.0-m[2])/mTotal/mTotal;
            d2xCaTs[2][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xJd[2][j]   = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xDi[3][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xHd[3][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xEn[3][j]   = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xCaTs[3][j] = (j == 3) ? -1.0/mTotal - (1.0-m[3])/mTotal/mTotal : -(1.0-m[3])/mTotal/mTotal;
            d2xJd[3][j]   = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xDi[4][j]   = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xHd[4][j]   = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xEn[4][j]   = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xCaTs[4][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xJd[4][j]   = (j == 4) ? -1.0/mTotal - (1.0-m[4])/mTotal/mTotal : -(1.0-m[4])/mTotal/mTotal;
        }
        double denom = xDi + xHd + 2.0*xEn;;
        double dinvDenom[5], d2invDenom[5][5];
        for (NSUInteger i=0; i<5; i++) {
            dinvDenom[i] = -(dxDi[i] + dxHd[i] + 2.0*dxEn[i])/denom/denom;
            for (NSUInteger j=0; j<5; j++) {
                d2invDenom[i][j] = -(d2xDi[i][j] + d2xHd[i][j] + 2.0*d2xEn[i][j])/denom/denom
                + 2.0*(dxDi[i] + dxHd[i] + 2.0*dxEn[i])*(dxDi[j] + dxHd[j] + 2.0*dxEn[j])/denom/denom/denom;
                // xJd
                d2x[0][i][j] = d2xJd[i][j];
                // xDi + xHd + xCaTs
                d2x[1][i][j] = d2xDi[i][j] + d2xHd[i][j] + d2xCaTs[i][j];
                // xEn*(1.0-xHd/denom)
                d2x[2][i][j] = d2xEn[i][j]*(1.0-xHd/denom) - dxEn[i]*dxHd[j]/denom - dxEn[i]*xHd*dinvDenom[j]
                - dxEn[j]*dxHd[i]/denom - xEn*d2xHd[i][j]/denom - xEn*dxHd[i]*dinvDenom[j]
                - dxEn[j]*xHd*dinvDenom[i] - xEn*dxHd[j]*dinvDenom[i] - xEn*xHd*d2invDenom[i][j];
                // xEn*xHd/denom
                d2x[3][i][j] = d2xEn[i][j]*xHd/denom + dxEn[i]*dxHd[j]/denom + dxEn[i]*xHd*dinvDenom[j]
                + dxEn[j]*dxHd[i]/denom + xEn*d2xHd[i][j]/denom + xEn*dxHd[i]*dinvDenom[j]
                + dxEn[j]*xHd*dinvDenom[i] + xEn*dxHd[j]*dinvDenom[i] + xEn*xHd*d2invDenom[i][j];
                // (xDi + xHd + xEn)*(1.0-xHd/denom)
                d2x[4][i][j] = (d2xDi[i][j] + d2xHd[i][j] + d2xEn[i][j])*(1.0-xHd/denom)
                - (dxDi[i] + dxHd[i] + dxEn[i])*dxHd[j]/denom - (dxDi[i] + dxHd[i] + dxEn[i])*xHd*dinvDenom[j]
                - (dxDi[j] + dxHd[j] + dxEn[j])*dxHd[i]/denom - (xDi + xHd + xEn)*d2xHd[i][j]/denom - (xDi + xHd + xEn)*dxHd[i]*dinvDenom[j]
                - (dxDi[j] + dxHd[j] + dxEn[j])*xHd*dinvDenom[i] - (xDi + xHd + xEn)*dxHd[j]*dinvDenom[i] - (xDi + xHd + xEn)*xHd*d2invDenom[i][j];
                // (xDi + xHd + xEn)*xHd/denom
                d2x[5][i][j] = (d2xDi[i][j] + d2xHd[i][j] + d2xEn[i][j])*xHd/denom + (dxDi[i] + dxHd[i] + dxEn[i])*dxHd[j]/denom
                + (dxDi[i] + dxHd[i] + dxEn[i])*xHd*dinvDenom[j]
                + (dxDi[j] + dxHd[j] + dxEn[j])*dxHd[i]/denom + (xDi + xHd + xEn)*d2xHd[i][j]/denom + (xDi + xHd + xEn)*dxHd[i]*dinvDenom[j]
                + (dxDi[j] + dxHd[j] + dxEn[j])*xHd*dinvDenom[i] + (xDi + xHd + xEn)*dxHd[j]*dinvDenom[i] + (xDi + xHd + xEn)*xHd*d2invDenom[i][j];
                // xCaTs + xJd
                d2x[6][i][j] = d2xCaTs[i][j] + d2xJd[i][j];
            }
        }
    };

    if ((self = [super initWithPhaseName:@"Clinopyroxene"
                             withSpecies:[NSArray arrayWithObjects:[[DiopsideStixrude alloc] init], // CaMgSi2O6
                                          [[HedenbergiteStixrude alloc] init],                      // CaFeSi2O6
                                          [[ClinoenstatiteStixrude alloc] init],                    // Mg2Si2O6
                                          [[CaTschermaksStixrude alloc] init],                      // CaAlSiAlO6
                                          [[JadeiteStixrude alloc] init],                           // NaAlSi2O6
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:3.5],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // CaFeSi2O6  - CaMgSi2O6
                                          [NSNumber numberWithDouble:24700.0],                      // Mg2Si2O6   - CaMgSi2O6
                                          [NSNumber numberWithDouble:24700.0],                      // Mg2Si2O6   - CaFeSi2O6
                                          [NSNumber numberWithDouble:26000.0],                      // CaAlSiAlO6 - CaMgSi2O6
                                          [NSNumber numberWithDouble:0.0],                          // CaAlSiAlO6 - CaFeSi2O6
                                          [NSNumber numberWithDouble:60600.0],                      // CaAlSiAlO6 - Mg2Si2O6
                                          [NSNumber numberWithDouble:24300.0],                      // NaAlSi2O6  - CaMgSi2O6
                                          [NSNumber numberWithDouble:0.0],                          // NaAlSi2O6  - CaFeSi2O6
                                          [NSNumber numberWithDouble:0.0],                          // NaAlSi2O6  - Mg2Si2O6
                                          [NSNumber numberWithDouble:10000.0],                      // NaAlSi2O6  - CaAlSiAlO6
                                          nil]
                               withSites:[NSArray arrayWithObjects: // XNaM2, XCaM2, XMgM2, XFeM2, XMgM1, XFeM1, XAlM1
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Na on M2 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Ca on M2 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on M2 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on M2 site
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:0.0], nil],
                                          nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on M1 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on M1 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on M1 site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],  // Fe M2 has been added; Al,Si on tetrahedral has been eliminated
                                          nil]
                              WithNAtoms:10
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 andWithDependentSpecies:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:@"Fe2Si2O6",
                                           [NSNumber numberWithDouble:-2.0],
                                           [NSNumber numberWithDouble:2.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:0.0], // Excess enthalpy
                                           nil],
                                          [NSArray arrayWithObjects:@"MgAlSiAlO6",
                                           [NSNumber numberWithDouble:-1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:8320.634920646902], // Excess enthalpy
                                           nil],
                                          [NSArray arrayWithObjects:@"FeAlSiAlO6",
                                           [NSNumber numberWithDouble:-2.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:-3234.920634894166], // Excess enthalpy
                                           nil],
                                          nil]
                   andWithSiteFunctions:siteFunction
                 andWithSiteFirstDerivative:DsiteFunction
                 andWithSiteSecondDerivative:D2siteFunction
				])) {
		// specific no initializers
	}
	return self;
}
/*
 double xCaM2 = m[0]/mTotal + m[1]/mTotal + m[3]/mTotal;
 double xMgM2 = m[2]/mTotal;
 double xNaM2 = m[4]/mTotal;
 double xMgM1 = m[0]/mTotal + m[2]/mTotal;
 double xFeM1 = m[1]/mTotal;
 double xAlM1 = m[3]/mTotal + m[4]/mTotal;
 double xAlTet = m[3]/2.0/mTotal;
 double xSiTet = m[0]/mTotal + m[1]/mTotal + m[2]/mTotal + m[3]/2.0/mTotal + m[4]/mTotal;

 result &= ((xCaM2  >= 0.0) && (xCaM2  <= 1.0));
 result &= ((xMgM2  >= 0.0) && (xMgM2  <= 1.0));
 result &= ((xNaM2  >= 0.0) && (xNaM2  <= 1.0));
 result &= ((xMgM1  >= 0.0) && (xMgM1  <= 1.0));
 result &= ((xFeM1  >= 0.0) && (xFeM1  <= 1.0));
 result &= ((xAlM1  >= 0.0) && (xAlM1  <= 1.0));
 result &= ((xAlTet >= 0.0) && (xAlTet <= 1.0));
 result &= ((xSiTet >= 0.0) && (xSiTet <= 1.0));
 NSLog(@"cpx permissible test, xAlM1 = %g. m[3] = %g. m[4] = %g", xAlM1, m[3], m[4]);
 */


// --> Override the SolutionPhaseProtocol function
-(BOOL)testPermissibleValuesOfComponents:(double *)m {
    // [0] CaMgSi2O6 [1] CaFeSi2O6 [2] Mg2Si2O6 [3] CaAl2SiO6 [4] NaAlSi2O6

    double mTotal = m[0] + m[1] + m[2] + m[3] + m[4];
    if (mTotal  < 0.0) return NO;
    if (mTotal == 0.0) mTotal = 1.0;

    double mNa = m[4]/mTotal;                             // 0 <-> 1
    double mMg = m[0]/mTotal + 2.0*m[2]/mTotal;           // 0 <-> 2
    double mFe = m[1]/mTotal;                             // 0 <-> 2
    double mAl = 2.0*m[3]/mTotal + m[4]/mTotal;           // 0 <-> 2
    double mCa = m[0]/mTotal + m[1]/mTotal + m[3]/mTotal; // 0 <-> 1
    double mSi = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + 2.0*m[2]/mTotal + m[3]/mTotal + 2.0*m[4]/mTotal; // 1 <-> 2

	BOOL result = YES;
	result &= ((mNa >= 0.0) && (mNa <= 1.0));
	result &= ((mMg >= 0.0) && (mMg <= 2.0));
	result &= ((mFe >= 0.0) && (mFe <= 2.0));
	result &= ((mAl >= 0.0) && (mAl <= 2.0));
	result &= ((mCa >= 0.0) && (mCa <= 1.0));
	result &= ((mSi >= 1.0) && (mSi <= 2.0));

	return result;
}

@end

@implementation HPClinopyroxeneStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12]/2.0;
		m[1] = e[26]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = 2.0*m[0]/mTotal;
		double mFe = 2.0*m[1]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fSi2O6", mMg, mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"HP-Clinopyroxene"
                             withSpecies:[NSArray arrayWithObjects:[[MgHpCpxStixrude alloc] init], [[FeHpCpxStixrude alloc] init], nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Mg on Oct sites
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil], nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Fe on large site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil], nil],
                                          nil]
                              WithNAtoms:10
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation AkimotoiteStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12];
		m[1] = e[26];
		m[2] = e[13]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg =     m[0]/mTotal;
		double mFe =     m[1]/mTotal;
		double mAl = 2.0*m[2]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fAl%4.2fSi%4.2fO3", mMg, mFe, mAl, mMg+mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<3; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Akimotoite"
                             withSpecies:[NSArray arrayWithObjects:[[MgAkimotoiteStixrude alloc] init], // MgSiO3
                                          [[FeAkimotoiteStixrude alloc] init],                          // FeSiO3
                                          [[AlAkimotoiteStixrude alloc] init],                          // AlAlO3
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // MgSiO3 - FeSiO3
                                          [NSNumber numberWithDouble:66000.0],                      // MgSiO3 - AlAlO3
                                          [NSNumber numberWithDouble:0.0],                          // FeSiO3 - AlAlO3
                                          nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Al on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Si on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          nil]
                              WithNAtoms:5
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
				])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation GarnetStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = (e[13] - 2.0*e[26]/3.0 - 2.0*e[20]/3.0 - 2.0*e[11]/2.0)/2.0;
		m[1] = e[26]/3.0;
		m[2] = e[20]/3.0;
		m[3] = (e[12] - 3.0*m[0])/4.0;
		m[4] = e[11]/2.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2] + m[3] + m[4];
		if (mTotal == 0.0) mTotal = 1.0;
		double mNa = 2.0*m[4]/mTotal;
		double mMg = 3.0*m[0]/mTotal + 4.0*m[3]/mTotal;
		double mFe = 3.0*m[1];
		double mAl = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + 2.0*m[2]/mTotal + 2.0*m[4]/mTotal;
		double mCa = 3.0*m[2]/mTotal;
		double mSi = 3.0 + m[3]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Na%4.2fCa%4.2fMg%4.2fFe%4.2fAl%4.2fSi%4.2fO12", mNa, mCa, mMg, mFe, mAl, mSi];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<6; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

    TypeSiteFunction siteFunction = ^(double m[NA], double x[NS]) {
        double mTotal = m[0] + m[1] + m[2] + m[3] + m[4]; // order Py, Al, Gr, Mj, Na
        if (mTotal == 0.0) mTotal = 1.0;
        double xPy = m[0]/mTotal;
        double xAl = m[1]/mTotal;
        double xGr = m[2]/mTotal;
        double xMj = m[3]/mTotal;
        double xNa = m[4]/mTotal;
        double denom = 1.0 - xGr - xNa + xMj/3.0;
        // sites XNaDo, XCaDo, XMgDo, XFeDo, XAlDo, XMgOc, XFeOc, XAlOc, XSiOc
        x[0] = 2.0*xNa/3.0;
        x[1] = xGr;
        x[2] = (1.0-xGr-xNa)*(1.0-xAl/denom);
        x[3] = xAl*(1.0-xGr-xNa)/denom;
        x[4] = xNa/3.0;
        x[5] = xMj*(1.0-xAl/denom)/2.0;
        x[6] = xMj*xAl/denom/2.0;
        x[7] = xPy + xAl + xGr + xNa/2.0;
        x[8] = xMj/2.0 + xNa/2.0;
    };

    TypeDsiteFunction DsiteFunction = ^(double m[NA], double dx[NS][NA]) {
        double mTotal = m[0] + m[1] + m[2] + m[3] + m[4]; // order Py, Al, Gr, Mj, Na
        if (mTotal == 0.0) mTotal = 1.0;
        // double xPy = m[0]/mTotal;
        double xAl = m[1]/mTotal;
        double xGr = m[2]/mTotal;
        double xMj = m[3]/mTotal;
        double xNa = m[4]/mTotal;
        double dxPy[5], dxAl[5], dxGr[5], dxMj[5], dxNa[5];
        for (NSUInteger i=0; i<5;i++) {
            dxPy[i] = (i == 0) ? (1.0-m[0])/mTotal : -m[0]/mTotal;
            dxAl[i] = (i == 1) ? (1.0-m[1])/mTotal : -m[1]/mTotal;
            dxGr[i] = (i == 2) ? (1.0-m[2])/mTotal : -m[2]/mTotal;
            dxMj[i] = (i == 3) ? (1.0-m[3])/mTotal : -m[3]/mTotal;
            dxNa[i] = (i == 4) ? (1.0-m[4])/mTotal : -m[4]/mTotal;
        }
        double denom = 1.0 - xGr - xNa + xMj/3.0;
        double dinvDenom[5];
        for (NSUInteger i=0; i<5; i++) {
            dinvDenom[i] = -(-dxGr[i] - dxNa[i] + dxMj[i]/3.0)/denom/denom;
            // 2.0*xNa/3.0
            dx[0][i] = 2.0*dxNa[i]/3.0;
            // xGr
            dx[1][i] = dxGr[i];
            // (1.0-xGr-xNa)*(1.0-xAl/denom)
            dx[2][i] = (-dxGr[i]-dxNa[i])*(1.0-xAl/denom) - (1.0-xGr-xNa)*dxAl[i]/denom - (1.0-xGr-xNa)*xAl*dinvDenom[i];
            // xAl*(1.0-xGr-xNa)/denom
            dx[3][i] = dxAl[i]*(1.0-xGr-xNa)/denom + xAl*(-dxGr[i]-dxNa[i])/denom + xAl*(1.0-xGr-xNa)*dinvDenom[i];
            // xNa/3.0
            dx[4][i] = dxNa[i]/3.0;
            // xMj*(1.0-xAl/denom)/2.0
            dx[5][i] = dxMj[i]*(1.0-xAl/denom)/2.0 - xMj*dxAl[i]/denom/2.0 - xMj*xAl*dinvDenom[i]/2.0;
            // xMj*xAl/denom/2.0
            dx[6][i] = dxMj[i]*xAl/denom/2.0 + xMj*dxAl[i]/denom/2.0 + xMj*xAl*dinvDenom[i]/2.0;
            // xPy + xAl + xGr + xNa/2.0
            dx[7][i] = dxPy[i] + dxAl[i] + dxGr[i] + dxNa[i]/2.0;
            // xMj/2.0 + xNa/2.0
            dx[8][i] = dxMj[i]/2.0 + dxNa[i]/2.0;
        }
    };

    TypeD2siteFunction D2siteFunction = ^(double m[NA], double d2x[NS][NA][NA]) {
        double mTotal = m[0] + m[1] + m[2] + m[3] + m[4]; // order Py, Al, Gr, Mj, Na
        if (mTotal == 0.0) mTotal = 1.0;
        // double xPy = m[0]/mTotal;
        double xAl = m[1]/mTotal;
        double xGr = m[2]/mTotal;
        double xMj = m[3]/mTotal;
        double xNa = m[4]/mTotal;
        double dxPy[5], dxAl[5], dxGr[5], dxMj[5], dxNa[5], d2xPy[5][5], d2xAl[5][5], d2xGr[5][5], d2xMj[5][5], d2xNa[5][5];
        for (NSUInteger i=0; i<5; i++) {
            dxPy[i] = (i == 0) ? (1.0-m[0])/mTotal : -m[0]/mTotal;
            dxAl[i] = (i == 1) ? (1.0-m[1])/mTotal : -m[1]/mTotal;
            dxGr[i] = (i == 2) ? (1.0-m[2])/mTotal : -m[2]/mTotal;
            dxMj[i] = (i == 3) ? (1.0-m[3])/mTotal : -m[3]/mTotal;
            dxNa[i] = (i == 4) ? (1.0-m[4])/mTotal : -m[4]/mTotal;
        }
        for (NSUInteger j=0; j<5; j++) {
            d2xPy[0][j] = (j == 0) ? -1.0/mTotal - (1.0-m[0])/mTotal/mTotal : -(1.0-m[0])/mTotal/mTotal;
            d2xAl[0][j] = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xGr[0][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xMj[0][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xNa[0][j] = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xPy[1][j] = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xAl[1][j] = (j == 1) ? -1.0/mTotal - (1.0-m[1])/mTotal/mTotal : -(1.0-m[1])/mTotal/mTotal;
            d2xGr[1][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xMj[1][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xNa[1][j] = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xPy[2][j] = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xAl[2][j] = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xGr[2][j] = (j == 2) ? -1.0/mTotal - (1.0-m[2])/mTotal/mTotal : -(1.0-m[2])/mTotal/mTotal;
            d2xMj[2][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xNa[2][j] = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xPy[3][j] = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xAl[3][j] = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xGr[3][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xMj[3][j] = (j == 3) ? -1.0/mTotal - (1.0-m[3])/mTotal/mTotal : -(1.0-m[3])/mTotal/mTotal;
            d2xNa[3][j] = (j == 4) ? -1.0/mTotal + m[4]/mTotal/mTotal : m[4]/mTotal/mTotal;
            d2xPy[4][j] = (j == 0) ? -1.0/mTotal + m[0]/mTotal/mTotal : m[0]/mTotal/mTotal;
            d2xAl[4][j] = (j == 1) ? -1.0/mTotal + m[1]/mTotal/mTotal : m[1]/mTotal/mTotal;
            d2xGr[4][j] = (j == 2) ? -1.0/mTotal + m[2]/mTotal/mTotal : m[2]/mTotal/mTotal;
            d2xMj[4][j] = (j == 3) ? -1.0/mTotal + m[3]/mTotal/mTotal : m[3]/mTotal/mTotal;
            d2xNa[4][j] = (j == 4) ? -1.0/mTotal - (1.0-m[4])/mTotal/mTotal : -(1.0-m[4])/mTotal/mTotal;
        }
        double denom = 1.0 - xGr - xNa + xMj/3.0;
        double dinvDenom[5], d2invDenom[5][5];
        for (NSUInteger i=0; i<5; i++) {
            dinvDenom[i] = -(-dxGr[i] - dxNa[i] + dxMj[i]/3.0)/denom/denom;
            for (NSUInteger j=0; j<5; j++) {
                d2invDenom[i][j] = -(-d2xGr[i][j] - d2xNa[i][j] + d2xMj[i][j]/3.0)/denom/denom
                +2.0*(-dxGr[i] - dxNa[i] + dxMj[i]/3.0)*(-dxGr[j] - dxNa[j] + dxMj[j]/3.0)/denom/denom/denom;
                // 2.0*xNa/3.0
                d2x[0][i][j] = 2.0*d2xNa[i][j]/3.0;
                // xGr
                d2x[1][i][j] = d2xGr[i][j];
                // (1.0-xGr-xNa)*(1.0-xAl/denom)
                d2x[2][i][j] = (-d2xGr[i][j]-d2xNa[i][j])*(1.0-xAl/denom) - (-dxGr[i]-dxNa[i])*dxAl[j]/denom
                - (-dxGr[i]-dxNa[i])*xAl*dinvDenom[j]
                - (-dxGr[j]-dxNa[j])*dxAl[i]/denom - (1.0-xGr-xNa)*d2xAl[i][j]/denom - (1.0-xGr-xNa)*dxAl[i]*dinvDenom[j]
                - (-dxGr[j]-dxNa[j])*xAl*dinvDenom[i] - (1.0-xGr-xNa)*dxAl[j]*dinvDenom[i] - (1.0-xGr-xNa)*xAl*d2invDenom[i][j];
                // xAl*(1.0-xGr-xNa)/denom
                d2x[3][i][j] = d2xAl[i][j]*(1.0-xGr-xNa)/denom + dxAl[i]*(-dxGr[j]-dxNa[j])/denom + dxAl[i]*(1.0-xGr-xNa)*dinvDenom[j]
                + dxAl[j]*(-dxGr[i]-dxNa[i])/denom + xAl*(-d2xGr[i][j]-d2xNa[i][j])/denom + xAl*(-dxGr[i]-dxNa[i])*dinvDenom[j]
                + dxAl[j]*(1.0-xGr-xNa)*dinvDenom[i] + xAl*(-dxGr[j]-dxNa[j])*dinvDenom[i] + xAl*(1.0-xGr-xNa)*d2invDenom[i][j];
                // xNa/3.0
                d2x[4][i][j] = d2xNa[i][j]/3.0;
                // xMj*(1.0-xAl/denom)/2.0
                d2x[5][i][j] = d2xMj[i][j]*(1.0-xAl/denom)/2.0 - dxMj[i]*dxAl[j]/denom/2.0 - dxMj[i]*xAl*dinvDenom[j]/2.0
                - dxMj[j]*dxAl[i]/denom/2.0 - xMj*d2xAl[i][j]/denom/2.0 - xMj*dxAl[i]*dinvDenom[j]/2.0
                - dxMj[j]*xAl*dinvDenom[i]/2.0 - xMj*dxAl[j]*dinvDenom[i]/2.0 - xMj*xAl*d2invDenom[i][j]/2.0;
                // xMj*xAl/denom/2.0
                d2x[6][i][j] = d2xMj[i][j]*xAl/denom/2.0 +  dxMj[i]*dxAl[j]/denom/2.0 +  dxMj[i]*xAl*dinvDenom[j]/2.0
                + dxMj[j]*dxAl[i]/denom/2.0 + xMj*d2xAl[i][j]/denom/2.0 + xMj*dxAl[i]*dinvDenom[j]/2.0
                + dxMj[j]*xAl*dinvDenom[i]/2.0 + xMj*dxAl[j]*dinvDenom[i]/2.0 + xMj*xAl*d2invDenom[i][j]/2.0;
                // xPy + xAl + xGr + xNa/2.0
                d2x[7][i][j] = d2xPy[i][j] + d2xAl[i][j] + d2xGr[i][j] + d2xNa[i][j]/2.0;
                // xMj/2.0 + xNa/2.0
                d2x[8][i][j] = d2xMj[i][j]/2.0 + d2xNa[i][j]/2.0;
            }
        }
    };

	if ((self = [super initWithPhaseName:@"Garnet"
                             withSpecies:[NSArray arrayWithObjects:[[PyropeStixrude alloc] init], // Mg3   AlAl Si3O12
                                          [[AlmandineStixrude alloc] init],                       // Fe3   AlAl Si3O12
                                          [[GrossularStixrude alloc] init],                       // Ca3   AlAl Si3O12
                                          [[MajoriteStixrude alloc] init],                        // Mg3   MgSi Si3O12
                                          [[NaMajoriteStixrude alloc] init],                      // Na2Al AlSi Si3O12
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // Fe3Al2Si3O12 - Mg3Al2Si3O12
                                          [NSNumber numberWithDouble:30000.0],                      // Ca3Al2Si3O12 - Mg3Al2Si3O12
                                          [NSNumber numberWithDouble:0.0],                          // Ca3Al2Si3O12 - Fe3Al2Si3O12
                                          [NSNumber numberWithDouble:21300.0],                      // Mg4Si4O12    - Mg3Al2Si3O12
                                          [NSNumber numberWithDouble:0.0],                          // Mg4Si4O12    - Fe3Al2Si3O12
                                          [NSNumber numberWithDouble:58000.0],                      // Mg4Si4O12    - Ca3Al2Si3O12
                                          [NSNumber numberWithDouble:0.0],                          // Na2Al2Si4O12 - Mg3Al2Si3O12
                                          [NSNumber numberWithDouble:0.0],                          // Na2Al2Si4O12 - Fe3Al2Si3O12
                                          [NSNumber numberWithDouble:0.0],                          // Na2Al2Si4O12 - Ca3Al2Si3O12
                                          [NSNumber numberWithDouble:0.0],                          // Na2Al2Si4O12 - Mg4Si4O12
                                          nil]
                               withSites:[NSArray arrayWithObjects: // XNaDo, XCaDo, XMgDo, XFeDo, XAlDo, XMgOc, XFeOc, XAlOc, XSiOc
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:3.0], // Na on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.666666666], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:3.0], // Ca on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:3.0], // Mg on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:3.0], // Fe on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:3.0], // Al on dodecahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.333333333], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Mg on first octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Fe on first octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Al on first octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Si on second octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          nil] // Removed octahedral site splitting, added Fe to octahedral site
                              WithNAtoms:20
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 andWithDependentSpecies:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:@"Fe3SiFeSi3O12",
                                           [NSNumber numberWithDouble:-4.0/3.0],
                                           [NSNumber numberWithDouble:4.0/3.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:1.0],
                                           [NSNumber numberWithDouble:0.0],
                                           [NSNumber numberWithDouble:-28400.0], // Excess enthalpy
                                           nil],
                                          nil]
                   andWithSiteFunctions:siteFunction
                 andWithSiteFirstDerivative:DsiteFunction
                 andWithSiteSecondDerivative:D2siteFunction
				])) {
		[self setMinimumAffinityDiffereneceToAllowConvergenceInSaturationMethod:1.0];
	}
	return self;
}

// --> Override the SolutionPhaseProtocol function
-(BOOL)testPermissibleValuesOfComponents:(double *)m {
    // [0] Mg3Al2Si3O12 [1] Fe3Al2Si3O12 [2] Ca3Al2Si3O12 [3] Mg3MgSiSi3O12 [4] (Na2Al)AlSiSi3O12
    double mTotal = m[0] + m[1] + m[2] + m[3] + m[4];
    if (mTotal  < 0.0) return NO;
    if (mTotal == 0.0) mTotal = 1.0;
    double mNa = 2.0*m[4]/mTotal;
    double mMg = 3.0*m[0]/mTotal + 4.0*m[3]/mTotal;
    double mFe = 3.0*m[1];
    double mAl = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + 2.0*m[2]/mTotal + 2.0*m[4]/mTotal;
    double mCa = 3.0*m[2]/mTotal;
    double mSi = 3.0 + m[3]/mTotal;

	BOOL result = YES;
    result &= ((mNa >= 0.0) && (mNa <= 2.0));
	result &= ((mMg >= 0.0) && (mMg <= 4.0));
	result &= ((mFe >= 0.0) && (mFe <= 4.0));
	result &= ((mAl >= 0.0) && (mAl <= 2.0));
	result &= ((mCa >= 0.0) && (mCa <= 3.0));
	result &= ((mSi >= 3.0) && (mSi <= 4.0));
	return result;
}

@end

@implementation FerropericlaseStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12];
		m[1] = e[26];
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = m[0]/mTotal;
		double mFe = m[1]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fO", mMg, mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Ferropericlase"
                             withSpecies:[NSArray arrayWithObjects:[[PericlaseStixrude alloc] init], [[WuestiteStixrude alloc] init], nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:13000.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil], nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil], nil],
                                          nil]
                              WithNAtoms:2
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation CaFerritePhaseStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12];
		m[1] = e[26];
		m[2] = e[11];
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1] + m[2];
		if (mTotal == 0.0) mTotal = 1.0;
		double mNa = m[2]/mTotal;
		double mMg = m[0]/mTotal;
		double mFe = m[1]/mTotal;
		double mAl = 2.0*m[0]/mTotal + 2.0*m[1]/mTotal + m[2]/mTotal;
		double mSi = m[2]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Na%4.2fMg%4.2fFe%4.2fAl%4.2fSi%4.2fO4", mNa, mMg, mFe, mAl, mSi];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<3; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"CaFerritePhase"
                             withSpecies:[NSArray arrayWithObjects:[[MgCFStixrude alloc] init], // Mg Al Al O4
                                          [[FeCFStixrude alloc] init],                          // Fe Al Al O4
                                          [[NaCFStixrude alloc] init],                          // Na Al Si O4
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], // MgSiO3 - FeSiO3
                                          [NSNumber numberWithDouble:0.0],                          // MgSiO3 - AlAlO3
                                          [NSNumber numberWithDouble:0.0],                          // FeSiO3 - AlAlO3
                                          nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Mg on first site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Fe on first site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], // Na on first site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:1.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Al on second site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:1.0],
                                            [NSNumber numberWithDouble:0.5], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:2.0], // Si on second site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.5], nil],
                                           nil],
                                          nil]
                              WithNAtoms:7
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end

@implementation SpinelStixrude
-(id)init {
	Convertor convertor = ^(double *e, double *m) {
		m[0] = e[12]/4.0;
		m[1] = e[26]/4.0;
	};

	Display display = ^(double *m) {
		double mTotal = m[0] + m[1];
		if (mTotal == 0.0) mTotal = 1.0;
		double mMg = 4.0*m[0]/mTotal;
		double mFe = 4.0*m[1]/mTotal;
		return (NSString *) [NSString stringWithFormat:@"Mg%4.2fFe%4.2fAl8O16", mMg, mFe];
	};

    Adjuster adjuster = ^(double *gamma, double *x)  {
        for (NSUInteger i=0, j=0; i<2; i++) if (x[i] != 0.0) {
            if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
            else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
            j++;
        }
    };

	if ((self = [super initWithPhaseName:@"Spinel"
                             withSpecies:[NSArray arrayWithObjects:[[MgSpinelStixrude alloc] init], // Mg3Al Al7Mg O16
                                          [[HercyniteStixrude alloc] init],                         // Fe3Al Al7Fe O16
                                          nil]
                      withSpeciesWeights:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                      withSpeciesXfactor:[NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0],
                                          [NSNumber numberWithDouble:1.0],
                                          nil]
                              withWarray:[NSArray arrayWithObjects:[NSNumber numberWithDouble:5000.0], nil]
                               withSites:[NSArray arrayWithObjects:
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:4.0], // Mg on tetrahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.75],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:4.0], // Fe on tetrahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.75], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:4.0], // Al on tetrahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.25],
                                            [NSNumber numberWithDouble:0.25], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:8.0], // Mg on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.125],
                                            [NSNumber numberWithDouble:0.0], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:8.0], // Fe on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0],
                                            [NSNumber numberWithDouble:0.125], nil],
                                           nil],
                                          [NSArray arrayWithObjects:[NSNumber numberWithDouble:8.0], // Al on octahedral site
                                           [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.25],
                                            [NSNumber numberWithDouble:0.875], nil],
                                           nil],
                                          nil]
                              WithNAtoms:28
                    withElementConvertor:convertor
                   andWithFormulaDisplay:display
          andWithSaturationStateAdjuster:adjuster
                 ])) {
		// specific no initializers
	}
	return self;
}
@end
