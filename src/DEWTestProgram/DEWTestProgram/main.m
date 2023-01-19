//
//  main.m
//  DEWTestProgram
//
//  Created by Mark Ghiorso on 11/7/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "DEWFluid.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

#define EPS 1.0e-6

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        DEWFluid *dewFluid;
        double t, p;
        double moles[17];

        dewFluid = [[DEWFluid alloc] init];
        t = 1273.15; // K
        p = 4000.0;  // bars

        moles[ 0] = 55.51; //    X0   H2O
        moles[ 1] =  0.0;  //    X1   CO2     r0
        moles[ 2] =  0.0;  //    X2   O2      r1
        moles[ 3] =  0.0;  //    X3   HF      r2
        moles[ 4] =  5.0;  //    X4   NaOH    r3 0.1
        moles[ 5] =  0.0;  //    X5   Mg(OH)2 r4
        moles[ 6] =  0.0;  //    X6   HAlO2   r5
        moles[ 7] =  0.0;  //    X7   SiO2    r6
        moles[ 8] =  0.0;  //    X8   H3PO4   r7
        moles[ 9] =  0.0;  //    X9   SO2     r8
        moles[10] =  5.1;  //    X10  HCl     r9 0.1
        moles[11] =  0.0;  //    X11  KOH     r10
        moles[12] =  0.0;  //    X12  Ca(OH)2 r11
        moles[13] =  0.0;  //    X13  H2CrO4  r12
        moles[14] =  0.0;  //    X14  Mn(OH)2 r13
        moles[15] =  0.0;  //    X15  Fe(OH)2 r14
        moles[16] =  0.0;  //    X16  Co(OH)2 r1

        NSDictionary *result = [dewFluid getSpeciesMoleFractionsForBulkComposition:moles aT:t andP:p];
        for (NSString *key in [result.allKeys sortedArrayUsingSelector:@selector(caseInsensitiveCompare:)]) {
            double value = [[result objectForKey:key] doubleValue];
            if (value != 0.0) {
                NSLog(@"Concentration of species %@ = %e", [key stringByPaddingToLength:20 withString:@" " startingAtIndex:0], value);
            }
        }

        [[NSUserDefaults standardUserDefaults] setBool:NO forKey:@"DEBUG.SATLOOP.VERBOSE"];
        [[NSUserDefaults standardUserDefaults] setBool:NO forKey:@"DEBUG.SATLOOP.SIMPLE"];

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Call to getActivities ... ");
        DoubleVector *activityWrapper = [dewFluid getActivityFromMolesOfComponents:moles andT:t andP:p];
        double *activity = [activityWrapper pointerToDouble];
        DoubleVector *chemicalPotentialWrapper = [dewFluid getChemicalPotentialFromMolesOfComponents:moles andT:t andP:p];
        double *chemicalPotential = [chemicalPotentialWrapper pointerToDouble];
        for (NSUInteger i=0; i<17; i++) {
            NSString *name = [dewFluid nameOfSolutionSpeciesAtIndex:i];
            double moleFrac = [[result objectForKey:name] doubleValue];
            if (moleFrac != 0.0) NSLog(@"Species: %@ X = %e mu = %e a = %e", [name stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                                       moleFrac, chemicalPotential[i], activity[i]);
        }

        NSUInteger NE = [dewFluid numberOfSolutionSpecies];
        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Call to getChemicalPotentialsOfSpeciesFromMolesOfComponents ... ");
        DoubleVector *chemicalPotentialsOfSpeciesWrapper = [dewFluid chemicalPotentialsOfSpeciesFromMolesOfComponents:moles andT:t andP:p];
        double *chemicalPotentialsOfSpecies = [chemicalPotentialsOfSpeciesWrapper pointerToDouble];
        DoubleVector *activitiesOfSpeciesWrapper = [dewFluid activitiesOfSpeciesFromMolesOfComponents:moles andT:t andP:p];
        double *activitiesOfSpecies = [activitiesOfSpeciesWrapper pointerToDouble];
        for (NSUInteger i=0; i<NE; i++) {
            NSString *name = [dewFluid nameOfSolutionSpeciesAtIndex:i];
            double moleFrac = [[result objectForKey:name] doubleValue];
            if (moleFrac != 0.0) NSLog(@"Species: %@ X = %e mu = %e a = %e", [name stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                                       moleFrac, chemicalPotentialsOfSpecies[i], activitiesOfSpecies[i]);
        }

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Call to getGibbs ...");
        double g = [dewFluid getGibbsFreeEnergyFromMolesOfComponents:moles andT:t andP:p];
        NSLog(@"... value %g", g);

        double gForward = 0.0, molesForward[17];
        for (NSUInteger i=0; i<17; i++) molesForward[i] = moles[i];

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Check chemical potentials, DgDm ...");
        DoubleVector *DgDmWrapper = [dewFluid getDgDmFromMolesOfComponents:moles andT:t andP:p];
        double *DgDm = [DgDmWrapper pointerToDouble];
        for (NSUInteger i=0; i<17; i++) if (moles[i] != 0.0) {
            molesForward[i] = moles[i]*(1.0+EPS);
            gForward = [dewFluid getGibbsFreeEnergyFromMolesOfComponents:molesForward andT:t andP:p];
            molesForward[i] = moles[i];
            NSString *name = [dewFluid nameOfSolutionSpeciesAtIndex:i];
            NSLog(@"... %@ value %g %g by diff %g from species %g", [name stringByPaddingToLength:12 withString:@" " startingAtIndex:0],
                  chemicalPotential[i], DgDm[i], (gForward-g)/moles[i]/EPS, chemicalPotentialsOfSpecies[i]);
        }

        double molesTotal = 0.0, r[16];
        for (NSUInteger i=0; i<17; i++) molesTotal += moles[i];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testD2GDR2withT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testPrimD2GDR2withT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testPrimD2GDS2withT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testPrimD2GDRDSwithT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testAzerowithT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testD2DHDR2withT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testD2DHDS2withT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testFillnSpecieswithT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testFillxSpecieswithT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testD2DHDRDSwithT:t p:p r:r];
        for (NSUInteger i=0; i<16; i++) r[i] = moles[i+1]/molesTotal;
        [dewFluid testDRDSwithT:t p:p r:r];

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Check D2gDm2 ...");
        DoubleMatrix *D2gDm2Wrapper = [dewFluid getD2gDm2FromMolesOfComponents:moles andT:t andP:p];
        double **D2gDm2 = [D2gDm2Wrapper pointerToPointerToDouble];
        for (NSUInteger i=0; i<17; i++) if (moles[i] != 0.0) {
            NSString *nameOfI = [dewFluid nameOfSolutionSpeciesAtIndex:i];

            molesForward[i] = moles[i]*(1.0+EPS);
            DoubleVector *DgDmForwardWrapper = [dewFluid getDgDmFromMolesOfComponents:molesForward andT:t andP:p];
            double *DgDmForward = [DgDmForwardWrapper pointerToDouble];
            molesForward[i] = moles[i]*(1.0-EPS);
            DoubleVector *DgDmBackwardWrapper = [dewFluid getDgDmFromMolesOfComponents:molesForward andT:t andP:p];
            double *DgDmBackward = [DgDmBackwardWrapper pointerToDouble];
            molesForward[i] = moles[i];

            for (NSUInteger j=0; j<17; j++) if (moles[j] != 0.0) {
                NSString *nameOfJ = [dewFluid nameOfSolutionSpeciesAtIndex:j];
                NSLog(@"... %@ %@ value %g by diff %g", [nameOfI stringByPaddingToLength:12 withString:@" " startingAtIndex:0],
                      [nameOfJ stringByPaddingToLength:12 withString:@" " startingAtIndex:0],
                      D2gDm2[i][j], (DgDmForward[j]-DgDmBackward[j])/moles[i]/EPS/2.0);
            }
        }

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Call to getEntropy ...");
        double s = [dewFluid getEntropyFromMolesOfComponents:moles andT:t andP:p];
        gForward = [dewFluid getGibbsFreeEnergyFromMolesOfComponents:moles andT:t*(1.0+EPS) andP:p];
        NSLog(@"... value %g by diff %g", s, -(gForward-g)/t/EPS);

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Call to getVolume ...");
        double v = [dewFluid getVolumeFromMolesOfComponents:moles andT:t andP:p];
        gForward = [dewFluid getGibbsFreeEnergyFromMolesOfComponents:moles andT:t andP:p*(1.0+EPS)];
        NSLog(@"... value %g by diff %g", v, (gForward-g)/p/EPS);

        double sForward = 0.0;

        NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
        NSLog(@"Call to getHeatCapacity ...");
        double cp = [dewFluid getHeatCapacityFromMolesOfComponents:moles andT:t andP:p];
        sForward = [dewFluid getEntropyFromMolesOfComponents:moles andT:t*(1.0+EPS) andP:p];
        NSLog(@"... value %g by diff %g", cp, t*(sForward-s)/t/EPS);

    }
    return 0;
}
