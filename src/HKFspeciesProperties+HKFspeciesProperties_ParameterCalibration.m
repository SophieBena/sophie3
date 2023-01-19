//
//  HKFspeciesProperties+HKFspeciesProperties_ParameterCalibration.m
//  PhaseDEWobjC
//
//  Created by Mark Ghiorso on 4/12/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import "HKFspeciesProperties+HKFspeciesProperties_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation HKFspeciesProperties (HKFspeciesProperties_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"deltaGibbsFreeEnergyOfFormationInTheReferenceState",
                         @"deltaEnthalpyOfFormationInTheReferenceState",
                         @"entropyInTheReferenceState",
                         @"volumeInTheReferenceState",
                         @"heatCapacityInTheReferenceState",
                         @"a1HKF",
                         @"a2HKF",
                         @"a3HKF",
                         @"a4HKF",
                         @"c1HKF",
                         @"c2HKF",
                         @"omegaHKF",
                         nil];
}

-(void)loadParameterUnitsArray {
    unitsOfParameters = [NSArray arrayWithObjects:
                         @"joules/mol",
                         @"joules/mol",
                         @"joules/K/mol",
                         @"joules/bar/mol",
                         @"joules/K/mol",
                         @"joules/bar/mol",
                         @"joules/mol",
                         @"joules-K/bar/mol",
                         @"joules-K/mol",
                         @"joules/K/mol",
                         @"joules-K/mol",
                         @"joules/mol",
                         nil];
}

#pragma mark -
#pragma mark Parameter Calibration Protocol methods

-(NSUInteger)getNumberOfFreeParameters {
    if (!namesOfParameters) [self loadParameterNameArray];
    return [namesOfParameters count];
}

-(NSArray *)getArrayOfNamesOfFreeParameters {
    if (!namesOfParameters) [self loadParameterNameArray];
    return namesOfParameters;
}

-(NSString *)getUnitsForParameterName:(NSString *)name {
    if (!namesOfParameters) [self loadParameterNameArray];
    if (!unitsOfParameters) [self loadParameterUnitsArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return @"";
    else return [unitsOfParameters objectAtIndex:index];
}

-(double)getValueForParameterName:(NSString *)name {
    if (!namesOfParameters) [self loadParameterNameArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    double value = PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
    switch (index) {
        case 0:
            value = self.deltaGibbsFreeEnergyOfFormationInTheReferenceState;
            break;
        case 1:
            value = self.deltaEnthalpyOfFormationInTheReferenceState;
            break;
        case 2:
            value = self.entropyInTheReferenceState;
            break;
        case 3:
            value = self.volumeInTheReferenceState;
            break;
        case 4:
            value = self.heatCapacityInTheReferenceState;
            break;
        case 5:
            value = self.a1HKF;
            break;
        case 6:
            value = self.a2HKF;
            break;
        case 7:
            value = self.a3HKF;
            break;
        case 8:
            value = self.a4HKF;
            break;
        case 9:
            value = self.c1HKF;
            break;
        case 10:
            value = self.c2HKF;
            break;
        case 11:
            value = self.omegaHKF;
            break;
        case NSNotFound:
            value = PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
            break;
        default:
            break;
    }
    return value;
}

-(BOOL)setParameterName:(NSString *)name tovalue:(double)value {
    if (!namesOfParameters) [self loadParameterNameArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    switch (index) {
        case 0:
            self.deltaGibbsFreeEnergyOfFormationInTheReferenceState = value;
            break;
        case 1:
            self.deltaEnthalpyOfFormationInTheReferenceState = value;
            break;
        case 2:
            self.entropyInTheReferenceState = value;
            break;
        case 3:
            self.volumeInTheReferenceState = value;
            break;
        case 4:
            self.heatCapacityInTheReferenceState = value;
            break;
        case 5:
            self.a1HKF = value;
            break;
        case 6:
            self.a2HKF = value;
            break;
        case 7:
            self.a3HKF = value;
            break;
        case 8:
            self.a4HKF = value;
            break;
        case 9:
            self.c1HKF = value;
            break;
        case 10:
            self.c2HKF = value;
            break;
        case 11:
            self.omegaHKF = value;
            break;
        case NSNotFound:
            return NO;
            break;
        default:
            break;
    }
    return YES;
}

-(DoubleVector *)getDgDwFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    return nil;
}

-(DoubleVector *)getChemicalPotentialDerivativesForParameter:(NSString *)name usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double value;
    if ((value = [self getValueForParameterName:name]) == PARAMETER_CALIBRATION_VALUE_NOT_FOUND) return nil;
    double chemicalPotential = [self getGibbsFreeEnergyFromT:t andP:p];

    double forwardValue = (value != 0.0) ? value*(1.0+sqrt(DBL_EPSILON)) : sqrt(DBL_EPSILON);
    [self setParameterName:name tovalue:forwardValue];
    double forwardChemicalPotential = [self getGibbsFreeEnergyFromT:t andP:p];

    DoubleVector *derivativeWrapper = [[DoubleVector alloc] initWithSize:1 andInitialValue:0.0];
    double *derivative = [derivativeWrapper pointerToDouble];
    derivative[0] = (forwardChemicalPotential - chemicalPotential)/(forwardValue - value);

    [self setParameterName:name tovalue:value];
    return derivativeWrapper;
}

-(DoubleMatrix *)getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    if (!array || ([array count] == 0)) return nil;
    NSUInteger lengthOfArray = [array count];
    DoubleMatrix *derivativeMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:lengthOfArray andWithColumnSize:1 andInitialValue:0.0];
    double **derivativeMatrix = [derivativeMatrixWrapper pointerToPointerToDouble];

    double chemicalPotential = [self getGibbsFreeEnergyFromT:t andP:p];

    for (NSUInteger j=0; j<lengthOfArray; j++) {
        double value;
        if ((value = [self getValueForParameterName:[array objectAtIndex:j]]) != PARAMETER_CALIBRATION_VALUE_NOT_FOUND) {
            double forwardValue = (value != 0.0) ? value*(1.0+sqrt(DBL_EPSILON)) : sqrt(DBL_EPSILON);
            [self setParameterName:[array objectAtIndex:j] tovalue:forwardValue];
            double forwardChemicalPotential = [self getGibbsFreeEnergyFromT:t andP:p];

            derivativeMatrix[j][0] = (forwardChemicalPotential - chemicalPotential)/(forwardValue - value);
            [self setParameterName:[array objectAtIndex:j] tovalue:value];
        }
    }

    return derivativeMatrixWrapper;
}

-(BOOL)supportsParameterCalibration {
    return YES;
}

@end
