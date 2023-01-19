//
//  FeldsparBerman+FeldsparBerman_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 1/2/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "FeldsparBerman+FeldsparBerman_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation FeldsparBerman (FeldsparBerman_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"whabor",
                         @"wsabor",
                         @"wvabor",
                         @"whorab",
                         @"wsorab",
                         @"wvorab",
                         @"whaban",
                         @"whanab",
                         @"whoran",
                         @"whanor",
                         @"wvanor",
                         @"whabanor",
                         @"wvabanor",
                         nil];
}

-(void)loadParameterUnitsArray {
    unitsOfParameters = [NSArray arrayWithObjects:
                         @"joules",
                         @"joules/K",
                         @"joules/bar",
                         @"joules",
                         @"joules/K",
                         @"joules/bar",
                         @"joules",
                         @"joules",
                         @"joules",
                         @"joules",
                         @"joules/bar",
                         @"joules",
                         @"joules/bar",
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
    if (!unitsOfParameters) [self loadParameterUnitsArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return @"";
    else return [unitsOfParameters objectAtIndex:index];
}

-(double)getValueForParameterName:(NSString *)name {
    NSUInteger index = [namesOfParameters indexOfObject:name];
    double value = PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
    switch (index) {
        case 0:
            value = self.whabor;
            break;
        case 1:
            value = self.wsabor;
            break;
        case 2:
            value = self.wvabor;
            break;
        case 3:
            value = self.whorab;
            break;
        case 4:
            value = self.wsorab;
            break;
        case 5:
            value = self.wvorab;
            break;
        case 6:
            value = self.whaban;
            break;
        case 7:
            value = self.whanab;
            break;
        case 8:
            value = self.whoran;
            break;
        case 9:
            value = self.whanor;
            break;
        case 10:
            value = self.wvanor;
            break;
        case 11:
            value = self.whabanor;
            break;
        case 12:
            value = self.wvabanor;
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
    NSUInteger index = [namesOfParameters indexOfObject:name];
    switch (index) {
        case 0:
            self.whabor = value;
            break;
        case 1:
            self.wsabor = value;
            break;
        case 2:
            self.wvabor = value;
            break;
        case 3:
            self.whorab = value;
            break;
        case 4:
            self.wsorab = value;
            break;
        case 5:
            self.wvorab = value;
            break;
        case 6:
            self.whaban = value;
            break;
        case 7:
            self.whanab = value;
            break;
        case 8:
            self.whoran = value;
            break;
        case 9:
            self.whanor = value;
            break;
        case 10:
            self.wvanor = value;
            break;
        case 11:
            self.whabanor = value;
            break;
        case 12:
            self.wvabanor = value;
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
    DoubleVector *moleFractionWrapper = [self convertMolesToMoleFractions:m];
    double *moleFraction = [moleFractionWrapper pointerToDouble];
    double xab = moleFraction[0];
    double xan = moleFraction[1];
    double xor = moleFraction[2];

    DoubleVector *dgdwWrapper = [[DoubleVector alloc] initWithSize:13];
    double *dgdw = [dgdwWrapper pointerToDouble];

    dgdw[ 0] =         xab*xor*(xor+xan/2.0); // whabor
    dgdw[ 1] =     - t*xab*xor*(xor+xan/2.0); // wsabor
    dgdw[ 2] = (p-1.0)*xab*xor*(xor+xan/2.0); // wvabor
    dgdw[ 3] =         xab*xor*(xab+xan/2.0); // whorab
    dgdw[ 4] =     - t*xab*xor*(xab+xan/2.0); // wsorab
    dgdw[ 5] = (p-1.0)*xab*xor*(xab+xan/2.0); // wvorab
    dgdw[ 6] =         xab*xan*(xan+xor/2.0); // whaban
    dgdw[ 7] =         xab*xan*(xab+xor/2.0); // whanab
    dgdw[ 8] =         xan*xor*(xan+xab/2.0); // whoran
    dgdw[ 9] =         xan*xor*(xor+xab/2.0); // whanor
    dgdw[10] = (p-1.0)*xan*xor*(xor+xab/2.0); // wvanor
    dgdw[11] =         xab*xan*xor;           // whabanor
    dgdw[12] = (p-1.0)*xab*xan*xor;           // wvabanor

    double total = [self totalMolesFromMolesOfComponents:m];
    for (NSUInteger i=0; i<12; i++) dgdw[i] *= total;

    return dgdwWrapper;
}

-(DoubleVector *)getChemicalPotentialDerivativesForParameter:(NSString *)name usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return nil;
    else {
        NSUInteger na = [self numberOfSolutionComponents];
        DoubleVector *chemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];

        double value = [self getValueForParameterName:name];
        double forwardValue = (value != 0.0) ? value*(1.0+sqrt(DBL_EPSILON)) : sqrt(DBL_EPSILON);
        [self setParameterName:name tovalue:forwardValue];
        DoubleVector *forwardChemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];

        DoubleVector *derivativeWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
        double *derivative = [derivativeWrapper pointerToDouble];
        double *chemicalPotentials = [chemicalPotentialsWrapper pointerToDouble];
        double *forwardChemicalPotentials = [forwardChemicalPotentialsWrapper pointerToDouble];
        for (NSUInteger i=0; i<na; i++) derivative[i] = (forwardChemicalPotentials[i] - chemicalPotentials[i])/(forwardValue - value);

        [self setParameterName:name tovalue:value];

        return derivativeWrapper;
    }
}

-(DoubleMatrix *)getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    if (!array || ([array count] == 0)) return nil;
    NSUInteger lengthOfArray = [array count];
    NSUInteger na = [self numberOfSolutionComponents];
    DoubleMatrix *derivativeMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:lengthOfArray andWithColumnSize:na andInitialValue:0.0];
    double **derivativeMatrix = [derivativeMatrixWrapper pointerToPointerToDouble];

    DoubleVector *chemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
    double *chemicalPotentials = [chemicalPotentialsWrapper pointerToDouble];

    for (NSUInteger j=0; j<lengthOfArray; j++) {
        NSUInteger index = [namesOfParameters indexOfObject:[array objectAtIndex:j]];
        if (index != NSNotFound) {
            double value = [self getValueForParameterName:[array objectAtIndex:j]];
            double forwardValue = (value != 0.0) ? value*(1.0+sqrt(DBL_EPSILON)) : sqrt(DBL_EPSILON);
            [self setParameterName:[array objectAtIndex:j] tovalue:forwardValue];
            DoubleVector *forwardChemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
            double *forwardChemicalPotentials = [forwardChemicalPotentialsWrapper pointerToDouble];

            for (NSUInteger i=0; i<na; i++) derivativeMatrix[j][i] = (forwardChemicalPotentials[i] - chemicalPotentials[i])/(forwardValue - value);
            [self setParameterName:[array objectAtIndex:j] tovalue:value];
        }
    }

    return derivativeMatrixWrapper;
}

-(BOOL)supportsParameterCalibration {
    return YES;
}

@end
