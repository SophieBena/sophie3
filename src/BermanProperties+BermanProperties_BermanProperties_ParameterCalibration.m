//
//  BermanProperties+BermanProperties_BermanProperties_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 1/20/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "BermanProperties+BermanProperties_BermanProperties_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation BermanProperties (BermanProperties_BermanProperties_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"delta H",
                         @"S",
                         @"Cp k0",
                         @"Cp k1",
                         @"Cp k2",
                         @"Cp k3",
                         @"Cp l1",
                         @"Cp l2",
                         @"Cp Tt",
                         @"Cp Ht",
                         @"V",
                         @"EOS v1",
                         @"EOS v2",
                         @"EOS v3",
                         @"EOS v4",
                         nil];
}

-(void)loadParameterUnitsArray {
    unitsOfParameters = [NSArray arrayWithObjects:
                         @"joules",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"K",
                         @"joules",
                         @"joules/bar",
                         @"bar^-1",
                         @"bar^-2",
                         @"K^-1",
                         @"K^-2",
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
            value = self.h;
            break;
        case 1:
            value = self.s;
            break;
        case 2:
            value = self.k0;
            break;
        case 3:
            value = self.k1;
            break;
        case 4:
            value = self.k2;
            break;
        case 5:
            value = self.k3;
            break;
        case 6:
            value = self.l1;
            break;
        case 7:
            value = self.l2;
            break;
        case 8:
            value = self.Tt;
            break;
        case 9:
            value = self.deltaH;
            break;
        case 10:
            value = self.v0;
            break;
        case 11:
            value = self.v1;
            break;
        case 12:
            value = self.v2;
            break;
        case 13:
            value = self.v3;
            break;
        case 14:
            value = self.v4;
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
            self.h = value;
            break;
        case 1:
            self.s = value;
            break;
        case 2:
            self.k0 = value;
            break;
        case 3:
            self.k1 = value;
            break;
        case 4:
            self.k2 = value;
            break;
        case 5:
            self.k3 = value;
            break;
        case 6:
            self.l1 = value;
            break;
        case 7:
            self.l2 = value;
            break;
        case 8:
            self.Tt = value;
            break;
        case 9:
            self.deltaH = value;
            break;
        case 10:
            self.v0 = value;
            break;
        case 11:
            self.v1 = value;
            break;
        case 12:
            self.v2 = value;
            break;
        case 13:
            self.v3 = value;
            break;
        case 14:
            self.v4 = value;
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
