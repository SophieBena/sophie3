//
//  TridymiteBerman+TridymiteBerman_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/16/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "TridymiteBerman+TridymiteBerman_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation TridymiteBerman (TridymiteBerman_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"delta H alpha",
                         @"S alpha",
                         @"Cp k0 alpha",
                         @"Cp k1 alpha",
                         @"Cp k2 alpha",
                         @"Cp k3 alpha",
                         @"Cp l1 alpha",
                         @"Cp l2 alpha",
                         @"Cp Tt alpha",
                         @"Cp Ht alpha",
                         @"V alpha",
                         @"EOS v1 alpha",
                         @"EOS v2 alpha",
                         @"EOS v3 alpha",
                         @"EOS v4 alpha",
                         @"delta H beta",
                         @"S beta",
                         @"Cp k0 beta",
                         @"Cp k1 beta",
                         @"Cp k2 beta",
                         @"Cp k3 beta",
                         @"V beta",
                         @"EOS v1 beta",
                         @"EOS v2 beta",
                         @"EOS v3 beta",
                         @"EOS v4 beta",
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
                         @"joules",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K",
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
            value = alphaTridymite.h;
            break;
        case 1:
            value = alphaTridymite.s;
            break;
        case 2:
            value = alphaTridymite.k0;
            break;
        case 3:
            value = alphaTridymite.k1;
            break;
        case 4:
            value = alphaTridymite.k2;
            break;
        case 5:
            value = alphaTridymite.k3;
            break;
        case 6:
            value = alphaTridymite.l1;
            break;
        case 7:
            value = alphaTridymite.l2;
            break;
        case 8:
            value = alphaTridymite.Tt;
            break;
        case 9:
            value = alphaTridymite.deltaH;
            break;
        case 10:
            value = alphaTridymite.v0;
            break;
        case 11:
            value = alphaTridymite.v1;
            break;
        case 12:
            value = alphaTridymite.v2;
            break;
        case 13:
            value = alphaTridymite.v3;
            break;
        case 14:
            value = alphaTridymite.v4;
            break;
        case 15:
            value = betaTridymite.h;
            break;
        case 16:
            value = betaTridymite.s;
            break;
        case 17:
            value = betaTridymite.k0;
            break;
        case 18:
            value = betaTridymite.k1;
            break;
        case 19:
            value = betaTridymite.k2;
            break;
        case 20:
            value = betaTridymite.k3;
            break;
        case 21:
            value = betaTridymite.v0;
            break;
        case 22:
            value = betaTridymite.v1;
            break;
        case 23:
            value = betaTridymite.v2;
            break;
        case 24:
            value = betaTridymite.v3;
            break;
        case 25:
            value = betaTridymite.v4;
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
        {
            double correction = value - alphaTridymite.h;
            alphaTridymite.h = value;
            betaTridymite.h += correction;
            break;
        }
        case 1:
        {
            double correction = value - alphaTridymite.s;
            alphaTridymite.s = value;
            betaTridymite.s += correction;
            break;
        }
        case 2:
        {
            double correction = value - alphaTridymite.k0;
            alphaTridymite.k0 = value;
            betaTridymite.k0 += correction;
            break;
        }
        case 3:
        {
            double correction = value - alphaTridymite.k1;
            alphaTridymite.k1 = value;
            betaTridymite.k1 += correction;
            break;
        }
        case 4:
        {
            double correction = value - alphaTridymite.k2;
            alphaTridymite.k2 = value;
            betaTridymite.k2 += correction;
            break;
        }
        case 5:
        {
            double correction = value - alphaTridymite.k3;
            alphaTridymite.k3 = value;
            betaTridymite.k3 += correction;
            break;
        }
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
        {
            double correction = value - alphaTridymite.v0;
            alphaTridymite.v0 = value;
            betaTridymite.v0 += correction;
            break;
        }
        case 11:
        {
            double correction = value - alphaTridymite.v1;
            alphaTridymite.v1 = value;
            betaTridymite.v1 += correction;
            break;
        }
        case 12:
        {
            double correction = value - alphaTridymite.v2;
            alphaTridymite.v2 = value;
            betaTridymite.v2 += correction;
            break;
        }
        case 13:
        {
            double correction = value - alphaTridymite.v3;
            alphaTridymite.v3 = value;
            betaTridymite.v3 += correction;
            break;
        }
        case 14:
        {
            double correction = value - alphaTridymite.v4;
            alphaTridymite.v4 = value;
            betaTridymite.v4 += correction;
            break;
        }
        case 15:
        {
            double correction = value - betaTridymite.h;
            betaTridymite.h = value;
            alphaTridymite.h += correction;
            break;
        }
        case 16:
        {
            double correction = value - betaTridymite.s;
            betaTridymite.s = value;
            alphaTridymite.s += correction;
            break;
        }
        case 17:
        {
            double correction = value - betaTridymite.k0;
            betaTridymite.k0 = value;
            alphaTridymite.k0 += correction;
            break;
        }
        case 18:
        {
            double correction = value - betaTridymite.k1;
            betaTridymite.k1 = value;
            alphaTridymite.k1 += correction;
            break;
        }
        case 19:
        {
            double correction = value - betaTridymite.k2;
            betaTridymite.k2 = value;
            alphaTridymite.k2 += correction;
            break;
        }
        case 20:
        {
            double correction = value - betaTridymite.k3;
            betaTridymite.k3 = value;
            alphaTridymite.k3 += correction;
            break;
        }
        case 21:
        {
            double correction = value - betaTridymite.v0;
            betaTridymite.v0 = value;
            alphaTridymite.v0 += correction;
            break;
        }
        case 22:
        {
            double correction = value - betaTridymite.v1;
            betaTridymite.v1 = value;
            alphaTridymite.v1 += correction;
            break;
        }
        case 23:
        {
            double correction = value - betaTridymite.v2;
            betaTridymite.v2 = value;
            alphaTridymite.v2 += correction;
            break;
        }
        case 24:
        {
            double correction = value - betaTridymite.v3;
            betaTridymite.v3 = value;
            alphaTridymite.v3 += correction;
            break;
        }
        case 25:
        {
            double correction = value - betaTridymite.v4;
            betaTridymite.v4 = value;
            alphaTridymite.v4 += correction;
            break;
        }
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
