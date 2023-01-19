//
//  CristobaliteBerman+CristobaliteBerman_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/16/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "CristobaliteBerman+CristobaliteBerman_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation CristobaliteBerman (CristobaliteBerman_ParameterCalibration)

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
                         @"dTtdp",
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
                         @"K/bar",
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
            value = alphaCristobalite.h;
            break;
        case 1:
            value = alphaCristobalite.s;
            break;
        case 2:
            value = alphaCristobalite.k0;
            break;
        case 3:
            value = alphaCristobalite.k1;
            break;
        case 4:
            value = alphaCristobalite.k2;
            break;
        case 5:
            value = alphaCristobalite.k3;
            break;
        case 6:
            value = alphaCristobalite.l1;
            break;
        case 7:
            value = alphaCristobalite.l2;
            break;
        case 8:
            value = alphaCristobalite.Tt;
            break;
        case 9:
            value = alphaCristobalite.deltaH;
            break;
        case 10:
            value = alphaCristobalite.v0;
            break;
        case 11:
            value = alphaCristobalite.v1;
            break;
        case 12:
            value = alphaCristobalite.v2;
            break;
        case 13:
            value = alphaCristobalite.v3;
            break;
        case 14:
            value = alphaCristobalite.v4;
            break;
        case 15:
            value = betaCristobalite.h;
            break;
        case 16:
            value = betaCristobalite.s;
            break;
        case 17:
            value = betaCristobalite.k0;
            break;
        case 18:
            value = betaCristobalite.k1;
            break;
        case 19:
            value = betaCristobalite.k2;
            break;
        case 20:
            value = betaCristobalite.k3;
            break;
        case 21:
            value = betaCristobalite.v0;
            break;
        case 22:
            value = betaCristobalite.v1;
            break;
        case 23:
            value = betaCristobalite.v2;
            break;
        case 24:
            value = betaCristobalite.v3;
            break;
        case 25:
            value = betaCristobalite.v4;
            break;
        case 26:
            value = dTtdp;
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
            double correction = value - alphaCristobalite.h;
            alphaCristobalite.h = value;
            betaCristobalite.h += correction;
            break;
        }
        case 1:
        {
            double correction = value - alphaCristobalite.s;
            alphaCristobalite.s = value;
            betaCristobalite.s += correction;
            break;
        }
        case 2:
        {
            double correction = value - alphaCristobalite.k0;
            alphaCristobalite.k0 = value;
            betaCristobalite.k0 += correction;
            break;
        }
        case 3:
        {
            double correction = value - alphaCristobalite.k1;
            alphaCristobalite.k1 = value;
            betaCristobalite.k1 += correction;
            break;
        }
        case 4:
        {
            double correction = value - alphaCristobalite.k2;
            alphaCristobalite.k2 = value;
            betaCristobalite.k2 += correction;
            break;
        }
        case 5:
        {
            double correction = value - alphaCristobalite.k3;
            alphaCristobalite.k3 = value;
            betaCristobalite.k3 += correction;
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
            double correction = value - alphaCristobalite.v0;
            alphaCristobalite.v0 = value;
            betaCristobalite.v0 += correction;
            break;
        }
        case 11:
        {
            double correction = value - alphaCristobalite.v1;
            alphaCristobalite.v1 = value;
            betaCristobalite.v1 += correction;
            break;
        }
        case 12:
        {
            double correction = value - alphaCristobalite.v2;
            alphaCristobalite.v2 = value;
            betaCristobalite.v2 += correction;
            break;
        }
        case 13:
        {
            double correction = value - alphaCristobalite.v3;
            alphaCristobalite.v3 = value;
            betaCristobalite.v3 += correction;
            break;
        }
        case 14:
        {
            double correction = value - alphaCristobalite.v4;
            alphaCristobalite.v4 = value;
            betaCristobalite.v4 += correction;
            break;
        }
        case 15:
        {
            double correction = value - betaCristobalite.h;
            betaCristobalite.h = value;
            alphaCristobalite.h += correction;
            break;
        }
        case 16:
        {
            double correction = value - betaCristobalite.s;
            betaCristobalite.s = value;
            alphaCristobalite.s += correction;
            break;
        }
        case 17:
        {
            double correction = value - betaCristobalite.k0;
            betaCristobalite.k0 = value;
            alphaCristobalite.k0 += correction;
            break;
        }
        case 18:
        {
            double correction = value - betaCristobalite.k1;
            betaCristobalite.k1 = value;
            alphaCristobalite.k1 += correction;
            break;
        }
        case 19:
        {
            double correction = value - betaCristobalite.k2;
            betaCristobalite.k2 = value;
            alphaCristobalite.k2 += correction;
            break;
        }
        case 20:
        {
            double correction = value - betaCristobalite.k3;
            betaCristobalite.k3 = value;
            alphaCristobalite.k3 += correction;
            break;
        }
        case 21:
        {
            double correction = value - betaCristobalite.v0;
            betaCristobalite.v0 = value;
            alphaCristobalite.v0 += correction;
            break;
        }
        case 22:
        {
            double correction = value - betaCristobalite.v1;
            betaCristobalite.v1 = value;
            alphaCristobalite.v1 += correction;
            break;
        }
        case 23:
        {
            double correction = value - betaCristobalite.v2;
            betaCristobalite.v2 = value;
            alphaCristobalite.v2 += correction;
            break;
        }
        case 24:
        {
            double correction = value - betaCristobalite.v3;
            betaCristobalite.v3 = value;
            alphaCristobalite.v3 += correction;
            break;
        }
        case 25:
        {
            double correction = value - betaCristobalite.v4;
            betaCristobalite.v4 = value;
            alphaCristobalite.v4 += correction;
            break;
        }
        case 26:
        {
            dTtdp = value;
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
