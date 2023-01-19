//
//  QuartzBerman+QuartzBerman_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/16/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "QuartzBerman+QuartzBerman_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation QuartzBerman (QuartzBerman_ParameterCalibration)

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
            value = alphaQuartz.h;
            break;
        case 1:
            value = alphaQuartz.s;
            break;
        case 2:
            value = alphaQuartz.k0;
            break;
        case 3:
            value = alphaQuartz.k1;
            break;
        case 4:
            value = alphaQuartz.k2;
            break;
        case 5:
            value = alphaQuartz.k3;
            break;
        case 6:
            value = alphaQuartz.l1;
            break;
        case 7:
            value = alphaQuartz.l2;
            break;
        case 8:
            value = alphaQuartz.Tt;
            break;
        case 9:
            value = alphaQuartz.deltaH;
            break;
        case 10:
            value = alphaQuartz.v0;
            break;
        case 11:
            value = alphaQuartz.v1;
            break;
        case 12:
            value = alphaQuartz.v2;
            break;
        case 13:
            value = alphaQuartz.v3;
            break;
        case 14:
            value = alphaQuartz.v4;
            break;
        case 15:
            value = betaQuartz.h;
            break;
        case 16:
            value = betaQuartz.s;
            break;
        case 17:
            value = betaQuartz.k0;
            break;
        case 18:
            value = betaQuartz.k1;
            break;
        case 19:
            value = betaQuartz.k2;
            break;
        case 20:
            value = betaQuartz.k3;
            break;
        case 21:
            value = betaQuartz.v0;
            break;
        case 22:
            value = betaQuartz.v1;
            break;
        case 23:
            value = betaQuartz.v2;
            break;
        case 24:
            value = betaQuartz.v3;
            break;
        case 25:
            value = betaQuartz.v4;
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
            double correction = value - alphaQuartz.h;
            alphaQuartz.h = value;
            betaQuartz.h += correction;
            break;
        }
        case 1:
        {
            double correction = value - alphaQuartz.s;
            alphaQuartz.s = value;
            betaQuartz.s += correction;
            break;
        }
        case 2:
        {
            double correction = value - alphaQuartz.k0;
            alphaQuartz.k0 = value;
            betaQuartz.k0 += correction;
            break;
        }
        case 3:
        {
            double correction = value - alphaQuartz.k1;
            alphaQuartz.k1 = value;
            betaQuartz.k1 += correction;
            break;
        }
        case 4:
        {
            double correction = value - alphaQuartz.k2;
            alphaQuartz.k2 = value;
            betaQuartz.k2 += correction;
            break;
        }
        case 5:
        {
            double correction = value - alphaQuartz.k3;
            alphaQuartz.k3 = value;
            betaQuartz.k3 += correction;
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
            double correction = value - alphaQuartz.v0;
            alphaQuartz.v0 = value;
            betaQuartz.v0 += correction;
            break;
        }
        case 11:
        {
            double correction = value - alphaQuartz.v1;
            alphaQuartz.v1 = value;
            betaQuartz.v1 += correction;
            break;
        }
        case 12:
        {
            double correction = value - alphaQuartz.v2;
            alphaQuartz.v2 = value;
            betaQuartz.v2 += correction;
            break;
        }
        case 13:
        {
            double correction = value - alphaQuartz.v3;
            alphaQuartz.v3 = value;
            betaQuartz.v3 += correction;
            break;
        }
        case 14:
        {
            double correction = value - alphaQuartz.v4;
            alphaQuartz.v4 = value;
            betaQuartz.v4 += correction;
            break;
        }
        case 15:
        {
            double correction = value - betaQuartz.h;
            betaQuartz.h = value;
            alphaQuartz.h += correction;
            break;
        }
        case 16:
        {
            double correction = value - betaQuartz.s;
            betaQuartz.s = value;
            alphaQuartz.s += correction;
            break;
        }
        case 17:
        {
            double correction = value - betaQuartz.k0;
            betaQuartz.k0 = value;
            alphaQuartz.k0 += correction;
            break;
        }
        case 18:
        {
            double correction = value - betaQuartz.k1;
            betaQuartz.k1 = value;
            alphaQuartz.k1 += correction;
            break;
        }
        case 19:
        {
            double correction = value - betaQuartz.k2;
            betaQuartz.k2 = value;
            alphaQuartz.k2 += correction;
            break;
        }
        case 20:
        {
            double correction = value - betaQuartz.k3;
            betaQuartz.k3 = value;
            alphaQuartz.k3 += correction;
            break;
        }
        case 21:
        {
            double correction = value - betaQuartz.v0;
            betaQuartz.v0 = value;
            alphaQuartz.v0 += correction;
            break;
        }
        case 22:
        {
            double correction = value - betaQuartz.v1;
            betaQuartz.v1 = value;
            alphaQuartz.v1 += correction;
            break;
        }
        case 23:
        {
            double correction = value - betaQuartz.v2;
            betaQuartz.v2 = value;
            alphaQuartz.v2 += correction;
            break;
        }
        case 24:
        {
            double correction = value - betaQuartz.v3;
            betaQuartz.v3 = value;
            alphaQuartz.v3 += correction;
            break;
        }
        case 25:
        {
            double correction = value - betaQuartz.v4;
            betaQuartz.v4 = value;
            alphaQuartz.v4 += correction;
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
