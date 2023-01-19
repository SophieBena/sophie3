//
//  LiquidMeltsCO2+LiquidMeltsCO2_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/16/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "LiquidMeltsCO2+LiquidMeltsCO2_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation LiquidMeltsCO2 (LiquidMeltsCO2_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"enthalpyCorrectionForCO2",
                         @"entropyCorrectionForCO2",
                         @"volumeForCO2",
                         @"dvdtForCO2",
                         @"dvdpForCO2",
                         @"d2vdtdpForCO2",
                         @"d2vdp2ForCO2",
                         nil];
}

-(void)loadParameterUnitsArray {
    unitsOfParameters = [NSArray arrayWithObjects:
                         @"joules",
                         @"joules/K",
                         @"joules/bar",
                         @"joules/bar/K",
                         @"joules/bar/bar",
                         @"joules/bar/bar/K",
                         @"joules/bar/bar/bar",
                         nil];
}

#pragma mark -
#pragma mark Parameter Calibration Protocol methods

- (NSUInteger)getNumberOfFreeParameters {
    if (!namesOfParameters) [self loadParameterNameArray];
    return [namesOfParameters count];
}

- (NSArray *)getArrayOfNamesOfFreeParameters {
    if (!namesOfParameters) [self loadParameterNameArray];
    return namesOfParameters;
}

- (double)getValueForParameterName:(NSString *)name {
    if (!namesOfParameters) [self loadParameterNameArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    double value = PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
    switch (index) {
        case 0:
            value = self.enthalpyCorrectionForCO2;
            break;
        case 1:
            value = self.entropyCorrectionForCO2;
            break;
        case 2:
            value = self.volumeForCO2;
            break;
        case 3:
            value = self.dvdtForCO2;
            break;
        case 4:
            value = self.dvdpForCO2;
            break;
        case 5:
            value = self.d2vdtdpForCO2;
            break;
        case 6:
            value = self.d2vdp2ForCO2;
            break;
        case NSNotFound:
            value = PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
            break;
        default:
            break;
    }
    return value;
}

- (BOOL)setParameterName:(NSString *)name tovalue:(double)value {
    if (!namesOfParameters) [self loadParameterNameArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    switch (index) {
        case 0:
            self.enthalpyCorrectionForCO2 = value;
            break;
        case 1:
            self.entropyCorrectionForCO2 = value;
            break;
        case 2:
            self.volumeForCO2 = value;
            break;
        case 3:
            self.dvdtForCO2 = value;
            break;
        case 4:
            self.dvdpForCO2 = value;
            break;
        case 5:
            self.d2vdtdpForCO2 = value;
            break;
        case 6:
            self.d2vdp2ForCO2 = value;
            break;
        case NSNotFound:
            return NO;
            break;
        default:
            break;
    }
    return YES;
}

- (NSString *)getUnitsForParameterName:(NSString *)name {
    if (!namesOfParameters) [self loadParameterNameArray];
    if (!unitsOfParameters) [self loadParameterUnitsArray];
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return @"";
    else return [unitsOfParameters objectAtIndex:index];
}

- (DoubleVector *)getDgDwFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    return nil;
}

- (DoubleVector *)getChemicalPotentialDerivativesForParameter:(NSString *)name usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
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

- (DoubleMatrix *)getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
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

- (BOOL)supportsParameterCalibration {
    return YES;
}

@end
