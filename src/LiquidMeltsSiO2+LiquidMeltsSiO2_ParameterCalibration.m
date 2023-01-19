//
//  LiquidMeltsSiO2+LiquidMeltsSiO2_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/16/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "LiquidMeltsSiO2+LiquidMeltsSiO2_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation LiquidMeltsSiO2 (LiquidMeltsSiO2_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"tr",
                         @"pr",
                         @"trl",
                         @"h0_sio2",
                         @"s0_sio2",
                         @"al_sio2",
                         @"bl_sio2",
                         @"cl_sio2",
                         @"dl_sio2",
                         @"tg_sio2",
                         @"cp_sio2",
                         @"vLiq",
                         @"dvdtLiq",
                         @"dvdpLiq",
                         @"d2vdtdpLiq",
                         @"d2vdp2Liq",
                         nil];
}

-(void)loadParameterUnitsArray {
    unitsOfParameters = [NSArray arrayWithObjects:
                         @"K",
                         @"bar",
                         @"K",
                         @"joules",
                         @"joules/K",
                         @"joules/K",
                         @"joules/K/K",
                         @"joules-K",
                         @"joules/K^1/2",
                         @"K",
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
            value = self.tr;
            break;
        case 1:
            value = self.pr;
            break;
        case 2:
            value = self.trl;
            break;
        case 3:
            value = self.h0_sio2;
            break;
        case 4:
            value = self.s0_sio2;
            break;
        case 5:
            value = self.al_sio2;
            break;
        case 6:
            value = self.bl_sio2;
            break;
        case 7:
            value = self.cl_sio2;
            break;
        case 8:
            value = self.dl_sio2;
            break;
        case 9:
            value = self.tg_sio2;
            break;
        case 10:
            value = self.cp_sio2;
            break;
        case 11:
            value = self.vLiq;
            break;
        case 12:
            value = self.dvdtLiq;
            break;
        case 13:
            value = self.dvdpLiq;
            break;
        case 14:
            value = self.d2vdtdpLiq;
            break;
        case 15:
            value = self.d2vdp2Liq;
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
            self.tr = value;
            break;
        case 1:
            self.pr = value;
            break;
        case 2:
            self.trl = value;
            break;
        case 3:
            self.h0_sio2 = value;
            break;
        case 4:
            self.s0_sio2 = value;
            break;
        case 5:
            self.al_sio2 = value;
            break;
        case 6:
            self.bl_sio2 = value;
            break;
        case 7:
            self.cl_sio2 = value;
            break;
        case 8:
            self.dl_sio2 = value;
            break;
        case 9:
            self.tg_sio2 = value;
            break;
        case 10:
            self.cp_sio2 = value;
            break;
        case 11:
            self.vLiq = value;
            break;
        case 12:
            self.dvdtLiq = value;
            break;
        case 13:
            self.dvdpLiq = value;
            break;
        case 14:
            self.d2vdtdpLiq = value;
            break;
        case 15:
            self.d2vdp2Liq = value;
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
