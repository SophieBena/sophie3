//
//  ParameterCalibrationProtocol.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 1/2/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>

@class DoubleVector;
@class DoubleMatrix;
@class DoubleTensor;

static double PARAMETER_CALIBRATION_VALUE_NOT_FOUND = 1234567890.0987654321;

@protocol ParameterCalibrationProtocol <NSObject>

-(NSUInteger)getNumberOfFreeParameters;
-(NSArray *)getArrayOfNamesOfFreeParameters;
-(double)getValueForParameterName:(NSString *)name;
-(BOOL)setParameterName:(NSString *)name tovalue:(double)value;
-(NSString *)getUnitsForParameterName:(NSString *)name;

-(DoubleVector *)getDgDwFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
-(DoubleVector *)getChemicalPotentialDerivativesForParameter:(NSString *)name usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
-(DoubleMatrix *)getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p;

-(BOOL)supportsParameterCalibration;

@end
