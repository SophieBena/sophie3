//
//  LiquidMelts+LiquidMelts_LiquidMelts_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 1/19/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "LiquidMelts+LiquidMelts_LiquidMelts_ParameterCalibration.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

@implementation LiquidMelts (LiquidMelts_LiquidMelts_ParameterCalibration)

static NSArray *namesOfParameters;
static NSArray *unitsOfParameters;

#pragma mark -
#pragma mark private methods

-(void)loadParameterNameArray {
    namesOfParameters = [NSArray arrayWithObjects:
                         @"W(TiO2      ,SiO2)",
                         @"W(Al2O3     ,SiO2)",
                         @"W(Fe2O3     ,SiO2)",
                         @"W(MgCr2O4   ,SiO2)",
                         @"W(Fe2SiO4   ,SiO2)",
                         @"W(MnSi0.5O2 ,SiO2)",
                         @"W(Mg2SiO4   ,SiO2)",
                         @"W(NiSi0.5O2 ,SiO2)",
                         @"W(CoSi0.5O2 ,SiO2)",
                         @"W(CaSiO3    ,SiO2)",
                         @"W(Na2SiO3   ,SiO2)",
                         @"W(KAlSiO4   ,SiO2)",
                         @"W(Ca3(PO4)2 ,SiO2)",
                         @"W(H2O       ,SiO2)",

                         @"W(Al2O3     ,TiO2)",
                         @"W(Fe2O3     ,TiO2)",
                         @"W(MgCr2O4   ,TiO2)",
                         @"W(Fe2SiO4   ,TiO2)",
                         @"W(MnSi0.5O2 ,TiO2)",
                         @"W(Mg2SiO4   ,TiO2)",
                         @"W(NiSi0.5O2 ,TiO2)",
                         @"W(CoSi0.5O2 ,TiO2)",
                         @"W(CaSiO3    ,TiO2)",
                         @"W(Na2SiO3   ,TiO2)",
                         @"W(KAlSiO4   ,TiO2)",
                         @"W(Ca3(PO4)2 ,TiO2)",
                         @"W(H2O       ,TiO2)",

                         @"W(Fe2O3     ,Al2O3)",
                         @"W(MgCr2O4   ,Al2O3)",
                         @"W(Fe2SiO4   ,Al2O3)",
                         @"W(MnSi0.5O2 ,Al2O3)",
                         @"W(Mg2SiO4   ,Al2O3)",
                         @"W(NiSi0.5O2 ,Al2O3)",
                         @"W(CoSi0.5O2 ,Al2O3)",
                         @"W(CaSiO3    ,Al2O3)",
                         @"W(Na2SiO3   ,Al2O3)",
                         @"W(KAlSiO4   ,Al2O3)",
                         @"W(Ca3(PO4)2 ,Al2O3)",
                         @"W(H2O       ,Al2O3)",

                         @"W(MgCr2O4   ,Fe2O3)",
                         @"W(Fe2SiO4   ,Fe2O3)",
                         @"W(MnSi0.5O2 ,Fe2O3)",
                         @"W(Mg2SiO4   ,Fe2O3)",
                         @"W(NiSi0.5O2 ,Fe2O3)",
                         @"W(CoSi0.5O2 ,Fe2O3)",
                         @"W(CaSiO3    ,Fe2O3)",
                         @"W(Na2SiO3   ,Fe2O3)",
                         @"W(KAlSiO4   ,Fe2O3)",
                         @"W(Ca3(PO4)2 ,Fe2O3)",
                         @"W(H2O       ,Fe2O3)",

                         @"W(Fe2SiO4   ,MgCr2O4)",
                         @"W(MnSi0.5O2 ,MgCr2O4)",
                         @"W(Mg2SiO4   ,MgCr2O4)",
                         @"W(NiSi0.5O2 ,MgCr2O4)",
                         @"W(CoSi0.5O2 ,MgCr2O4)",
                         @"W(CaSiO3    ,MgCr2O4)",
                         @"W(Na2SiO3   ,MgCr2O4)",
                         @"W(KAlSiO4   ,MgCr2O4)",
                         @"W(Ca3(PO4)2 ,MgCr2O4)",
                         @"W(H2O       ,MgCr2O4)",

                         @"W(MnSi0.5O2 ,Fe2SiO4)",
                         @"W(Mg2SiO4   ,Fe2SiO4)",
                         @"W(NiSi0.5O2 ,Fe2SiO4)",
                         @"W(CoSi0.5O2 ,Fe2SiO4)",
                         @"W(CaSiO3    ,Fe2SiO4)",
                         @"W(Na2SiO3   ,Fe2SiO4)",
                         @"W(KAlSiO4   ,Fe2SiO4)",
                         @"W(Ca3(PO4)2 ,Fe2SiO4)",
                         @"W(H2O       ,Fe2SiO4)",

                         @"W(Mg2SiO4   ,MnSi0.5O2)",
                         @"W(NiSi0.5O2 ,MnSi0.5O2)",
                         @"W(CoSi0.5O2 ,MnSi0.5O2)",
                         @"W(CaSiO3    ,MnSi0.5O2)",
                         @"W(Na2SiO3   ,MnSi0.5O2)",
                         @"W(KAlSiO4   ,MnSi0.5O2)",
                         @"W(Ca3(PO4)2 ,MnSi0.5O2)",
                         @"W(H2O       ,MnSi0.5O2)",

                         @"W(NiSi0.5O2 ,Mg2SiO4)",
                         @"W(CoSi0.5O2 ,Mg2SiO4)",
                         @"W(CaSiO3    ,Mg2SiO4)",
                         @"W(Na2SiO3   ,Mg2SiO4)",
                         @"W(KAlSiO4   ,Mg2SiO4)",
                         @"W(Ca3(PO4)2 ,Mg2SiO4)",
                         @"W(H2O       ,Mg2SiO4)",

                         @"W(CoSi0.5O2 ,NiSi0.5O2)",
                         @"W(CaSiO3    ,NiSi0.5O2)",
                         @"W(Na2SiO3   ,NiSi0.5O2)",
                         @"W(KAlSiO4   ,NiSi0.5O2)",
                         @"W(Ca3(PO4)2 ,NiSi0.5O2)",
                         @"W(H2O       ,NiSi0.5O2)",

                         @"W(CaSiO3    ,CoSi0.5O2)",
                         @"W(Na2SiO3   ,CoSi0.5O2)",
                         @"W(KAlSiO4   ,CoSi0.5O2)",
                         @"W(Ca3(PO4)2 ,CoSi0.5O2)",
                         @"W(H2O       ,CoSi0.5O2)",

                         @"W(Na2SiO3   ,CaSiO3)",
                         @"W(KAlSiO4   ,CaSiO3)",
                         @"W(Ca3(PO4)2 ,CaSiO3)",
                         @"W(H2O       ,CaSiO3)",

                         @"W(KAlSiO4   ,Na2SiO3)",
                         @"W(Ca3(PO4)2 ,Na2SiO3)",
                         @"W(H2O       ,Na2SiO3)",

                         @"W(Ca3(PO4)2 ,KAlSiO4)",
                         @"W(H2O       ,KAlSiO4)",

                         @"W(H2O       ,Ca3(PO4)2)",

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

-(double)getValueForParameterName:(NSString *)name {
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
    else {
        double *modelParameters = [[self modelParametersWrapper] pointerToDouble];
        return modelParameters[index];
    }
}

-(BOOL)setParameterName:(NSString *)name tovalue:(double)value {
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return NO;
    else {
        double *modelParameters = [[self modelParametersWrapper] pointerToDouble];
        modelParameters[index] = value;
        return YES;
    }
}

-(NSString *)getUnitsForParameterName:(NSString *)name {
    return @"joules";
}

-(DoubleVector *)getDgDwFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    return nil;
}

-(DoubleVector *)getChemicalPotentialDerivativesForParameter:(NSString *)name usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    NSUInteger index = [namesOfParameters indexOfObject:name];
    if (index == NSNotFound) return nil;
    else {
        NSUInteger na = [self numberOfSolutionComponents];
        double *modelParameters = [[self modelParametersWrapper] pointerToDouble];
        DoubleVector *chemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];

        double value = modelParameters[index];
        double forwardValue = (value != 0.0) ? value*(1.0+sqrt(DBL_EPSILON)) : sqrt(DBL_EPSILON);
        modelParameters[index] = forwardValue;
        DoubleVector *forwardChemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];

        DoubleVector *derivativeWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
        double *derivative = [derivativeWrapper pointerToDouble];
        double *chemicalPotentials = [chemicalPotentialsWrapper pointerToDouble];
        double *forwardChemicalPotentials = [forwardChemicalPotentialsWrapper pointerToDouble];
        for (NSUInteger i=0; i<na; i++) derivative[i] = (forwardChemicalPotentials[i] - chemicalPotentials[i])/(forwardValue - value);

        modelParameters[index] = value;

        return derivativeWrapper;
    }
}

-(DoubleMatrix *)getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    if (!array || ([array count] == 0)) return nil;
    NSUInteger lengthOfArray = [array count];
    NSUInteger na = [self numberOfSolutionComponents];
    DoubleMatrix *derivativeMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:lengthOfArray andWithColumnSize:na andInitialValue:0.0];
    double **derivativeMatrix = [derivativeMatrixWrapper pointerToPointerToDouble];

    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];
    DoubleVector *chemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
    double *chemicalPotentials = [chemicalPotentialsWrapper pointerToDouble];

    for (NSUInteger j=0; j<lengthOfArray; j++) {
        NSUInteger index = [namesOfParameters indexOfObject:[array objectAtIndex:j]];
        if (index != NSNotFound) {
            double value = modelParameters[index];
            double forwardValue = (value != 0.0) ? value*(1.0+sqrt(DBL_EPSILON)) : sqrt(DBL_EPSILON);
            modelParameters[index] = forwardValue;
            DoubleVector *forwardChemicalPotentialsWrapper = [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
            double *forwardChemicalPotentials = [forwardChemicalPotentialsWrapper pointerToDouble];

            for (NSUInteger i=0; i<na; i++) derivativeMatrix[j][i] = (forwardChemicalPotentials[i] - chemicalPotentials[i])/(forwardValue - value);
            modelParameters[index] = value;
        }
    }

    return derivativeMatrixWrapper;
}

-(BOOL)supportsParameterCalibration {
    return YES;
}

@end
