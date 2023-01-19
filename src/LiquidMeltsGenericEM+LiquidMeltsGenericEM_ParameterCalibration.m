//
//  LiquidMeltsGenericEM+LiquidMeltsGenericEM_ParameterCalibration.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/16/15.
//  Copyright (c) 2015 Mark Ghiorso. All rights reserved.
//

#import "LiquidMeltsGenericEM+LiquidMeltsGenericEM_ParameterCalibration.h"
#import "BermanProperties+BermanProperties_BermanProperties_ParameterCalibration.h"

@implementation LiquidMeltsGenericEM (LiquidMeltsGenericEM_ParameterCalibration)

static NSArray *namesOfAdditionalParameters;
static NSArray *unitsOfAdditionalParameters;

#pragma mark -
#pragma mark private methods

-(void)loadAdditionalParameterNameArray {
    namesOfAdditionalParameters = [NSArray arrayWithObjects:
                         @"vLiq",
                         @"dvdtLiq",
                         @"dvdpLiq",
                         @"d2vdtdpLiq",
                         @"d2vdp2Liq",
                         @"tFusion",
                         @"sFusion",
                         @"cpLiq",
                         nil];
}

-(void)loadAdditionalParameterUnitsArray {
    unitsOfAdditionalParameters = [NSArray arrayWithObjects:
                         @"joules/bar",
                         @"joules/bar/K",
                         @"joules/bar/bar",
                         @"joules/bar/bar/K",
                         @"joules/bar/bar/bar",
                         @"K",
                         @"joules/K",
                         @"joules/K",
                         nil];
}

#pragma mark -
#pragma mark Parameter Calibration Protocol methods

//
// Additional required protocol methods are inherited from BermanProperties category extension.  These methods call the
// methods overridden below:

- (NSUInteger)getNumberOfFreeParameters {
    if (!namesOfAdditionalParameters) [self loadAdditionalParameterNameArray];
    return [super getNumberOfFreeParameters] + namesOfAdditionalParameters.count;
}

- (NSArray *)getArrayOfNamesOfFreeParameters {
    if (!namesOfAdditionalParameters) [self loadAdditionalParameterNameArray];
    return [[super getArrayOfNamesOfFreeParameters] arrayByAddingObjectsFromArray:namesOfAdditionalParameters];
}

-(double)getValueForParameterName:(NSString *)name {
    if (!namesOfAdditionalParameters) [self loadAdditionalParameterNameArray];
    NSUInteger index = [namesOfAdditionalParameters indexOfObject:name];
    if (index == NSNotFound) return [super getValueForParameterName:name];
    double value = PARAMETER_CALIBRATION_VALUE_NOT_FOUND;
    switch (index) {
        case 0:
            value = self.vLiq;
            break;
        case 1:
            value = self.dvdtLiq;
            break;
        case 2:
            value = self.dvdpLiq;
            break;
        case 3:
            value = self.d2vdtdpLiq;
            break;
        case 4:
            value = self.d2vdp2Liq;
            break;
        case 5:
            value = self.tFusion;
            break;
        case 6:
            value = self.sFusion;
            break;
        case 7:
            value = self.cpLiq;
            break;
        default:
            break;
    }
    return value;
}

-(BOOL)setParameterName:(NSString *)name tovalue:(double)value {
    if (!namesOfAdditionalParameters) [self loadAdditionalParameterNameArray];
    NSUInteger index = [namesOfAdditionalParameters indexOfObject:name];
    if (index == NSNotFound) return [super setParameterName:name tovalue:value];
    switch (index) {
        case 0:
            self.vLiq = value;
            break;
        case 1:
            self.dvdtLiq = value;
            break;
        case 2:
            self.dvdpLiq = value;
            break;
        case 3:
            self.d2vdtdpLiq = value;
            break;
        case 4:
            self.d2vdp2Liq = value;
            break;
        case 5:
            self.tFusion = value;
            break;
        case 6:
            self.sFusion = value;
            break;
        case 7:
            self.cpLiq = value;
            break;
        default:
            break;
    }
    return YES;
}

-(NSString *)getUnitsForParameterName:(NSString *)name {
    if (!namesOfAdditionalParameters) [self loadAdditionalParameterNameArray];
    if (!unitsOfAdditionalParameters) [self loadAdditionalParameterUnitsArray];
    NSUInteger index = [namesOfAdditionalParameters indexOfObject:name];
    if (index == NSNotFound) return [super getUnitsForParameterName:name];
    else return [unitsOfAdditionalParameters objectAtIndex:index];
}

@end
