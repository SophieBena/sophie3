//
//  HoltenJPCRD2014.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 2/29/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "PhaseBase.h"
#import "StoichiometricPhaseProtocol.h"

@interface HoltenJPCRD2014 : PhaseBase <StoichiometricPhaseProtocol>

- (double)homogeneousIceNucleationTemperatureForPressureInBars:(double)p;
- (double)homogeneousIceNucleationPressureForTemperatureInK:(double)t;

@end
