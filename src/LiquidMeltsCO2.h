//
//  LiquidMeltsCO2.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/26/12.
//  Copyright 2012 OFM Research Inc.. All rights reserved.
//

//#import <Foundation/Foundation.h>
#import "PhaseBase.h"
#import "StoichiometricPhaseProtocol.h"
#import "FluidDuan.h"

@interface LiquidMeltsCO2 : PhaseBase <StoichiometricPhaseProtocol> {
    DuanCO2 *duanCO2;
}

@property double enthalpyCorrectionForCO2;
@property double entropyCorrectionForCO2;
@property double volumeForCO2;
@property double dvdtForCO2;
@property double dvdpForCO2;
@property double d2vdtdpForCO2;
@property double d2vdp2ForCO2;

@property (atomic, assign, getter = isGibbsFreeEnergyReferenceStateUsed) BOOL gibbsFreeEnergyReferenceStateUsed;

@end
