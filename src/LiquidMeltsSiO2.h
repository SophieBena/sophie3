//
//  LiquidMeltsSiO2.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/18/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "PhaseBase.h"
#import "StoichiometricPhaseProtocol.h"

@interface LiquidMeltsSiO2 : PhaseBase <StoichiometricPhaseProtocol>

@property (readwrite) double tr;
@property (readwrite) double pr;
@property (readwrite) double trl;
@property (readwrite) double h0_sio2;
@property (readwrite) double s0_sio2;
@property (readwrite) double al_sio2;
@property (readwrite) double bl_sio2;
@property (readwrite) double cl_sio2;
@property (readwrite) double dl_sio2;
@property (readwrite) double tg_sio2;
@property (readwrite) double cp_sio2;
@property (readwrite) double vLiq;
@property (readwrite) double dvdtLiq;
@property (readwrite) double dvdpLiq;
@property (readwrite) double d2vdtdpLiq;
@property (readwrite) double d2vdp2Liq;

@property (atomic, assign, getter = isGibbsFreeEnergyReferenceStateUsed) BOOL gibbsFreeEnergyReferenceStateUsed;

@end
