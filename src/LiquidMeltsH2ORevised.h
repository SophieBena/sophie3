//
//  LiquidMeltsH2ORevised.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 9/12/13.
//
//

#import <Foundation/Foundation.h>
#import "PhaseBase.h"
#import "StoichiometricPhaseProtocol.h"

@interface LiquidMeltsH2ORevised : PhaseBase <StoichiometricPhaseProtocol> {

}

@property double enthalpyCorrectionForH2O;
@property double entropyCorrectionForH2O;
@property double volumeForH2O;
@property double dvdtForH2O;
@property double dvdpForH2O;

@property (atomic, assign, getter = isGibbsFreeEnergyReferenceStateUsed) BOOL gibbsFreeEnergyReferenceStateUsed;

@end
