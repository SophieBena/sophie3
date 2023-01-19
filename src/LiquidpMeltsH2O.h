//
//  LiquidpMeltsH2O.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/16/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "LiquidpMeltsGenericEM.h"


@interface LiquidpMeltsH2O : PhaseBase <StoichiometricPhaseProtocol> {
	LiquidpMeltsGenericEM *correction;
}

@end

@interface WaterpMelts : PhaseBase <StoichiometricPhaseProtocol> {
}

@end
