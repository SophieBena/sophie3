//
//  HornblendeBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/10/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         2    // Independent composition variables
#define NA         3    // Endmember components

@interface HornblendeBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
}

#undef NR
#undef NA

@end
