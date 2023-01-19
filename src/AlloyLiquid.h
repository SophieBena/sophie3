//
//  AlloyLiquid.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/11/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         1    // Independent composition variables
#define NA         2    // Endmember components

@interface AlloyLiquid : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
}

#undef NR
#undef NA

@end
