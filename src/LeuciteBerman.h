//
//  LeuciteBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/5/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         2    // Two independent composition variables
#define NA         3    // Three endmember components

@interface LeuciteBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
}

#undef NR
#undef NA

@end
