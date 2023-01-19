//
//  LiquidpMelts.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/16/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>


#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR 14    /* Eighteen independent composition variables */
#define NA 15    /* Nineteen endmember compositions            */

@interface LiquidpMelts : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
}

#undef NR
#undef NA

@end
