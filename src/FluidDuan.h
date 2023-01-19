//
//  FluidDuan.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/25/12.
//  Copyright 2012 OFM Research Inc.. All rights reserved.
//

//#import <Foundation/Foundation.h>


#import "SolutionPhaseProtocol.h"
#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"

#define NR 1    /* One independent composition variables */
#define NA 2    /* Two endmember compositions            */

@interface FluidDuan : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
}

#undef NR
#undef NA

@end

@interface DuanH2O : PhaseBase <StoichiometricPhaseProtocol> {
}

@end

@interface DuanCO2 : PhaseBase <StoichiometricPhaseProtocol> {
}

@end
