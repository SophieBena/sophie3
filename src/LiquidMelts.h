//
//  LiquidMelts.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/18/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>


#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

@class DoubleVector;

@interface LiquidMelts : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
}

@property (readonly) DoubleVector *modelParametersWrapper;
@property (readonly) NSArray *endmembers;

@end
