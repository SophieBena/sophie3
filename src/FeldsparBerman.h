//
//  FeldsparBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/30/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         2    // Two independent composition variables
#define NA         3    // Three endmember components

@interface FeldsparBerman : PhaseBase <SolutionPhaseProtocol> {
	NSString *operationParent;

@private
	BOOL computeMixingQuantities;
    BOOL isMELTS;
}

@property (readwrite,copy) NSString *operationParent;

@property (readwrite, nonatomic) double whabor;
@property (readwrite, nonatomic) double wsabor;
@property (readwrite, nonatomic) double wvabor;
@property (readwrite, nonatomic) double whorab;
@property (readwrite, nonatomic) double wsorab;
@property (readwrite, nonatomic) double wvorab;
@property (readwrite, nonatomic) double whaban;
@property (readwrite, nonatomic) double whanab;
@property (readwrite, nonatomic) double whoran;
@property (readwrite, nonatomic) double whanor;
@property (readwrite, nonatomic) double wvanor;
@property (readwrite, nonatomic) double whabanor;
@property (readwrite, nonatomic) double wvabanor;

@property (readonly) NSArray *endmembers;

-(id)initWithCompositionConstraint:(NSString *)name;
-(void)incrementInstanceCountOfPhase;
-(void)decrementInstanceCountOfPhase;
-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p;

-(NSString *)nameOfPhaseWithComposition:(double *)refMoles;

#undef NR
#undef NA

@end
