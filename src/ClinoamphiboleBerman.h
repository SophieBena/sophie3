//
//  ClinoamphiboleBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/8/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         2    /* Two independent composition variables */
#define NS         2    /* Two ordering parameters               */
#define NA         3    /* Three endmember compositions          */

@interface ClinoamphiboleBerman : PhaseBase <SolutionPhaseProtocol> {
	NSString *operationParent;

@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS];
	double xfe2m13, xmg2m13, xfe2m2, xmg2m2, xfe2m4, xmg2m4, xca2m4;
}

@property (readwrite,copy) NSString *operationParent;

-(id)initWithCompositionConstraint:(NSString *)name;
-(void)incrementInstanceCountOfPhase;
-(void)decrementInstanceCountOfPhase;
-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p;

#undef NR
#undef NS
#undef NA

@end
