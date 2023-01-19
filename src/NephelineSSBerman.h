//
//  NephelineBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/7/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"
#import "BermanProperties.h"

#define NR         3    /* Independent composition variables */
#define NS         1    /* Ordering parameters               */
#define NA         4    /* Endmember compositions            */

@interface NephelineSSBerman : PhaseBase <SolutionPhaseProtocol> {
	NSString *operationParent;

@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS];
	double xkls, xvcls, xnals, xkss, xcass, xnass;
}

@property (readwrite,copy) NSString *operationParent;

-(id)initWithCompositionConstraint:(NSString *)name;
-(void)incrementInstanceCountOfPhase;
-(void)decrementInstanceCountOfPhase;
-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p;

-(NSString *)nameOfPhaseWithComposition:(double *)refMoles;

#undef NR
#undef NS
#undef NA

@end

@interface VacancyNepheline : BermanProperties  {
@private
	BermanProperties *highAlbite;
	BermanProperties *betaNepheline;
}

@end

@interface CalciumNepheline : BermanProperties  {
@private
	BermanProperties *anorthite;
	BermanProperties *betaNepheline;
}

@end
