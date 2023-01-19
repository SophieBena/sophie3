//
//  RhombohedralBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/4/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         4    /* Five independent composition variables */
#define NS         3    /* Four ordering parameters               */
#define NA         5    /* Six endmember compositions             */

@interface RhombohedralBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS], tOldPure, pOldPure, sOldPure[NS], d2gds2Pure[NS];
	double xmg2a, xfe2a, xmn2a, xti4a, xfe3a, xal3a;
	double xmg2b, xfe2b, xmn2b, xti4b, xfe3b, xal3b;
	double xmg2ID, xfe2ID, xmn2ID, xti4ID, xfe3ID, xal3ID;
}

#undef NR
#undef NS
#undef NA

@end
