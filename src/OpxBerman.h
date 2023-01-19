//
//  OpxBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 8/11/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         6    /* Five independent composition variables */
#define NS         2    /* Four ordering parameters               */
#define NA         7    /* Six endmember compositions             */

@interface OpxBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS], tOldPure, pOldPure, sOldPure, d2gds2Pure;
	double xal3m1, xfe2m1, xfe3m1, xmg2m1, xti4m1, xca2m2, xfe2m2, xmg2m2, xna1m2, xal3tet, xfe3tet, xsi4tet;
}

#undef NR
#undef NS
#undef NA


@end
