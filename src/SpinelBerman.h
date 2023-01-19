//
//  SpinelBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 8/10/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         4    /* Five independent composition variables */
#define NS         3    /* Four ordering parameters               */
#define NA         5    /* Six endmember compositions             */

@interface SpinelBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS], tOldPure, pOldPure, sOldPure[NS], d2gds2Pure[NS];
	double xmg2tet, xfe2tet, xal3tet, xfe3tet, xmg2oct, xfe2oct, xal3oct, xfe3oct, xcr3oct, xti4oct;
}

#undef NR
#undef NS
#undef NA

@end
