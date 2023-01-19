//
//  OrthoOxideBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/8/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"
#import "BermanProperties.h"

#define NR         2    /* Independent composition variables */
#define NS         3    /* Ordering parameters               */
#define NA         3    /* Endmember compositions            */

@interface OrthoOxideBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS], tOldPure, pOldPure, sOldPure[NS], d2gds2Pure[NS];
	double xmg2M1, xfe2M1, xti4M1, xfe3M1, xmg2M2, xfe2M2, xti4M2, xfe3M2;
}

#undef NR
#undef NS
#undef NA

@end
