//
//  OlivineBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR         5    /* Five independent composition variables */
#define NS         4    /* Four ordering parameters               */
#define NA         6    /* Six endmember compositions             */

@interface OlivineBerman : PhaseBase <SolutionPhaseProtocol> {
	@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS];
	double xm1mg, xm1mn, xm1fe, xm1co, xm1ni, xm2mg, xm2mn, xm2fe, xm2co, xm2ni, xm2ca;
}

#undef NR
#undef NS
#undef NA

@end
