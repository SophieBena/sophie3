//
//  MeliliteBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/7/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"
#import "BermanProperties.h"

#define NR         3    /* Independent composition variables */
#define NS         1    /* Ordering parameters               */
#define NA         4    /* Endmember compositions            */

@interface MeliliteBerman : PhaseBase <SolutionPhaseProtocol> {
@private
	BOOL computeMixingQuantities;
	double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS], tOldPure, pOldPure, sOldPure[NS], d2gds2Pure[NS];
	double xCaOct, xNaOct, xMg2T1, xFe2T1, xAl3T1, xSi4T1, xAl3T2, xSi4T2;
}

#undef NR
#undef NS
#undef NA

@end

@interface Gehlenite : BermanProperties {
@private
	double d0, d1, d2, d3, d4, d5;
}
@end


@interface IronAkermanite : BermanProperties  {
@private
	BermanProperties *akermanite;
	BermanProperties *fayalite;
	BermanProperties *forsterite;
}

@end

@interface SodaMelilite : BermanProperties  {
@private
	BermanProperties *gehlenite;
	BermanProperties *monalbite;
	BermanProperties *anorthite;
}

@end
