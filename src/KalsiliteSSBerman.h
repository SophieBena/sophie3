//
//  KalsiliteSSBerman.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/10/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"
#import "BermanProperties.h"

#define NR         3    /* Independent composition variables */
#define NS         0    /* Ordering parameters               */
#define NA         4    /* Endmember compositions            */

@interface KalsiliteSSBerman : PhaseBase <SolutionPhaseProtocol> {
	@private
	BOOL computeMixingQuantities;
}

#undef NR
#undef NS
#undef NA

@end

@interface VacancyKalsilite : BermanProperties  {
@private
	BermanProperties *highAlbite;
	BermanProperties *betaNepheline;
}

@end

@interface CalciumKalsilite : BermanProperties  {
@private
	BermanProperties *anorthite;
	BermanProperties *betaNepheline;
}

@end
