//
//  LiquidMeltsPlusOldH2OandNewCO2.h
//  PhaseObjC
//
//  Created by Mark Ghiorso on 7/6/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NR 15    /* Fifteen independent composition variables */
#define NS  1    /* Three ordering parameters                 */
#define NA 16    /* Sixteen endmember compositions            */

@interface LiquidMeltsPlusOldH2OandNewCO2 : PhaseBase <SolutionPhaseProtocol> {
@private
    BOOL computeMixingQuantities;
    double tOld, pOld, rOld[NR], sOld[NS], invd2gds2[NS][NS];
    double deltaHsp[NS], deltaSsp[NS], deltaVsp[NS];
    BOOL modelParametersHaveBeenAlteredSinceLastOrderingDetermination;
}

@property (readonly) DoubleVector *modelParametersWrapper;
@property (readonly) NSArray *endmembers;

-(NSDictionary *)getSpeciesMoleFractionsForBulkComposition:(double *)m aT:(double)t andP:(double)p;

#undef NR
#undef NS
#undef NA

@end
