//
//  PseudoPhase.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 8/20/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"

@interface PseudoPhase : PhaseBase <StoichiometricPhaseProtocol> {
    double g, h, s, cp, dcpdt, v, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2;
    double gLast;
}

@property double g, h, s, cp, dcpdt, v, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2;
@property double gLast;

@end
