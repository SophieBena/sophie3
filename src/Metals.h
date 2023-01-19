//
//  Metals.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/11/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "PhaseBase.h"
#import "StoichiometricPhaseProtocol.h"

@interface FeSolid : PhaseBase <StoichiometricPhaseProtocol> {
@private
	double tr, pr, k0, k1, k2, k3, ct, v1, v2, v3, v4, h0, s0, v0;
}
@end

@interface FeLiquid : PhaseBase <StoichiometricPhaseProtocol> {
@private
	double tr, tr1, pr, h0, s0, v0, v1, v2, v3, v4, cp;
}
@end

@interface NiSolid : PhaseBase <StoichiometricPhaseProtocol> {
@private
	double tr, tr1, pr, h0, s0, v0, vni1, vni3a, vni3b, cpa, cpb, cpc;
}
@end

@interface NiLiquid : PhaseBase <StoichiometricPhaseProtocol> {
@private
	double tr, tr1, pr, h0, s0, v0, v1, v2, v3, v4, cp;
}
@end
