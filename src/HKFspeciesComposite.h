//
//  HKFspeciesComposite.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 1/19/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"

@class HKFspeciesProperties;
@class DEWspecies;

@interface HKFspeciesComposite : PhaseBase <StoichiometricPhaseProtocol, NSSecureCoding>

@property (weak)DEWspecies *dewSpecies;
@property (strong)NSMutableArray <NSString *> *hkfSpecies;
@property (strong)NSMutableArray <NSNumber *> *hkfStoichiometry;
@property (assign)NSInteger count;
@property (assign)double correctionToGibbsEnergy;

- (instancetype)initWithDEWspeciesInstance:(DEWspecies *)dewSpecies
                 andWithReactionDictionary:(NSDictionary <NSString *, NSNumber *> *)nameKeysAndStoichiometryValues;

- (double)gSolventAtT:(double)t andP:(double)p;
- (double)DgSolventDtAtT:(double)t andP:(double)p;
- (double)DgSolventDpAtT:(double)t andP:(double)p;
- (double)D2gSolventDt2AtT:(double)t andP:(double)p;
- (double)D2gSolventDtDpAtT:(double)t andP:(double)p;
- (double)D2gSolventDp2AtT:(double)t andP:(double)p;

- (double)omegaAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z;
- (double)DomegaDtAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z;
- (double)DomegaDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z;
- (double)D2omegaDt2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z;
- (double)D2omegaDtDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z;
- (double)D2omegaDp2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z;

@end
