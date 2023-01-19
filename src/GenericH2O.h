//
//  GenericH2O.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 3/1/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "PhaseBase.h"
#import "StoichiometricPhaseProtocol.h"

@interface GenericH2O : PhaseBase <StoichiometricPhaseProtocol>

- (void)forceModeChoiceAutomatic;
- (void)forceModeChoiceTo:(NSString *)region;
- (void)setModeToDuanAndZhang2006;
- (void)setModeToZhangAndDuan2005;
- (void)setModeToHoltenEtAl2014;
- (void)setModeToWagnerEtAl2002;

-(void)useZhangAndDuan2005forDEW;
-(void)useDuanAndZhang2006forDEW;
-(void)useZhangAndDuan2009forDEW;

- (NSString *)eosForRegionAtT:(double)t andP:(double)p;
- (double)lowerTemperatureLimitAtPinBars:(double)p;

- (double)tTransHoltenWagner;
- (double)tTransWagnerDZ2006;
- (double)tTransHoltenZD2005;
- (double)pTransWagnerZD2005;
- (double)pTransDZ2006ZD2005;

@end
