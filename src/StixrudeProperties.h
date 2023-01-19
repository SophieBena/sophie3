//
//  StixrudeProperties.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/10/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"


/**
 Computes thermodynamic properties using a Berman-type database model.
 */

@interface StixrudeProperties : PhaseBase <StoichiometricPhaseProtocol> {
@protected
	NSArray *parameters;
	double currentT, currentP;
	NSDictionary *currentResults;
	double TC0, VD, SD;
	double numberOfFeAtoms;
}

-(id)initWithParameters:(NSArray *)paramsIn;
-(id)initWithParameters:(NSArray *)paramsIn andLandauTerms:(NSArray *)landauIn;

@end
