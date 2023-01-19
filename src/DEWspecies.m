//
//  DEWspecies.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 11/17/15.
//  Copyright © 2015 Mark Ghiorso. All rights reserved.
//

#import "DEWspecies.h"
#import "HKFspeciesProperties.h"

// Load latest DEW model parameters
#import "DEWmodelParameters.h"

@implementation DEWspecies

static NSArray <NSString *> *_chargeDesignators;

#pragma mark -
#pragma mark Class functions

+ (void)initialize {
    _chargeDesignators = [NSArray arrayWithObjects:@"(0)", @"(-)", @"(+)",
                          @"(-1)", @"(-2)", @"(-3)", @"(-4)", @"(-5)",
                          @"(+1)", @"(+2)", @"(+3)", @"(+4)", @"(+5)",
                          nil];
}

#pragma mark -
#pragma mark Initializers

- (NSString *)formulaFromDEWsymbol:(NSString *)symbol {
    NSRange range = NSMakeRange(NSNotFound, 0);
    for (NSString *suffix in _chargeDesignators) {
        if ((range = [symbol rangeOfString:suffix]).location != NSNotFound) break;
    }
    if (range.location != NSNotFound) return [symbol substringToIndex:range.location];
    else return symbol;
}

#define CALtoJOULES   4.184

- (instancetype)init {
    self = [super init];
    if (self) {
        _dictionaryOfSpecies = [NSMutableDictionary dictionaryWithCapacity:nDEW];
        _setOfCationSpecies  = [NSMutableSet setWithCapacity:nDEW];
        _setOfAnionSpecies   = [NSMutableSet setWithCapacity:nDEW];
        _setOfNeutralSpecies = [NSMutableSet setWithCapacity:nDEW];

        for (NSUInteger i=0; i<nDEW; i++) {
            NSString *name   = [NSString stringWithCString:dewTable[i].name   encoding:NSASCIIStringEncoding];
            NSString *symbol = [NSString stringWithCString:dewTable[i].symbol encoding:NSASCIIStringEncoding];
            double g     = atof(dewTable[i].g)*CALtoJOULES;
            double h     = atof(dewTable[i].h)*CALtoJOULES;
            double s     = atof(dewTable[i].s)*CALtoJOULES;
            double v     = atof(dewTable[i].v)/10.0; // cm^3 to J/bar
            double cp    = atof(dewTable[i].cp)*CALtoJOULES;
            double omega = atof(dewTable[i].omega)*CALtoJOULES*SCALEforOmega;
            double z     = atof(dewTable[i].z);
            double a1    = atof(dewTable[i].a1)*CALtoJOULES*SCALEforA1;
            double a2    = atof(dewTable[i].a2)*CALtoJOULES*SCALEforA2;
            double a3    = atof(dewTable[i].a3)*CALtoJOULES;
            double a4    = atof(dewTable[i].a4)*CALtoJOULES*SCALEforA4;
            double c1    = atof(dewTable[i].c1)*CALtoJOULES;
            double c2    = atof(dewTable[i].c2)*CALtoJOULES*SCALEforC2;
            NSString *comment = [NSString stringWithCString:dewTable[i].comment encoding:NSUTF8StringEncoding];

            HKFspeciesProperties *hkf = [[HKFspeciesProperties alloc] initWithGibbsFreeEnergyOfFormation:g
                                                                              andWithEnthalpyOfFormation:h
                                                                                 andWithReferenceEntropy:s
                                                                                  andWithReferenceVolume:v
                                                                            andWithReferenceHeatCapacity:cp
                                                                                   andWithHKFa1Parameter:a1
                                                                                   andWithHKFa2Parameter:a2
                                                                                   andWithHKFa3Parameter:a3
                                                                                   andWithHKFa4Parameter:a4
                                                                                   andWithHKFc1Parameter:c1
                                                                                   andWithHKFc2Parameter:c2
                                                                                andWithHKFomegaParameter:omega
                                                                                           andWithCharge:z];
            hkf.phaseName = name;
            hkf.phaseFormula = [self formulaFromDEWsymbol:symbol];

            // type HKFSpeciesArray enumeration
            NSArray *entry = [NSArray arrayWithObjects:hkf, symbol, comment, nil];
            [_dictionaryOfSpecies setObject:entry forKey:name];

            if (z < 0.0) [_setOfAnionSpecies addObject:name];
            else if (z > 0.0) [_setOfCationSpecies addObject:name];
            else [_setOfNeutralSpecies addObject:name];
        }
    }
    return self;
}

#pragma mark -
#pragma mark API

- (BOOL)namedSpeciesInCollection:(NSString *)name {
    if ([_dictionaryOfSpecies objectForKey:name]) return YES;
    else return NO;
}

- (HKFspeciesProperties *)HKFclassForSpeciesNamed:(NSString *)name {
    return [[_dictionaryOfSpecies objectForKey:name] objectAtIndex:HKFSpeciesOBJECT];
}

- (NSString *)HKFsymbolForSpeciesNamed:(NSString *)name {
    return [[_dictionaryOfSpecies objectForKey:name] objectAtIndex:HKFSpeciesSYMBOL];
}

- (NSString *)HKFcommentForSpeciesNamed:(NSString *)name {
    return [[_dictionaryOfSpecies objectForKey:name] objectAtIndex:HKFspeciesCOMMENTS];
}

- (NSUInteger)numberofCationSpecies {
    return [_setOfCationSpecies count];
}

- (NSUInteger)numberOfAnionSpecies {
    return [_setOfAnionSpecies count];
}

- (NSUInteger)numberOfNeutralSpecies {
    return [_setOfNeutralSpecies count];
}

- (NSUInteger)numberOfSpecies {
    return [_dictionaryOfSpecies count];
}

@end
