//
//  HKFspeciesComposite.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 1/19/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "HKFspeciesComposite.h"
#import "HKFspeciesProperties.h"
#import "DEWspecies.h"

@implementation HKFspeciesComposite

- (instancetype)initWithDEWspeciesInstance:(DEWspecies *)dewSpecies
                 andWithReactionDictionary:(NSDictionary <NSString *, NSNumber *> *)nameKeysAndStoichiometryValues {
    self = [super init];
    if (self) {
        if (!dewSpecies || !nameKeysAndStoichiometryValues) return nil;
        if (nameKeysAndStoichiometryValues.count == 0) return nil;

        _hkfSpecies = [NSMutableArray arrayWithCapacity:2];
        _hkfStoichiometry = [NSMutableArray arrayWithCapacity:2];
        _dewSpecies = dewSpecies;

        _count = 0;
        for (NSString *name in nameKeysAndStoichiometryValues.allKeys) {
            NSNumber *stoichiometry = [nameKeysAndStoichiometryValues objectForKey:name];
            [_hkfSpecies addObject:name];
            [_hkfStoichiometry addObject:stoichiometry];
            _count++;
        }

        _correctionToGibbsEnergy = 0.0;
    }
    return self;
}

- (double)gSolventAtT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf gSolventAtT:t andP:p];
    }
    return result;
}

- (double)DgSolventDtAtT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf DgSolventDtAtT:t andP:p];
    }
    return result;
}

- (double)DgSolventDpAtT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf DgSolventDpAtT:t andP:p];
    }
    return result;
}

- (double)D2gSolventDt2AtT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf D2gSolventDt2AtT:t andP:p];
    }
    return result;
}

- (double)D2gSolventDtDpAtT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf D2gSolventDtDpAtT:t andP:p];
    }
    return result;
}

- (double)D2gSolventDp2AtT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf D2gSolventDp2AtT:t andP:p];
    }
    return result;
}

- (double)omegaAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf omegaAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
    }
    return result;
}

- (double)DomegaDtAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf DomegaDtAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
    }
    return result;
}

- (double)DomegaDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf DomegaDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
    }
    return result;
}

- (double)D2omegaDt2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf D2omegaDt2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
    }
    return result;
}

- (double)D2omegaDtDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf D2omegaDtDpAtT:t andP:p withOmegaRef:omegaRef andCharge:z];
    }
    return result;
}

- (double)D2omegaDp2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf D2omegaDp2AtT:t andP:p withOmegaRef:omegaRef andCharge:z];
    }
    return result;
}

#pragma mark -
#pragma mark StoichiometricPhaseProtocol protocol

- (double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getGibbsFreeEnergyFromT:t andP:p];
    }
    return result + self.correctionToGibbsEnergy;
}
- (double)getEnthalpyFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getEnthalpyFromT:t andP:p];
    }
    return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getEntropyFromT:t andP:p];
    }
    return result;
}

- (double)getHeatCapacityFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getHeatCapacityFromT:t andP:p];
    }
    return result;
}

- (double)getDcpDtFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getDcpDtFromT:t andP:p];
    }
    return result;
}

- (double)getVolumeFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getVolumeFromT:t andP:p];
    }
    return result;
}

- (double)getDvDtFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getDvDtFromT:t andP:p];
    }
    return result;
}

- (double)getDvDpFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getDvDpFromT:t andP:p];
    }
    return result;
}

- (double)getD2vDt2FromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getD2vDt2FromT:t andP:p];
    }
    return result;
}

- (double)getD2vDtDpFromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getD2vDtDpFromT:t andP:p];
    }
    return result;
}

- (double)getD2vDp2FromT:(double)t andP:(double)p {
    double result = 0.0;
    for (NSUInteger i=0; i<self.count; i++) {
        HKFspeciesProperties *hkf = [self.dewSpecies HKFclassForSpeciesNamed:[self.hkfSpecies objectAtIndex:i]];
        double stoichiometry = [[self.hkfStoichiometry objectAtIndex:i] doubleValue];
        result += stoichiometry*[hkf getD2vDp2FromT:t andP:p];
    }
    return result;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *khkfSpecies              = @"hkfSpecies";
static NSString *khkfStoichiomtry         = @"hkfStoichiometry";
static NSString *kcount                   = @"count";
static NSString *kcorrectionToGibbsEnergy = @"correctionToGibbsEnergy";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super init])) {
#ifdef __APPLE__
        _hkfSpecies = (NSMutableArray *) [aDecoder decodeObjectOfClass:[NSArray class] forKey:khkfSpecies];
        _hkfStoichiometry = (NSMutableArray *) [aDecoder decodeObjectOfClass:[NSArray class] forKey:khkfStoichiomtry];
#else
        _hkfSpecies = (NSMutableArray *) [aDecoder decodeObjectForKey:khkfSpecies];
        _hkfStoichiometry = (NSMutableArray *) [aDecoder decodeObjectForKey:khkfStoichiomtry];
#endif
        _count = [aDecoder decodeIntegerForKey:kcount];
        _correctionToGibbsEnergy = [aDecoder decodeDoubleForKey:kcorrectionToGibbsEnergy];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeObject:_hkfSpecies forKey:khkfSpecies];
        [aCoder encodeObject:_hkfStoichiometry forKey:khkfStoichiomtry];
        [aCoder encodeInteger:_count forKey:kcount];
        [aCoder encodeDouble:_correctionToGibbsEnergy forKey:kcorrectionToGibbsEnergy];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class HKFspeciesProperties only supports NSKeyedArchiver coders."];
}

+ (BOOL)supportsSecureCoding {
    return YES;
}

@end
