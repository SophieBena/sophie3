//
//  DEWspecies.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 11/17/15.
//  Copyright © 2015 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>

@class HKFspeciesProperties;

typedef NS_ENUM(NSInteger, HKFSpeciesARRAY) {
    HKFSpeciesOBJECT = 0,
    HKFSpeciesSYMBOL,
    HKFspeciesCOMMENTS
};

@interface DEWspecies : NSObject {
    NSMutableDictionary <NSString *, NSArray *> *_dictionaryOfSpecies;
    NSMutableSet <NSString *> *_setOfCationSpecies;
    NSMutableSet <NSString *> *_setOfAnionSpecies;
    NSMutableSet <NSString *> *_setOfNeutralSpecies;
}

@property (strong, readwrite)NSMutableDictionary *dictionaryOfSpecies;
@property (strong, readonly)NSMutableSet *setOfCationSpecies;
@property (strong, readonly)NSMutableSet *setOfAnionSpecies;
@property (strong, readonly)NSMutableSet *setOfNeutralSpecies;

- (HKFspeciesProperties *)HKFclassForSpeciesNamed:(NSString *)name;
- (NSString *)HKFsymbolForSpeciesNamed:(NSString *)name;
- (NSString *)HKFcommentForSpeciesNamed:(NSString *)name;
- (BOOL)namedSpeciesInCollection:(NSString *)name;

- (NSUInteger)numberofCationSpecies;
- (NSUInteger)numberOfAnionSpecies;
- (NSUInteger)numberOfNeutralSpecies;
- (NSUInteger)numberOfSpecies;

@end
