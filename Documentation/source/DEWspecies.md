## DEWspecies Class  
### Class Inheritance  
NSObject ▶️ DEWspecies  

### Protocols Implemented  
None  

### Properties  

```
@property (strong, readonly)NSMutableDictionary *dictionaryOfSpecies
```

```
@property (strong, readonly)NSMutableSet *setOfCationSpecies
```

```
@property (strong, readonly)NSMutableSet *setOfAnionSpecies
```

```
@property (strong, readonly)NSMutableSet *setOfNeutralSpecies
```

### Class Methods  
None  

### Instance Methods  

```
- (HKFspeciesProperties *)HKFclassForSpeciesNamed:(NSString *)name
```

```
- (NSString *)HKFsymbolForSpeciesNamed:(NSString *)name
```

```
- (NSString *)HKFcommentForSpeciesNamed:(NSString *)name
```

```
- (BOOL)namedSpeciesInCollection:(NSString *)name
```

```

- (NSUInteger)numberofCationSpecies
```

```
- (NSUInteger)numberOfAnionSpecies
```

```
- (NSUInteger)numberOfNeutralSpecies
```

```
- (NSUInteger)==numberOfSpecies==
```