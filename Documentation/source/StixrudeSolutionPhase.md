## StixrudeSolutionPhase Class  
### Class Inheritance  
NSObject ▶️  [PhaseBase](PhaseBase.html) ▶️  StixrudeSolutionPhase  

### Protocols Implemented  
[SolutionPhaseProtocol](SolutionPhaseProtocol.html)    

### Properties  
```
@property NSUInteger numberOfIterationsAllowedInSaturationMethod
```
```
@property double minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod
```

### Class Methods  
None  

### Instance Methods  

```
-(id)initWithPhaseName:(NSString *)name 
		   withSpecies:(NSArray *)speciesProperties
	withSpeciesWeights:(NSArray *)speciesWeights
	withSpeciesXfactor:(NSArray *)speciesXfactor
			withWarray:(NSArray *)wIn 
			 withSites:(NSArray *)sites 
			WithNAtoms:(NSUInteger)nAtomIn
  withElementConvertor:(void (^)(double *, double *))convertor
 andWithFormulaDisplay:(NSString *(^)(double *))display
andWithSaturationStateAdjuster:(void (^)(double *, double *))adjuster

```