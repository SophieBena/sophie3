## Equilibrate Class  
### Class Inheritance  
NSObject ▶️ Equilibrate  

### Protocols Implemented  
None  

### Properties  

```
@property (strong, readonly) NSArray <NSString *> *systemOxides
```
```
@property (strong, readonly) NSArray <PhaseBase *> *bulkComposition  
```
```
@property (strong, readonly) NSArray <NSNumber *> *elementsKnownToSystem 
```
```
@property (strong, readonly) NSDictionary <NSString *, PhaseBase *> *dictionaryOfPhasesInSystem  
```
```
@property (assign, readwrite) BOOL debugV 
```
```
@property (assign, readwrite) BOOL debugS  
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_DISTRIBUTION 
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_CHEM_POTENTIAL  
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_RANK  
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_ZERO_LINEAR  
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_QUAD_ITERS  
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_FORCE_COMPONENTS  
```
```
@property (strong, readwrite) NSNumber *PPOPTIONS_DETECT_MINIMAL  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_CHEM_POTENTIAL  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_RANK  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_QUAD_OPTIMAL  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_QUAD_SUBOPTIMAL  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_QUAD_MAX_ITERS  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_MAX_ITERS  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_RELTEST  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_MIN_STEPLENGTH  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_COMPONENT_MINIMUM  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_MASSOUT  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_MOLESIN  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_MINIMAL  
```
```
@property (strong, readwrite) NSNumber *PPPARAMETERS_CYCLES_MAX  
```

### Class Methods  
None  

### Instance Methods  

```
-(NSMutableArray *) oxideNames
```
 🔹 Provides an array of NSString objects that correspond to the names of phases known to this class.
 
 @return  
 An array of Cocoa strings


```
-(NSArray *) phaseNames
```
   🔹 Converts moles of elements into moles of oxides.
 
 @param  
 e vector of length 107 containing moles of elements
 
 @return  
 Double vector of length 15 containing moles of oxides


```
-(DoubleVector *)molesOfOxidesFromMolesOfElements:(double *)e
```
   🔹 Converts moles of oxides into moles of elements.
 
 @param  
 e vector of length 15 containing moles of oxides
 
 @return  
 Double vector of length 107 containing moles of elements



```
-(DoubleVector *)gramsOfOxidesFromMolesOfElements:(double *)e
```

 🔹 Provides a test on the feasibility of an input composition.
 
@param  
compositionInWtPercentOxides array of NSNumber objects containing wt% values of system oxides  
   
@return  
YES = feasible, NO = unfeasible

```
-(BOOL)compositionIsFeasible:(NSArray *)compositionInWtPercentOxides forSolution:(id <SolutionPhaseProtocol>)omniComponentPhase
```
```
-(instancetype)init
```
```
-(instancetype)initWithSystemOxides:(NSArray <NSString *> *)systemOxides
      andDictionaryOfPhasesInSystem:(NSDictionary <NSString *, PhaseBase *> *)dictionaryOfPhasesInSystem
```
  

```
-(instancetype)initWithInformationalDebug:(BOOL)debugS
                          andVerboseDebug:(BOOL)debugV
                          andSystemOxides:(NSArray <NSString *> *)systemOxides
            andDictionaryOfPhasesInSystem:(NSDictionary <NSString *, PhaseBase *> *)dictionaryOfPhasesInSystem
```
 🔹 Determines if a phase is supersaturated in the system realative to a liquid phase or a collection of solid phases.   
   
 The most supersaturated phase is transferred from the potentialPhaseList dictionary to the phasesInSystem dictionary.  
 
 @return   
   A YES if supersaturation is detected. A NO if the current system is stable.  

```
-(BOOL)evaluateSaturationState:(NSString *)lastPhaseRemovedFromSystem
```
  🔹 Sets the temperature.  
  
 @param  
   t in K  

```
-(void)setTemperature:(double)t
```
  🔹 Increments the temperature (assumes a previously successful execute step).  
   
 @param  
 t in K  


```
-(void)incrementTemperature:(double)t
```
  🔹 Sets the pressure.  
 
 @param  
   p in bars  


```
-(void)setPressure:(double)p
```
  🔹 Increments the pressure (assumes a previously successful execute step).  
 
 @param  
 p in bars  

```
-(void)incrementPressure:(double)p
```
 🔹 Sets the entropy.  
   
 @param  
 s in J/K  

```
-(void)setEntropy:(double)s
```
  🔹 Increments the entropy (assumes a previously successful execute step).  
 
 @param  
 s in J/K  

```
-(void)incrementEntropy:(double)s
```
 🔹 Sets the volume.  
 
 @param  
 v in J/bar  


```
-(void)setVolume:(double)v
```
 🔹 Increments the volume (assumes a previously successful execute step).  
   
 @param  
 v in J/bar  

```
-(void)incrementVolume:(double)v
```
 🔹 Sets the bulk composition of the system.  
   
 Wipes the system clean of all phases and allocates an instance of a liquid phase (which may be metastable).  
   
 @param  
   wtOxides NSArray of length 15 containing objects of class NSNumber. Each object holds a doubleValue of grams of an oxide.  



```
-(void)setComposition:(double *)wtOxides
```
 🔹 Sets the permissable phase state.  
   
 Turns on or off entries in the permissablePhasesState array that determine inclusion of phases in the system.  
   
 @param  
 arrayOfStateValues NSArray of containing objects of class NSNumber.  Each object holds a BOOL value. The order is the same as returned by the method phaseNames.  
  

```
-(void)setPermissablePhasesState:(NSArray *)arrayOfStateValues
```
 🔹 Set the options for the Phase equilibration.  
   
 @param  
   dictionaryOfCalculationOptions NSDictionary containing options corresponding to known keys:  
   - value of the key @"ordinate" is either @"Temperature", @"Enthalpy", or @"Entropy"  
   - value of the key @"abscissa" is either @"Pressure" or @"Volume"  
   - value of the key @"buffer" is either @"HM", @"IW", @"NNO", or @"QFM"  
   - value of the key @"imposeBuffer" is either the boolean YES or NO contained in an NSNumber object  
   - value of the key @"bufferOffset" is a double contained in an NSNumber object  

```
-(void)setCalculationOptions:(NSDictionary *)dictionaryOfCalculationOptions
```
  🔹 Executes an equilibration calculation.  
   
 @return  
   An NSDictionary of status and output items:  
     - status key:  A flag indicating success or failure of the calculation  
     - results key: A EquilibrateState object  


```
-(NSDictionary *)execute
```
 🔹 Computes chemical potentials of elements from current system phase assemblage.  
   
 The solution will be unique if the system is in equilibrium and the norm of the residuals will be zero.   
 A non-zero residual norm indicates a disequilibrium state.  

 @param  
 elementMu of moles of elements on length 107. Solution is returned in this vector  
   
 @return  
 Residual norm.  


```
-(double)chemicalPotentialsOfElements:(double *)elementMu
```
 🔹 
```
-(NSString *)equilibrateResultsAsXML
```