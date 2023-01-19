## FeldsparBerman Class
### Class Inheritance
NSObject ▶️ [PhaseBase](PhaseBase.html) ▶️ FeldsparBerman
### Protocols Implemented
[SolutionPhaseProtocol](SolutionPhaseProtocol.html)<br />    
[ParameterCalibrationProtocol](ParameterCalibrationProtocol.html)    

### Properties
🔹 Name of phase assigned at initialization (provided so that multiple instances of the class can be easily identified)
```
@property (readwrite,copy) NSString *operationParent
```
🔹 Parameters of the model
```
@property (readwrite, nonatomic) double whabor
@property (readwrite, nonatomic) double wsabor
@property (readwrite, nonatomic) double wvabor
@property (readwrite, nonatomic) double whorab
@property (readwrite, nonatomic) double wsorab
@property (readwrite, nonatomic) double wvorab
@property (readwrite, nonatomic) double whaban
@property (readwrite, nonatomic) double whanab
@property (readwrite, nonatomic) double whoran
@property (readwrite, nonatomic) double whanor
@property (readwrite, nonatomic) double wvanor
@property (readwrite, nonatomic) double whabanor
@property (readwrite, nonatomic) double wvabanor
```
🔹 An array of class instances of thermodynamic properties of endmember components (These classes generally are rooted in PhaseBase and implement the Stoichiometric Phase Protocol.)
```
@property (readonly) NSArray *endmembers
```
### Class Methods
None  

### Instance Methods
🔹 Initialize instance of class with name for phase
```
-(id)initWithCompositionConstraint:(NSString *)name
```
🔹 Increase or decrease reference count of the number of instances of this class
```
-(void)incrementInstanceCountOfPhase
-(void)decrementInstanceCountOfPhase
```
🔹 Determinine if an instance of this phase with specified composition at specified T and P is table vis-a-vis phase separation
```
-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p
```
🔹 Classify a phase with specifed composition, e.g., Sanidine, Plaguioclase, Anorthoclase
```
-(NSString *)nameOfPhaseWithComposition:(double *)refMoles
```