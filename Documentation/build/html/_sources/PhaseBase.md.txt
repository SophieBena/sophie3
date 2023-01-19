## PhaseBase Class

### Class Inheritance
NSObject ▶️ PhaseBase

### Protocols Implemented
NSSecureCoding

### Properties
🔹 Chemical formula of the phase. A custom setter routine is defined to assign the property. This class automatically initializes the formula conversion array and molecular weight.

```
@property (nonatomic, readwrite, copy) NSString phaseFormula; 
```
🔹 Name of the phase
```
@property (readwrite, copy) NSString phaseName
```
🔹 Molecular weight of the phase. Computed when formula is set. Reeadonly.
```
@property (readwrite, assign) double mw
```
🔹 Number of moles of each element calculated from the phase formula. Computed when formula is set. Reeadonly.
```
@property (readonly) DoubleVector formulaAsElementArray
```

🔹  Entropy of the elements calculated from phase stoichiometry using values tabulated by
 Robie, Hemingway and Fisher (1979) USGS Bull 1452
 ```
 @property (readwrite, assign) double entropyFromRobieEtAl1979
 ```
 
### Class Methods
🔹  Name of element from element index
```
+(NSString *)elementNameFromAtomicNumber:(NSUInteger)atomicNumber
```
### Instance Methods
🔹  Moles of elements (standard order) => Moles of phase
```
-(double)convertElementsToMolesOfPhase:(double )e
```
🔹  Moles of elements (standard order) => Total mass of the phase (g)
```
-(double)convertElementsToMassOfPhase:(double )e
```