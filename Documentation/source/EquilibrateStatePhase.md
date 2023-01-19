## EquilibrateStatePhase Class  
### Class Inheritance  
NSObject ▶️ EquilibrateStatePhase  

### Protocols Implemented  
None  

### Properties  

```
@property (readwrite, strong) id phaseClassInstance
```

```
@property (readwrite, assign) double mass
```

```
@property (readwrite, assign) DoubleVector *bulkCompositionInElements
```

```
@property (readwrite, assign) double affinity
```

### Class Methods  
None  

### Instance Methods 
  
🔹  Returns the molar sum of the bulkCompositionInElements array.

```
-(double)sumOfMolesOfElements
```

🔹 Determines the indices of non-zero entries for moles of components in the phase.  
 
 @return  
 An array of indices of non-zero components.  
```
-(NSArray *)indicesOfNonZeroMolesOfEndmemberComponents
```

🔹   Determines non-zero entries for moles of components and composition/potential energy arrays for the phase.  
 
 **Procedure**:    
1. Determine if the phase is a solution or a stoichiometric phase.
1. For each component of the solution or for the stoichiometric phase, report:  
   - Total moles.  
   - A pointer to a vector that represents the stoichiometry in terms of elements.
   - Chemical potential.  
   - Hessian (row/column), if the phase is a solution. Otherwise, the entry is missing: An array that holds the derivatives of the chemical potential with repect to all solution components whose molar abundances are non-zero.  
3. Repeat the above as required for all the endmembers of the solution (See return section for format).  
 
 @param constraint  
   Specifies constraints for calculation of potential energy functions:  
   - independentT & independentP =>> derivatives of the Gibbs free energy  
   - independentT & independentV =>> derivatives of the Helmholtz energy  
   - independentS & independentP =>> derivatives of the enthalpy  
   - independentS & independentV =>> derinatives of the internal energy  
   - independentH & independentP =>> derivatives of the negative entropy  
 
 @param t  
   Temperature in K.  
 
 @param p  
   Pressure in bars.  
 
 @return  
   An NSArray object with the following structure:  

  NSArray object: One entry for each non-zero component or one entry for the stoichiometric phase:  
   - NSNumber object: Contains a double value.  This is the number of moles of the component/phase.  
    - NSValue object: Contains a pointer to double value. This is an array of length 107 containing stoiciometric abundances of the elements in one mole of the component/phase.  
    - NSNumber object: Contains a double value.  This is the first derivative of the potential function with respect to the component/phase.  
   - (Only for components) NSArray object: Contains NSNumber objects that hold a double value. These values are the elements of the Hessian (second derivative matrix of the potential function) row for the non-zero components in the solution.  
 
NOTE: Additional entries are included to handle constraint vectors. See comments in code.  

 The number of objects in the array is one for a stoichiometric phase and equal to the number of components with non-zero molar abundances for a solution.  

```
-(NSArray *)nonZeroMolesOfEndmemberComponents:(NSUInteger)constraint andT:(double)t andP:(double)p
```

🔹 Computes the potential function for the phase.  
 
 @param constraint  
 Specifies constraints for calculation of potential energy functions:  
 - independentT & independentP =>> derivatives of the Gibbs free energy  
 - independentT & independentV =>> derivatives of the Helmholtz energy  
 - independentS & independentP =>> derivatives of the enthalpy  
 - independentS & independentV =>> derinatives of the internal energy  
 - independentH & independentP =>> derivatives of the negative entropy  
   
 @param t  
 Temperature in K.  
 
 @param p  
 Pressure in bars.  
   
 @param deltaMoles  
 Corrections to non-zero molar abundances of components (or the moles of the phase). Usually determined by a quadratic minimization solution returned from HFTI.   
 
 @param scaleFractor  
 A proposed scaler multiplier for computation of the molar abundances of each component or the phase, as: mNew = mOld + scaleFactor * deltaMoles.  
   
 @return  
 The value of the potential or the number compositionIsUnacceptable if the computed composition is infeasible.  

```
-(double)potentialFunctionFor:(NSUInteger)constraint andT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```

```
-(double)entropyFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```

```
-(double)volumeFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```

```
-(double)enthalpyFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```

```
-(double)heatCapacityFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```

```
-(double)dvdpFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```

```
-(double)dvdtFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor

```
🔹 Resets molar abundances of phases and components in phases by updating the bulkCompositionInElements array.
 
 @param deltaMoles  
 Corrections to non-zero molar abundances of components (or the moles of the phase). Usually determined by a quadratic minimization solution returned from HFTI.   
 
 @param scaleFractor  
 A proposed scaler multiplier for computation of the molar abundances of each component or the phase, as: mNew = mOld + scaleFactor * deltaMoles.  
```
-(void)addDeltaMolesToNonZeroMolesOfEndmemberComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor
```