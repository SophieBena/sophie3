# Parameter Calibration Protocol
🔹 Retrieves the number of parameters in the model
```
-(NSUInteger)getNumberOfFreeParameters
```
🔹 Retrieves an array of parameter names
```
-(NSArray *)getArrayOfNamesOfFreeParameters
```
🔹 Retrieves the value of a model paramter of the given name
```
-(double)getValueForParameterName:(NSString *)name
```
🔹 Sets the value of a model parameter of the given name
```
-(BOOL)setParameterName:(NSString *)name tovalue:(double)value
```
🔹 Retrieves the units of a model parameter of the given name
```
-(NSString *)getUnitsForParameterName:(NSString *)name
```
🔹 Retrieve a vector of derivatives of the Gibbs Free Energy with respect to the model parameters
```
-(DoubleVector *)getDgDwFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Retrieve a vector of derivatives of the chemical potentials of all solution components with respect to the named model parameter
```
-(DoubleVector *)getChemicalPotentialDerivativesForParameter:(NSString *)name usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Retrieve a matrix of derivatives of the chemical potentials of all solution components (rows) with respect to all model parameters (columns)
```
-(DoubleMatrix *)getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Returns YES if the calibration propotocol is supported by this phase (testing for the existence of this method should be done first).
```
-(BOOL)supportsParameterCalibration
```