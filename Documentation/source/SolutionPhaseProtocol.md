## Solution Phase Protocol  

<p style="border-style: solid; padding:3px">
\({a_i}\) = activity of the i<sup>th</sup> endmember component<br />  
\({n_j}\) = moles of the j<sup>th</sup> endmember component<br />  
<em>G</em> = Gibbs free energy<br />
<em>S</em> = entropy<br />
<em>G</em> = volume<br />
\({C_P}\) = isobaric heat capacity
</p>

🔹 By default all thermodynamic quantities are calculated to include endmember properties (yesForMixing = NO). Calling this function with YES as an argument converts output to mixing quantities.
```
-(void)setResultsToMixingQuantities:(BOOL)yesForMixing
```
🔹 Retrieves the number of endmember components in the system
```
-(NSUInteger)numberOfSolutionComponents
```
🔹 Retrieves superclass instance of PhaseBase object for component at specified index
```
-(id)componentAtIndex:(NSUInteger)index
```
🔹 Moles of endmember components => validity of input values
```
-(BOOL)testPermissibleValuesOfComponents:(double *)m
```
🔹 Moles of elements (standard order) => Moles of end member components of the phase
```
-(DoubleVector )convertElementsToMoles:(double *)e
```
🔹 Moles of elements (standard order) => Total moles of end member components of the phase
```
-(double)convertElementsToTotalMoles:(double )e
```
🔹 Moles of elements (standard order) => Total mass of the phase (g)
```
-(double)convertElementsToTotalMass:(double *)e
```
🔹 Moles of endmember components => Mole fractions of endmember components
```
-(DoubleVector )convertMolesToMoleFractions:(double *)m
```
🔹 Moles of endmember components => Moles of elements (standard order)
```
-(DoubleVector )convertMolesToElements:(double )m
```
🔹 Moles of endmember components => Molar sum
```
-(double)totalMolesFromMolesOfComponents:(double *)m
```
🔹 Moles of components, T (K), P (bars) => activities of endmember components
```
-(DoubleVector )getActivityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => chemical potentials of endmember components (J)
```
-(DoubleVector )getChemicalPotentialFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial G}}{{\partial {n_i}}}\\)
```
-(DoubleMatrix )getDaDmFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => Gibbs free energy (J)
```
-(double)getGibbsFreeEnergyFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial G}}{{\partial {n_i}}}\\) (J)
```
-(DoubleVector )getDgDmFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^2}G}}{{\partial {n_i}\partial {n_j}}}\\) (J)
```
-(DoubleMatrix )getD2gDm2FromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^3}G}}{{\partial {n_i}\partial {n_j}\partial {n_k}}}\\) (J)
```
-(DoubleTensor )getD3gDm3FromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => enthalpy (J)
```
-(double)getEnthalpyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => entropy (J/K)
```
-(double)getEntropyFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial S}}{{\partial {n_i}}}\\) (J/K)
```
-(DoubleVector )getDsDmFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^2}S}}{{\partial {n_i}\partial {n_j}}}\\) (J/K)
```
-(DoubleMatrix )getD2sDm2FromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => isobaric heat capacity (J/K)
```
-(double)getHeatCapacityFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial {C_P}}}{{\partial T}}\\) (J/K<sup>2</sup>)
```
-(double)getDcpDtFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial {C_P}}}{{\partial {n_i}}}\\) (J/K)
```
-(DoubleVector )getDCpDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => volume (J/bar)
```
-(double)getVolumeFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial V}}{{\partial {n_i}}}\\) (J/bar)
```
-(DoubleVector )getDvDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial {n_i}\partial {n_j}}}\\) (J/bar)
```
-(DoubleMatrix )getD2vDm2FromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial V}}{{\partial T}}\\)
```
-(double)getDvDtFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial V}}{{\partial P}}\\) (J/bar<sup>2</sup>)
```
-(double)getDvDpFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial {T^2}}}\\) (J/bar-K<sup>2</sup>)
```
-(double)getD2vDt2FromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial T\partial P}}\\) (J/bar<sup>2</sup>-K)
```
-(double)getD2vDtDpFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial {P^2}}}\\) (J/bar<sup>3</sup>)
```
-(double)getD2vDp2FromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial V}}{{\partial {n_i}\partial T}}\\) (J/bar-K)
```
-(DoubleVector )getD2vDmDtFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => \\(\frac{{\partial V}}{{\partial {n_i}\partial P}}\\) (J/bar<sup>2</sup>)
```
-(DoubleVector )getD2vDmDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p
```
🔹 Moles of components, T (K), P (bars) => formulae as an NSString object
```
-(NSString )getFormulaFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Retrieves the number of species (dependent endmembers with positive mole fractions) in the soluton
```
-(NSUInteger)numberOfSolutionSpecies
```
🔹 Retrieves the name of the solution species at the specified index
```
-(NSString )nameOfSolutionSpeciesAtIndex:(NSUInteger)index
```
🔹 Moles of solution species => moles of endmember components
```
-(DoubleVector )convertMolesOfSpeciesToMolesOfComponents:(double )mSpecies
```
🔹 Retrieves an elemental stoichiometry vector for the species at the specified index
```
-(DoubleVector )elementalCompositionOfSpeciesAtIndex:(NSUInteger)index
```
🔹 Moles of components, T (K), P (bars) => chemical potentials of solution species (J)
```
-(DoubleVector )chemicalPotentialsOfSpeciesFromMolesOfComponents:(double )m andT:(double)t andP:(double)p
```
🔹 Chemical potentials of endmembers, T (K), P (bars) => chemical affinity and composition of the phase, etc.<br />

@param chemicalPotentials <br />
Chemical potentials of endmember components in the solution. A zero entry indicates that a component is absent. <br /> 

@param t<br />  
  temperature in K <br /> 
  
@param p<br />  
  pressure in bars<br />  
  
@return <br /> 
NSArray output structure:<br />  
\- (0): NSNumber object wrapping a double - chemical affinity<br />  
\- (1...NA): NSnumber object wrapping a double - X[0] - X[NA-1], mole fraction compositional variables<br />  
\- (NA+1): NSNumber object wapping a BOOL - convergence flag<br />  
\- (NA+2): NSNumber object wrapping an NSUInteger - iteration count<br />  
\- (NA+3): NSNumber object wrapping a double - number of atoms in formula unit to scale affinity<br />  
\- (NA+4): NSNumber object wrapping a double - approximate error in calculated affinity<br />

```
-(NSArray )affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p
```
🔹 
```
-(void)incrementInstanceCountOfPhase
```
🔹 
```
-(void)decrementInstanceCountOfPhase
```
🔹 
```
-(NSString *)nameOfPhaseWithComposition:(double *)refMoles
```
🔹 
```
-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p
```
