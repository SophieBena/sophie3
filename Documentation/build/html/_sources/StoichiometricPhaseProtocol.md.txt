## Stoichiometric Phase Protocol


<p style="border-style: solid; padding:3px">
\({{C_P}}\) = isobaric heat capacity<br />
<em>V</em> = volume
</p>

🔹  T (K), P (bars) => Gibbs free energy (J)
```
-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => enthalpy (J)
```
-(double)getEnthalpyFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => entropy (J/K)
```
-(double)getEntropyFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\({{C_P}}\\) (J/K)
```
-(double)getHeatCapacityFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\(\frac{{\partial {C_P}}}{{\partial T}}\\) (J/K<sup>2</sup>)
```
-(double)getDcpDtFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => *V* (J/bar)
```
-(double)getVolumeFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\(\frac{{\partial V}}{{\partial T}}\\) (J/bar-K)
```
-(double)getDvDtFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\(\frac{{\partial V}}{{\partial P}}\\) (J/bar<sup>2</sup>)<br />
```
-(double)getDvDpFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial {T^2}}}\\) (J/bar-K<sup>2</sup>)<br />
```
-(double)getD2vDt2FromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial T\partial P}}\\) (J/bar<sup>2</sup>-K)<br />  
```
-(double)getD2vDtDpFromT:(double)t andP:(double)p
```
🔹  T (K), P (bars) => \\(\frac{{{\partial ^2}V}}{{\partial {P^2}}}\\) (J/bar<sup>3</sup>)<br />
```
-(double)getD2vDp2FromT:(double)t andP:(double)p
```
🔹  optional: T (K), P (bars) => chemical potential (J)
```
-(double)getChemicalPotentialFromT:(double)t andP:(double)p
```
🔹  optional: formulae of the phase
```
-(NSString )getFormulaFromInternalVariables
```