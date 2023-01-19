## DEWDielectricConstant Class  
### Class Inheritance  
NSObject ▶️ DEWDielectricConstant  

### Protocols Implemented  
NSSecureCoding  

### Properties  

```
@property (assign) double tLast  
```
```
@property (assign) double pLast  
```
```
@property (assign) double rhoLast  
```
```
@property (assign) double DrhoDtLast  
```
```
@property (assign) double DrhoDpLast  
```
```
@property (assign) double D2rhoDt2Last  
```
```
@property (assign) double D2rhoDtDpLast  
```
```
@property (assign) double D2rhoDp2Last   
```
```
@property (strong) GenericH2O *genericH2O  
```

### Class Methods  
None  

### Instance Methods  

🔹 
```
- (void)loadDensityPropertiesWithT:(double)t andWithP:(double)p
```

🔹  T (K), P (bars) => \\(\varepsilon\\)
 ```
- (double)epsilonFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\(\frac{{\partial \varepsilon }}{{\partial T}}\\) ( 1/K)
```
- (double)dEpsilonDtFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\(\frac{{\partial \varepsilon }}{{\partial P}}\\) ( 1/bar)
```
- (double)dEpsilonDpFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\(\frac{{{\partial ^2}\varepsilon }}{{\partial {T^2}}}\\) ( 1/K/K)
```
- (double)d2EpsilonDt2FromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\(\frac{{{\partial ^2}\varepsilon }}{{\partial T\partial P}}\\) ( 1/K/bar)
```
- (double)d2EpsilonDtDpFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\(\frac{{{\partial ^2}\varepsilon }}{{\partial {P^2}}}\\) ( 1/bar/bar)
```
- (double)d2EpsilonDp2FromT:(double)t andP:(double)p
``` 



🔹 T (K), P (bars) => Q Born parameter
```
- (double)QfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => N Born parameter
```
- (double)NfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => U Born parameter
```
- (double)UfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => Y Born parameter
```
- (double)YfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => X Born parameter
```
- (double)XfromT:(double)t andP:(double)p
```

🔹 Numerical derivatives <br/> 

- \\(\frac{{\partial U}}{{\partial T}}\\)  
- \\(\frac{{\partial U}}{{\partial P}}\\)  
- \\(\frac{{\partial N}}{{\partial T}}\\)  
- \\(\frac{{\partial N}}{{\partial P}}\\)  
- \\(\frac{{\partial X}}{{\partial T}}\\)  


```
- (double)dUdTfromT:(double)t andP:(double)p
- (double)dUdPfromT:(double)t andP:(double)p
- (double)dNdTfromT:(double)t andP:(double)p
- (double)dNdPfromT:(double)t andP:(double)p
- (double)dXdTfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\({A_\gamma }\\) Debye-Huckel parameter
```
- (double)AgammaFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\({B_\gamma }\\) Debye-Huckel parameter
```
- (double)BgammaFromT:(double)t andP:(double)p
```

🔹  T (K), P (bars) => *A<sub>G</sub>* Debye-Huckel parameter (J)
```
- (double)AsubGfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => *A<sub>H</sub>* Debye-Huckel parameter (J)
```
- (double)AsubHfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => *A<sub>J</sub>* Debye-Huckel parameter (J)
```
- (double)AsubJfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => *A<sub>V</sub>* Debye-Huckel parameter (J)

```
- (double)AsubVfromT:(double)t andP:(double)p
```

🔹  T (K), P (bars) => \\({A_\kappa }\\) Debye-Huckel parameter (J)
```
- (double)AsubKappaFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => *A<sub>Ex</sub>* Debye-Huckel parameter (J)
```
- (double)AsubExFromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => *B<sub>G</sub>* Debye-Huckel parameter (J)
```
- (double)BsubGfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => *B<sub>H</sub>* Debye-Huckel parameter (J)
```
- (double)BsubHfromT:(double)t andP:(double)p
```

🔹  T (K), P (bars) => *B<sub>J</sub>* Debye-Huckel parameter (J)
 ```
- (double)BsubJfromT:(double)t andP:(double)p
```

 🔹 T (K), P (bars) => *B<sub>V</sub>* Debye-Huckel parameter (J)
```
- (double)BsubVfromT:(double)t andP:(double)p
```

🔹 T (K), P (bars) => \\({B_\kappa }\\) Debye-Huckel parameter (J)
```
- (double)BsubKappaFromT:(double)t andP:(double)p
```

🔹  T (K), P (bars) => *B<sub>Ex</sub>* Debye-Huckel parameter (J)
```
- (double)BsubExFromT:(double)t andP:(double)p

```