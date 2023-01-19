## HKFspeciesProperties Class  
### Class Inheritance  
NSObject ▶️ [PhaseBase](PhaseBase.html) ▶️ HKFspeciesProperties  

### Protocols Implemented  
NSSecureCoding<br />  
[StoichiometricPhaseProtocol](StoichiometricPhaseProtocol.html)  

### Properties  

```
@property (assign) double deltaGibbsFreeEnergyOfFormationInTheReferenceState

```

```
@property (assign) double deltaEnthalpyOfFormationInTheReferenceState

```

```
@property (assign) double entropyInTheReferenceState

```

```
@property (assign) double volumeInTheReferenceState

```

```
@property (assign) double heatCapacityInTheReferenceState

```

```
@property (assign) double a1HKF

```

```
@property (assign) double a2HKF

```

```
@property (assign) double a3HKF

```

```
@property (assign) double a4HKF

```

```
@property (assign) double c1HKF

```

```
@property (assign) double c2HKF

```

```
@property (assign) double omegaHKF

```

```
@property (assign) double charge

```

### Class Methods  
None  

### Instance Methods  

```

- (instancetype)initWithGibbsFreeEnergyOfFormation:(double)gIn
                        andWithEnthalpyOfFormation:(double)hIn
                           andWithReferenceEntropy:(double)sIn
                            andWithReferenceVolume:(double)vIn
                      andWithReferenceHeatCapacity:(double)cpIn
                             andWithHKFa1Parameter:(double)a1In
                             andWithHKFa2Parameter:(double)a2In
                             andWithHKFa3Parameter:(double)a3In
                             andWithHKFa4Parameter:(double)a4In
                             andWithHKFc1Parameter:(double)c1In
                             andWithHKFc2Parameter:(double)c2In
                          andWithHKFomegaParameter:(double)omegaIn
                                     andWithCharge:(double)chargeIn
```
```
- (void)updateDensityPropertiesWithT:(double)t andP:(double)p
```
```
- (double)gSolventAtT:(double)t andP:(double)p
```
```
- (double)DgSolventDtAtT:(double)t andP:(double)p
```
```
- (double)DgSolventDpAtT:(double)t andP:(double)p
```
```
- (double)D2gSolventDt2AtT:(double)t andP:(double)p
```
```
- (double)D2gSolventDtDpAtT:(double)t andP:(double)p
```
```
- (double)D2gSolventDp2AtT:(double)t andP:(double)p
```
```

- (double)omegaAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z
```
```
- (double)DomegaDtAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z
```
```
- (double)DomegaDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z
```
```
- (double)D2omegaDt2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z
```
```
- (double)D2omegaDtDpAtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z
```
```
- (double)D2omegaDp2AtT:(double)t andP:(double)p withOmegaRef:(double)omegaRef andCharge:(double)z

```