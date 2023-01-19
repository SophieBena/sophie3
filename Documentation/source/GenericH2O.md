## GenericH2O Class  
### Class Inheritance  
NSObject ▶️ [PhaseBase](PhaseBase.html) ▶️ GenericH2O  

### Protocols Implemented  
[StoichiometricPhaseProtocol](StoichiometricPhaseProtocol.html)     

### Properties  
None  

### Class Methods  
None  

### Instance Methods  

```
- (void)forceModeChoiceAutomatic
```
```
- (void)forceModeChoiceTo:(NSString *)region
```
```
- (void)setModeToDuanAndZhang2006
```
```
- (void)setModeToZhangAndDuan2005
```
```
- (void)setModeToHoltenEtAl2014
```
```
- (void)setModeToWagnerEtAl2002
```
```
-(void)useZhangAndDuan2005forDEW
```
```
-(void)useDuanAndZhang2006forDEW
```
```
-(void)useZhangAndDuan2009forDEW
```
```
- (NSString *)eosForRegionAtT:(double)t andP:(double)p
```
```
- (double)lowerTemperatureLimitAtPinBars:(double)p
```
```
- (double)tTransHoltenWagner
```
```
- (double)tTransWagnerDZ2006
```
```
- (double)tTransHoltenZD2005
```
```
- (double)pTransWagnerZD2005
```
```
- (double)pTransDZ2006ZD2005
```
