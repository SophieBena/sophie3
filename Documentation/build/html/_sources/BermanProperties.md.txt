## BermanProperties Class  
### Class Inheritance  
NSObject ▶️  [PhaseBase](PhaseBase.html) ▶️  BermanProperties  

### Protocols Implemented  
NSSecureCoding  

[StoichiometricPhaseProtocol](StoichiometricPhaseProtocol.html)  

### Properties  

```
@property (readwrite, nonatomic) double h
```

```
@property (readwrite, nonatomic) double s
```

```
@property (readwrite, nonatomic) double k0
```

```
@property (readwrite, nonatomic) double k1
```

```
@property (readwrite, nonatomic) double k2
```

```
@property (readwrite, nonatomic) double k3
```

```
@property (readwrite, nonatomic) double l1
```

```
@property (readwrite, nonatomic) double l2
```

```
@property (readwrite, nonatomic) double Tt
```

```
@property (readwrite, nonatomic) double deltaH
```

```
@property (readwrite, nonatomic) double v0
```

```
@property (readwrite, nonatomic) double v1
```

```
@property (readwrite, nonatomic) double v2
```

```
@property (readwrite, nonatomic) double v3
```

```
@property (readwrite, nonatomic) double v4
```


### Class Methods  
 🔹 These class functions set the behavior of all instances of the class, including those previously instantiated.
 Enables or Disables the convention of using the Gibbs free energy of formation from the elements as the reference state. The default is to use the enthalpy of formation from the elements as the reference state.  
 
- Default, Disabled: 

$$\Delta {G_{T,P}} = \Delta {H_{{T_r},{P_r}}} + \int_{{T_r}}^T {{C_P}} dT - T{S_{{T_r},{P_r}}} - T\int_{{T_r}}^T {\frac{{{C_P}}}{T}} dT + \int_{{P_r}}^P {VdP}$$

- Enabled: 
$$\Delta {G_{T,P}} = \Delta {G_{{T_r},{P_r}}} + \int_{{T_r}}^T {{C_P}} dT - \left( {T - {T_r}} \right){S_{{T_r},{P_r}}} - T\int_{{T_r}}^T {\frac{{{C_P}}}{T}} dT + \int_{{P_r}}^P {VdP}$$

where \\(\Delta {G_{{T_r},{P_r}}} = \Delta {H_{{T_r},{P_r}}} - {T_r}{S_{{T_r},{P_r}}} + {T_r}S_{{T_r},{P_r}}^{elements}\\)

or, equivalently:

- Default, Disabled:  

$$\Delta {G_{T,P}} = \Delta {H_{{T_r},{P_r}}} + \int_{{T_r}}^T {{C_P}} dT - T{S_{{T_r},{P_r}}} - T\int_{{T_r}}^T {\frac{{{C_P}}}{T}} dT + \int_{{P_r}}^P {VdP}$$


- Enabled: 

$$\Delta {G_{T,P}} = \Delta {H_{{T_r},{P_r}}} + \int_{{T_r}}^T {{C_P}} dT - T{S_{{T_r},{P_r}}} - T\int_{{T_r}}^T {\frac{{{C_P}}}{T}} dT + \int_{{P_r}}^P {VdP}  + {T_r}S_{{T_r},{P_r}}^{elements}$$

so, the only difference is the constant \\({T_r}S_{{T_r},{P_r}}^{elements}\\)
 
 ```
 +(void)enableGibbsFreeEnergyReferenceStateUsed

 ```
 ```
 +(void)disableGibbsFreeEnergyReferenceStateUsed
 ```

### Instance Methods  

```
-(id)initWithH:(double)hIn 
			 S:(double)sIn 
			k0:(double)k0In 
			k1:(double)k1In 
			k2:(double)k2In 
			k3:(double)k3In 
			l1:(double)l1In  
			l2:(double)l2In 
			Tt:(double)TtIn 
		deltaH:(double)deltaHIn
			v0:(double)v0In 
			v1:(double)v1In 
			v2:(double)v2In 
			v3:(double)v3In 
			v4:(double)v4In  

```
```
-(id)initWithH:(double)hIn 
			 S:(double)sIn 
			k0:(double)k0In 
			k1:(double)k1In 
			k2:(double)k2In 
			k3:(double)k3In 
			v0:(double)v0In 
			v1:(double)v1In 
			v2:(double)v2In 
			v3:(double)v3In 
			v4:(double)v4In  

```
```
-(id)initWithH:(double)hIn 
			 S:(double)sIn 
			k0:(double)k0In 
			k1:(double)k1In 
			k2:(double)k2In 
			k3:(double)k3In  

```
🔹 Set the reference temperature for the Cp integration (default 298.15 K)

```
-(void)setTr:(double)trIn  

```
🔹  Set the reference pressure (default 1.0 bar)

```
-(void)setPr:(double)prIn

```
🔹  Set the reference temperature for the lambda Cp correction (default 298.15 K)

```
-(void)setTrl:(double)trlIn 
```