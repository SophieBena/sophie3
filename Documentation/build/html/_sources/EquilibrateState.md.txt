## EquilibrateState Class  
### Class Inheritance  
NSObject ▶️  EquilibrateState  

### Protocols Implemented  
None  

### Properties  

```
@property (readwrite, assign) double mass
```
```
@property (readwrite, assign) DoubleVector *bulkCompositionInOxides
```
```
@property (readwrite, assign) DoubleVector *bulkCompositionInElements
```
```
@property (readwrite, strong) NSMutableDictionary *phasesInSystem
```
```
@property (readwrite, strong) NSMutableDictionary *potentialPhaseList
```
```
@property (readwrite, assign) double T
```
```
@property (readwrite, assign) double P
```
```
@property (readwrite, assign) double referenceEntropyOfSystem
```
```
@property (readwrite, assign) double referenceVolumeOfSystem
```
```
@property (readwrite, assign) double referenceEnthalpyOfSystem
```
```
@property (readwrite, assign) double correctionToReferenceEntropyOfSystem
```
```
@property (readwrite, assign) double correctionToReferenceVolumeOfSystem
```
```
@property (readwrite, assign) double fo2
```
```
@property (readwrite, assign) NSUInteger fo2Path
```
```
@property (readwrite, assign) double fo2Delta
```
```
@property (readwrite, assign) double referenceOxygen
```

### Class Methods  
🔹 elements[] += scaler * reference[]. 
 
 Takes an input vector of lenth 107 containig mole numbers of elements indexed on atomic number and adds to that vector a scaler multple of another vector of length 107 which represents a phase stoichiometry.
 
 @param elements  
   double vector of length 107. Input and output.  
 
 @param scaler  
   input only  
 
 @param reference  
   double vector of length 107. Input only.

 ```
+(void)accumulateIntoVectorOfElements:(double *)elements aScaler:(double)scaler timesElementReferenceVector:(double *)reference
``` 
🔹 elements[] *= scaler.  
 
 Takes an input vector of lenth 107 containig mole numbers of elements indexed on atomic number and multiplies each element of that vector by a scaler.  
 
 @param elements  
 input/output double vector of length 107. Input and output.  
 
 @param scaler  
 input only  

 ```
+(void)scaleVectorOfElements:(double *)elements aScaler:(double)scaler
 ```

### Instance Methods  

🔹  Initialize with the number of oxides (system components)  
 
 @param na  
   Number of components in the system (in the usual case the number of oxides)  


```
-(id)initWithComponentNumber:(NSUInteger)na

```

🔹  Determines the number of non-zero entries in the element array instance variable bulkCompositionInElements  
 
 @return  
   number of non-zero entries  

 
 ```
-(NSUInteger)numberOfNonZeroElements
```

 🔹 Generates and returns a table of element indices for non-zero entries in the element array instance variable bulkCompositionInElements  
 
 @return  
   An NSArray of NSNumber objects each containing a NSUInteger index   
```
-(NSArray *)hashTableOfNonZeroElementEntries
```

 🔹  Formats an instance of the EquilibrateState class as an XML document  
   
 @return  
 Pointer to an NSXMLDocument with root element EquilibrateState that contains the current instance   
 
 ```
-(NSXMLElement *)getEquilibrateStateAsXMLElement:(id)delegate selector:(SEL)convertMolesOfElementsToGramsOfOxides

```
