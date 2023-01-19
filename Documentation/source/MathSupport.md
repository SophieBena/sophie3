## MathSupport Class  
### Class Inheritance  
NSObject ▶️ MathSupport  

### Protocols Implemented  
None  

### Properties  
None  

### Class Methods  
```
+(void)hfti:(double **)a m:(NSInteger)m n:(NSInteger)n b:(double **)b nb:(NSInteger)nb tau:(double)tau k:(NSInteger *)k
	  rnorm:(double *)rnorm h:(double *)h g:(double *)g p:(NSInteger *)p 
```
```
+(void)householderColCol:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vCol:(NSInteger)vCol 
					   h:(double *)h c:(double **)c cColStart:(NSInteger)cColStart cColEnd:(NSInteger)cColEnd

```
```
+(void)householderColRow:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vCol:(NSInteger)vCol 
					   h:(double *)h c:(double **)c cRowStart:(NSInteger)cRowStart cRowEnd:(NSInteger)cRowEnd
```
```
+(void)householderRowCol:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vRow:(NSInteger)vRow 
					   h:(double *)h c:(double **)c cColStart:(NSInteger)cColStart cColEnd:(NSInteger)cColEnd
```
```
+(void)householderRowRow:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vRow:(NSInteger)vRow 
					   h:(double *)h c:(double **)c cRowStart:(NSInteger)cRowStart cRowEnd:(NSInteger)cRowEnd
```
```
+(NSInteger)min1d:(double *)bb st:(double *)st reltest:(double)reltest ifn:(NSInteger *)ifn fnminval:(double *)fnminval fn1d:(double (^)(double bb, NSInteger *notcomp))fn1d
```
```
+(NSInteger)modmrt:(NSInteger)n m:(NSInteger)m Bvec:(double *)Bvec Fmin:(double *)Fmin reltest:(double)reltest maxIter:(NSInteger)maxIter 
	   nlres:(double (^)(NSInteger i, NSInteger n, double *Bvec, NSInteger *notcomp))nlres nljac:(void (^)(NSInteger i, NSInteger n, double *Bvec, double *X))nljac
```
```	   
+(NSInteger)rqmcg:(NSInteger)n A:(double **)A B:(double **)B X:(double *)X ipr:(NSInteger *)ipr rq:(double *)rq
```
```
+(NSInteger)vmmin:(NSInteger)n Bvec:(double *)Bvec Fmin:(double *)Fmin reltest:(double)reltest 
	 fminfn:(double (^)(NSInteger n, double *Bvec, NSInteger *notcomp))fminfn fmingr:(void (^)(NSInteger n, double *Bvec, double *g))fmingr
```
```
+(void)G1withA:(double)a andB:(double)b andCterm:(double *)cterm andSterm:(double *)sterm andSig:(double *)sig
```
```
+(void)G2withCterm:(double)cterm andSterm:(double)sterm andX:(double *)x andY:(double *)y
```
```
+(BOOL)nnlsWithConstraintMatrix:(double **)a 
				andRowDimension:(NSInteger)m
			 andColumnDimension:(NSInteger)n
			andConstraintVector:(double *)b
			  andSolutionVector:(double *)x
				andResidualNorm:(double *)rnorm
		  andDualSolutionVector:(double *)w
		 andDoubleWorkingVector:(double *)zz
	   andUIntegerWorkingVector:(NSInteger *)index
	   
```
```	   
+(void)simplx:(double **)a m:(NSInteger)m n:(NSInteger)n m1:(NSInteger)m1 m2:(NSInteger)m2 m3:(NSInteger)m3 icase:(NSInteger *)icase izrov:(NSInteger *)izrov iposv:(NSInteger *)iposv
```
```
+(BOOL)choldc:(double **)a dimension:(NSUInteger)n
```
```
+(void)polint:(double *)xa ya:(double *)ya n:(NSInteger)n x:(double)x y:(double *)y dy:(double *)dy
```
```
+(void)polin2:(double *)x1a x2a:(double *)x2a ya:(double **)ya m:(NSInteger)m n:(NSInteger)n x1:(double)x1 x2:(double)x2 y:(double *)y dy:(double *)dy
```
```
+(BOOL)ldpWithConstraintMatrix:(double **)g
               andRowDimension:(NSInteger)m
            andColumnDimension:(NSInteger)n
           andConstraintVector:(double *)h
             andSolutionVector:(double *)x
         andSolutionVectorNorm:(double *)xnorm
        andDoubleWorkingVector:(double *)w
        andDoubleWorkingMatrix:(double **)e
       andIntegerWorkingVector:(NSInteger *)index
	   
```

### Instance Methods  
None  