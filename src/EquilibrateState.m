//
//  EquilibrateState.m
//  PhaseObjC ENKI package (derived from PhasePlot package)
//
//  Originally created by Mark Ghiorso on 6/20/10.
//  Modified by Mark S. Ghiorso on 2/14/17.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "EquilibrateState.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "SolutionPhaseProtocol.h"
#import "StoichiometricPhaseProtocol.h"
#import "PhaseBase.h"
#import "BermanProperties.h"

@implementation EquilibrateStatePhase

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
        bulkCompositionInElements = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		mass = 0.0;
		affinity = 0.0;
	}
	return self;
}

-(double)sumOfMolesOfElements {
	double result = 0.0;
	double *e = [[self bulkCompositionInElements] pointerToDouble];
	for (NSUInteger i=1; i<107; i++) result += e[i];
	return result;
}

-(NSArray *)indicesOfNonZeroMolesOfEndmemberComponents {
	NSMutableArray *result = [NSMutableArray arrayWithCapacity:1];
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		for (NSUInteger i=0; i<na; i++) if (m[i] != 0.0) [result addObject:[NSNumber numberWithUnsignedInteger:i]];
	} else [result addObject:[NSNumber numberWithUnsignedInteger:0]];
	return [NSArray arrayWithArray:result];
}

-(NSArray *)nonZeroMolesOfEndmemberComponents:(NSUInteger)constraint andT:(double)t andP:(double)p {
	NSMutableArray *result = [NSMutableArray arrayWithCapacity:0];
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];

		if ((constraint & independentT) && (constraint & independentP)) {
            DoubleVector *dpWrapper = [phaseClassInstance getDgDmFromMolesOfComponents:m andT:t andP:p];
			double *dp  = [dpWrapper pointerToDouble];
            DoubleMatrix *d2pWrapper = [phaseClassInstance getD2gDm2FromMolesOfComponents:m andT:t andP:p];
			double **d2p = [d2pWrapper pointerToPointerToDouble];
            for (NSUInteger i=0; i<na; i++) {
                if (m[i] != 0.0) {
                    NSMutableArray *hess = [NSMutableArray arrayWithCapacity:0];
                    for (NSUInteger j=0; j<na; j++) if (m[j] != 0.0) [hess addObject:[NSNumber numberWithDouble:d2p[i][j]]];
                    NSArray *object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                       [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                       [NSNumber numberWithDouble:dp[i]],
                                       [NSArray arrayWithArray:hess],
                                       nil];
                    [result addObject:object];
                }
            }
        } else if ((constraint & independentS) && (constraint & independentP)) {
            DoubleVector *dpWrapperG = [phaseClassInstance getDgDmFromMolesOfComponents:m andT:t andP:p];
            double *dpG = [dpWrapperG pointerToDouble];
            DoubleVector *dpWrapperS = [phaseClassInstance getDsDmFromMolesOfComponents:m andT:t andP:p];
            double *dpS = [dpWrapperS pointerToDouble];
            DoubleMatrix *d2pWrapperG = [phaseClassInstance getD2gDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pG = [d2pWrapperG pointerToPointerToDouble];
            DoubleMatrix *d2pWrapperS = [phaseClassInstance getD2sDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pS = [d2pWrapperS pointerToPointerToDouble];

            DoubleVector *dsdmWrapper = [phaseClassInstance getDsDmFromMolesOfComponents:m andT:t andP:p];
            double *dsdm = [dsdmWrapper pointerToDouble];
            DoubleVector *dcpdmWrapper = [phaseClassInstance getDCpDmFromMolesOfComponents:m andT:t andP:p];
            double *dcpdm = [dcpdmWrapper pointerToDouble];
            double cp = [phaseClassInstance getHeatCapacityFromMolesOfComponents:m andT:t andP:p];
            double dcpdt = [phaseClassInstance getDcpDtFromMolesOfComponents:m andT:t andP:p];

            BOOL firstEntry = YES;
            for (NSUInteger i=0; i<na; i++) {
                if (m[i] != 0.0) {
                    NSMutableArray *hess = [NSMutableArray arrayWithCapacity:0];
                    NSMutableArray *hessLambda = [NSMutableArray arrayWithCapacity:0];

                    for (NSUInteger j=0; j<na; j++) if (m[j] != 0.0) {
                        [hess addObject:[NSNumber numberWithDouble:(d2pG[i][j]+t*d2pS[i][j])]];
                        [hessLambda addObject:[NSNumber numberWithDouble:d2pS[i][j]]];
                    }

                    NSArray *object = nil;
                    if (firstEntry) {
                        object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                  [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                  [NSNumber numberWithDouble:(dpG[i]+t*dpS[i])], // 2. Gradient of [i] term
                                  [NSArray arrayWithArray:hess],                 // 3. Hessian row/column of [i] vector
                                  [NSNumber numberWithDouble:dcpdm[i]],          // 4. Hessian dmdp derivative
                                  [NSNumber numberWithDouble:dsdm[i]],           // 5. Constraint Jacobian entry for [i] term
                                  [NSArray arrayWithArray:hessLambda],           // 6. Constraint Hession row/column of [i] vector
                                  [NSNumber numberWithDouble:dcpdm[i]/t],        // 7. Constraint Hession dmdp derivative
                                  [NSNumber numberWithDouble:cp],                // 8. Last row of the gradient vector
                                  [NSNumber numberWithDouble:cp/t],              // 9. Last row/column entry of the constraint Jacobian
                                  [NSNumber numberWithDouble:dcpdt],             // 10. Last row/column entry of the Hessian
                                  [NSNumber numberWithDouble:-cp/t/t+dcpdt/t],   // 11. Last row/column entry of the constraint Hessian
                                  nil];
                        firstEntry = NO;
                    } else {
                        object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                  [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                  [NSNumber numberWithDouble:(dpG[i]+t*dpS[i])], // 2. Gradient of [i] term
                                  [NSArray arrayWithArray:hess],                 // 3. Hessian row/column of [i] vector
                                  [NSNumber numberWithDouble:dcpdm[i]],          // 4. Hessian dmdp derivative
                                  [NSNumber numberWithDouble:dsdm[i]],           // 5. Constraint Jacobian entry for [i] term
                                  [NSArray arrayWithArray:hessLambda],           // 6. Constraint Hession row/column of [i] vector
                                  [NSNumber numberWithDouble:dcpdm[i]/t],        // 7. Constraint Hession dmdp derivative in slot na
                                  nil];
                    }
                    [result addObject:object];
                }
            }
        } else if ((constraint & independentT) && (constraint & independentV)) {
            DoubleVector *dpWrapperG = [phaseClassInstance getDgDmFromMolesOfComponents:m andT:t andP:p];
            double *dpG = [dpWrapperG pointerToDouble];
            DoubleVector *dpWrapperV = [phaseClassInstance getDvDmFromMolesOfComponents:m andT:t andP:p];
            double *dpV = [dpWrapperV pointerToDouble];
            DoubleMatrix *d2pWrapperG = [phaseClassInstance getD2gDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pG = [d2pWrapperG pointerToPointerToDouble];
            DoubleMatrix *d2pWrapperV = [phaseClassInstance getD2vDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pV = [d2pWrapperV pointerToPointerToDouble];

            DoubleVector *dvdmWrapper = [phaseClassInstance getDvDmFromMolesOfComponents:m andT:t andP:p];
            double *dvdm = [dvdmWrapper pointerToDouble];
            DoubleVector *d2vdmdpWrapper = [phaseClassInstance getD2vDmDpFromMolesOfComponents:m andT:t andP:p];
            double *d2vdmdp = [d2vdmdpWrapper pointerToDouble];
            double dvdp = [phaseClassInstance getDvDpFromMolesOfComponents:m andT:t andP:p];
            double d2vdp2 = [phaseClassInstance getD2vDp2FromMolesOfComponents:m andT:t andP:p];

            BOOL firstEntry = YES;
            for (NSUInteger i=0; i<na; i++) {
                if (m[i] != 0.0) {
                    NSMutableArray *hess = [NSMutableArray arrayWithCapacity:0];
                    NSMutableArray *hessLambda = [NSMutableArray arrayWithCapacity:0];

                    for (NSUInteger j=0; j<na; j++) if (m[j] != 0.0) {
                        [hess addObject:[NSNumber numberWithDouble:(d2pG[i][j]-p*d2pV[i][j])]];
                        [hessLambda addObject:[NSNumber numberWithDouble:d2pV[i][j]]];
                    }

                    NSArray *object = nil;
                    if (firstEntry) {
                        object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                  [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                  [NSNumber numberWithDouble:(dpG[i]-p*dpV[i])], // 2. Gradient of [i] term
                                  [NSArray arrayWithArray:hess],                 // 3. Hessian row/column of [i] vector
                                  [NSNumber numberWithDouble:-p*d2vdmdp[i]],     // 4. Hessian dmdp derivative
                                  [NSNumber numberWithDouble:dvdm[i]],           // 5. Constraint Jacobian entry for [i] term
                                  [NSArray arrayWithArray:hessLambda],           // 6. Constraint Hession row/column of [i] vector
                                  [NSNumber numberWithDouble:d2vdmdp[i]],        // 7. Constraint Hession dmdp derivative
                                  [NSNumber numberWithDouble:-p*dvdp],           // 8. Last row of the gradient vector
                                  [NSNumber numberWithDouble:dvdp],              // 9. Last row/column entry of the constraint Jacobian
                                  [NSNumber numberWithDouble:-dvdp-p*d2vdp2],    // 10. Last row/column entry of the Hessian
                                  [NSNumber numberWithDouble:d2vdp2],            // 11. Last row/column entry of the constraint Hessian
                                  nil];
                        firstEntry = NO;
                    } else {
                        object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                  [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                  [NSNumber numberWithDouble:(dpG[i]-p*dpV[i])], // 2. Gradient of [i] term
                                  [NSArray arrayWithArray:hess],                 // 3. Hessian row/column of [i] vector
                                  [NSNumber numberWithDouble:-p*d2vdmdp[i]],     // 4. Hessian dmdp derivative
                                  [NSNumber numberWithDouble:dvdm[i]],           // 5. Constraint Jacobian entry for [i] term
                                  [NSArray arrayWithArray:hessLambda],           // 6. Constraint Hession row/column of [i] vector
                                  [NSNumber numberWithDouble:d2vdmdp[i]],        // 7. Constraint Hession dmdp derivative in slot na
                                  nil];
                    }
                    [result addObject:object];
                }
            }
        } else if ((constraint & independentS) && (constraint & independentV)) {
            DoubleVector *dpWrapperG = [phaseClassInstance getDgDmFromMolesOfComponents:m andT:t andP:p];
            double *dpG = [dpWrapperG pointerToDouble];
            DoubleVector *dpWrapperS = [phaseClassInstance getDsDmFromMolesOfComponents:m andT:t andP:p];
            double *dpS = [dpWrapperS pointerToDouble];
            DoubleVector *dpWrapperV = [phaseClassInstance getDvDmFromMolesOfComponents:m andT:t andP:p];
            double *dpV = [dpWrapperV pointerToDouble];
            DoubleMatrix *d2pWrapperG = [phaseClassInstance getD2gDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pG = [d2pWrapperG pointerToPointerToDouble];
            DoubleMatrix *d2pWrapperS = [phaseClassInstance getD2sDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pS = [d2pWrapperS pointerToPointerToDouble];
            DoubleMatrix *d2pWrapperV = [phaseClassInstance getD2vDm2FromMolesOfComponents:m andT:t andP:p];
            double **d2pV = [d2pWrapperV pointerToPointerToDouble];

            DoubleVector *dsdmWrapper = [phaseClassInstance getDsDmFromMolesOfComponents:m andT:t andP:p];
            double *dsdm = [dsdmWrapper pointerToDouble];
            DoubleVector *dcpdmWrapper = [phaseClassInstance getDCpDmFromMolesOfComponents:m andT:t andP:p];
            double *dcpdm = [dcpdmWrapper pointerToDouble];
            DoubleVector *dvdmWrapper = [phaseClassInstance getDvDmFromMolesOfComponents:m andT:t andP:p];
            double *dvdm = [dvdmWrapper pointerToDouble];
            DoubleVector *d2vdmdtWrapper = [phaseClassInstance getD2vDmDtFromMolesOfComponents:m andT:t andP:p];
            double *d2vdmdt = [d2vdmdtWrapper pointerToDouble];
            DoubleVector *d2vdmdpWrapper = [phaseClassInstance getD2vDmDpFromMolesOfComponents:m andT:t andP:p];
            double *d2vdmdp = [d2vdmdpWrapper pointerToDouble];
            double cp      = [phaseClassInstance getHeatCapacityFromMolesOfComponents:m andT:t andP:p];
            double dcpdt   = [phaseClassInstance getDcpDtFromMolesOfComponents:m andT:t andP:p];
            double dvdt    = [phaseClassInstance getDvDtFromMolesOfComponents:m andT:t andP:p];
            double dvdp    = [phaseClassInstance getDvDpFromMolesOfComponents:m andT:t andP:p];
            double d2vdt2  = [phaseClassInstance getD2vDt2FromMolesOfComponents:m andT:t andP:p];
            double d2vdtdp = [phaseClassInstance getD2vDtDpFromMolesOfComponents:m andT:t andP:p];
            double d2vdp2  = [phaseClassInstance getD2vDp2FromMolesOfComponents:m andT:t andP:p];

            BOOL firstEntry = YES;
            for (NSUInteger i=0; i<na; i++) {
                if (m[i] != 0.0) {
                    NSMutableArray *hess = [NSMutableArray arrayWithCapacity:0];
                    NSMutableArray *hessLambdaS = [NSMutableArray arrayWithCapacity:0];
                    NSMutableArray *hessLambdaV = [NSMutableArray arrayWithCapacity:0];

                    for (NSUInteger j=0; j<na; j++) if (m[j] != 0.0) {
                        [hess addObject:[NSNumber numberWithDouble:(d2pG[i][j]+t*d2pS[i][j]-p*d2pV[i][j])]];
                        [hessLambdaS addObject:[NSNumber numberWithDouble:d2pS[i][j]]];
                        [hessLambdaV addObject:[NSNumber numberWithDouble:d2pV[i][j]]];
                    }

                    NSArray *object = nil;
                    if (firstEntry) {
                        object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                  [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                  [NSNumber numberWithDouble:dpG[i]+t*dpS[i]-p*dpV[i]],    // 2. Gradient of [i] term
                                  [NSArray arrayWithArray:hess],                           // 3. Hessian row/column of [i] vector
                                  [NSArray arrayWithObjects:                               // 4. Hessian dmd* derivatives
                                   [NSNumber numberWithDouble:dcpdm[i]-p*d2vdmdt[i]],      //    dmdt derivative
                                   [NSNumber numberWithDouble:-t*d2vdmdt[i]-p*d2vdmdp[i]], //    dmdp derivative
                                   nil],
                                  [NSArray arrayWithObjects:                               // 5. Constraint Jacobian entry for [i] term
                                   [NSNumber numberWithDouble:dsdm[i]],                    //    ds
                                   [NSNumber numberWithDouble:dvdm[i]],                    //    dv
                                   nil],
                                  [NSArray arrayWithObjects:                               // 6. Constraint Hession row/column of [i] vector
                                   [NSArray arrayWithArray:hessLambdaS],                   //    S
                                   [NSArray arrayWithArray:hessLambdaV],                   //    V
                                   nil],
                                  [NSArray arrayWithObjects:                               // 7. Constraint Hession dmd* derivative
                                   [NSNumber numberWithDouble:dcpdm[i]/t],                 //    S:dmdt
                                   [NSNumber numberWithDouble:-d2vdmdt[i]],                //    S:dmdp
                                   [NSNumber numberWithDouble:d2vdmdt[i]],                 //    V:dmdt
                                   [NSNumber numberWithDouble:d2vdmdp[i]],                 //    V:dmdp
                                   nil],
                                  [NSArray arrayWithObjects:                               // 8. Last two by two rows of the gradient vector
                                   [NSNumber numberWithDouble:cp-p*dvdt],                  //    dt
                                   [NSNumber numberWithDouble:-t*dvdt-p*dvdp],             //    dp
                                   nil],
                                  [NSArray arrayWithObjects:                               // 9. Last two by two row/column of the constraint Jacobian
                                   [NSNumber numberWithDouble:cp/t],                       //    dsdt
                                   [NSNumber numberWithDouble:-dvdt],                      //    dsdp
                                   [NSNumber numberWithDouble:dvdt],                       //    dvdt
                                   [NSNumber numberWithDouble:dvdp],                       //    dvdp
                                   nil],
                                  [NSArray arrayWithObjects:                               // 10. Last two by two row/column of the Hessian
                                   [NSNumber numberWithDouble:dcpdt-p*d2vdt2],             //     dt2
                                   [NSNumber numberWithDouble:-dvdt-t*d2vdt2-p*d2vdtdp],   //     dtdp
                                   [NSNumber numberWithDouble:-dvdp-p*d2vdp2-t*d2vdtdp],   //     dp2
                                   nil],
                                  [NSArray arrayWithObjects:                               // 11. Last two by two row/column of the constraint Hessian
                                   [NSNumber numberWithDouble:-cp/t/t+dcpdt/t],            //     S:dt2
                                   [NSNumber numberWithDouble:-d2vdt2],                    //     S:dtdp
                                   [NSNumber numberWithDouble:-d2vdtdp],                   //     S:dp2
                                   [NSNumber numberWithDouble:d2vdt2],                     //     V:dt2
                                   [NSNumber numberWithDouble:d2vdtdp],                    //     V:dtdp
                                   [NSNumber numberWithDouble:d2vdp2],                     //     V:dp2
                                   nil],
                                  nil];
                        firstEntry = NO;
                    } else {
                        object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:m[i]],
                                  [NSValue valueWithPointer:[[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble]],
                                  [NSNumber numberWithDouble:dpG[i]+t*dpS[i]-p*dpV[i]],    // 2. Gradient of [i] term
                                  [NSArray arrayWithArray:hess],                           // 3. Hessian row/column of [i] vector
                                  [NSArray arrayWithObjects:                               // 4. Hessian dmd* derivatives
                                   [NSNumber numberWithDouble:dcpdm[i]-p*d2vdmdt[i]],      //    dmdt derivative
                                   [NSNumber numberWithDouble:-t*d2vdmdt[i]-p*d2vdmdp[i]], //    dmdp derivative
                                   nil],
                                  [NSArray arrayWithObjects:                               // 5. Constraint Jacobian entry for [i] term
                                   [NSNumber numberWithDouble:dsdm[i]],                    //    ds
                                   [NSNumber numberWithDouble:dvdm[i]],                    //    dv
                                   nil],
                                  [NSArray arrayWithObjects:                               // 6. Constraint Hession row/column of [i] vector
                                   [NSArray arrayWithArray:hessLambdaS],                   //    S
                                   [NSArray arrayWithArray:hessLambdaV],                   //    V
                                   nil],
                                  [NSArray arrayWithObjects:                               // 7. Constraint Hession dmd* derivative
                                   [NSNumber numberWithDouble:dcpdm[i]/t],                 //    S:dmdt
                                   [NSNumber numberWithDouble:-d2vdmdt[i]],                //    S:dmdp
                                   [NSNumber numberWithDouble:d2vdmdt[i]],                 //    V:dmdt
                                   [NSNumber numberWithDouble:d2vdmdp[i]],                 //    V:dmdp
                                   nil],
                                  nil];
                    }
                    [result addObject:object];
                }
            }
        }
	} else {
		double *ePure = [[(PhaseBase *)phaseClassInstance formulaAsElementArray] pointerToDouble];
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		double dp = 0.0;
        NSArray *object = nil;

		if ((constraint & independentT) && (constraint & independentP)) {
            dp = [(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p];
            object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:moles], [NSValue valueWithPointer:ePure], [NSNumber numberWithDouble:dp], nil];
            [result addObject:object];

        } else if ((constraint & independentS) && (constraint & independentP)) {
            double ds = [(BermanProperties *)phaseClassInstance getEntropyFromT:t andP:p];
            dp = [(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] + t*ds;
            double cp = [(BermanProperties *)phaseClassInstance getHeatCapacityFromT:t andP:p];
            double dcpdt = [(BermanProperties *)phaseClassInstance getDcpDtFromT:t andP:p];
            object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:moles], [NSValue valueWithPointer:ePure],
                      [NSNumber numberWithDouble:dp],                           // 2. Gradient of [i] term
                      [NSNumber numberWithDouble:0.0],                          // 3. Hessian row/column of [i] vector
                      [NSNumber numberWithDouble:cp],                           // 4. Hessian dmdp derivative
                      [NSNumber numberWithDouble:ds],                           // 5. Constraint Jacobian entry for [i] term
                      [NSNumber numberWithDouble:0.0],                          // 6. Constraint Hession row/column of [i] vector
                      [NSNumber numberWithDouble:cp/t],                         // 7. Constraint Hession dmdp derivative
                      [NSNumber numberWithDouble:moles*cp],                     // 8. Last row of the gradient vector
                      [NSNumber numberWithDouble:moles*cp/t],                   // 9. Last row/column entry of the constraint Jacobian
                      [NSNumber numberWithDouble:moles*dcpdt],                  // 10. Last row/column entry of the Hessian
                      [NSNumber numberWithDouble:-moles*cp/t/t+moles*dcpdt/t],  // 11. Last row/column entry of the constraint Hessian
                      nil];
            [result addObject:object];

        } else if ((constraint & independentT) && (constraint & independentV)) {
            double dv = [(BermanProperties *)phaseClassInstance getVolumeFromT:t andP:p];
            dp = [(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] - p*dv;
            double dvdp = [(BermanProperties *)phaseClassInstance getDvDpFromT:t andP:p];
            double d2vdp2 = [(BermanProperties *)phaseClassInstance getD2vDp2FromT:t andP:p];
            object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:moles], [NSValue valueWithPointer:ePure],
                      [NSNumber numberWithDouble:dp],                           // 2. gradient, first row
                      [NSNumber numberWithDouble:0.0],                          // 3. Hessian row/column of [i] vector
                      [NSNumber numberWithDouble:-p*dvdp],                      // 4. Hessian dmdp derivative
                      [NSNumber numberWithDouble:dv],                           // 5. Jacobian, last row first column
                      [NSNumber numberWithDouble:0.0],                          // 6. Constraint Hessian, row/column of [i] vector
                      [NSNumber numberWithDouble:dvdp],                         // 7. Constraint Hession dmdp derivative
                      [NSNumber numberWithDouble:-p*moles*dvdp],                // 8. Last row of the gradient vector
                      [NSNumber numberWithDouble:moles*dvdp],                   // 9. Last row/column entry of the constraint Jacobian
                      [NSNumber numberWithDouble:-moles*dvdp - p*moles*d2vdp2], // 10. Last row/column entry of the Hessian
                      [NSNumber numberWithDouble:moles*d2vdp2],                 // 11. Last row/column entry of the constraint Hessian
                      nil];
            [result addObject:object];

        } else if ((constraint & independentS) && (constraint & independentV)) {
            double ds = [(BermanProperties *)phaseClassInstance getEntropyFromT:t andP:p];
            double dv = [(BermanProperties *)phaseClassInstance getVolumeFromT:t andP:p];
            dp = [(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] + t*ds - p*dv;
            double cp      = [(BermanProperties *)phaseClassInstance getHeatCapacityFromT:t andP:p];
            double dcpdt   = [(BermanProperties *)phaseClassInstance getDcpDtFromT:t andP:p];
            double dvdt    = [(BermanProperties *)phaseClassInstance getDvDtFromT:t andP:p];
            double dvdp    = [(BermanProperties *)phaseClassInstance getDvDpFromT:t andP:p];
            double d2vdt2  = [(BermanProperties *)phaseClassInstance getD2vDt2FromT:t andP:p];
            double d2vdtdp = [(BermanProperties *)phaseClassInstance getD2vDtDpFromT:t andP:p];
            double d2vdp2  = [(BermanProperties *)phaseClassInstance getD2vDp2FromT:t andP:p];
            object = [NSArray arrayWithObjects:[NSNumber numberWithDouble:moles], [NSValue valueWithPointer:ePure],
                      [NSNumber numberWithDouble:dp],                          // 2. Gradient of [i] term
                      [NSNumber numberWithDouble:0.0],                         // 3. Hessian row/column of [i] vector
                      [NSArray arrayWithObjects:                               // 4. Hessian dmd* derivatives
                       [NSNumber numberWithDouble:cp-p*dvdt],                  //    dmdt derivative
                       [NSNumber numberWithDouble:-t*dvdt-p*dvdp],             //    dmdp derivative
                       nil],
                      [NSArray arrayWithObjects:                               // 5. Constraint Jacobian entry for [i] term
                       [NSNumber numberWithDouble:ds],                         //    ds
                       [NSNumber numberWithDouble:dv],                         //    dv
                       nil],
                      [NSNumber numberWithDouble:0.0],                         // 6. Constraint Hession row/column of [i] vector
                      [NSArray arrayWithObjects:                               // 7. Constraint Hession dmd* derivative
                       [NSNumber numberWithDouble:cp/t],                       //    S:dmdt
                       [NSNumber numberWithDouble:-dvdt],                      //    S:dmdp
                       [NSNumber numberWithDouble:dvdt],                       //    V:dmdt
                       [NSNumber numberWithDouble:dvdp],                       //    V:dmdp
                       nil],
                      [NSArray arrayWithObjects:                               // 8. Last two by two rows of the gradient vector
                       [NSNumber numberWithDouble:moles*cp-p*moles*dvdt],      //    dt
                       [NSNumber numberWithDouble:-t*moles*dvdt-p*moles*dvdp], //    dp
                       nil],
                      [NSArray arrayWithObjects:                               // 9. Last two by two row/column of the constraint Jacobian
                       [NSNumber numberWithDouble:cp*moles/t],                 //    dsdt
                       [NSNumber numberWithDouble:-dvdt*moles],                //    dsdp
                       [NSNumber numberWithDouble:dvdt*moles],                 //    dvdt
                       [NSNumber numberWithDouble:dvdp*moles],                 //    dvdp
                       nil],
                      [NSArray arrayWithObjects:                               // 10. Last two by two row/column of the Hessian
                       [NSNumber numberWithDouble:moles*dcpdt-p*moles*d2vdt2],                   //     dt2
                       [NSNumber numberWithDouble:-moles*dvdt-t*moles*d2vdt2-p*moles*d2vdtdp],   //     dtdp
                       [NSNumber numberWithDouble:-moles*dvdp-p*moles*d2vdp2-t*moles*d2vdtdp],   //     dp2
                       nil],
                      [NSArray arrayWithObjects:                               // 11. Last two by two row/column of the constraint Hessian
                       [NSNumber numberWithDouble:-moles*cp/t/t+moles*dcpdt/t],//     S:dt2
                       [NSNumber numberWithDouble:-moles*d2vdt2],              //     S:dtdp
                       [NSNumber numberWithDouble:-moles*d2vdtdp],             //     S:dp2
                       [NSNumber numberWithDouble:moles*d2vdt2],               //     V:dt2
                       [NSNumber numberWithDouble:moles*d2vdtdp],              //     V:dtdp
                       [NSNumber numberWithDouble:moles*d2vdp2],               //     V:dp2
                       nil],
                      nil];
            [result addObject:object];
        }

	}
	return [NSArray arrayWithArray:result];
}

-(void)addDeltaMolesToNonZeroMolesOfEndmemberComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		for (NSUInteger i=0, j=0; i<na; i++) {
			if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
			if (![[NSUserDefaults standardUserDefaults] boolForKey:@"PPOPTIONS.FORCE.COMPONENTS"] &&
				(fabs(m[i]) < [[NSUserDefaults standardUserDefaults] doubleForKey:@"PPPARAMETERS.COMPONENT.MINIMUM"])) m[i] = 0.0;
		}
		[self setBulkCompositionInElements:[phaseClassInstance convertMolesToElements:m]];
		[self setMass:[phaseClassInstance convertElementsToTotalMass:[[self bulkCompositionInElements] pointerToDouble]]];
	} else {
        double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		double *eActual = [[self bulkCompositionInElements] pointerToDouble];
		[EquilibrateState scaleVectorOfElements:eActual aScaler:(moles+scaleFractor*deltaMoles[0])/moles];
		[self setMass:[phaseClassInstance convertElementsToMassOfPhase:eActual]];
	}
}

-(double)potentialFunctionFor:(NSUInteger)constraint andT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles
	andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

		if ((constraint & independentT) && (constraint & independentP)) {
			double potential = [phaseClassInstance getGibbsFreeEnergyFromMolesOfComponents:m andT:t andP:p];
			if (isnan(potential)) {
				NSLog(@"nan potential(G) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potential = compositionIsUnacceptable;
			}
			return potential;

		} else if ((constraint & independentS) && (constraint & independentP)) {
            double potentialG = [phaseClassInstance getGibbsFreeEnergyFromMolesOfComponents:m andT:t andP:p];
			if (isnan(potentialG)) {
				NSLog(@"nan potential(G) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialG = compositionIsUnacceptable;
			}
            double potentialS = [phaseClassInstance getEntropyFromMolesOfComponents:m andT:t andP:p];
            if (isnan(potentialS)) {
				NSLog(@"nan potential(S) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialS = compositionIsUnacceptable;
			}
			if ((potentialG != compositionIsUnacceptable) && (potentialS != compositionIsUnacceptable)) return potentialG + t*potentialS;
            else return compositionIsUnacceptable;

        } else if ((constraint & independentT) && (constraint & independentV)) {
            double potentialG = [phaseClassInstance getGibbsFreeEnergyFromMolesOfComponents:m andT:t andP:p];
            if (isnan(potentialG)) {
				NSLog(@"nan potential(G) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialG = compositionIsUnacceptable;
			}
            double potentialV = [phaseClassInstance getVolumeFromMolesOfComponents:m andT:t andP:p];
            if (isnan(potentialV)) {
				NSLog(@"nan potential(V) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialV = compositionIsUnacceptable;
			}
            if ((potentialG != compositionIsUnacceptable) && (potentialV != compositionIsUnacceptable)) return potentialG - p*potentialV;
            else return compositionIsUnacceptable;

        } else if ((constraint & independentS) && (constraint & independentV)) {
            double potentialG = [phaseClassInstance getGibbsFreeEnergyFromMolesOfComponents:m andT:t andP:p];
            if (isnan(potentialG)) {
				NSLog(@"nan potential(G) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialG = compositionIsUnacceptable;
			}
            double potentialS = [phaseClassInstance getEntropyFromMolesOfComponents:m andT:t andP:p];
            if (isnan(potentialS)) {
				NSLog(@"nan potential(S) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialS = compositionIsUnacceptable;
			}
            double potentialV = [phaseClassInstance getVolumeFromMolesOfComponents:m andT:t andP:p];
            if (isnan(potentialV)) {
				NSLog(@"nan potential(V) for %@ detected at t = %f °C p = %f MPa and scaleLength = %g", [[self phaseClassInstance] phaseName], t-273.15, p/10.0, scaleFractor);
				NSLog(@"... Composition:");
				for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
													  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
													  m[i]);
				NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
				potentialV = compositionIsUnacceptable;
			}
            if ((potentialG != compositionIsUnacceptable) && (potentialV != compositionIsUnacceptable)) return potentialG + t*potentialS - p*potentialV;
            else return compositionIsUnacceptable;
        }
	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]] + scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		if      ((constraint & independentT) && (constraint & independentP)) return moles*[ (BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p];
        else if ((constraint & independentS) && (constraint & independentP)) return moles*([(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] + t*[(BermanProperties *)phaseClassInstance getEntropyFromT:t andP:p]);
        else if ((constraint & independentT) && (constraint & independentV)) return moles*([(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] - p*[(BermanProperties *)phaseClassInstance getVolumeFromT:t andP:p]);
        else if ((constraint & independentS) && (constraint & independentV))
            return moles*([(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] + t*[(BermanProperties *)phaseClassInstance getEntropyFromT:t andP:p] - p*[(BermanProperties *)phaseClassInstance getVolumeFromT:t andP:p]);
	}
	return compositionIsUnacceptable;
}

-(double)entropyFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
        if (scaleFractor != 0.0) for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

        double potential = [phaseClassInstance getEntropyFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potential)) {
            NSLog(@"nan entropy for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potential = compositionIsUnacceptable;
        }
        return potential;

	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
        if (scaleFractor != 0.0) moles += scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		return moles*[(BermanProperties *)phaseClassInstance getEntropyFromT:t andP:p];
	}
	return compositionIsUnacceptable;
}

-(double)volumeFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		if (scaleFractor != 0.0) for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

        double potential = [phaseClassInstance getVolumeFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potential)) {
            NSLog(@"nan volume for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potential = compositionIsUnacceptable;
        }
        return potential;

	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		if (scaleFractor != 0.0) moles += scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		return moles*[(BermanProperties *)phaseClassInstance getVolumeFromT:t andP:p];
	}
	return compositionIsUnacceptable;
}

-(double)enthalpyFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		if (scaleFractor != 0.0) for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

        double potentialG = [phaseClassInstance getGibbsFreeEnergyFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potentialG)) {
            NSLog(@"nan Gibbs free energy for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potentialG = compositionIsUnacceptable;
        }
        double potentialS = [phaseClassInstance getEntropyFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potentialS)) {
            NSLog(@"nan entropy for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potentialS = compositionIsUnacceptable;
        }
        if ((potentialG != compositionIsUnacceptable) && (potentialS != compositionIsUnacceptable)) return potentialG + t*potentialS;
        else return compositionIsUnacceptable;

	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		if (scaleFractor != 0.0) moles += scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		return moles*([(BermanProperties *)phaseClassInstance getGibbsFreeEnergyFromT:t andP:p] + t*[(BermanProperties *)phaseClassInstance getEntropyFromT:t andP:p]);
	}
	return compositionIsUnacceptable;
}

-(double)heatCapacityFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		if (scaleFractor != 0.0) for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

        double potential = [phaseClassInstance getHeatCapacityFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potential)) {
            NSLog(@"nan heat capacity for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potential = compositionIsUnacceptable;
        }
        return potential;

	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		if (scaleFractor != 0.0) moles += scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		return moles*[(BermanProperties *)phaseClassInstance getHeatCapacityFromT:t andP:p];
	}
	return compositionIsUnacceptable;
}

-(double)dvdpFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		if (scaleFractor != 0.0) for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

        double potential = [phaseClassInstance getDvDpFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potential)) {
            NSLog(@"nan heat capacity for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potential = compositionIsUnacceptable;
        }
        return potential;

	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		if (scaleFractor != 0.0) moles += scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		return moles*[(BermanProperties *)phaseClassInstance getDvDpFromT:t andP:p];
	}
	return compositionIsUnacceptable;
}

-(double)dvdtFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor {
	if ([[self phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
		NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
        DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[self bulkCompositionInElements] pointerToDouble]];
		double *m = [mWrapper pointerToDouble];
		if (scaleFractor != 0.0) for (NSUInteger i=0, j=0; i<na; i++) if (m[i] != 0.0) m[i] += scaleFractor*deltaMoles[j++];
		if (![phaseClassInstance testPermissibleValuesOfComponents:m]) return compositionIsUnacceptable;

        double potential = [phaseClassInstance getDvDtFromMolesOfComponents:m andT:t andP:p];
        if (isnan(potential)) {
            NSLog(@"nan heat capacity for %@ detected at t = %f °C p = %f MPa", [[self phaseClassInstance] phaseName], t-273.15, p/10.0);
            NSLog(@"... Composition:");
            for (NSUInteger i=0; i<na; i++) NSLog(@"... %@ %g",
                                                  [[[[self phaseClassInstance] componentAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                                                  m[i]);
            NSLog(@"... %@", [[self phaseClassInstance] getFormulaFromMolesOfComponents:m andT:t andP:p]);
            potential = compositionIsUnacceptable;
        }
        return potential;

	} else {
		double moles = [phaseClassInstance convertElementsToMolesOfPhase:[[self bulkCompositionInElements] pointerToDouble]];
		if (scaleFractor != 0.0) moles += scaleFractor*deltaMoles[0];
		if (moles < 0.0) return compositionIsUnacceptable;
		return moles*[(BermanProperties *)phaseClassInstance getDvDtFromT:t andP:p];
	}
	return compositionIsUnacceptable;
}

#pragma mark -
#pragma mark pointer to double getter/setter methods

-(DoubleVector *)bulkCompositionInElements {
    return bulkCompositionInElements;
}

-(void)setBulkCompositionInElements:(DoubleVector *)bulkCompositionInElementsWrapper {
    double *bulkCompositionInElementsIn = [bulkCompositionInElementsWrapper pointerToDouble];
    double *e = [bulkCompositionInElements pointerToDouble];
    for (NSUInteger i=1; i<107; i++) e[i] = bulkCompositionInElementsIn[i];
}

#pragma mark -
#pragma mark synthesized methods

@synthesize phaseClassInstance, mass, affinity;

@end

#pragma mark -

@implementation EquilibrateState

#pragma mark -
#pragma mark class methods

+(void)accumulateIntoVectorOfElements:(double *)elements aScaler:(double)scaler timesElementReferenceVector:(double *)reference {
	for (NSUInteger i=1; i<107; i++) elements[i] += scaler*reference[i];
}

+(void)scaleVectorOfElements:(double *)elements aScaler:(double)scaler {
	for (NSUInteger i=1; i<107; i++) elements[i] *= scaler;
}

#pragma mark -
#pragma mark instance methods

-(id)initWithComponentNumber:(NSUInteger)na {
	if ((self = [super init])) {
        bulkCompositionInOxides = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
        bulkCompositionInElements = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		phasesInSystem = [NSMutableDictionary dictionaryWithCapacity:1];
        nOxides = na;

		T = 0.0;
		P = 0.0;

		referenceEnthalpyOfSystem = 0.0;
		referenceEntropyOfSystem = 0.0;
		referenceVolumeOfSystem = 0.0;

        correctionToReferenceEntropyOfSystem = 0.0;
        correctionToReferenceVolumeOfSystem = 0.0;

		fo2 = 0.0;
		fo2Path = FO2_NONE;
		fo2Delta = 0.0;
		referenceOxygen = 0.0;

		isIsenthalpic = NO;
		isIsentropic = NO;
		isIsochoric = NO;
	}
	return self;
}

-(NSUInteger)numberOfNonZeroElements {
	NSUInteger result = 0;
	double *e = [[self bulkCompositionInElements] pointerToDouble];
	for (NSUInteger i=0; i<107; i++) if (e[i] != 0.0) result++;
	return result;
}

-(NSArray *)hashTableOfNonZeroElementEntries {
	NSMutableArray *result = [NSMutableArray arrayWithCapacity:0];
	double *e = [[self bulkCompositionInElements] pointerToDouble];
	for (NSUInteger i=0; i<107; i++) if (e[i] != 0.0) [result addObject:[NSNumber numberWithUnsignedInteger:i]];
	return [NSArray arrayWithArray:result];
}

-(NSXMLElement *)getEquilibrateStateAsXMLElement:(id)delegate selector:(SEL)convertMolesOfElementsToGramsOfOxides {
	NSXMLElement *root = (NSXMLElement *)[NSXMLNode elementWithName:@"EquilibrateState"];
	[root addChild:[NSXMLNode commentWithStringValue:@"EquilibrateState class object"]];
	[root addChild:[NSXMLNode elementWithName:@"Temperature" stringValue:[NSString localizedStringWithFormat:@"%.2f",    [self T]-273.15]]];
	[root addChild:[NSXMLNode elementWithName:@"Pressure"    stringValue:[NSString localizedStringWithFormat:@"%.4f",    [self P]/10000.0]]];
	[root addChild:[NSXMLNode elementWithName:@"Mass"        stringValue:[NSString localizedStringWithFormat:@"%.20f", [self mass]]]];

	NSXMLElement *bcElements = (NSXMLElement *)[NSXMLNode elementWithName:@"Composition"];
	             [bcElements addChild:[NSXMLNode commentWithStringValue:@"Composition of the system in moles of elements."]];
	             [bcElements addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:@"Elements"]];
	double *moles = [[self bulkCompositionInElements] pointerToDouble];
	NSArray *nonZeroElements = [self hashTableOfNonZeroElementEntries];
	for (NSNumber *n in nonZeroElements) {
		NSXMLElement *entry = (NSXMLElement *)[NSXMLNode elementWithName:@"Element" stringValue:[NSString localizedStringWithFormat:@"%.20f", moles[[n unsignedIntValue]]]];
		             [entry addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:[PhaseBase elementNameFromAtomicNumber:[n unsignedIntValue]]]];
		[bcElements addChild:entry];
	}
	[root addChild:bcElements];

	NSMutableArray *molesTemp = [NSMutableArray arrayWithCapacity:106];
	for (NSUInteger i=1; i<107; i++) [molesTemp addObject:[NSNumber numberWithDouble:moles[i]]];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warc-performSelector-leaks"
	NSArray *gramsOfOxides = (NSArray *)[delegate performSelector:convertMolesOfElementsToGramsOfOxides withObject:molesTemp];
#pragma clang diagnostic pop

	NSXMLElement *bcOxides = (NSXMLElement *)[NSXMLNode elementWithName:@"Composition"];
	[bcOxides addChild:[NSXMLNode commentWithStringValue:@"Composition of the system in grams of oxides."]];
	[bcOxides addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:@"Oxides"]];

	NSMutableSet *namesOfOxidesInSystem = [NSMutableSet setWithCapacity:0];

	for (NSArray *object in gramsOfOxides) {
		NSString *key = [object objectAtIndex:0];
		double grams = [[object objectAtIndex:1] doubleValue];
	    if (grams != 0.0) {
			NSXMLElement *entry = (NSXMLElement *)[NSXMLNode elementWithName:@"Oxide" stringValue:[NSString localizedStringWithFormat:@"%.20f", grams]];
			[entry addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:key]];
			[namesOfOxidesInSystem addObject:key];
			[bcOxides addChild:entry];
	    }
	}
	[root addChild:bcOxides];

	NSXMLElement *system = (NSXMLElement *)[NSXMLNode elementWithName:@"System"];
	             [system addChild:[NSXMLNode commentWithStringValue:@"Phases present in the system."]];
	for (NSString *key in [[self phasesInSystem] allKeys]) {
		EquilibrateStatePhase *phaseWrapper = [[self phasesInSystem] objectForKey:key];
		id phase = [phaseWrapper phaseClassInstance];
		NSXMLElement *entry = (NSXMLElement *)[NSXMLNode elementWithName:@"Phase"];
		             [entry addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:[phase phaseName]]];
		[entry addChild:[NSXMLNode elementWithName:@"Mass" stringValue:[NSString localizedStringWithFormat:@"%.20f", [phaseWrapper mass]]]];

		if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
			// endmember component mole fractions
			NSUInteger na = [phase numberOfSolutionComponents];
			moles = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
            DoubleVector *mWrapper = [phase convertElementsToMoles:moles];
			double *m = [mWrapper pointerToDouble];
			double mTotal = [phase convertElementsToTotalMoles:moles];

			[entry addChild:[NSXMLNode elementWithName:@"Formula" stringValue:[phase getFormulaFromMolesOfComponents:m andT:[self T] andP:[self P]]]];
            [entry addChild:[NSXMLNode elementWithName:@"GibbsFreeEnergy"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getGibbsFreeEnergyFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"Enthalpy"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getEnthalpyFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"Entropy"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getEntropyFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"HeatCapacity"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getHeatCapacityFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"DcpDt"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getDcpDtFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"Volume"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getVolumeFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"DvDt"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getDvDtFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"DvDp"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getDvDpFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"D2vDt2"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getD2vDt2FromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"D2vDtDp"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getD2vDtDpFromMolesOfComponents:m andT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"D2vDp2"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        [phase getD2vDp2FromMolesOfComponents:m andT:[self T] andP:[self P]]]]];

            DoubleVector *muWrapper = [phase chemicalPotentialsOfSpeciesFromMolesOfComponents:m andT:[self T] andP:[self P]];
            double *mu = [muWrapper pointerToDouble];
			for (NSUInteger i=0; i<na; i++) {
				NSString *componentName = [[phase componentAtIndex:i] phaseName];
				double X = (mTotal != 0.0) ? m[i]/mTotal : 0.0;
				if (X != 0.0) {
					NSXMLElement *component = (NSXMLElement *)[NSXMLNode elementWithName:@"Component" stringValue:[NSString localizedStringWithFormat:@"%.20f", X]];
								 [component addAttribute:[NSXMLNode attributeWithName:@"Name" stringValue:componentName]];
					[entry addChild:component];
                    NSXMLElement *potential = (NSXMLElement *)[NSXMLNode elementWithName:@"ChemicalPotential" stringValue:[NSString localizedStringWithFormat:@"%.20f", mu[i]]];
                    [potential addAttribute:[NSXMLNode attributeWithName:@"Name" stringValue:componentName]];
                    [entry addChild:potential];
                    id endmember = [phase componentAtIndex:i];
                    double gEndmember = [endmember getGibbsFreeEnergyFromT:[self T] andP:[self P]];
                    NSXMLElement *excess = (NSXMLElement *)[NSXMLNode elementWithName:@"ExcessChemicalPotential" stringValue:[NSString localizedStringWithFormat:@"%.20f", mu[i]-gEndmember]];
                    [excess addAttribute:[NSXMLNode attributeWithName:@"Name" stringValue:componentName]];
                    [entry addChild:excess];
				}
			}

			[molesTemp removeAllObjects];
			for (NSUInteger i=1; i<107; i++) [molesTemp addObject:[NSNumber numberWithDouble:moles[i]]];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warc-performSelector-leaks"
			gramsOfOxides = (NSArray *)[delegate performSelector:convertMolesOfElementsToGramsOfOxides withObject:molesTemp];
#pragma clang diagnostic pop
			double total = 0.0;
			for (NSArray *object in gramsOfOxides) total += [[object objectAtIndex:1] doubleValue];
			if (total == 0.0) total = 100.0;
			for (NSArray *object in gramsOfOxides) {
				NSString *key = [object objectAtIndex:0];
				double grams = [[object objectAtIndex:1] doubleValue];
				if ([namesOfOxidesInSystem member:key]) {
					NSXMLElement *oxide = (NSXMLElement *)[NSXMLNode elementWithName:@"Oxide" stringValue:[NSString localizedStringWithFormat:@"%.20f", 100.0*grams/total]];
					[oxide addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:key]];
					[entry addChild:oxide];
				}
			}

		} else {
			moles = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
            double molesOfPhase = [phase convertElementsToMolesOfPhase:moles];
			[entry addChild:[NSXMLNode elementWithName:@"Formula" stringValue:[phase phaseFormula]]];
            [entry addChild:[NSXMLNode elementWithName:@"GibbsFreeEnergy"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getGibbsFreeEnergyFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"Enthalpy"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getEnthalpyFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"Entropy"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getEntropyFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"HeatCapacity"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getHeatCapacityFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"DcpDt"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getDcpDtFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"Volume"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getVolumeFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"DvDt"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getDvDtFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"DvDp"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getDvDpFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"D2vDt2"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getD2vDt2FromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"D2vDtDp"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getD2vDtDpFromT:[self T] andP:[self P]]]]];
            [entry addChild:[NSXMLNode elementWithName:@"D2vDp2"
                                           stringValue:[NSString localizedStringWithFormat:@"%.20f",
                                                        molesOfPhase*[phase getD2vDp2FromT:[self T] andP:[self P]]]]];
			[molesTemp removeAllObjects];
			for (NSUInteger i=1; i<107; i++) [molesTemp addObject:[NSNumber numberWithDouble:moles[i]]];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warc-performSelector-leaks"
			gramsOfOxides = (NSArray *)[delegate performSelector:convertMolesOfElementsToGramsOfOxides withObject:molesTemp];
#pragma clang diagnostic pop
			double total = 0.0;
			for (NSArray *object in gramsOfOxides) total += [[object objectAtIndex:1] doubleValue];
			if (total == 0.0) total = 100.0;
			for (NSArray *object in gramsOfOxides) {
				NSString *key = [object objectAtIndex:0];
				double grams = [[object objectAtIndex:1] doubleValue];
				if ([namesOfOxidesInSystem member:key]) {
					NSXMLElement *oxide = (NSXMLElement *)[NSXMLNode elementWithName:@"Oxide" stringValue:[NSString localizedStringWithFormat:@"%.20f", 100.0*grams/total]];
					[oxide addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:key]];
					[entry addChild:oxide];
				}
			}
		}


		[system addChild:entry];
	}
	[root addChild:system];

    NSXMLElement *potential = (NSXMLElement *)[NSXMLNode elementWithName:@"Potential"];
    [potential addChild:[NSXMLNode commentWithStringValue:@"Phases not present in the system."]];
    for (EquilibrateStatePhase *phaseWrapper in [[self potentialPhaseList] allValues]) {
        id phase = [phaseWrapper phaseClassInstance];
        NSXMLElement *entry = (NSXMLElement *)[NSXMLNode elementWithName:@"Phase"];
        [entry addAttribute:[NSXMLNode attributeWithName:@"Type" stringValue:[phase phaseName]]];

        [entry addChild:[NSXMLNode elementWithName:@"Affinity"
                                       stringValue:[NSString localizedStringWithFormat:@"%.20f", [phaseWrapper affinity]]]];

        if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
            double *m = [[phase convertElementsToMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]] pointerToDouble];
            [entry addChild:[NSXMLNode elementWithName:@"Formula"
                                           stringValue:[phase getFormulaFromMolesOfComponents:m
                                                                                         andT:[self T]
                                                                                         andP:[self P]]]];
        } else {
            [entry addChild:[NSXMLNode elementWithName:@"Formula" stringValue:[phase phaseFormula]]];
        }

        [potential addChild:entry];
    }
    [root addChild:potential];

	return root;
}

#ifdef NEVER
// Assemble table data - all entries are filled (either by phases in or potentially in the system)
tableRecords = [NSMutableArray arrayWithCapacity:numberOfPhases];
for (NSUInteger i=0; i<numberOfPhases; i++) [tableRecords addObject:[NSNull	null]];

for (MeltsStatePhase *phaseWrapper in [[meltsState phasesInSystem] allValues]) {
    NSUInteger index = [phaseNames indexOfObject:[[phaseWrapper phaseClassInstance] phaseName]];
    id phase = [phaseWrapper phaseClassInstance];

    NSMutableDictionary *record = [NSMutableDictionary dictionaryWithCapacity:3];
    [record setObject:[phase phaseName] forKey:@"phase"];
    [record setObject:[NSString localizedStringWithFormat:@"%.2f",[phaseWrapper mass]] forKey:@"massAndAffinity"];

    if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
        [record setObject:[phase getFormulaFromMolesOfComponents:[[phase convertElementsToMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]] pointerToDouble]
                                                            andT:[meltsState T]
                                                            andP:[meltsState P]] forKey:@"formula"];
    } else {
        [record setObject:[phase phaseFormula] forKey:@"formula"];
    }

    [tableRecords replaceObjectAtIndex:index withObject:record];
}

for (MeltsStatePhase *phaseWrapper in [[meltsState potentialPhaseList] allValues]) {
    NSUInteger index = [phaseNames indexOfObject:[[phaseWrapper phaseClassInstance] phaseName]];
    id phase = [phaseWrapper phaseClassInstance];

    NSMutableDictionary *record = [NSMutableDictionary dictionaryWithCapacity:3];
    [record setObject:[phase phaseName] forKey:@"phase"];
    [record setObject:[NSString localizedStringWithFormat:@"(%.2f)",[phaseWrapper affinity]] forKey:@"massAndAffinity"];

    if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
        [record setObject:[phase getFormulaFromMolesOfComponents:[[phase convertElementsToMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]] pointerToDouble]
                                                            andT:[meltsState T]
                                                            andP:[meltsState P]] forKey:@"formula"];
    } else {
        [record setObject:[phase phaseFormula] forKey:@"formula"];
    }

    [tableRecords replaceObjectAtIndex:index withObject:record];
}

// Assemble browser data
browserTree = [[PopoverBrowserDataTree alloc] initWithDisplayName:@"System"];
NSMutableArray *children = [NSMutableArray arrayWithCapacity:0];

[self addNodeToChildArray:children withName:@"Temperature" andValue:[meltsState T]-273.15 andFormat:@"%.2f (°C)"];
[self addNodeToChildArray:children withName:@"Pressure"    andValue:[meltsState P]/10.0   andFormat:@"%.1f (MPa)"];
[self addNodeToChildArray:children withName:@"Mass"        andValue:[meltsState mass]     andFormat:@"%.2f (g)"];
for (NSString *key in orderedListOfPhasesInSystem) {
    MeltsStatePhase *phaseWrapper = [[meltsState phasesInSystem] objectForKey:key];
    id phase = [phaseWrapper phaseClassInstance];

    PopoverBrowserDataTree *phaseNode = [[PopoverBrowserDataTree alloc] initWithDisplayName:[phase phaseName]];
    [children addObject:phaseNode];

    NSMutableArray *phaseChildren = [NSMutableArray arrayWithCapacity:0];
    [self addNodeToChildArray:phaseChildren withName:@"Mass" andValue:[phaseWrapper mass] andFormat:@"%.2f (g)"];

    if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
        // endmember component mole fractions
        NSUInteger na = [phase numberOfSolutionComponents];
        DoubleVector *mWrapper = [phase convertElementsToMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]];
        double *m = [mWrapper pointerToDouble];
        double mTotal = [phase convertElementsToTotalMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]];

        NSMutableArray *componentChildren = [NSMutableArray arrayWithCapacity:0];
        PopoverBrowserDataTree *componentNode = [[PopoverBrowserDataTree alloc] initWithDisplayName:@"Component mole fractions"];
        [phaseChildren addObject:componentNode];

        for (NSUInteger i=0; i<na; i++) {
            NSString *componentName = [[phase componentAtIndex:i] phaseName];
            double X = (mTotal != 0.0) ? m[i]/mTotal : 0.0;
            if (X != 0.0) [self addNodeToChildArray:componentChildren withName:componentName andValue:X andFormat:@"%.6f"];
        }

        [componentNode setChildren:componentChildren];

        // oxide wt %
        NSArray *oxideNames = [Melts oxideNames];
        na = [oxideNames count];
        DoubleVector *gWrapper = [Melts gramsOfOxidesFromMolesOfElements:[[phaseWrapper bulkCompositionInElements] pointerToDouble]];
        double *g = [gWrapper pointerToDouble];
        double mass = [phaseWrapper mass];

        NSMutableArray *oxideChildren = [NSMutableArray arrayWithCapacity:na];
        PopoverBrowserDataTree *oxideNode = [[PopoverBrowserDataTree alloc] initWithDisplayName:@"Oxide Wt%"];
        [phaseChildren addObject:oxideNode];

        for (NSUInteger i=0; i<na; i++) {
            double value = (mass != 0.0) ? 100.0*g[i]/mass : 0.0;
            if (value != 0.0) [self addNodeToChildArray:oxideChildren withName:[oxideNames objectAtIndex:i] andValue:value andFormat:@"%.3f"];
        }

        [oxideNode setChildren:oxideChildren];
    }

    [phaseNode setChildren:phaseChildren];
}
#endif

#pragma mark -
#pragma mark pointer to double getter/setter methods

-(DoubleVector *)bulkCompositionInElements {
    return bulkCompositionInElements;
}

-(void)setBulkCompositionInElements:(DoubleVector *)bulkCompositionInElementsWrapper {
    double *bulkCompositionInElementsIn = [bulkCompositionInElementsWrapper pointerToDouble];
    double *e = [bulkCompositionInElements pointerToDouble];
    for (NSUInteger i=1; i<107; i++) e[i] = bulkCompositionInElementsIn[i];
}

-(DoubleVector *)bulkCompositionInOxides {
    return bulkCompositionInOxides;
}

-(void)setBulkCompositionInOxides:(DoubleVector *)bulkCompositionInOxidesWrapper {
    double *bulkCompositionInOxidesIn = [bulkCompositionInOxidesWrapper pointerToDouble];
    double *m = [bulkCompositionInOxides pointerToDouble];
    for (NSUInteger i=0; i<nOxides; i++) m[i] = bulkCompositionInOxidesIn[i];
}

#pragma mark -
#pragma mark synthesized methods

@synthesize mass, phasesInSystem, potentialPhaseList, T, P;
@synthesize referenceEntropyOfSystem, referenceVolumeOfSystem, referenceEnthalpyOfSystem, referenceOxygen;
@synthesize fo2, fo2Path, fo2Delta;
@synthesize correctionToReferenceEntropyOfSystem, correctionToReferenceVolumeOfSystem;

@end
