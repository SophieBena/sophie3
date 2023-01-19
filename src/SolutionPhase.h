/*
 *  SolutionPhase.h
 *  PhasePlot
 *
 *  Created by Mark Ghiorso on 6/30/10.
 *  Copyright 2010 OFM Research Inc.. All rights reserved.
 *
 */

/**
 Support functions for implementation of the SolutionPhase Protocol

 NA must be defined as a macro.  It represents the total number of endmember components.
 NR must be defined as a macro.  It represents the total number of independent internal composition variables.
 */

// Support functions

-(void)convertIntensiveToExtensiveGradients:(double)pMix          // Input: Mixing energy
									  dpMix:(double [NR])dpMix    // Input: Mixing energy gradient (internal composition variables)
                                     mToTal:(double)mTotal        // Input: Total moles
                                       drdm:(double [NR][NA])drdm // Input: Conversion matrix, dr[NR]/dm[NA]
                                         dp:(double [NA])dp       // Output: Mixing energy gradient (endmember moles)
{
	for (NSUInteger j=0; j<NA; j++) {
		dp[j] = 0.0;
		for (NSUInteger i=0; i<NR; i++) dp[j] += dpMix[i]*drdm[i][j];
			dp[j] = mTotal*dp[j] + pMix;
			}
}

-(void)convertIntensiveToExtensiveHessian:(double)pMix
									dpMix:(double [NR])dpMix
                                   d2pMix:(double [NR][NR])d2pMix
                                   mTotal:(double)mTotal
                                     drdm:(double [NR][NA])drdm
                                   d2rdm2:(double [NR][NA][NA])d2rdm2
                                      d2p:(double **)d2p
{
	for (NSUInteger j=0; j<NA; j++) {
		for (NSUInteger l=0; l<NA; l++) {
			d2p[j][l] = 0.0;
			for (NSUInteger i=0; i<NR; i++) {
				double temp = 0.0;
				for (NSUInteger k=0; k<NR; k++) temp += d2pMix[i][k]*drdm[k][l];
					d2p[j][l] += dpMix[i]*(drdm[i][j] + drdm[i][l] + mTotal*d2rdm2[i][j][l]) + mTotal*drdm[i][j]*temp;
					}
		}
	}
}

-(void)convertIntensiveToExtensiveTensor:(double)pMix
								   dpMix:(double [NR])dpMix
                                  d2pMix:(double [NR][NR])d2pMix
                                  d3pMix:(double [NR][NR][NR])d3pMix
                                  mTotal:(double)mTotal
                                    drdm:(double [NR][NA])drdm
                                  d2rdm2:(double [NR][NA][NA])d2rdm2
                                  d3rdm3:(double [NR][NA][NA][NA])d3rdm3
                                     d3p:(double ***)d3p
{
	for (NSUInteger j=0; j<NA; j++) {
		for (NSUInteger l=0; l<NA; l++) {
			for (NSUInteger v=0; v<NA; v++) {
				d3p[j][l][v] = 0.0;

				for (NSUInteger i=0; i<NR; i++) {
					d3p[j][l][v] = dpMix[i]*(d2rdm2[i][l][v]+d2rdm2[i][j][l]+d2rdm2[i][j][v]+mTotal*d3rdm3[i][j][l][v]);
					for (NSUInteger k=0; k<NR; k++) {
						d3p[j][l][v] += d2pMix[i][k]*(drdm[i][l]*drdm[k][v] + drdm[i][j]*drdm[k][v] + drdm[i][j]*drdm[k][l]
													  + mTotal*(d2rdm2[i][j][l]*drdm[k][v] + d2rdm2[i][j][v]*drdm[k][l] + d2rdm2[k][l][v]*drdm[i][j]));
						for (NSUInteger w=0; w<NR; w++) d3p[j][l][v] += mTotal*d3pMix[i][k][w]*drdm[i][j]*drdm[w][v]*drdm[k][l];
							}
				}

			}
		}
	}
}

// Begin objective-C <SolutionPhaseProtocol> functions

// --> SolutionPhaseProtocol public function
-(void)setResultsToMixingQuantities:(BOOL)yesForMixing {
	computeMixingQuantities = yesForMixing;
}

// --> SolutionPhaseProtocol public function
-(NSUInteger)numberOfSolutionComponents {
	return NA;
}

// --> SolutionPhaseProtocol public function
-(id)componentAtIndex:(NSUInteger)index {
	return [endmembers objectAtIndex:index];
}

// Test routine for permissible component numbers

// --> SolutionPhaseProtocol public function
-(BOOL)testPermissibleValuesOfInternalVariables:(double	*)r {
	return [self test:FIFTH t:0.0 p:0.0 na:0 nr:0 names:NULL formulas:NULL r:r m:NULL];
}

// --> SolutionPhaseProtocol public function
-(BOOL)testPermissibleValuesOfComponents:(double *)m {
	return [self test:SIXTH t:0.0 p:0.0 na:0 nr:0 names:NULL formulas:NULL r:NULL m:m];
}

// Component conversion routines

// --> SolutionPhaseProtocol public function
-(DoubleVector *)convertElementsToMoles:(double *)e {
    DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *m = [mWrapper pointerToDouble];
	[self convert:FIRST outMask:SECOND t:0.0 p:0.0 e:e m:m r:NULL x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	return mWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)convertElementsToTotalMoles:(double *)e {
	DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:NA];
    double *m = [mWrapper pointerToDouble];
	[self convert:FIRST outMask:SECOND t:0.0 p:0.0 e:e m:m r:NULL x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	double result = 0.0;
	for (NSUInteger i=0; i<NA; i++) result += m[i];
    return result;
}

// --> SolutionPhaseProtocol public function
-(double)convertElementsToTotalMass:(double *)e {
	DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:NA];
    double *m = [mWrapper pointerToDouble];
	[self convert:FIRST outMask:SECOND t:0.0 p:0.0 e:e m:m r:NULL x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	double result = 0.0;
	for (NSUInteger i=0; i<NA; i++) result += m[i]*[[endmembers objectAtIndex:i] mw];
    return result;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)convertMolesToMoleFractions:(double *)m {
    DoubleVector *xWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *x = [xWrapper pointerToDouble];
	[self convert:SECOND outMask:FOURTH t:0.0 p:0.0 e:NULL m:m r:NULL x:x dm:NULL d2m:NULL dr:NULL d3m:NULL];
	return xWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)convertMolesToElements:(double *)m {
    DoubleVector *eWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
	double *e = [eWrapper pointerToDouble];
    for (NSUInteger i=0; i<NA; i++) {
		double *componentToElements = [[(BermanProperties *)[endmembers objectAtIndex:i] formulaAsElementArray] pointerToDouble];
		for (NSUInteger j=1; j<107; j++) e[j] += m[i]*componentToElements[j];
    }
	return eWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)totalMolesFromMolesOfComponents:(double *)m {
	double result = 0.0;
	for (int i=0; i<NA; i++) result += m[i];
		return result;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getActivityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *activityWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *activity = [activityWrapper pointerToDouble];
	double r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self activity:FIRST t:t p:p r:r a:activity mu:NULL dx:NULL];
	return activityWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getChemicalPotentialFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *muWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *mu = [muWrapper pointerToDouble];
	double r[NA];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self activity:SECOND t:t p:p r:r a:NULL mu:mu dx:NULL];

	if (computeMixingQuantities) return muWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		mu[i] += [component getGibbsFreeEnergyFromT:t andP:p];
	}

	return muWrapper;
}

// --> SolutionPhaseProtocol public function
//-(DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
//    DoubleMatrix *dadmWrapper = [[DoubleMatrix alloc] initWithRowSize:NA andWithColumnSize:NA];
//	double **dadm = [dadmWrapper pointerToPointerToDouble];
//	double activity[NA], dadr[NA][NR], r[NR], drdm[NR][NA];
//	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
//	[self activity:FIRST | THIRD t:t p:p r:r a:activity mu:NULL dx:dadr];
//	for (NSUInteger i=0; i<NA; i++) {
//		[self convertIntensiveToExtensiveGradients:activity[i] dpMix:dadr[i] mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dadm[i]];
//	}
//	return dadmWrapper;
//}
-(DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *dadmWrapper = [[DoubleMatrix alloc] initWithRowSize:NA andWithColumnSize:NA];
    double **dadm = [dadmWrapper pointerToPointerToDouble];
    double activity[NA], dadr[NA][NR], r[NR], drdm[NR][NA];
    [self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
    [self activity:FIRST | THIRD t:t p:p r:r a:activity mu:NULL dx:dadr];
    for (NSUInteger i=0; i<NA; i++) {
        for (NSUInteger j=0; j<NA; j++) {
            dadm[i][j] = 0.0;
            for (NSUInteger k=0; k<NR; k++) {
                dadm[i][j] += dadr[i][k]*drdm[k][j];
            }
        }
    }
    return dadmWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getGibbsFreeEnergyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double g, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self gmix:FIRST t:t p:p r:r gmix:&g dx:NULL dx2:NULL dx3:NULL];
	g *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return g;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		g += m[i]*[component getGibbsFreeEnergyFromT:t andP:p];
	}
	return g;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDgDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dgdmWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *dgdm = [dgdmWrapper pointerToDouble];
	double g, dgdr[NR], r[NR], drdm[NR][NA];
	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
	[self gmix:FIRST | SECOND t:t p:p r:r gmix:&g dx:dgdr dx2:NULL dx3:NULL];
	[self convertIntensiveToExtensiveGradients:g dpMix:dgdr mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dgdm];

	if (computeMixingQuantities) return dgdmWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dgdm[i] += [component getGibbsFreeEnergyFromT:t andP:p];
	}
	return dgdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getD2gDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *d2gdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:NA andWithColumnSize:NA];
	double **d2gdm2 = [d2gdm2Wrapper pointerToPointerToDouble];
	double g, dgdr[NR], d2gdr2[NR][NR], r[NR], drdm[NR][NA], d2rdm2[NR][NA][NA];
	[self convert:SECOND outMask:THIRD | FIFTH | SIXTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:d2rdm2 dr:NULL d3m:NULL];
	[self gmix:FIRST | SECOND | THIRD t:t p:p r:r gmix:&g dx:dgdr dx2:d2gdr2 dx3:NULL];
	[self convertIntensiveToExtensiveHessian:g dpMix:dgdr d2pMix:d2gdr2 mTotal:[self totalMolesFromMolesOfComponents:m] drdm:drdm d2rdm2:d2rdm2 d2p:d2gdm2];

	return d2gdm2Wrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleTensor *)getD3gDm3FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleTensor *d3gdm3Wrapper = [[DoubleTensor alloc] initWithFirstSize:NA andWithSecondSize:NA andWithThirdSize:NA];
	double ***d3gdm3 = [d3gdm3Wrapper pointerToPointerToPointerToDouble];
	double g, dgdr[NR], d2gdr2[NR][NR], d3gdr3[NR][NR][NR], r[NR], drdm[NR][NA], d2rdm2[NR][NA][NA], d3rdm3[NR][NA][NA][NA];
	[self convert:SECOND outMask:THIRD | FIFTH | SIXTH | EIGHTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:d2rdm2 dr:NULL d3m:d3rdm3];
	[self gmix:FIRST | SECOND | THIRD | FOURTH t:t p:p r:r gmix:&g dx:dgdr dx2:d2gdr2 dx3:d3gdr3];
	[self convertIntensiveToExtensiveTensor:g dpMix:dgdr d2pMix:d2gdr2 d3pMix:d3gdr3 mTotal:[self totalMolesFromMolesOfComponents:m] drdm:drdm d2rdm2:d2rdm2 d3rdm3:d3rdm3 d3p:d3gdm3];

	return d3gdm3Wrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getEnthalpyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double h, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self hmix:FIRST t:t p:p r:r hmix:&h];
	h *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return h;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		h += m[i]*[component getEnthalpyFromT:t andP:p];
	}
	return h;
}

// --> SolutionPhaseProtocol public function
-(double)getEntropyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double s, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self smix:FIRST t:t p:p r:r smix:&s dx:NULL dx2:NULL];
	s *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return s;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		s += m[i]*[component getEntropyFromT:t andP:p];
	}
	return s;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDsDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dsdmWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *dsdm = [dsdmWrapper pointerToDouble];
	double s, dsdr[NR], r[NR], drdm[NR][NA];
	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
	[self smix:FIRST | SECOND t:t p:p r:r smix:&s dx:dsdr dx2:NULL];
	[self convertIntensiveToExtensiveGradients:s dpMix:dsdr mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dsdm];

	if (computeMixingQuantities) return dsdmWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dsdm[i] += [component getEntropyFromT:t andP:p];
	}
	return dsdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getD2sDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *d2sdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:NA andWithColumnSize:NA];
	double **d2sdm2 = [d2sdm2Wrapper pointerToPointerToDouble];
	double s, dsdr[NR], d2sdr2[NR][NR], r[NR], drdm[NR][NA], d2rdm2[NR][NA][NA];
	[self convert:SECOND outMask:THIRD | FIFTH | SIXTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:d2rdm2 dr:NULL d3m:NULL];
	[self smix:FIRST | SECOND | THIRD t:t p:p r:r smix:&s dx:dsdr dx2:d2sdr2];
	[self convertIntensiveToExtensiveHessian:s dpMix:dsdr d2pMix:d2sdr2 mTotal:[self totalMolesFromMolesOfComponents:m] drdm:drdm d2rdm2:d2rdm2 d2p:d2sdm2];

	return d2sdm2Wrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getHeatCapacityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double cp, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self cpmix:FIRST t:t p:p r:r cpmix:&cp dt:NULL dx:NULL];
	cp *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return cp;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		cp += m[i]*[component getHeatCapacityFromT:t andP:p];
	}
	return cp;
}

// --> SolutionPhaseProtocol public function
-(double)getDcpDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double dcpdt, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self cpmix:SECOND t:t p:p r:r cpmix:NULL dt:&dcpdt dx:NULL];
	dcpdt *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return dcpdt;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dcpdt += m[i]*[component getDcpDtFromT:t andP:p];
	}
	return dcpdt;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDCpDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dcpdmWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *dcpdm = [dcpdmWrapper pointerToDouble];
	double cp, dcpdr[NR], r[NR], drdm[NR][NA];
	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
	[self cpmix:FIRST | THIRD t:t p:p r:r cpmix:&cp dt:NULL dx:dcpdr];
	[self convertIntensiveToExtensiveGradients:cp dpMix:dcpdr mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dcpdm];

	if (computeMixingQuantities) return dcpdmWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dcpdm[i] += [component getHeatCapacityFromT:t andP:p];
	}
	return dcpdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getVolumeFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double v, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self vmix:FIRST t:t p:p r:r vmix:&v dx:NULL dx2:NULL dt:NULL dp:NULL dt2:NULL dtdp:NULL dp2:NULL dxdt:NULL dxdp:NULL];
	v *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return v;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		v += m[i]*[component getVolumeFromT:t andP:p];
	}
	return v;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getDvDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dvdmWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *dvdm = [dvdmWrapper pointerToDouble];
	double v, dvdr[NR], r[NR], drdm[NR][NA];
	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
	[self vmix:FIRST | SECOND t:t p:p r:r vmix:&v dx:dvdr dx2:NULL dt:NULL dp:NULL dt2:NULL dtdp:NULL dp2:NULL dxdt:NULL dxdp:NULL];
	[self convertIntensiveToExtensiveGradients:v dpMix:dvdr mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dvdm];

	if (computeMixingQuantities) return dvdmWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdm[i] += [component getVolumeFromT:t andP:p];
	}
	return dvdmWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getD2vDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *d2vdm2Wrapper = [[DoubleMatrix alloc] initWithRowSize:NA andWithColumnSize:NA];
	double **d2vdm2 = [d2vdm2Wrapper pointerToPointerToDouble];
	double v, dvdr[NR], d2vdr2[NR][NR], r[NR], drdm[NR][NA], d2rdm2[NR][NA][NA];
	[self convert:SECOND outMask:THIRD | FIFTH | SIXTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:d2rdm2 dr:NULL d3m:NULL];
	[self vmix:FIRST | SECOND | THIRD t:t p:p r:r vmix:&v dx:dvdr dx2:d2vdr2 dt:NULL dp:NULL dt2:NULL dtdp:NULL dp2:NULL dxdt:NULL dxdp:NULL];
	[self convertIntensiveToExtensiveHessian:v dpMix:dvdr d2pMix:d2vdr2 mTotal:[self totalMolesFromMolesOfComponents:m] drdm:drdm d2rdm2:d2rdm2 d2p:d2vdm2];

	return d2vdm2Wrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getDvDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double dvdt, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self vmix:FOURTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:&dvdt dp:NULL dt2:NULL dtdp:NULL dp2:NULL dxdt:NULL dxdp:NULL];
	dvdt *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return dvdt;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdt += m[i]*[component getDvDtFromT:t andP:p];
	}
	return dvdt;
}

// --> SolutionPhaseProtocol public function
-(double)getDvDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double dvdp, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self vmix:FIFTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:NULL dp:&dvdp dt2:NULL dtdp:NULL dp2:NULL dxdt:NULL dxdp:NULL];
	dvdp *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return dvdp;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdp += m[i]*[component getDvDpFromT:t andP:p];
	}
	return dvdp;
}

// --> SolutionPhaseProtocol public function
-(double)getD2vDt2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double d2vdt2, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self vmix:SIXTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:NULL dp:NULL dt2:&d2vdt2 dtdp:NULL dp2:NULL dxdt:NULL dxdp:NULL];
	d2vdt2 *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return d2vdt2;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		d2vdt2 += m[i]*[component getD2vDt2FromT:t andP:p];
	}
	return d2vdt2;
}

// --> SolutionPhaseProtocol public function
-(double)getD2vDtDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double d2vdtdp, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self vmix:SEVENTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:NULL dp:NULL dt2:NULL dtdp:&d2vdtdp dp2:NULL dxdt:NULL dxdp:NULL];
	d2vdtdp *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return d2vdtdp;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		d2vdtdp += m[i]*[component getD2vDtDpFromT:t andP:p];
	}
	return d2vdtdp;
}

// --> SolutionPhaseProtocol public function
-(double)getD2vDp2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double d2vdp2, r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	[self vmix:EIGHTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:NULL dp:NULL dt2:NULL dtdp:NULL dp2:&d2vdp2 dxdt:NULL dxdp:NULL];
	d2vdp2 *= [self totalMolesFromMolesOfComponents:m];

	if (computeMixingQuantities) return d2vdp2;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		d2vdp2 += m[i]*[component getD2vDp2FromT:t andP:p];
	}
	return d2vdp2;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getD2vDmDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dvdmdtWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *dvdmdt = [dvdmdtWrapper pointerToDouble];
	double dvdt, dvdrdt[NR], r[NR], drdm[NR][NA];
	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
	[self vmix:FOURTH | NINTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:&dvdt dp:NULL dt2:NULL dtdp:NULL dp2:NULL dxdt:dvdrdt dxdp:NULL];
	[self convertIntensiveToExtensiveGradients:dvdt dpMix:dvdrdt mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dvdmdt];

	if (computeMixingQuantities) return dvdmdtWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdmdt[i] += [component getDvDtFromT:t andP:p];
	}
	return dvdmdtWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleVector *)getD2vDmDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *dvdmdpWrapper = [[DoubleVector alloc] initWithSize:NA];
	double *dvdmdp = [dvdmdpWrapper pointerToDouble];
	double dvdp, dvdrdp[NR], r[NR], drdm[NR][NA];
	[self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
	[self vmix:FIFTH | TENTH t:t p:p r:r vmix:NULL dx:NULL dx2:NULL dt:NULL dp:&dvdp dt2:NULL dtdp:NULL dp2:NULL dxdt:NULL dxdp:dvdrdp];
	[self convertIntensiveToExtensiveGradients:dvdp dpMix:dvdrdp mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dvdmdp];

	if (computeMixingQuantities) return dvdmdpWrapper;

	for (NSUInteger i=0; i<NA; i++) {
		id component = [endmembers objectAtIndex:i];
		dvdmdp[i] += [component getDvDpFromT:t andP:p];
	}
	return dvdmdpWrapper;
}

// --> SolutionPhaseProtocol public function
-(NSString *)getFormulaFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	double r[NR];
	[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
	return [self displayFormula:t p:p r:r];
}
