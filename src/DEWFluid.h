//
//  DEWFluid.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 12/5/15.
//  Copyright © 2015 Mark Ghiorso. All rights reserved.
//

#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

@interface DEWFluid : PhaseBase <SolutionPhaseProtocol>

@property (readonly) NSMutableArray *endmembers;

@property (assign, readwrite) BOOL debugV;
@property (assign, readwrite) BOOL debugS;

- (instancetype)initWithDuanCO2:(BOOL)useDuanCO2;
- (NSDictionary *)getSpeciesMoleFractionsForBulkComposition:(double *)m aT:(double)t andP:(double)p;
- (DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;
- (DoubleVector *)activitiesOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p;

-(void)testD2GDR2withT:(double)t p:(double)p r:(double [16])r;
-(void)testPrimD2GDR2withT:(double)t p:(double)p r:(double [16])r;
-(void)testPrimD2GDS2withT:(double)t p:(double)p r:(double [16])r;
-(void)testPrimD2GDRDSwithT:(double)t p:(double)p r:(double [16])r;
-(void)testAzerowithT:(double)t p:(double)p r:(double [16])r;
-(void)testD2DHDR2withT:(double)t p:(double)p r:(double [16])r;
-(void)testD2DHDRDSwithT:(double)t p:(double)p r:(double [16])r;
-(void)testD2DHDS2withT:(double)t p:(double)p r:(double [16])r;
-(void)testDRDSwithT:(double)t p:(double)p r:(double [16])r;

-(void)testFillnSpecieswithT:(double)t p:(double)p r:(double [16])r;
-(void)testFillxSpecieswithT:(double)t p:(double)p r:(double [16])r;

@end
