//
//  DEWFluid.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 12/5/15.
//  Copyright © 2015 Mark Ghiorso. All rights reserved.
//

#import "DEWFluid.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

#import "FluidDuan.h"
#import "GenericH2O.h"
#import "DEWspecies.h"
#import "HKFspeciesProperties.h"
#import "HKFspeciesComposite.h"
#import "DEWDielectricConstant.h"
#import "StoichiometricPhaseProtocol.h"

#ifdef __APPLE__
#import <Accelerate/Accelerate.h>
#else
#include "clapack.h"
#endif

#define NR      16      // Number of independent mole fraction variables
#define NS      75      // Three ordering parameters
#define NA      17      // Number of components
#define NE      (NA+NS) // Total number of species

#define mIndH2O      0
#define mIndHplus   17
#define mIndOHminus 18

#pragma mark -
#pragma mark Private instance variables

@interface DEWFluid () {
    BOOL computeMixingQuantities;
    // persistent variables
    double tOld, pOld, rOld[NR], sOld[NS];
    // global variables = f(r, s, t, p)
    double DcoeffDs[NS];
    double xSpecies[NE], dxSpecies[NE][NS], dxSpeciesDr[NE][NR];
    double nSpecies, dnSpeciesds[NS], d2nSpeciesds2[NS][NS];
    double DcbDs[NS];
    // local dielectric function for Debye-Hückel terms
    DEWDielectricConstant *dielectricConstant;
    // Instance variable to hold DewSpecies
    DEWspecies *dewSpecies;
    BOOL lowConcFlagForS[NS];

    // matrix utilized by the ordering routines for higher order derivatives
    double invd2gds2[NS][NS];

    // These variables are needed for computation of Debye-Hückel Excess terms
    BOOL initializeDHexcessTermsZerothOrder;
    double aZero, aZeroSum, mStar, AsubG, Bgamma, CapGamSubG, sqrtOfI, CapLambda;
    BOOL initializeDHexcessTermsFirstOrderDs;
    double dCapGamSubGds[NS], DsqrtOfIds[NS], dCapLambdaDs[NS], DmStarDs[NS];
    BOOL initializeDHexcessTermsFirstOrderDr;
    double dCapGamSubGdr[NR], DsqrtOfIdr[NR], dCapLambdaDr[NR], DmStarDr[NR];
    BOOL initializeDHexcessTermsSecondOrderDs2;
    double d2CapGamSubGds2[NS][NS], D2sqrtOfIds2[NS][NS], d2CapLambdaDs2[NS][NS], D2mStarDs2[NS][NS];
    BOOL initializeDHexcessTermsFirstOrderDt;
    double dAsubGDt, dBgammaDt, dCapGamSubGDt, dCapLambdaDt;
    BOOL initializeDHexcessTermsZerothOrderDp;
    double dAsubGDp, dBgammaDp, dCapLambdaDp;
    BOOL initializeDHexcessTermsSecondOrderDr2;
    double d2CapGamSubGdr2[NR][NR], D2sqrtOfIdr2[NR][NR], d2CapLambdaDr2[NR][NR], D2mStarDr2[NR][NR];
    BOOL initializeDHexcessTermsSecondOrderDrDs;
    double d2CapGamSubGdrds[NR][NS], D2sqrtOfIdrds[NR][NS], d2CapLambdaDrDs[NR][NS], D2mStarDrDs[NR][NS];

    //
    // These are instance and working storage for the LAPACK routine dspsvx (Diagonal pivoting factorization and solution of linear
    //   systems for a packed symmetric positive definite matrix
    //
    char dspsvx_fact;                          // equilibrate the system; 'F' denotes matrix is factored already and stored in AF
    char dspsvx_uplo;                          // 'U' stored in upper triangle of A; 'L' would be lower triangle
    __CLPK_integer dspsvx_n;                   // number of linear equations
    __CLPK_integer dspsvx_nrhs;                // number of right-hand sides
    __CLPK_doublereal dspsvx_ap[NS*(NS-1)/2];  // content depends on fact, equed, and uplo
    __CLPK_doublereal dspsvx_afp[NS*(NS-1)/2]; // input/output Factorization of A
    __CLPK_integer dspsvx_ipiv[NS];            // input/output Pivot elements for diagonal interchanges
    __CLPK_doublereal dspsvx_b[NS*NR];         // input/output right-hand side vector
    __CLPK_integer dspsvx_ldb;                 // leading dimension of B
    __CLPK_doublereal dspsvx_x[NS*NR];         // output solution vector; if equed='y', then solution is x_i/s_i
    __CLPK_integer dspsvx_ldx;                 // leading dimension of X
    __CLPK_doublereal dspsvx_rcond;            // output condition number; should be > machine precision
    __CLPK_doublereal dspsvx_ferr[NR];         // output relative forward error
    __CLPK_doublereal dspsvx_berr[NR];         // output relative backward error
    __CLPK_doublereal dspsvx_work[3*NS];       // working array
    __CLPK_integer dspsvx_iwork[NS];           // working array
    __CLPK_integer dspsvx_info;                // output;  0 == successful exit
                                               //         -i == ith argument is bad
                                               //          i == ith diagonal element is exactly zero
                                               //       NS+1 == rcond < machine precision, no solution or errors reported

    //
    // These are instance and working storage for the LAPACK routine dgglse (Linear Equality Constrained Least Sqaures Problem)
    //
    __CLPK_integer    dgglse_m;            // rows of matrix A
    __CLPK_integer    dgglse_n;            // columns of matrices A and B
    __CLPK_integer    dgglse_p;            // rows of matrix B
    __CLPK_doublereal dgglse_a[NS*NS];     // matrix A(lda, n), least squares solution matrix
    __CLPK_integer    dgglse_lda;          // row dimension of A, lda > m
    __CLPK_doublereal dgglse_b[NS];        // matrix B(lba, n), equality constraints
    __CLPK_integer    dgglse_ldb;          // row dimension of B, ldb > p
    __CLPK_doublereal dgglse_c[NS];        // vector C(m), rhs of least sqaures problem.  Residual sum of squares in n-p+1 to m
    __CLPK_doublereal dgglse_d[1];         // vector D(p), constraint constant
    __CLPK_doublereal dgglse_x[NS];        // solution vector (n)
    __CLPK_doublereal dgglse_work[2*NS+1]; // work array, minimal size: n+p+max(n,m,p)
    __CLPK_integer    dgglse_lwork;        // dimension of work
    __CLPK_integer    dgglse_info;         // = 0 on successful exit
                                           // < 0 if -i, the ith argument had an illegal value

    // Indexes for inversion routines to accomodate zero values of xSpecies
    BOOL incOfxSpecies[NS];
}

@end

@implementation DEWFluid

@synthesize endmembers;

#pragma mark -
#pragma mark instance methods

- (instancetype)initWithDuanCO2:(BOOL)useDuanCO2 {
    if ((self = [super init])) {
        NSString *orderedNamesOfSpecies[] = {
            @"H2O", @"CO2,aq",
            @"O2,aq", @"HF,aq", @"NaOH,aq", @"Mg(OH)2,aq", @"HAlO2,aq", @"SiO2,aq", @"H3PO4,aq", @"SO2,aq", @"HCl,aq",
            @"KOH,aq", @"Ca(OH)2,aq", @"H2CrO4,aq", @"Mn(OH)2,aq", @"Fe(OH)2,aq", @"Co(OH)2,aq",
            @"H+", @"OH-", @"H2,aq", @"CO3-2", @"HCO3-", @"CO,aq", @"F-",
            @"NaCl,aq", @"Na+", @"NaCO3-", @"NaHCO3,aq", @"NaHSiO3,aq",
            @"MgCO3,aq", @"Mg(HSiO3)+", @"Mg(HCO3)+", @"Mg+2", @"MgCl+", @"MgOH+", @"MgSO4,aq",
            @"Al+3", @"AlO2-",
            @"HSiO3-", @"Si2O4,aq",
            @"H2PO4-", @"HPO4-2", @"PO4-3", @"H3P2O7-", @"H2P2O7-2",
            @"H2S,aq", @"HS-", @"S2-2", @"S2O3-2", @"S2O4-2", @"S2O5-2", @"S2O6-2", @"S2O8-2",
            @"S3-2", @"S3O6-2", @"S4-2", @"S4O6-2", @"S5-2", @"S5O6-2", @"SO3-2", @"HSO3-",
            @"SO4-2", @"HSO4-", @"HSO5-",
            @"Cl-",
            @"K+", @"KCl,aq", @"KSO4-",
            @"CaCO3,aq", @"Ca(HCO3)+", @"Ca(OH)+", @"Ca+2", @"CaCl+", @"CaCl2,aq", @"CaSO4,aq",
            @"Cr+2", @"Cr+3", @"Cr2O7-2", @"CrO4-2", @"HCrO4-",
            @"Mn+2", @"MnCl+", @"MnO4-", @"MnO4-2", @"MnSO4,aq",
            @"Fe+2", @"Fe+3", @"FeCl+", @"FeCl+2", @"FeCl2,aq",
            @"Co+2", @"Co+3"
        };

        NSUInteger check = sizeof(orderedNamesOfSpecies)/sizeof(NSString *);
        NSAssert(check == NE, @"Inconsistency in ordered species namelist and species total.");

        dewSpecies = [[DEWspecies alloc] init];
        endmembers = [NSMutableArray arrayWithCapacity:[dewSpecies numberOfSpecies]];

        for (NSUInteger i=0; i<NE; i++) {
            NSString *speciesName = orderedNamesOfSpecies[i];
            if ([speciesName isEqualToString:@"H2O"]) {
                [endmembers addObject:[[GenericH2O alloc] init]];  // X0 - H2O is a special case taken from SWIM
            } else if ([speciesName isEqualToString:@"CO2,aq"] && useDuanCO2) {
                [endmembers addObject:[[DuanCO2 alloc] init]];     // X1 - CO2 is a special case taken from Duan and Zhang, 2006
            } else if ([speciesName isEqualToString:@"HF,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"H+":@1.0, @"F-":@1.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"HF,aq";
                hkfComposite.phaseFormula = @"HF";
                [endmembers addObject:hkfComposite];
            } else if ([speciesName isEqualToString:@"Mg(OH)2,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"Mg+2":@1.0, @"OH-":@2.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"Mg(OH)2,aq";
                hkfComposite.phaseFormula = @"Mg(OH)2";
                [endmembers addObject:hkfComposite];
            } else if ([speciesName isEqualToString:@"Ca(OH)2,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"Ca+2":@1.0, @"OH-":@2.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"Ca(OH)2,aq";
                hkfComposite.phaseFormula = @"Ca(OH)2";
                [endmembers addObject:hkfComposite];
            } else if ([speciesName isEqualToString:@"H2CrO4,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"H+":@2.0, @"CrO4-2":@1.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"H2CrO4,aq";
                hkfComposite.phaseFormula = @"H2CrO4";
                [endmembers addObject:hkfComposite];
            } else if ([speciesName isEqualToString:@"Mn(OH)2,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"Mn+2":@1.0, @"OH-":@2.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"Mn(OH)2,aq";
                hkfComposite.phaseFormula = @"Mn(OH)2";
                [endmembers addObject:hkfComposite];
            } else if ([speciesName isEqualToString:@"Fe(OH)2,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"Fe+2":@1.0, @"OH-":@2.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"Fe(OH)2,aq";
                hkfComposite.phaseFormula = @"Fe(OH)2";
                [endmembers addObject:hkfComposite];
            } else if ([speciesName isEqualToString:@"Co(OH)2,aq"]) {
                HKFspeciesComposite *hkfComposite =
                [[HKFspeciesComposite alloc] initWithDEWspeciesInstance:dewSpecies
                                              andWithReactionDictionary:@{@"Co+2":@1.0, @"OH-":@2.0}];
                [hkfComposite setCorrectionToGibbsEnergy:20000.0]; // arbitrary, set to destablize the species
                hkfComposite.phaseName = @"Co(OH)2,aq";
                hkfComposite.phaseFormula = @"Co(OH)2";
                [endmembers addObject:hkfComposite];
            } else {
#ifdef __APPLE__
                NSAssert([dewSpecies namedSpeciesInCollection:speciesName], @"Species %@ not found in DEW collection", speciesName);
#else
                NSAssert([dewSpecies namedSpeciesInCollection:speciesName], @"Species not found in DEW collection");
#endif
                [endmembers addObject:[dewSpecies HKFclassForSpeciesNamed:speciesName]];
            }
        }

        dielectricConstant = [[DEWDielectricConstant alloc] init];
        initializeDHexcessTermsZerothOrder     = YES;
        initializeDHexcessTermsFirstOrderDs    = YES;
        initializeDHexcessTermsFirstOrderDr    = YES;
        initializeDHexcessTermsSecondOrderDs2  = YES;
        initializeDHexcessTermsFirstOrderDt    = YES;
        initializeDHexcessTermsZerothOrderDp   = YES;
        initializeDHexcessTermsSecondOrderDr2  = YES;
        initializeDHexcessTermsSecondOrderDrDs = YES;

        computeMixingQuantities = NO;
        tOld = -9999.0;
        pOld = -9999.0;
        for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
        for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
        [self setPhaseName:@"DEWFluid"];
        _debugS = YES;
        _debugV = YES;
    }
    return self;
}

- (instancetype)init {
    return [self initWithDuanCO2:NO];
}

#pragma mark -
#pragma mark NSSecureCoding protocol methods

static NSString *kendmembers = @"endmembers";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        computeMixingQuantities = NO;
        tOld = -9999.0;
        pOld = -9999.0;
        for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
        for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
#ifdef __APPLE__
        endmembers = (NSMutableArray *) [aDecoder decodeObjectOfClass:[NSArray class] forKey:kendmembers];
#else
        endmembers = (NSMutableArray *) [aDecoder decodeObjectForKey:kendmembers];
#endif
        dielectricConstant = [[DEWDielectricConstant alloc] init];
        initializeDHexcessTermsZerothOrder     = YES;
        initializeDHexcessTermsFirstOrderDs    = YES;
        initializeDHexcessTermsFirstOrderDr    = YES;
        initializeDHexcessTermsSecondOrderDs2  = YES;
        initializeDHexcessTermsFirstOrderDt    = YES;
        initializeDHexcessTermsZerothOrderDp   = YES;
        initializeDHexcessTermsSecondOrderDr2  = YES;
        initializeDHexcessTermsSecondOrderDrDs = YES;
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeObject:endmembers forKey:kendmembers];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", [self className]];
}

#pragma mark -
#pragma mark Solution properties

#define R      8.3143
#define CapGam (1000.0/18.01528)

#define FIRST       00000001 /* octal for binary 00000000000000000001 */
#define SECOND      00000002 /* octal for binary 00000000000000000010 */
#define THIRD       00000004 /* octal for binary 00000000000000000100 */
#define FOURTH      00000010 /* octal for binary 00000000000000001000 */
#define FIFTH       00000020 /* octal for binary 00000000000000010000 */
#define SIXTH       00000040 /* octal for binary 00000000000000100000 */
#define SEVENTH     00000100 /* octal for binary 00000000000001000000 */
#define EIGHTH      00000200 /* octal for binary 00000000000010000000 */
#define NINTH       00000400 /* octal for binary 00000000000100000000 */
#define TENTH       00001000 /* octal for binary 00000000001000000000 */
#define ELEVENTH    00002000 /* octal for binary 00000000010000000000 */
#define TWELFTH     00004000 /* octal for binary 00000000100000000000 */
#define THIRTEENTH  00010000 /* octal for binary 00000001000000000000 */
#define FOURTEENTH  00020000 /* octal for binary 00000010000000000000 */
#define FIFTEENTH   00040000 /* octal for binary 00000100000000000000 */
#define SIXTEENTH   00100000 /* octal for binary 00001000000000000000 */
#define SEVENTEENTH 00200000 /* octal for binary 00010000000000000000 */
#define EIGHTEENTH  00400000 /* octal for binary 00100000000000000000 */
#define NINETEENTH  01000000 /* octal for binary 01000000000000000000 */
#define TWENTIETH   02000000 /* octal for binary 10000000000000000000 */

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))


/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives
 */

/*
 X0   H2O
 X1   CO2     r0
 X2   O2      r1
 X3   HF      r2
 X4   NaOH    r3
 X5   Mg(OH)2 r4
 X6   HAlO2   r5
 X7   SiO2    r6
 X8   H3PO4   r7
 X9   SO2     r8
 X10  HCl     r9
 X11  KOH     r10
 X12  Ca(OH)2 r11
 X13  H2CrO4  r12
 X14  Mn(OH)2 r13
 X15  Fe(OH)2 r14
 X16  Co(OH)2 r15

 X17  H+       s0
 X18  OH-      s1
 X19  H2       s2
 X20  CO3=     s3
 X21  HCO3-    s4
 X22  CO       s5
 X23  F-       s6
 X24  NaCl     s7
 X25  Na+      s8
 X26  NaCO3-   s9
 X27  NaHCO3   s10
 X28  NaHSiO3  s11
 X29  MgCO3    s12
 X30  MgHSiO3+ s13
 X31  MgHCO3+  s14
 X32  Mg+2     s15
 X33  MgCl+    s16
 X34  MgOH+    s17
 X35  MgSO4    s18
 X36  Al+3     s19
 X37  AlO2-    s20
 X38  HSiO3-   s21
 X39  Si2O4    s22
 X40  H2PO4-   s23
 X41  HPO4-2   s24
 X42  PO4-3    s25
 X43  H3P2O7-  s26
 X44  H2P2O7-2 s27
 X45  H2S      s28
 X46  HS-      s29
 X47  S2-2     s30
 X48  S2O3-2   s31
 X49  S2O4-2   s32
 X50  S2O5-2   s33
 X51  S2O6-2   s34
 X52  S2O8-2   s35
 X53  S3-2     s36
 X54  S3O6-2   s37
 X55  S4-2     s38
 X56  S4O6-2   s39
 X57  S5-2     s40
 X58  S5O6-2   s41
 X59  SO3-2    s42
 X60  HSO3-    s43
 X61  SO4-2    s44
 X62  HSO4-    s45
 X63  HSO5-    s46
 X64  Cl-      s47
 X65  K+       s48
 X66  KCl      s49
 X67  KSO4-    s50
 X68  CaCO3    s51
 X69  CaHCO3+  s52
 X70  CaOH+    s53
 X71  Ca+2     s54
 X72  CaCl+    s55
 X73  CaCl2    s56
 X74  CaSO4    s57
 X75  Cr+2     s58
 X76  Cr+3     s59
 X77  Cr2O7-2  s60
 X78  CrO4-2   s61
 X79  HCrO4-   s62
 X80  Mn+2     s63
 X81  MnCl+    s64
 X82  MnO4-    s65
 X83  MnO4-2   s66
 X84  MnSO4    s67
 X85  Fe2+     s68
 X86  Fe3+     s69
 X87  FeCl+    s70
 X88  FeCl+2   s71
 X89  FeCl2    s72
 X90  Co+2     s73
 X91  Co+3     s74
 */

static const NSUInteger indexOfHydrogenIon = 0;

// only called by - (BOOL)fillxSpeciesWith:(double [NR])r andWith:(double [NS])s
- (double)speciationCoefficientWith:(double [NS])s {
    double coeff = 1.0 - 3.0*s[0]/4.0 - s[1]/4.0 - s[2]/2.0 + 3.0*s[4]/4.0  + s[10] + s[11] + s[13]/4.0 + s[14]/4.0
    - 3.0*s[17]/4.0 + 3.0*s[21]/4.0 - s[23]/4.0 - s[24]/2.0 - s[26]/4.0 - s[27]/2.0 - s[28]/2.0 - 3.0*s[29]/4.0
    + 3.0*s[43]/4.0 + 5.0*s[45]/4.0 + 7.0*s[46]/4.0 + s[52]/4.0 - 3.0*s[53]/4.0 - s[62]/4.0 + s[3]/2.0 - s[5]/2.0
    + 3.0*s[9]/4.0 + s[18]/2.0 - s[20]/4.0 + s[22] - 3.0*s[25]/4.0 + s[31]/2.0 + s[32] + 3.0*s[33]/2.0 + 2.0*s[34]
    + 3.0*s[35] + 2.0*s[37] + 2.0*s[39] + 2.0*s[41] + s[42]/2.0 + s[44] + 5.0*s[50]/4.0 + s[57]/2.0 - s[60]/2.0
    - s[61]/2.0 + s[65]/2.0 + s[66]/2.0 + s[67]/2.0 - s[6]/4.0 - 3.0*s[8]/4.0 - 3.0*s[15]/2.0 - 3.0*s[16]/4.0
    - 5.0*s[19]/4.0 - s[30] - s[36] - s[38] - s[40] - s[47]/4.0 - 3.0*s[55]/4.0 - 3.0*s[64]/4.0 - 3.0*s[70]/4.0
    - 3.0*s[71]/4.0 - 3.0*s[48]/4.0 - 3.0*s[54]/2.0 - 5.0*s[58]/2.0 - 5.0*s[59]/2.0 - 3.0*s[63]/2.0
    - 3.0*s[68]/2.0 - 3.0*s[69]/2.0 - 3.0*s[73]/2.0 - 3.0*s[74]/2.0;
    return coeff;
}

// only called by - (void)DxSpeciesDsWith:(double [NR])r andWith:(double [NS])s
- (void)DspeciationCoefficientDsWith:(double [NS])s {
    DcoeffDs[ 0] = - 3.0/4.0;
    DcoeffDs[ 1] = - 1.0/4.0;
    DcoeffDs[ 2] = - 1.0/2.0;
    DcoeffDs[ 3] =   1.0/2.0;
    DcoeffDs[ 4] =   3.0/4.0;
    DcoeffDs[ 5] = - 1.0/2.0;
    DcoeffDs[ 6] = - 1.0/4.0;
    DcoeffDs[ 7] =   0.0;
    DcoeffDs[ 8] = - 3.0/4.0;
    DcoeffDs[ 9] =   3.0/4.0;
    DcoeffDs[10] =   1.0;
    DcoeffDs[11] =   1.0;
    DcoeffDs[12] =   0.0;
    DcoeffDs[13] =   1.0/4.0;
    DcoeffDs[14] =   1.0/4.0;
    DcoeffDs[15] = - 3.0/2.0;
    DcoeffDs[16] = - 3.0/4.0;
    DcoeffDs[17] = - 3.0/4.0;
    DcoeffDs[18] =   1.0/2.0;
    DcoeffDs[19] = - 5.0/4.0;
    DcoeffDs[20] = - 1.0/4.0;
    DcoeffDs[21] =   3.0/4.0;
    DcoeffDs[22] =   1.0;
    DcoeffDs[23] = - 1.0/4.0;
    DcoeffDs[24] = - 1.0/2.0;
    DcoeffDs[25] = - 3.0/4.0;
    DcoeffDs[26] = - 1.0/4.0;
    DcoeffDs[27] = - 1.0/2.0;
    DcoeffDs[28] = - 1.0/2.0;
    DcoeffDs[29] = - 3.0/4.0;
    DcoeffDs[30] = - 1.0;
    DcoeffDs[31] =   1.0/2.0;
    DcoeffDs[32] =   1.0;
    DcoeffDs[33] =   3.0/2.0;
    DcoeffDs[34] =   2.0;
    DcoeffDs[35] =   3.0;
    DcoeffDs[36] = - 1.0;
    DcoeffDs[37] =   2.0;
    DcoeffDs[38] = - 1.0;
    DcoeffDs[39] =   2.0;
    DcoeffDs[40] = - 1.0;
    DcoeffDs[41] =   2.0;
    DcoeffDs[42] =   1.0/2.0;
    DcoeffDs[43] =   3.0/4.0;
    DcoeffDs[44] =   1.0;
    DcoeffDs[45] =   5.0/4.0;
    DcoeffDs[46] =   7.0/4.0;
    DcoeffDs[47] = - 1.0/4.0;
    DcoeffDs[48] = - 3.0/4.0;
    DcoeffDs[49] =   0.0;
    DcoeffDs[50] =   5.0/4.0;
    DcoeffDs[51] =   0.0;
    DcoeffDs[52] =   1.0/4.0;
    DcoeffDs[53] = - 3.0/4.0;
    DcoeffDs[54] = - 3.0/2.0;
    DcoeffDs[55] = - 3.0/4.0;
    DcoeffDs[56] =   0.0;
    DcoeffDs[57] =   1.0/2.0;
    DcoeffDs[58] = - 5.0/2.0;
    DcoeffDs[59] = - 5.0/2.0;
    DcoeffDs[60] = - 1.0/2.0;
    DcoeffDs[61] = - 1.0/2.0;
    DcoeffDs[62] = - 1.0/4.0;
    DcoeffDs[63] = - 3.0/2.0;
    DcoeffDs[64] = - 3.0/4.0;
    DcoeffDs[65] =   1.0/2.0;
    DcoeffDs[66] =   1.0/2.0;
    DcoeffDs[67] =   1.0/2.0;
    DcoeffDs[68] = - 3.0/2.0;
    DcoeffDs[69] = - 3.0/2.0;
    DcoeffDs[70] = - 3.0/4.0;
    DcoeffDs[71] = - 3.0/4.0;
    DcoeffDs[72] =   0.0;
    DcoeffDs[73] = - 3.0/2.0;
    DcoeffDs[74] = - 3.0/2.0;
}

- (BOOL)fillxSpeciesWith:(double [NR])r andWith:(double [NS])s {
    double rTotal = 0.0;
    for (NSUInteger i=0; i<NR; i++) rTotal += r[i];

    double coeff = [self speciationCoefficientWith:s];

    xSpecies[0] = (1.0-rTotal)*coeff + 2.0*(coeff-1.0)/3.0
                - s[1]/3.0 - 2.0*s[2]/3.0 - s[4] - 2.0*s[10]/3.0 - 2.0*s[11]/3.0 + s[13]/3.0 + s[14]/3.0 + s[17]
                - s[21] + 2.0*s[23]/3.0 + 4.0*s[24]/3.0 + 5.0*s[26]/3.0 + 7.0*s[27]/3.0 - 2.0*s[28]/3.0 - s[43]
                - 4.0*s[45]/3.0 - 5.0*s[46]/3.0 + s[52]/3.0 + s[53] + 2.0*s[62]/3.0- s[3]/3.0 + s[5]/3.0 + s[12]
                + s[51] + 2.0*s[18]/3.0 + 2.0*s[20]/3.0 - 2.0*s[22]/3.0 + 2.0*s[25] - s[31]/3.0 - 2.0*s[32]/3.0
                - s[33] - 4.0*s[34]/3.0 - 2.0*s[35] - 4.0*s[37]/3.0 - 4.0*s[39]/3.0 - 4.0*s[41]/3.0 - s[42]/3.0
                - 2.0*s[44]/3.0 - s[50]/3.0 + 2.0*s[57]/3.0 + 7.0*s[60]/3.0 + 4.0*s[61]/3.0 + 2.0*s[65]/3.0
                + 2.0*s[66]/3.0 + 2.0*s[67]/3.0 + 2.0*s[6]/3.0 + s[7] + s[8] + 2.0*s[15] + 2.0*s[16] + 4.0*s[19]/3.0
                + 2.0*s[30]/3.0 + 2.0*s[36]/3.0 + 2.0*s[38]/3.0 + 2.0*s[40]/3.0 + 2.0*s[47]/3.0 + s[49]
                + 2.0*s[55] + 2.0*s[56] + 2.0*s[64] + 2.0*s[70] + 2.0*s[71] + 2.0*s[72] + s[48] + 2.0*s[54]
                + 8.0*s[58]/3.0 + 8.0*s[59]/3.0 + 2.0*s[63] + 2.0*s[68] + 2.0*s[69] + 2.0*s[73] + 2.0*s[74];
    xSpecies[1] = (r[0] > 0.0) ? r[0]*coeff - s[4] - s[10] - s[14] - s[52] - s[3] - s[5] - s[9] - s[12] - s[51] : 0.0;
    xSpecies[2] = (r[1] > 0.0) ? r[1]*coeff + (1.0-coeff)/3.0 - s[1]/3.0 + s[2]/3.0 + s[10]/3.0 + s[11]/3.0 + s[13]/3.0 + s[14]/3.0
                - s[23]/3.0 - 2.0*s[24]/3.0 - s[26]/3.0 - 2.0*s[27]/3.0 + 4.0*s[28]/3.0 + s[29] - s[45]/3.0
                - 2.0*s[46]/3.0 + s[52]/3.0 - s[62]/3.0 - s[3]/3.0 + s[5]/3.0 - s[18]/3.0 - s[20]/3.0 + s[22]/3.0
                - s[25] + 2.0*s[31]/3.0 + s[32]/3.0 - s[34]/3.0 - s[35] + 2.0*s[37]/3.0 + 5.0*s[39]/3.0
                + 8.0*s[41]/3.0 - s[42]/3.0 - 2.0*s[44]/3.0 - s[50]/3.0 - s[57]/3.0 - 2.0*s[60]/3.0 - 2.0*s[61]/3.0
                - 4.0*s[65]/3.0 - 4.0*s[66]/3.0 - s[67]/3.0 - s[6]/3.0 + s[19]/3.0 + 5.0*s[30]/3.0 + 8.0*s[36]/3.0
                + 11.0*s[38]/3.0 + 14.0*s[40]/3.0 - s[47]/3.0 + 2.0*s[58]/3.0 + 2.0*s[59]/3.0 : 0.0;
    xSpecies[3] = (r[2] > 0.0) ? r[2]*coeff - s[6] : 0.0;
    xSpecies[4] = (r[3] > 0.0) ? r[3]*coeff - s[10] - s[11] - s[9] - s[7] - s[8] : 0.0;
    xSpecies[5] = (r[4] > 0.0) ? r[4]*coeff - s[13] - s[14] - s[17] - s[12] - s[18] - s[15] - s[16] : 0.0;
    xSpecies[6] = (r[5] > 0.0) ? r[5]*coeff - s[20] - s[19] : 0.0;
    xSpecies[7] = (r[6] > 0.0) ? r[6]*coeff - s[11] - s[13] - s[21] - 2.0*s[22] : 0.0;
    xSpecies[8] = (r[7] > 0.0) ? r[7]*coeff - s[23] - s[24] - 2.0*s[26] - 2.0*s[27] - s[25] : 0.0;
    xSpecies[9] = (r[8] > 0.0) ? r[8]*coeff - s[28] - s[29] - s[43] - s[45] - s[46] - s[18] - 2.0*s[31] - 2.0*s[32] - 2.0*s[33]
                - 2.0*s[34] - 2.0*s[35] - 3.0*s[37] - 4.0*s[39] - 5.0*s[41] - s[42] - s[44] - s[50] - s[57]
                - s[67] - 2.0*s[30] - 3.0*s[36] - 4.0*s[38] - 5.0*s[40] : 0.0;
    xSpecies[10] = (r[9] > 0.0) ? r[9]*coeff - s[7] - s[16] - s[47] - s[49] - s[55] - 2.0*s[56] - s[64] - s[70] - s[71] - 2.0*s[72] : 0.0;
    xSpecies[11] = (r[10] > 0.0) ? r[10]*coeff - s[50] - s[49] - s[48] : 0.0;
    xSpecies[12] = (r[11] > 0.0) ? r[11]*coeff - s[52] - s[53] - s[51] - s[57] - s[55] - s[54] - s[56] : 0.0;
    xSpecies[13] = (r[12] > 0.0) ? r[12]*coeff - s[62] - 2.0*s[60] - s[61] - s[58] - s[59] : 0.0;
    xSpecies[14] = (r[13] > 0.0) ? r[13]*coeff - s[65] - s[66] - s[67] - s[64] - s[63] : 0.0;
    xSpecies[15] = (r[14] > 0.0) ? r[14]*coeff - s[70] - s[71] - s[72] - s[68] - s[69] : 0.0;
    xSpecies[16] = (r[15] > 0.0) ? r[15]*coeff - s[73] - s[74] : 0.0;

    for (NSUInteger i=0; i<NS; i++) xSpecies[NA+i] = s[i];

    BOOL result = YES;
    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] < 0.0) result = NO;
    //if (xSpecies[0] < 0.9) result = NO;
    return result;
}

- (void)DxSpeciesDrWith:(double [NR])r andWith:(double [NS])s {
    memset(dxSpeciesDr, 0, sizeof(dxSpeciesDr));
    double coeff = [self speciationCoefficientWith:s];
    for (NSUInteger i=0; i<NR; i++) {
        dxSpeciesDr[0][i]   = -coeff;
        dxSpeciesDr[i+1][i] = coeff;
    }
}

- (void)DxSpeciesDsWith:(double [NR])r andWith:(double [NS])s {
    double rTotal = 0.0;
    for (NSUInteger i=0; i<NR; i++) rTotal += r[i];
    memset(dxSpecies, 0, sizeof(dxSpecies));

    [self DspeciationCoefficientDsWith:s];

    for (NSUInteger i=0; i<NS; i++) dxSpecies[0][i] = (1.0-rTotal)*DcoeffDs[i] + 2.0*DcoeffDs[i]/3.0;
    dxSpecies[0][ 1] += - 1.0/3.0;
    dxSpecies[0][ 2] += - 2.0/3.0;
    dxSpecies[0][ 4] += - 1.0;
    dxSpecies[0][10] += - 2.0/3.0;
    dxSpecies[0][11] += - 2.0/3.0;
    dxSpecies[0][13] +=   1.0/3.0;
    dxSpecies[0][14] +=   1.0/3.0;
    dxSpecies[0][17] +=   1.0;
    dxSpecies[0][21] += - 1.0;
    dxSpecies[0][23] +=   2.0/3.0;
    dxSpecies[0][24] +=   4.0/3.0;
    dxSpecies[0][26] +=   5.0/3.0;
    dxSpecies[0][27] +=   7.0/3.0;
    dxSpecies[0][28] += - 2.0/3.0;
    dxSpecies[0][43] += - 1.0;
    dxSpecies[0][45] += - 4.0/3.0;
    dxSpecies[0][46] += - 5.0/3.0;
    dxSpecies[0][52] +=   1.0/3.0;
    dxSpecies[0][53] +=   1.0;
    dxSpecies[0][62] +=   2.0/3.0;
    dxSpecies[0][ 3] += - 1.0/3.0;
    dxSpecies[0][ 5] +=   1.0/3.0;
    dxSpecies[0][12] +=   1.0;
    dxSpecies[0][51] +=   1.0;
    dxSpecies[0][18] +=   2.0/3.0;
    dxSpecies[0][20] +=   2.0/3.0;
    dxSpecies[0][22] += - 2.0/3.0;
    dxSpecies[0][25] +=   2.0;
    dxSpecies[0][31] += - 1.0/3.0;
    dxSpecies[0][32] += - 2.0/3.0;
    dxSpecies[0][33] += - 1.0;
    dxSpecies[0][34] += - 4.0/3.0;
    dxSpecies[0][35] += - 2.0;
    dxSpecies[0][37] += - 4.0/3.0;
    dxSpecies[0][39] += - 4.0/3.0;
    dxSpecies[0][41] += - 4.0/3.0;
    dxSpecies[0][42] += - 1.0/3.0;
    dxSpecies[0][44] += - 2.0/3.0;
    dxSpecies[0][50] += - 1.0/3.0;
    dxSpecies[0][57] +=   2.0/3.0;
    dxSpecies[0][60] +=   7.0/3.0;
    dxSpecies[0][61] +=   4.0/3.0;
    dxSpecies[0][65] +=   2.0/3.0;
    dxSpecies[0][66] +=   2.0/3.0;
    dxSpecies[0][67] +=   2.0/3.0;
    dxSpecies[0][ 6] +=   2.0/3.0;
    dxSpecies[0][ 7] +=   1.0;
    dxSpecies[0][ 8] +=   1.0;
    dxSpecies[0][15] +=   2.0;
    dxSpecies[0][16] +=   2.0;
    dxSpecies[0][19] +=   4.0/3.0;
    dxSpecies[0][30] +=   2.0/3.0;
    dxSpecies[0][36] +=   2.0/3.0;
    dxSpecies[0][38] +=   2.0/3.0;
    dxSpecies[0][40] +=   2.0/3.0;
    dxSpecies[0][47] +=   2.0/3.0;
    dxSpecies[0][49] +=   1.0;
    dxSpecies[0][55] +=   2.0;
    dxSpecies[0][56] +=   2.0;
    dxSpecies[0][64] +=   2.0;
    dxSpecies[0][70] +=   2.0;
    dxSpecies[0][71] +=   2.0;
    dxSpecies[0][72] +=   2.0;
    dxSpecies[0][48] +=   1.0;
    dxSpecies[0][54] +=   2.0;
    dxSpecies[0][58] +=   8.0/3.0;
    dxSpecies[0][59] +=   8.0/3.0;
    dxSpecies[0][63] +=   2.0;
    dxSpecies[0][68] +=   2.0;
    dxSpecies[0][69] +=   2.0;
    dxSpecies[0][73] +=   2.0;
    dxSpecies[0][74] +=   2.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[1][i] = r[0]*DcoeffDs[i];
    dxSpecies[1][ 4] += - 1.0;
    dxSpecies[1][10] += - 1.0;
    dxSpecies[1][14] += - 1.0;
    dxSpecies[1][52] += - 1.0;
    dxSpecies[1][ 3] += - 1.0;
    dxSpecies[1][ 5] += - 1.0;
    dxSpecies[1][ 9] += - 1.0;
    dxSpecies[1][12] += - 1.0;
    dxSpecies[1][51] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[2][i] = r[1]*DcoeffDs[i] - DcoeffDs[i]/3.0;
    dxSpecies[2][ 1] += - 1.0/3.0;
    dxSpecies[2][ 2] +=   1.0/3.0;
    dxSpecies[2][10] +=   1.0/3.0;
    dxSpecies[2][11] +=   1.0/3.0;
    dxSpecies[2][13] +=   1.0/3.0;
    dxSpecies[2][14] +=   1.0/3.0;
    dxSpecies[2][23] += - 1.0/3.0;
    dxSpecies[2][24] += - 2.0/3.0;
    dxSpecies[2][26] += - 1.0/3.0;
    dxSpecies[2][27] += - 2.0/3.0;
    dxSpecies[2][28] +=   4.0/3.0;
    dxSpecies[2][29] +=   1.0;
    dxSpecies[2][45] += - 1.0/3.0;
    dxSpecies[2][46] += - 2.0/3.0;
    dxSpecies[2][52] +=   1.0/3.0;
    dxSpecies[2][62] += - 1.0/3.0;
    dxSpecies[2][ 3] += - 1.0/3.0;
    dxSpecies[2][ 5] +=   1.0/3.0;
    dxSpecies[2][18] += - 1.0/3.0;
    dxSpecies[2][20] += - 1.0/3.0;
    dxSpecies[2][22] +=   1.0/3.0;
    dxSpecies[2][25] += - 1.0;
    dxSpecies[2][31] +=   2.0/3.0;
    dxSpecies[2][32] +=   1.0/3.0;
    dxSpecies[2][34] += - 1.0/3.0;
    dxSpecies[2][35] += - 1.0;
    dxSpecies[2][37] +=   2.0/3.0;
    dxSpecies[2][39] +=   5.0/3.0;
    dxSpecies[2][41] +=   8.0/3.0;
    dxSpecies[2][42] += - 1.0/3.0;
    dxSpecies[2][44] += - 2.0/3.0;
    dxSpecies[2][50] += - 1.0/3.0;
    dxSpecies[2][57] += - 1.0/3.0;
    dxSpecies[2][60] += - 2.0/3.0;
    dxSpecies[2][61] += - 2.0/3.0;
    dxSpecies[2][65] += - 4.0/3.0;
    dxSpecies[2][66] += - 4.0/3.0;
    dxSpecies[2][67] += - 1.0/3.0;
    dxSpecies[2][ 6] += - 1.0/3.0;
    dxSpecies[2][19] +=   1.0/3.0;
    dxSpecies[2][30] +=   5.0/3.0;
    dxSpecies[2][36] +=   8.0/3.0;
    dxSpecies[2][38] +=   11.0/3.0;
    dxSpecies[2][40] +=   14.0/3.0;
    dxSpecies[2][47] += - 1.0/3.0;
    dxSpecies[2][58] +=   2.0/3.0;
    dxSpecies[2][59] +=   2.0/3.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[3][i] = r[2]*DcoeffDs[i];
    dxSpecies[3][ 6] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[4][i] = r[3]*DcoeffDs[i];
    dxSpecies[4][10] += - 1.0;
    dxSpecies[4][11] += - 1.0;
    dxSpecies[4][ 9] += - 1.0;
    dxSpecies[4][ 7] += - 1.0;
    dxSpecies[4][ 8] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[5][i] = r[4]*DcoeffDs[i];
    dxSpecies[5][13] += - 1.0;
    dxSpecies[5][14] += - 1.0;
    dxSpecies[5][17] += - 1.0;
    dxSpecies[5][12] += - 1.0;
    dxSpecies[5][18] += - 1.0;
    dxSpecies[5][15] += - 1.0;
    dxSpecies[5][16] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[6][i] = r[5]*DcoeffDs[i];
    dxSpecies[6][20] += - 1.0;
    dxSpecies[6][19] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[7][i] = r[6]*DcoeffDs[i];
    dxSpecies[7][11] += - 1.0;
    dxSpecies[7][13] += - 1.0;
    dxSpecies[7][21] += - 1.0;
    dxSpecies[7][22] += - 2.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[8][i] = r[7]*DcoeffDs[i];
    dxSpecies[8][23] += - 1.0;
    dxSpecies[8][24] += - 1.0;
    dxSpecies[8][26] += - 2.0;
    dxSpecies[8][27] += - 2.0;
    dxSpecies[8][25] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[9][i] = r[8]*DcoeffDs[i];
    dxSpecies[9][28] += - 1.0;
    dxSpecies[9][29] += - 1.0;
    dxSpecies[9][43] += - 1.0;
    dxSpecies[9][45] += - 1.0;
    dxSpecies[9][46] += - 1.0;
    dxSpecies[9][18] += - 1.0;
    dxSpecies[9][31] += - 2.0;
    dxSpecies[9][32] += - 2.0;
    dxSpecies[9][33] += - 2.0;
    dxSpecies[9][34] += - 2.0;
    dxSpecies[9][35] += - 2.0;
    dxSpecies[9][37] += - 3.0;
    dxSpecies[9][39] += - 4.0;
    dxSpecies[9][41] += - 5.0;
    dxSpecies[9][42] += - 1.0;
    dxSpecies[9][44] += - 1.0;
    dxSpecies[9][50] += - 1.0;
    dxSpecies[9][57] += - 1.0;
    dxSpecies[9][67] += - 1.0;
    dxSpecies[9][30] += - 2.0;
    dxSpecies[9][36] += - 3.0;
    dxSpecies[9][38] += - 4.0;
    dxSpecies[9][40] += - 5.0;;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[10][i] = r[9]*DcoeffDs[i];
    dxSpecies[10][ 7] += - 1.0;
    dxSpecies[10][16] += - 1.0;
    dxSpecies[10][47] += - 1.0;
    dxSpecies[10][49] += - 1.0;
    dxSpecies[10][55] += - 1.0;
    dxSpecies[10][56] += - 2.0;
    dxSpecies[10][64] += - 1.0;
    dxSpecies[10][70] += - 1.0;
    dxSpecies[10][71] += - 1.0;
    dxSpecies[10][72] += - 2.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[11][i] = r[10]*DcoeffDs[i];
    dxSpecies[11][50] += - 1.0;
    dxSpecies[11][49] += - 1.0;
    dxSpecies[11][48] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[12][i] = r[11]*DcoeffDs[i];
    dxSpecies[12][52] += - 1.0;
    dxSpecies[12][53] += - 1.0;
    dxSpecies[12][51] += - 1.0;
    dxSpecies[12][57] += - 1.0;
    dxSpecies[12][55] += - 1.0;
    dxSpecies[12][54] += - 1.0;
    dxSpecies[12][56] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[13][i] = r[12]*DcoeffDs[i];
    dxSpecies[13][62] += - 1.0;
    dxSpecies[13][60] += - 2.0;
    dxSpecies[13][61] += - 1.0;
    dxSpecies[13][58] += - 1.0;
    dxSpecies[13][59] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[14][i] = r[13]*DcoeffDs[i];
    dxSpecies[14][65] += - 1.0;
    dxSpecies[14][66] += - 1.0;
    dxSpecies[14][67] += - 1.0;
    dxSpecies[14][64] += - 1.0;
    dxSpecies[14][63] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[15][i] = r[14]*DcoeffDs[i];
    dxSpecies[15][70] += - 1.0;
    dxSpecies[15][71] += - 1.0;
    dxSpecies[15][72] += - 1.0;
    dxSpecies[15][68] += - 1.0;
    dxSpecies[15][69] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[16][i] = r[15]*DcoeffDs[i];
    dxSpecies[16][73] += - 1.0;
    dxSpecies[16][74] += - 1.0;

    for (NSUInteger i=0; i<NS; i++) dxSpecies[NA+i][i] = 1.0;
}

- (double)D2xSpeciesDrDsForIndex:(NSUInteger)index andDr:(NSUInteger)i andDs:(NSUInteger)j {
    if (index == 0)          return -DcoeffDs[j];
    else if ((index-1) == i) return  DcoeffDs[j];
    else                     return 0.0;
}

- (void)fillnSpeciesWith:(double [NS])s {
    double coeff = [self speciationCoefficientWith:s];
    nSpecies = (coeff != 0.0) ? 1.0/coeff : 0.0;
}

- (void) DnSpeciesDsWith:(double [NS])s {
    double coeff = [self speciationCoefficientWith:s];
    if (coeff == 0.0) {
        memset(dnSpeciesds, 0, sizeof(dnSpeciesds));
    } else {
        [self DspeciationCoefficientDsWith:s];
        for (NSUInteger i=0; i<NS; i++) dnSpeciesds[i] = - DcoeffDs[i]/coeff/coeff;
    }
}

- (void) D2nSpeciesDs2With:(double [NS])s {
    double coeff = [self speciationCoefficientWith:s];
    if (coeff == 0.0) {
        memset(d2nSpeciesds2, 0, sizeof(d2nSpeciesds2));
    } else {
        [self DspeciationCoefficientDsWith:s];
        for (NSUInteger i=0; i<NS; i++) for (NSUInteger j=0; j<NS; j++)
            d2nSpeciesds2[i][j] = 2.0*DcoeffDs[i]*DcoeffDs[j]/coeff/coeff/coeff;
    }
}

#pragma mark -
#pragma mark Charge balance functions

- (double)chargeBalanceWith:(double [NS])s {
    double cations = s[0] + s[8] + s[13] + s[14] + 2.0*s[15] + s[16] + s[17] + 3.0*s[19] + s[48] + s[52] + s[53] + 2.0*s[54]
                   + s[55] + 2.0*s[58] + 3.0*s[59] + 2.0*s[63] + s[64] + 2.0*s[68] + 3.0*s[69] + s[70] + 2.0*s[71]
                   + 2.0*s[73] + 3.0*s[74];
    double anions = s[1] + 2.0*s[3] + s[4] + s[6] + s[9] + s[20] + s[21] + s[23] + 2.0*s[24] + 3.0*s[25] + s[26] + 2.0*s[27]
                  + s[29] + 2.0*s[30] + 2.0*s[31] + 2.0*s[32] + 2.0*s[33] + 2.0*s[34] + 2.0*s[35] + 2.0*s[36] + 2.0*s[37]
                  + 2.0*s[38] + 2.0*s[39] + 2.0*s[40] + 2.0*s[41] + 2.0*s[42] + s[43] + 2.0*s[44] + s[45] + s[46] + s[47]
                  + s[50] + 2.0*s[60] + 2.0*s[61] + s[62] + s[65] + 2.0*s[66];
    return (cations - anions)*nSpecies;
}

- (void)DchargeBalanceDsWith:(double [NS])s {
    double coeff = [self speciationCoefficientWith:s];
    memset(DcbDs, 0, sizeof(DcbDs));
    if (coeff == 0.0) return;

    double DcationsDs[NS], DanionsDs[NS];
    memset(DcationsDs, 0, sizeof(DcationsDs));
    memset(DanionsDs, 0, sizeof(DanionsDs));

    double cations = s[0] + s[8] + s[13] + s[14] + 2.0*s[15] + s[16] + s[17] + 3.0*s[19] + s[48] + s[52] + s[53] + 2.0*s[54]
    + s[55] + 2.0*s[58] + 3.0*s[59] + 2.0*s[63] + s[64] + 2.0*s[68] + 3.0*s[69] + s[70] + 2.0*s[71]
    + 2.0*s[73] + 3.0*s[74];
    DcationsDs[ 0] = 1.0;
    DcationsDs[ 8] = 1.0;
    DcationsDs[13] = 1.0;
    DcationsDs[14] = 1.0;
    DcationsDs[15] = 2.0;
    DcationsDs[16] = 1.0;
    DcationsDs[17] = 1.0;
    DcationsDs[19] = 3.0;
    DcationsDs[48] = 1.0;
    DcationsDs[52] = 1.0;
    DcationsDs[53] = 1.0;
    DcationsDs[54] = 2.0;
    DcationsDs[55] = 1.0;
    DcationsDs[58] = 2.0;
    DcationsDs[59] = 3.0;
    DcationsDs[63] = 2.0;
    DcationsDs[64] = 1.0;
    DcationsDs[68] = 2.0;
    DcationsDs[69] = 3.0;
    DcationsDs[70] = 1.0;
    DcationsDs[71] = 2.0;
    DcationsDs[73] = 2.0;
    DcationsDs[74] = 3.0;

    double anions = s[1] + 2.0*s[3] + s[4] + s[6] + s[9] + s[20] + s[21] + s[23] + 2.0*s[24] + 3.0*s[25] + s[26] + 2.0*s[27]
    + s[29] + 2.0*s[30] + 2.0*s[31] + 2.0*s[32] + 2.0*s[33] + 2.0*s[34] + 2.0*s[35] + 2.0*s[36] + 2.0*s[37]
    + 2.0*s[38] + 2.0*s[39] + 2.0*s[40] + 2.0*s[41] + 2.0*s[42] + s[43] + 2.0*s[44] + s[45] + s[46] + s[47]
    + s[50] + 2.0*s[60] + 2.0*s[61] + s[62] + s[65] + 2.0*s[66];
    DanionsDs[ 1] = 1.0;
    DanionsDs[ 3] = 2.0;
    DanionsDs[ 4] = 1.0;
    DanionsDs[ 6] = 1.0;
    DanionsDs[ 9] = 1.0;
    DanionsDs[20] = 1.0;
    DanionsDs[21] = 1.0;
    DanionsDs[23] = 1.0;
    DanionsDs[24] = 2.0;
    DanionsDs[25] = 3.0;
    DanionsDs[26] = 1.0;
    DanionsDs[27] = 2.0;
    DanionsDs[29] = 1.0;
    DanionsDs[30] = 2.0;
    DanionsDs[31] = 2.0;
    DanionsDs[32] = 2.0;
    DanionsDs[33] = 2.0;
    DanionsDs[34] = 2.0;
    DanionsDs[35] = 2.0;
    DanionsDs[36] = 2.0;
    DanionsDs[37] = 2.0;
    DanionsDs[38] = 2.0;
    DanionsDs[39] = 2.0;
    DanionsDs[40] = 2.0;
    DanionsDs[41] = 2.0;
    DanionsDs[42] = 2.0;
    DanionsDs[43] = 1.0;
    DanionsDs[44] = 2.0;
    DanionsDs[45] = 1.0;
    DanionsDs[46] = 1.0;
    DanionsDs[47] = 1.0;
    DanionsDs[50] = 1.0;
    DanionsDs[60] = 2.0;
    DanionsDs[61] = 2.0;
    DanionsDs[62] = 1.0;
    DanionsDs[65] = 1.0;
    DanionsDs[66] = 2.0;

    // (cations - anions)*nSpecies = (cations - anions)/coeff
    [self DspeciationCoefficientDsWith:s];

    for (NSUInteger i=0; i<NS; i++) DcbDs[i] = - (cations-anions)*DcoeffDs[i]/coeff/coeff + (DcationsDs[i]-DanionsDs[i])/coeff;
    return;
}

#pragma mark -
#pragma mark True ionic strength functions

- (double)trueIonicStrength {
    if (xSpecies[mIndH2O] == 0.0) return 0.0;
    double result = 0.0;
    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        id object = [endmembers objectAtIndex:i];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            HKFspeciesProperties *species = (HKFspeciesProperties *)object;
            double charge = species.charge;
            if (charge != 0.0) {
                double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
                result += molality*charge*charge;
            }
        }
    }
    return result/2.0;
}

- (void)DtrueIonicStrengthDrIn:(double [NR])didr {
    memset(didr, 0, sizeof(didr[0])*NR);
    if (xSpecies[mIndH2O] == 0.0) return;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        id object = [endmembers objectAtIndex:i];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            HKFspeciesProperties *species = (HKFspeciesProperties *)object;
            double charge = species.charge;
            if (charge != 0.0) {
                for (NSUInteger j=0; j<NR; j++) {
                    double DmolalityDr = dxSpeciesDr[i][j]*CapGam/xSpecies[mIndH2O]
                    - xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O],2.0);
                    didr[j] += DmolalityDr*charge*charge;
                }
            }
        }
    }
    for (NSUInteger i=0; i<NR; i++) didr[i] /= 2.0;
}

- (void)DtrueIonicStrengthDsIn:(double [NS])dids {
    memset(dids, 0, sizeof(dids[0])*NS);
    if (xSpecies[mIndH2O] == 0.0) return;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        id object = [endmembers objectAtIndex:i];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            HKFspeciesProperties *species = (HKFspeciesProperties *)object;
            double charge = species.charge;
            if (charge != 0.0) {
                for (NSUInteger j=0; j<NS; j++) {
                    double DmolalityDs = dxSpecies[i][j]*CapGam/xSpecies[mIndH2O]
                    - xSpecies[i]*CapGam*dxSpecies[mIndH2O][j]/pow(xSpecies[mIndH2O],2.0);
                    dids[j] += DmolalityDs*charge*charge;
                }
            }
        }
    }
    for (NSUInteger i=0; i<NS; i++) dids[i] /= 2.0;
}

- (void)D2trueIonicStrengthDr2In:(double [NR][NR])d2idr2 {
    memset(d2idr2, 0, sizeof(d2idr2[0][0])*NR*NR);
    if (xSpecies[mIndH2O] == 0.0) return;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        id object = [endmembers objectAtIndex:i];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            HKFspeciesProperties *species = (HKFspeciesProperties *)object;
            double charge = species.charge;
            if (charge != 0.0) {
                for (NSUInteger j=0; j<NR; j++) {
                    for (NSUInteger k=j; k<NR; k++) {
                        double D2molalityDr2 = - dxSpeciesDr[i][j]*dxSpeciesDr[mIndH2O][k]*CapGam/pow(xSpecies[mIndH2O], 2.0)
                                               - dxSpeciesDr[i][k]*dxSpeciesDr[mIndH2O][j]*CapGam/pow(xSpecies[mIndH2O], 2.0)
                          + 2.0*xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]*dxSpeciesDr[mIndH2O][k]/pow(xSpecies[mIndH2O],3.0);
                        d2idr2[j][k] += D2molalityDr2*charge*charge;
                        d2idr2[k][j] = d2idr2[j][k];
                    }
                }
            }
        }
    }
    for (NSUInteger i=0; i<NR; i++) for (NSUInteger j=0; j<NR; j++) d2idr2[i][j] /= 2.0;
}

- (void)D2trueIonicStrengthDrDsIn:(double [NR][NS])d2idrds {
    memset(d2idrds, 0, sizeof(d2idrds[0][0])*NR*NS);
    if (xSpecies[mIndH2O] == 0.0) return;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        id object = [endmembers objectAtIndex:i];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            HKFspeciesProperties *species = (HKFspeciesProperties *)object;
            double charge = species.charge;
            if (charge != 0.0) {
                for (NSUInteger j=0; j<NR; j++) {
                    for (NSUInteger k=0; k<NS; k++) {
                        double d2xSpeciesDrDs = [self D2xSpeciesDrDsForIndex:i andDr:j andDs:k];
                        double d2xH2ODrDs     = [self D2xSpeciesDrDsForIndex:mIndH2O andDr:j andDs:k];

                        double D2molalityDrDs = d2xSpeciesDrDs*CapGam/xSpecies[mIndH2O]
                                              - dxSpeciesDr[i][j]*dxSpecies[mIndH2O][k]*CapGam/pow(xSpecies[mIndH2O],2.0)
                                              - dxSpecies[i][k]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O],2.0)
                                              - xSpecies[i]*CapGam*d2xH2ODrDs/pow(xSpecies[mIndH2O],2.0)
                                              + 2.0*xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O],3.0);
                        d2idrds[j][k] += D2molalityDrDs*charge*charge;
                    }
                }
            }
        }
    }
    for (NSUInteger i=0; i<NR; i++) for (NSUInteger k=0; k<NS; k++) d2idrds[i][k] /= 2.0;
}


- (void)D2trueIonicStrengthDs2In:(double [NS][NS])d2ids2 {
    memset(d2ids2, 0, sizeof(d2ids2[0][0])*NS*NS);
    if (xSpecies[mIndH2O] == 0.0) return;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        id object = [endmembers objectAtIndex:i];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            HKFspeciesProperties *species = (HKFspeciesProperties *)object;
            double charge = species.charge;
            if (charge != 0.0) {
                for (NSUInteger j=0; j<NS; j++) {
                    for (NSUInteger k=j; k<NS; k++) {
                        double D2molalityDs2 = - dxSpecies[i][j]*dxSpecies[mIndH2O][k]*CapGam/pow(xSpecies[mIndH2O], 2.0)
                                               - dxSpecies[i][k]*dxSpecies[mIndH2O][j]*CapGam/pow(xSpecies[mIndH2O], 2.0)
                          + 2.0*xSpecies[i]*CapGam*dxSpecies[mIndH2O][j]*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O],3.0);
                        d2ids2[j][k] += D2molalityDs2*charge*charge;
                        d2ids2[k][j] = d2ids2[j][k];
                    }
                }
            }
        }
    }
    for (NSUInteger i=0; i<NS; i++) for (NSUInteger j=0; j<NS; j++) d2ids2[i][j] /= 2.0;
}

#pragma mark -
#pragma mark Debye-Huckel terms

// azero is a function ONLY of bulk composition, NOT speciation and NOT t and p
// it is called only from the ordering routine if composition of the fluid changes
- (void)azeroWithR:(double [NR])r {
    //    X1   CO2     r0
    //    X2   O2      r1
    //    X3   HF      r2
    //    X4   NaOH    r3
    //    X5   Mg(OH)2 r4
    //    X6   HAlO2   r5
    //    X7   SiO2    r6 - ignored, assumed neutral
    //    X8   H3PO4   r7
    //    X9   SO2     r8
    //    X10  HCl     r9
    //    X11  KOH     r10
    //    X12  Ca(OH)2 r11
    //    X13  H2CrO4  r12
    //    X14  Mn(OH)2 r13
    //    X15  Fe(OH)2 r14
    //    X16  Co(OH)2 r15
    double CO3minus2  = r[0];
    double Fminus     = r[2];
    double Naplus     = r[3];
    double Mgplus2    = r[4];
    double Alplus3    = r[5];
    double PO4minus3  = r[7];
    double SO4minus2  = 3.0*r[8]/4.0;
    double HSminus    = r[8]/4.0;
    double Clminus    = r[9];
    double Kplus      = r[10];
    double Caplus2    = r[11];
    double CrO4minus2 = r[12];
    double Mnplus2    = r[13] + r[15]; // pretend Co++ is Mn++
    double Feplus2    = r[14];

    double Hplus   = 2.0*r[0] + r[2] + 3.0*r[7] + 3.0*r[8]/2.0 + r[8]/4.0 + r[9] + 2.0*r[12];
    double OHminus = r[3] + 2.0*r[4] + 3.0*r[5] + r[10] + 2.0*r[11] + 2.0*r[13] + 2.0*r[14] + 2.0*r[15];

    aZeroSum = CO3minus2 + Fminus + Naplus + Mgplus2 + Alplus3 + PO4minus3 + SO4minus2 + HSminus + Clminus + Kplus
             + Caplus2 + CrO4minus2 + Mnplus2 + Feplus2 + Hplus + OHminus;

    if (aZeroSum == 0.0) { aZero = 1.0e-8; return; }

    // copnstants from Table 7, Helgeson, Kirkham and Flowers (1976); PO4-3 set to ClO4-
    double result = CO3minus2*2.81 + Fminus*1.33 + Naplus*1.91 + Mgplus2*2.54 + Alplus3*3.33 + PO4minus3*3.59 + SO4minus2*3.15
                  + HSminus*1.84 + Clminus*1.81 + Kplus*2.27 + Caplus2*2.87 + CrO4minus2*3.59 + Mnplus2*2.68 + Feplus2*2.62
                  + Hplus*3.08 + OHminus*1.40;

    // eqn 125 Helgeson, Kirkham and Flowers (1980), p 1297
    result *= 2.0/aZeroSum;
    aZero = (result != 0.0) ? result*1.0e-8: 1.0e-8;
}

- (double)DazeroDrForIndex:(NSUInteger)i {
    if (aZero == 1.0e-8) return 0.0;

    double dOHminusDr = 0.0, dHplusDr = 0.0, DsumDr = 0.0, DresultDr = 0.0;
    switch (i) {
        case 0: // CO3minus2
        {
            dHplusDr   =  2.0;;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.81 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 2: // Fminus
        {
            dHplusDr   =  1.0;;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 1.33 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 3: // Naplus
        {
            dOHminusDr =  1.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 1.91 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 4: // Mgplus2
        {
            dOHminusDr =  2.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.54 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 5: // Alplus3
        {
            dOHminusDr =  3.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 3.33 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 7: // PO4minus3
        {
            dHplusDr   =  3.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 3.59 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 8: // SO4minus2 and HSminus
        {
            dHplusDr   =  2.0*3.0/4.0 + 1.0/4.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 3.15*3.0/4.0 + 1.84/4.0 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 9: // Clminus
        {
            dHplusDr   =  1.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 1.81 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 10: // Kplus
        {
            dOHminusDr =  1.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.27 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 11: // Caplus2
        {
            dOHminusDr =  2.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.87 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 12: // CrO4minus2
        {
            dHplusDr   =  2.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 3.59 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 13: // Mnplus2
        {
            dOHminusDr =  2.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.68 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 14: // Feplus2
        {
            dOHminusDr =  2.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.62 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        case 15: // Mnplus2
        {
            dOHminusDr =  2.0;
            DsumDr = 1.0 + dHplusDr + dOHminusDr;
            DresultDr = 2.68 + dHplusDr*3.08 + dOHminusDr*1.40;
            break;
        }
        default:
            break;
    }
    return DresultDr*2.0*1.0e-8/aZeroSum - aZero*DsumDr/aZeroSum;
}

- (double)D2azeroDr2ForIndex:(NSUInteger)i andIndex:(NSUInteger)j {
    if (aZero == 1.0) return 0.0;

    double dOHminusDrI = 0.0, dHplusDrI = 0.0, DsumDrI = 0.0, DresultDrI = 0.0;
    switch (i) {
        case 0:
        {
            dHplusDrI   =  2.0;;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.81 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 2:
        {
            dHplusDrI   =  1.0;;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 1.33 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 3:
        {
            dOHminusDrI =  1.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 1.91 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 4:
        {
            dOHminusDrI =  2.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.54 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 5:
        {
            dOHminusDrI =  3.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 3.33 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 7:
        {
            dHplusDrI   =  3.0;;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 3.59 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 8:
        {
            dHplusDrI   =  2.0*3.0/4.0 + 1.0/4.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 3.15*3.0/4.0 + 1.84/4.0 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 9:
        {
            dHplusDrI   =  1.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 1.81 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 10:
        {
            dOHminusDrI =  1.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.27 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 11:
        {
            dOHminusDrI =  2.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.87 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 12:
        {
            dHplusDrI   =  2.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 3.59 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 13:
        {
            dOHminusDrI =  2.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.68 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 14:
        {
            dOHminusDrI =  2.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.62 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        case 15:
        {
            dOHminusDrI =  2.0;
            DsumDrI = 1.0 + dHplusDrI + dOHminusDrI;
            DresultDrI = 2.68 + dHplusDrI*3.08 + dOHminusDrI*1.40;
            break;
        }
        default:
            break;
    }

    double dOHminusDrJ = 0.0, dHplusDrJ = 0.0, DsumDrJ = 0.0, DresultDrJ = 0.0;
    switch (j) {
        case 0:
        {
            dHplusDrJ   =  2.0;;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.81 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 2:
        {
            dHplusDrJ   =  1.0;;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 1.33 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 3:
        {
            dOHminusDrJ =  1.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 1.91 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 4:
        {
            dOHminusDrJ =  2.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.54 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 5:
        {
            dOHminusDrJ =  3.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 3.33 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 7:
        {
            dHplusDrJ   =  3.0;;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 3.59 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 8:
        {
            dHplusDrJ   =  2.0*3.0/4.0 + 1.0/4.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 3.15*3.0/4.0 + 1.84/4.0 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 9:
        {
            dHplusDrJ   =  1.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 1.81 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 10:
        {
            dOHminusDrJ =  1.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.27 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 11:
        {
            dOHminusDrJ =  2.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.87 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 12:
        {
            dHplusDrJ   =  2.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 3.59 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 13:
        {
            dOHminusDrJ =  2.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.68 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 14:
        {
            dOHminusDrJ =  2.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.62 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        case 15:
        {
            dOHminusDrJ =  2.0;
            DsumDrJ = 1.0 + dHplusDrJ + dOHminusDrJ;
            DresultDrJ = 2.68 + dHplusDrJ*3.08 + dOHminusDrJ*1.40;
            break;
        }
        default:
            break;
    }

    return - DresultDrI*2.0*DsumDrJ*1.0e-8/pow(aZeroSum, 2.0) - DresultDrJ*2.0*DsumDrI*1.0e-8/pow(aZeroSum, 2.0)
           + 2.0*aZero*DsumDrJ*DsumDrI/pow(aZeroSum, 2.0);
}

// molality of i = (nSpecies X[i] 1000 gm solvent) / (nSpecies X[H2O] MW_H2O)
// eq 120 (HKF 80)  mu - mu0 = (psi_i A_G sqrt (I)) / CapLambda + CapGam_G
// A_G = -2 (2.303) R T A_gamma
// CapGam_G = - R T ln (1 + sum_m/55.51)
// CapLambda = 1 + a^0 B_gamma sqrt(I)
// psi_i = z^2 / 2

- (double)fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:(NSUInteger)index
                                                   withT:(double)t
                                                andWithP:(double)p {
    NSAssert(index < NE, @"DEWFluid:fillDebyeHuckelExcessFreeEnergyForSpeciesIndex: with bad index.");
    double result = 0.0;
    if (initializeDHexcessTermsZerothOrder) {
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) mStar += xSpecies[i];
        if (mStar == 0.0) return 0.0;
        mStar *= CapGam/xSpecies[mIndH2O];

        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];

        CapGamSubG = - R*t*log(1.0+mStar/CapGam);
        sqrtOfI = sqrt([self trueIonicStrength]);
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;

        initializeDHexcessTermsZerothOrder = NO;
    }
    if (xSpecies[index] == 0.0) return result;

    result += CapGamSubG;
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        result += - AsubG*pow(sqrtOfI, 3.0)*sigma/3.0/CapGam;
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)object charge];
            if (z != 0.0) result += z*z*AsubG*sqrtOfI/2.0/CapLambda;
        }
    }

    return result;
}

- (void)fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:(NSUInteger)index
                                                       in:(double [NR])DdhMuExcessDr
                                                    withT:(double)t
                                                 andWithP:(double)p {
    NSAssert(index < NE, @"DEWFluid:fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex: with bad index.");
    memset(DdhMuExcessDr, 0, sizeof(DdhMuExcessDr[0])*NR);

    if (initializeDHexcessTermsFirstOrderDr) {
        memset(DmStarDr, 0, sizeof(DmStarDr));
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) {
            mStar += xSpecies[i];
            for (NSUInteger j=0; j<NR; j++) DmStarDr[j] += dxSpeciesDr[i][j];
        }
        if (mStar == 0.0) return;

        mStar *= CapGam/xSpecies[mIndH2O];
        for (NSUInteger j=0; j<NR; j++)
            DmStarDr[j] = DmStarDr[j]*CapGam/xSpecies[mIndH2O] - mStar*dxSpeciesDr[mIndH2O][j]/xSpecies[mIndH2O];

        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        sqrtOfI = sqrt([self trueIonicStrength]);
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;

        [self DtrueIonicStrengthDrIn:DsqrtOfIdr];
        for (NSUInteger j=0; j<NR; j++) {
            dCapGamSubGdr[j]  = - R*t*(DmStarDr[j]/CapGam)/(1.0+mStar/CapGam);
            DsqrtOfIdr[j]    /= 2.0*sqrtOfI;
            double DazeroDr   = [self DazeroDrForIndex:j];
            dCapLambdaDr[j]   = DazeroDr*Bgamma*sqrtOfI + aZero*Bgamma*DsqrtOfIdr[j];
        }

        initializeDHexcessTermsFirstOrderDr = NO;
    }
    if (xSpecies[index] == 0.0) return;

    for (NSUInteger j=0; j<NR; j++) DdhMuExcessDr[j] += dCapGamSubGdr[j];
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        for (NSUInteger j=0; j<NR; j++) {
            double DsigmaDr = 3.0*(dCapLambdaDr[j] + dCapLambdaDr[j]/pow(CapLambda,2.0)
                                   - 2.0*dCapLambdaDr[j]/CapLambda)/pow(CapLambda-1.0, 3.0)
                            - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDr[j]/pow(CapLambda-1.0, 4.0);
            DdhMuExcessDr[j] += - 3.0*AsubG*DsqrtOfIdr[j]*pow(sqrtOfI, 2.0)*sigma/3.0/CapGam
                              - AsubG*pow(sqrtOfI, 3.0)*DsigmaDr/3.0/CapGam;
        }
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0) for (NSUInteger j=0; j<NR; j++) {
                DdhMuExcessDr[j] += z*z*AsubG*DsqrtOfIdr[j]/2.0/CapLambda
                                  - z*z*AsubG*sqrtOfI*dCapLambdaDr[j]/2.0/pow(CapLambda, 2.0);
            }
        }
    }
}

- (void)fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:(NSUInteger)index
                                                       in:(double [NS])DdhMuExcessDs
                                                    withT:(double)t
                                                 andWithP:(double)p {
    NSAssert(index < NE, @"DEWFluid:fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex: with bad index.");
    memset(DdhMuExcessDs, 0, sizeof(DdhMuExcessDs[0])*NS);

    if (initializeDHexcessTermsFirstOrderDs) {
        memset(DmStarDs, 0, sizeof(DmStarDs));
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) {
            mStar += xSpecies[i];
            for (NSUInteger j=0; j<NS; j++) DmStarDs[j] += dxSpecies[i][j];
        }
        if (mStar == 0.0) return;

        mStar *= CapGam/xSpecies[mIndH2O];
        for (NSUInteger j=0; j<NS; j++)
            DmStarDs[j] = DmStarDs[j]*CapGam/xSpecies[mIndH2O] - mStar*dxSpecies[mIndH2O][j]/xSpecies[mIndH2O];

        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        sqrtOfI = sqrt([self trueIonicStrength]);
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;

        [self DtrueIonicStrengthDsIn:DsqrtOfIds];
        for (NSUInteger j=0; j<NS; j++) {
            dCapGamSubGds[j]  = - R*t*(DmStarDs[j]/CapGam)/(1.0+mStar/CapGam);
            DsqrtOfIds[j]    /= 2.0*sqrtOfI;
            dCapLambdaDs[j]   = aZero*Bgamma*DsqrtOfIds[j];
        }

        initializeDHexcessTermsFirstOrderDs = NO;
    }
    if (xSpecies[index] == 0.0) return;

    for (NSUInteger j=0; j<NS; j++) DdhMuExcessDs[j] += dCapGamSubGds[j];
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        for (NSUInteger j=0; j<NS; j++) {
            double DsigmaDs = 3.0*(dCapLambdaDs[j] + dCapLambdaDs[j]/pow(CapLambda,2.0)
                                   - 2.0*dCapLambdaDs[j]/CapLambda)/pow(CapLambda-1.0, 3.0)
                            - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDs[j]/pow(CapLambda-1.0, 4.0);
            DdhMuExcessDs[j] += - 3.0*AsubG*DsqrtOfIds[j]*pow(sqrtOfI, 2.0)*sigma/3.0/CapGam
                                - AsubG*pow(sqrtOfI, 3.0)*DsigmaDs/3.0/CapGam;
        }
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0) for (NSUInteger j=0; j<NS; j++) {
                DdhMuExcessDs[j] += z*z*AsubG*DsqrtOfIds[j]/2.0/CapLambda
                                  - z*z*AsubG*sqrtOfI*dCapLambdaDs[j]/2.0/pow(CapLambda, 2.0);
            }
        }
    }
}

- (double)fillDebyeHuckelExcessFreeEnergyDtForSpeciesIndex:(NSUInteger)index
                                                   withT:(double)t
                                                andWithP:(double)p {
    NSAssert(index < NE, @"DEWFluid:fillDebyeHuckelExcessFreeEnergyDtForSpeciesIndex: with bad index.");
    double result = 0.0;
    if (initializeDHexcessTermsFirstOrderDt) {
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) mStar += xSpecies[i];
        if (mStar == 0.0) return 0.0;
        mStar *= CapGam/xSpecies[mIndH2O];

        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        double AsubH = [dielectricConstant AsubHfromT:t andP:p];
        dAsubGDt = AsubG/t - AsubH/t;
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        dBgammaDt = [dielectricConstant BsubHfromT:t andP:p]/2.0/R/log(10.0)/t/t;

        CapGamSubG = - R*t*log(1.0+mStar/CapGam);
        dCapGamSubGDt = - R*log(1.0+mStar/CapGam);
        sqrtOfI = sqrt([self trueIonicStrength]);
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;
        dCapLambdaDt = aZero*dBgammaDt*sqrtOfI;

        initializeDHexcessTermsFirstOrderDt = NO;
    }
    if (xSpecies[index] == 0.0) return result;

    result += dCapGamSubGDt;
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        double DsigmaDt = 3.0*(dCapLambdaDt + dCapLambdaDt/CapLambda/CapLambda - 2.0*dCapLambdaDt/CapLambda)/pow(CapLambda-1.0, 3.0)
                        - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDt/pow(CapLambda-1.0, 4.0);
        result += - dAsubGDt*pow(sqrtOfI, 3.0)*sigma/3.0/CapGam
                  - AsubG*pow(sqrtOfI, 3.0)*DsigmaDt/3.0/CapGam;
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0) result += z*z*dAsubGDt*sqrtOfI/2.0/CapLambda
                                  - z*z*AsubG*sqrtOfI*dCapLambdaDt/2.0/CapLambda/CapLambda;
        }
    }

    return result;
}

- (double)fillDebyeHuckelExcessFreeEnergyDpForSpeciesIndex:(NSUInteger)index
                                                   withT:(double)t
                                                andWithP:(double)p {
    NSAssert(index < NE, @"DEWFluid:fillDebyeHuckelExcessFreeEnergyDpForSpeciesIndex: with bad index.");
    double result = 0.0;
    if (initializeDHexcessTermsZerothOrderDp) {
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) mStar += xSpecies[i];
        if (mStar == 0.0) return 0.0;
        mStar *= CapGam/xSpecies[mIndH2O];

        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        dAsubGDp = [dielectricConstant AsubVfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        dBgammaDp = -[dielectricConstant BsubVfromT:t andP:p]/2.0/log(10.0)/R/t;

        CapGamSubG = - R*t*log(1.0+mStar/CapGam);
        sqrtOfI = sqrt([self trueIonicStrength]);
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;
        dCapLambdaDp = aZero*dBgammaDp*sqrtOfI;

        initializeDHexcessTermsZerothOrderDp = NO;
    }
    if (xSpecies[index] == 0.0) return result;

    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        double DsigmaDp = 3.0*(dCapLambdaDp + dCapLambdaDp/CapLambda/CapLambda - 2.0*dCapLambdaDp/CapLambda)/pow(CapLambda-1.0, 3.0)
                        - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDp/pow(CapLambda-1.0, 4.0);
        result += - dAsubGDp*pow(sqrtOfI, 3.0)*sigma/3.0/CapGam - AsubG*pow(sqrtOfI, 3.0)*DsigmaDp/3.0/CapGam;
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0) result += z*z*dAsubGDp*sqrtOfI/2.0/CapLambda - z*z*AsubG*sqrtOfI*dCapLambdaDp/2.0/CapLambda/CapLambda;
        }
    }

    return result;
}

- (void)fillD2debyeHuckelExcessFreeEnergyDr2ForSpeciesIndex:(NSUInteger)index
                                                         in:(double [NR][NR])D2dhMuExcessDr2
                                                      withT:(double)t
                                                   andWithP:(double)p {
    memset(D2dhMuExcessDr2, 0, sizeof(D2dhMuExcessDr2[0][0])*NR*NR);

    if (initializeDHexcessTermsSecondOrderDr2) {
        memset(DmStarDr, 0, sizeof(DmStarDr));
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) {
            mStar += xSpecies[i];
            for (NSUInteger j=0; j<NR; j++) DmStarDr[j] += dxSpeciesDr[i][j];
        }
        if (mStar == 0.0) return;

        mStar *= CapGam/xSpecies[mIndH2O];
        for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=j; k<NR; k++) {
            D2mStarDr2[j][k] = - DmStarDr[j]*CapGam*dxSpeciesDr[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0)
                               - DmStarDr[k]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0)
                             + 2.0*mStar*dxSpeciesDr[mIndH2O][j]*dxSpeciesDr[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0);
            D2mStarDr2[k][j] = D2mStarDr2[j][k];
        }
        for (NSUInteger j=0; j<NR; j++)
            DmStarDr[j] = DmStarDr[j]*CapGam/xSpecies[mIndH2O] - mStar*dxSpeciesDr[mIndH2O][j]/xSpecies[mIndH2O];

        sqrtOfI = sqrt([self trueIonicStrength]);
        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;

        [self DtrueIonicStrengthDrIn:DsqrtOfIdr];
        for (NSUInteger j=0; j<NR; j++) {
            dCapGamSubGdr[j]  = - R*t*(DmStarDr[j]/CapGam)/(1.0+mStar/CapGam);
            DsqrtOfIdr[j]    /= 2.0*sqrtOfI;
            double DazeroDr   = [self DazeroDrForIndex:j];
            dCapLambdaDr[j]   = DazeroDr*Bgamma*sqrtOfI + aZero*Bgamma*DsqrtOfIdr[j];
        }

        [self D2trueIonicStrengthDr2In:D2sqrtOfIdr2];
        for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=j; k<NR; k++) {
            d2CapGamSubGdr2[j][k] = - R*t*(D2mStarDr2[j][k]/CapGam)/(1.0+mStar/CapGam)
                                  + R*t*(DmStarDr[j]/CapGam)*(DmStarDr[k]/CapGam)/pow(1.0+mStar/CapGam, 2.0);
            d2CapGamSubGdr2[k][j] = d2CapGamSubGdr2[j][k];
        }
        for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=j; k<NR; k++) {
            D2sqrtOfIdr2[j][k] = - DsqrtOfIdr[j]*DsqrtOfIdr[k]/sqrtOfI + D2sqrtOfIdr2[j][k]/2.0/sqrtOfI;
            D2sqrtOfIdr2[k][j] = D2sqrtOfIdr2[j][k];
        }
        for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=j; k<NR; k++) {
            double DazeroDrJ  = [self DazeroDrForIndex:j];
            double DazeroDrK  = [self DazeroDrForIndex:k];
            double D2azeroDr2 = [self D2azeroDr2ForIndex:j andIndex:k];
            d2CapLambdaDr2[j][k] = D2azeroDr2*Bgamma*sqrtOfI + DazeroDrJ*Bgamma*DsqrtOfIdr[k]
                                 + DazeroDrK*Bgamma*DsqrtOfIdr[j] + aZero*Bgamma*D2sqrtOfIdr2[j][k];
            d2CapLambdaDr2[k][j] = d2CapLambdaDr2[j][k];
        }

        initializeDHexcessTermsSecondOrderDr2 = NO;
    }
    if (xSpecies[index] == 0.0) return;

    for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=j; k<NR; k++) {
        D2dhMuExcessDr2[j][k] += d2CapGamSubGdr2[j][k];
        D2dhMuExcessDr2[k][j]  = D2dhMuExcessDr2[j][k];
    }
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        for (NSUInteger j=0; j<NR; j++) {
            double DsigmaDrj = 3.0*(dCapLambdaDr[j] + dCapLambdaDr[j]/pow(CapLambda,2.0)
                             - 2.0*dCapLambdaDr[j]/CapLambda)/pow(CapLambda-1.0, 3.0)
                             - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDr[j]
                                    /pow(CapLambda-1.0, 4.0);
            for (NSUInteger k=j; k<NR; k++) {
                double DsigmaDrk = 3.0*(dCapLambdaDr[k] + dCapLambdaDr[k]/pow(CapLambda,2.0)
                                 - 2.0*dCapLambdaDr[k]/CapLambda)/pow(CapLambda-1.0, 3.0)
                                 - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDr[k]
                                        /pow(CapLambda-1.0, 4.0);
                double D2sigmaDrjk = 3.0*(d2CapLambdaDr2[j][k] + d2CapLambdaDr2[j][k]/pow(CapLambda,2.0)
                                          - 2.0*dCapLambdaDr[j]*dCapLambdaDr[k]/pow(CapLambda,3.0)
                                          - 2.0*d2CapLambdaDr2[j][k]/CapLambda + 2.0*dCapLambdaDr[j]*dCapLambdaDr[k]/pow(CapLambda, 2.0))
                                          /pow(CapLambda-1.0, 3.0)
                                   - 9.0*(dCapLambdaDr[j] + dCapLambdaDr[j]/pow(CapLambda,2.0)
                                          - 2.0*dCapLambdaDr[j]/CapLambda)*dCapLambdaDr[k]/pow(CapLambda-1.0, 4.0)
                                   - 9.0*(dCapLambdaDr[k] + dCapLambdaDr[k]/pow(CapLambda,2.0)
                                          - 2.0*dCapLambdaDr[k]/CapLambda)*dCapLambdaDr[j]/pow(CapLambda-1.0, 4.0)
                                   - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*d2CapLambdaDr2[j][k]/pow(CapLambda-1.0, 4.0)
                                   + 36.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))
                                          *dCapLambdaDr[j]*dCapLambdaDr[k]/pow(CapLambda-1.0, 5.0);

                D2dhMuExcessDr2[j][k] += - 3.0*AsubG*D2sqrtOfIdr2[j][k]*pow(sqrtOfI, 2.0)*sigma/3.0/CapGam
                                         - 6.0*AsubG*DsqrtOfIdr[j]*DsqrtOfIdr[k]*sqrtOfI*sigma/3.0/CapGam
                                         - 3.0*AsubG*DsqrtOfIdr[j]*pow(sqrtOfI, 2.0)*DsigmaDrk/3.0/CapGam
                                         - 3.0*AsubG*DsqrtOfIdr[k]*pow(sqrtOfI, 2.0)*DsigmaDrj/3.0/CapGam
                                         - AsubG*pow(sqrtOfI, 3.0)*D2sigmaDrjk/3.0/CapGam;
                D2dhMuExcessDr2[k][j] = D2dhMuExcessDr2[j][k];
            }
        }
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0)
                for (NSUInteger j=0; j<NR; j++) {
                    for (NSUInteger k=j; k<NR; k++) {
                        D2dhMuExcessDr2[j][k] += z*z*AsubG*D2sqrtOfIdr2[j][k]/2.0/CapLambda
                                               - z*z*AsubG*DsqrtOfIdr[j]*dCapLambdaDr[k]/2.0/pow(CapLambda, 2.0)
                                               - z*z*AsubG*DsqrtOfIdr[k]*dCapLambdaDr[j]/2.0/pow(CapLambda, 2.0)
                                               - z*z*AsubG*sqrtOfI*d2CapLambdaDr2[j][k]/2.0/pow(CapLambda, 2.0)
                                               + 2.0*z*z*AsubG*sqrtOfI*dCapLambdaDr[j]*dCapLambdaDr[k]/2.0/pow(CapLambda, 3.0);
                        D2dhMuExcessDr2[k][j] = D2dhMuExcessDr2[j][k];
                    }
                }
        }
    }
}

- (void)fillD2debyeHuckelExcessFreeEnergyDrDsForSpeciesIndex:(NSUInteger)index
                                                       in:(double [NR][NS])D2dhMuExcessDrDs
                                                    withT:(double)t
                                                 andWithP:(double)p {
    NSAssert(index < NE, @"DEWFluid:fillD2debyeHuckelExcessFreeEnergyDrDsForSpeciesIndex: with bad index.");
    memset(D2dhMuExcessDrDs, 0, sizeof(D2dhMuExcessDrDs[0][0])*NR*NS);

    if (initializeDHexcessTermsSecondOrderDrDs) {
        memset(DmStarDr, 0, sizeof(DmStarDr));
        memset(DmStarDs, 0, sizeof(DmStarDs));
        memset(D2mStarDrDs, 0, sizeof(D2mStarDrDs));

        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) {
            mStar += xSpecies[i];
            for (NSUInteger j=0; j<NR; j++) DmStarDr[j] += dxSpeciesDr[i][j];
            for (NSUInteger j=0; j<NS; j++) DmStarDs[j] += dxSpecies[i][j];
        }
        if (mStar == 0.0) return;
        mStar *= CapGam/xSpecies[mIndH2O];
        for (NSUInteger j=0; j<NR; j++) {
            for (NSUInteger k=0; k<NS; k++) {
                for (NSUInteger i=1; i<NE; i++) D2mStarDrDs[j][k] += [self D2xSpeciesDrDsForIndex:i andDr:j andDs:k];
                double d2xH2ODrDs = [self D2xSpeciesDrDsForIndex:mIndH2O andDr:j andDs:k];
                D2mStarDrDs[j][k] = D2mStarDrDs[j][k]*CapGam/xSpecies[mIndH2O]
                                  - DmStarDr[j]*CapGam*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0)
                                  - DmStarDs[k]*dxSpeciesDr[mIndH2O][j]*CapGam/pow(xSpecies[mIndH2O],2.0)
                                  - mStar*d2xH2ODrDs/xSpecies[mIndH2O]
                                  + 2.0*mStar*dxSpeciesDr[mIndH2O][j]*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0);
            }
        }
        for (NSUInteger j=0; j<NR; j++)
            DmStarDr[j] = DmStarDr[j]*CapGam/xSpecies[mIndH2O] - mStar*dxSpeciesDr[mIndH2O][j]/xSpecies[mIndH2O];
        for (NSUInteger j=0; j<NS; j++)
            DmStarDs[j] = DmStarDs[j]*CapGam/xSpecies[mIndH2O] - mStar*dxSpecies[mIndH2O][j]/xSpecies[mIndH2O];

        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        sqrtOfI = sqrt([self trueIonicStrength]);
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;

        [self DtrueIonicStrengthDrIn:DsqrtOfIdr];
        [self DtrueIonicStrengthDsIn:DsqrtOfIds];
        for (NSUInteger j=0; j<NS; j++) DsqrtOfIds[j] /= 2.0*sqrtOfI; // must be done here!
        [self D2trueIonicStrengthDrDsIn:D2sqrtOfIdrds];

        for (NSUInteger j=0; j<NR; j++) {
            dCapGamSubGdr[j]  = - R*t*(DmStarDr[j]/CapGam)/(1.0+mStar/CapGam);
            DsqrtOfIdr[j]    /= 2.0*sqrtOfI;
            double DazeroDr   = [self DazeroDrForIndex:j];
            dCapLambdaDr[j]   = DazeroDr*Bgamma*sqrtOfI + aZero*Bgamma*DsqrtOfIdr[j];

            for (NSUInteger k=0; k<NS; k++) {
                d2CapGamSubGdrds[j][k] = - R*t*(D2mStarDrDs[j][k]/CapGam)/(1.0+mStar/CapGam)
                                       + R*t*(DmStarDr[j]/CapGam)*(DmStarDs[k]/CapGam)/pow(1.0+mStar/CapGam, 2.0);
                D2sqrtOfIdrds[j][k]    = - DsqrtOfIdr[j]*DsqrtOfIds[k]/sqrtOfI + D2sqrtOfIdrds[j][k]/2.0/sqrtOfI;
                d2CapLambdaDrDs[j][k]  = DazeroDr*Bgamma*DsqrtOfIds[k] + aZero*Bgamma*D2sqrtOfIdrds[j][k];
            }
        }

        for (NSUInteger j=0; j<NS; j++) {
            dCapGamSubGds[j]  = - R*t*(DmStarDs[j]/CapGam)/(1.0+mStar/CapGam);
            dCapLambdaDs[j]   = aZero*Bgamma*DsqrtOfIds[j];
        }

        initializeDHexcessTermsSecondOrderDrDs = NO;
    }
    if (xSpecies[index] == 0.0) return;

    for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=0; k<NS; k++) D2dhMuExcessDrDs[j][k] += d2CapGamSubGdrds[j][k];
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        for (NSUInteger j=0; j<NR; j++) {
            double DsigmaDr = 3.0*(dCapLambdaDr[j] + dCapLambdaDr[j]/pow(CapLambda,2.0)
                                   - 2.0*dCapLambdaDr[j]/CapLambda)/pow(CapLambda-1.0, 3.0)
            - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDr[j]/pow(CapLambda-1.0, 4.0);
            for (NSUInteger k=0; k<NS; k++) {
                double DsigmaDs = 3.0*(dCapLambdaDs[k] + dCapLambdaDs[k]/pow(CapLambda,2.0)
                                       - 2.0*dCapLambdaDs[k]/CapLambda)/pow(CapLambda-1.0, 3.0)
                       - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDs[k]/pow(CapLambda-1.0, 4.0);
                double D2sigmaDrDs =
                  3.0*(d2CapLambdaDrDs[j][k] + d2CapLambdaDrDs[j][k]/pow(CapLambda,2.0)
                       - 2.0*dCapLambdaDr[j]*dCapLambdaDs[k]/pow(CapLambda,3.0)
                       + 2.0*dCapLambdaDr[j]*dCapLambdaDs[k]/pow(CapLambda, 2.0) - 2.0*d2CapLambdaDrDs[j][k]/CapLambda)
                  /pow(CapLambda-1.0, 3.0)
                - 9.0*(dCapLambdaDr[j] + dCapLambdaDr[j]/pow(CapLambda,2.0) - 2.0*dCapLambdaDr[j]/CapLambda)
                  *dCapLambdaDs[k]/pow(CapLambda-1.0, 4.0)

                - 9.0*(dCapLambdaDs[k] + dCapLambdaDs[k]/pow(CapLambda,2.0)
                       - 2.0*dCapLambdaDs[k]/CapLambda)*dCapLambdaDr[j]/pow(CapLambda-1.0, 4.0)
                - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*d2CapLambdaDrDs[j][k]/pow(CapLambda-1.0, 4.0)
                + 36.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDr[j]*dCapLambdaDs[k]/pow(CapLambda-1.0, 5.0);
                D2dhMuExcessDrDs[j][k] +=
                - 3.0*AsubG*D2sqrtOfIdrds[j][k]*pow(sqrtOfI, 2.0)*sigma/3.0/CapGam
                - 6.0*AsubG*DsqrtOfIdr[j]*DsqrtOfIds[k]*sqrtOfI*sigma/3.0/CapGam
                - 3.0*AsubG*DsqrtOfIdr[j]*pow(sqrtOfI, 2.0)*DsigmaDs/3.0/CapGam
                - 3.0*AsubG*DsqrtOfIds[k]*pow(sqrtOfI, 2.0)*DsigmaDr/3.0/CapGam
                - AsubG*pow(sqrtOfI, 3.0)*D2sigmaDrDs/3.0/CapGam;
            }
        }
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0) for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=0; k<NS; k++) {
                D2dhMuExcessDrDs[j][k] += z*z*AsubG*D2sqrtOfIdrds[j][k]/2.0/CapLambda
                                        - z*z*AsubG*DsqrtOfIdr[j]*dCapLambdaDs[k]/2.0/pow(CapLambda, 2.0)
                                        - z*z*AsubG*DsqrtOfIds[k]*dCapLambdaDr[j]/2.0/pow(CapLambda, 2.0)
                                        - z*z*AsubG*sqrtOfI*d2CapLambdaDrDs[j][k]/2.0/pow(CapLambda, 2.0)
                                        + 2.0*z*z*AsubG*sqrtOfI*dCapLambdaDr[j]*dCapLambdaDs[k]/2.0/pow(CapLambda, 3.0);
            }
        }
    }
}


- (void)fillD2debyeHuckelExcessFreeEnergyDs2ForSpeciesIndex:(NSUInteger)index
                                                         in:(double [NS][NS])D2dhMuExcessDs2
                                                      withT:(double)t
                                                   andWithP:(double)p {
    memset(D2dhMuExcessDs2, 0, sizeof(D2dhMuExcessDs2[0][0])*NS*NS);

    if (initializeDHexcessTermsSecondOrderDs2) {
        memset(DmStarDs, 0, sizeof(DmStarDs));
        memset(D2mStarDs2, 0, sizeof(D2mStarDs2));
        mStar = 0.0;
        for (NSUInteger i=1; i<NE; i++) {
            mStar += xSpecies[i];
            for (NSUInteger j=0; j<NS; j++) DmStarDs[j] += dxSpecies[i][j];
        }
        if (mStar == 0.0) return;

        mStar *= CapGam/xSpecies[mIndH2O];
        for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
            D2mStarDs2[j][k] = - DmStarDs[j]*CapGam*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0)
                               - DmStarDs[k]*CapGam*dxSpecies[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0)
                             + 2.0*mStar*dxSpecies[mIndH2O][j]*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0);
            D2mStarDs2[k][j] = D2mStarDs2[j][k];
        }
        for (NSUInteger j=0; j<NS; j++)
            DmStarDs[j] = DmStarDs[j]*CapGam/xSpecies[mIndH2O] - mStar*dxSpecies[mIndH2O][j]/xSpecies[mIndH2O];

        sqrtOfI = sqrt([self trueIonicStrength]);
        AsubG = [dielectricConstant AsubGfromT:t andP:p];
        Bgamma = [dielectricConstant BgammaFromT:t andP:p];
        CapLambda = 1.0 + aZero*Bgamma*sqrtOfI;

        [self DtrueIonicStrengthDsIn:DsqrtOfIds];
        for (NSUInteger j=0; j<NS; j++) {
            dCapGamSubGds[j]  = - R*t*(DmStarDs[j]/CapGam)/(1.0+mStar/CapGam);
            DsqrtOfIds[j]    /= 2.0*sqrtOfI;
            dCapLambdaDs[j]   = aZero*Bgamma*DsqrtOfIds[j];
        }

        [self D2trueIonicStrengthDs2In:D2sqrtOfIds2];
        for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
            d2CapGamSubGds2[j][k] = - R*t*(D2mStarDs2[j][k]/CapGam)/(1.0+mStar/CapGam)
                                  + R*t*(DmStarDs[j]/CapGam)*(DmStarDs[k]/CapGam)/pow(1.0+mStar/CapGam, 2.0);
            d2CapGamSubGds2[k][j] = d2CapGamSubGds2[j][k];
        }
        for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
            D2sqrtOfIds2[j][k] = - DsqrtOfIds[j]*DsqrtOfIds[k]/sqrtOfI + D2sqrtOfIds2[j][k]/2.0/sqrtOfI;
            D2sqrtOfIds2[k][j] = D2sqrtOfIds2[j][k];
        }
        for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
            d2CapLambdaDs2[j][k] = aZero*Bgamma*D2sqrtOfIds2[j][k];
            d2CapLambdaDs2[k][j] = d2CapLambdaDs2[j][k];
        }

        initializeDHexcessTermsSecondOrderDs2 = NO;
    }
    if (xSpecies[index] == 0.0) return;

    for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
        D2dhMuExcessDs2[j][k] += d2CapGamSubGds2[j][k];
        D2dhMuExcessDs2[k][j]  = D2dhMuExcessDs2[j][k];
    }
    if (index == mIndH2O) {
        double sigma = 3.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))/pow(CapLambda-1.0, 3.0);
        for (NSUInteger j=0; j<NS; j++) {
            double DsigmaDsj = 3.0*(dCapLambdaDs[j] + dCapLambdaDs[j]/pow(CapLambda,2.0)
                                    - 2.0*dCapLambdaDs[j]/CapLambda)/pow(CapLambda-1.0, 3.0)
                             - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDs[j]
                                    /pow(CapLambda-1.0, 4.0);
            for (NSUInteger k=j; k<NS; k++) {
                double DsigmaDsk = 3.0*(dCapLambdaDs[k] + dCapLambdaDs[k]/pow(CapLambda,2.0)
                                        - 2.0*dCapLambdaDs[k]/CapLambda)/pow(CapLambda-1.0, 3.0)
                                 - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*dCapLambdaDs[k]
                                        /pow(CapLambda-1.0, 4.0);
                double D2sigmaDsjk = 3.0*(d2CapLambdaDs2[j][k] + d2CapLambdaDs2[j][k]/pow(CapLambda,2.0)
                                          - 2.0*dCapLambdaDs[j]*dCapLambdaDs[k]/pow(CapLambda,3.0)
                                          - 2.0*d2CapLambdaDs2[j][k]/CapLambda + 2.0*dCapLambdaDs[j]*dCapLambdaDs[k]/pow(CapLambda, 2.0))
                                     /pow(CapLambda-1.0, 3.0)
                                   - 9.0*(dCapLambdaDs[j] + dCapLambdaDs[j]/pow(CapLambda,2.0)
                                          - 2.0*dCapLambdaDs[j]/CapLambda)*dCapLambdaDs[k]/pow(CapLambda-1.0, 4.0)
                                   - 9.0*(dCapLambdaDs[k] + dCapLambdaDs[k]/pow(CapLambda,2.0)
                                          - 2.0*dCapLambdaDs[k]/CapLambda)*dCapLambdaDs[j]/pow(CapLambda-1.0, 4.0)
                                   - 9.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))*d2CapLambdaDs2[j][k]/pow(CapLambda-1.0, 4.0)
                                   + 36.0*(CapLambda - 1.0/CapLambda - 2.0*log(CapLambda))
                                         *dCapLambdaDs[j]*dCapLambdaDs[k]/pow(CapLambda-1.0, 5.0);

                D2dhMuExcessDs2[j][k] += - 3.0*AsubG*D2sqrtOfIds2[j][k]*pow(sqrtOfI, 2.0)*sigma/3.0/CapGam
                                       - 6.0*AsubG*DsqrtOfIds[j]*DsqrtOfIds[k]*sqrtOfI*sigma/3.0/CapGam
                                       - 3.0*AsubG*DsqrtOfIds[j]*pow(sqrtOfI, 2.0)*DsigmaDsk/3.0/CapGam
                                       - 3.0*AsubG*DsqrtOfIds[k]*pow(sqrtOfI, 2.0)*DsigmaDsj/3.0/CapGam
                                       - AsubG*pow(sqrtOfI, 3.0)*D2sigmaDsjk/3.0/CapGam;
                D2dhMuExcessDs2[k][j] = D2dhMuExcessDs2[j][k];
            }
        }
    } else {
        id object = [endmembers objectAtIndex:index];
        if ([object isKindOfClass:[HKFspeciesProperties class]]) {
            double z = [(HKFspeciesProperties *)[endmembers objectAtIndex:index] charge];
            if (z != 0.0)
                for (NSUInteger j=0; j<NS; j++) {
                    for (NSUInteger k=j; k<NS; k++) {
                        D2dhMuExcessDs2[j][k] += z*z*AsubG*D2sqrtOfIds2[j][k]/2.0/CapLambda
                                               - z*z*AsubG*DsqrtOfIds[j]*dCapLambdaDs[k]/2.0/pow(CapLambda, 2.0)
                                               - z*z*AsubG*DsqrtOfIds[k]*dCapLambdaDs[j]/2.0/pow(CapLambda, 2.0)
                                               - z*z*AsubG*sqrtOfI*d2CapLambdaDs2[j][k]/2.0/pow(CapLambda, 2.0)
                                               + 2.0*z*z*AsubG*sqrtOfI*dCapLambdaDs[j]*dCapLambdaDs[k]/2.0/pow(CapLambda, 3.0);
                        D2dhMuExcessDs2[k][j] = D2dhMuExcessDs2[j][k];
                    }
                }
        }
    }
}

#pragma mark -
#pragma mark Thermodynamic functions

- (double)fillGwithT:(double)t andwithP:(double)p {
    double g = 0.0;
    initializeDHexcessTermsZerothOrder = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        // standard state (solvent and solute)
        g += nSpecies*xSpecies[i]*[[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
        // mixing properties
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        if (i == mIndH2O) {
            // solvent - osmotic term
            g += nSpecies*xSpecies[mIndH2O]*dhMuExcess;
        } else {
            // solute - ideal term
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            g += nSpecies*xSpecies[i]*R*t*log(molality);
            // solute - actvity coefficient
            g += nSpecies*xSpecies[i]*dhMuExcess;
        }
    }
    return g;
}

- (double)fillSwithT:(double)t andWithP:(double)p {
    double dgdt = 0.0;
    initializeDHexcessTermsFirstOrderDt = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        // standard state (solvent and solute)
        dgdt += -nSpecies*xSpecies[i]*[[endmembers objectAtIndex:i] getEntropyFromT:t andP:p];
        // mixing properties
        double DdhMuExcessDt = [self fillDebyeHuckelExcessFreeEnergyDtForSpeciesIndex:i withT:t andWithP:p];
        if (i == mIndH2O) {
            // solvent - osmotic term
            dgdt += nSpecies*xSpecies[mIndH2O]*DdhMuExcessDt;
        } else {
            // solute - ideal term
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            dgdt += nSpecies*xSpecies[i]*R*log(molality);
            // solute - actvity coefficient
            dgdt += nSpecies*xSpecies[i]*DdhMuExcessDt;
        }
    }
    return -dgdt;
}

- (double)fillVwithT:(double)t andWithP:(double)p {
    double dgdp = 0.0;
    initializeDHexcessTermsZerothOrderDp = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        // standard state (solvent and solute)
        dgdp += nSpecies*xSpecies[i]*[[endmembers objectAtIndex:i] getVolumeFromT:t andP:p];
        // mixing properties
        double DdhMuExcessDp = [self fillDebyeHuckelExcessFreeEnergyDpForSpeciesIndex:i withT:t andWithP:p];
        if (i == mIndH2O) {
            // solvent - osmotic term
            dgdp += nSpecies*xSpecies[mIndH2O]*DdhMuExcessDp;
        } else {
            // solute - actvity coefficient
            dgdp += nSpecies*xSpecies[i]*DdhMuExcessDp;
        }
    }
    return dgdp;
}

- (void)fillDGDRin:(double [NR])dgdr withT:(double)t andWithP:(double)p {
    memset(dgdr, 0, sizeof(dgdr[0])*NR);
    initializeDHexcessTermsZerothOrder = YES;
    initializeDHexcessTermsFirstOrderDr = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        double DdhMuExcessDr[NR];
        double gSS = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:DdhMuExcessDr withT:t andWithP:p];

        for (NSUInteger j=0; j<NR; j++) {
            dgdr[j] += nSpecies*dxSpeciesDr[i][j]*gSS;
            if (i == mIndH2O) {
                dgdr[j] += nSpecies*dxSpeciesDr[mIndH2O][j]*dhMuExcess
                         + nSpecies*xSpecies[mIndH2O]*DdhMuExcessDr[j];
            } else {
                double molality    = xSpecies[i]*CapGam/xSpecies[mIndH2O];
                double DmolalityDr = dxSpeciesDr[i][j]*CapGam/xSpecies[mIndH2O]
                                   - xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0);
                dgdr[j] += nSpecies*dxSpeciesDr[i][j]*R*t*log(molality)
                         + nSpecies*xSpecies[i]*R*t*DmolalityDr/molality;
                dgdr[j] += nSpecies*dxSpeciesDr[i][j]*dhMuExcess
                         + nSpecies*xSpecies[i]*DdhMuExcessDr[j];
            }
        }
    }
}

- (void)fillDGDSin:(double [NS])dgds withT:(double)t andWithP:(double)p {
    memset(dgds, 0, sizeof(dgds[0])*NS);
    initializeDHexcessTermsZerothOrder = YES;
    initializeDHexcessTermsFirstOrderDs = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        double DdhMuExcessDs[NS];
        double gSS = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:DdhMuExcessDs withT:t andWithP:p];

        for (NSUInteger j=0; j<NS; j++) {
            dgds[j] += (dnSpeciesds[j]*xSpecies[i] + nSpecies*dxSpecies[i][j])*gSS;
            if (i == mIndH2O) {
                dgds[j] += dnSpeciesds[j]*xSpecies[mIndH2O]*dhMuExcess
                         + nSpecies*dxSpecies[mIndH2O][j]*dhMuExcess
                         + nSpecies*xSpecies[mIndH2O]*DdhMuExcessDs[j];
            } else {
                double molality    = xSpecies[i]*CapGam/xSpecies[mIndH2O];
                double DmolalityDs = dxSpecies[i][j]*CapGam/xSpecies[mIndH2O]
                                   - xSpecies[i]*CapGam*dxSpecies[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0);
                dgds[j] += dnSpeciesds[j]*xSpecies[i]*R*t*log(molality)
                         + nSpecies*dxSpecies[i][j]*R*t*log(molality)
                         + nSpecies*xSpecies[i]*R*t*DmolalityDs/molality;
                dgds[j] += dnSpeciesds[j]*xSpecies[i]*dhMuExcess
                         + nSpecies*dxSpecies[i][j]*dhMuExcess
                         + nSpecies*xSpecies[i]*DdhMuExcessDs[j];
            }
        }
    }
}

- (void)fillD2GDR2in:(double [NR][NR])d2gdr2 withT:(double)t andWithP:(double)p {
    memset(d2gdr2, 0, sizeof(d2gdr2[0][0])*NR*NR);
    initializeDHexcessTermsZerothOrder    = YES;
    initializeDHexcessTermsFirstOrderDr   = YES;
    initializeDHexcessTermsSecondOrderDr2 = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        double DdhMuExcessDr[NR], D2dhMuExcessDr2[NR][NR];
        [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:DdhMuExcessDr withT:t andWithP:p];
        [self fillD2debyeHuckelExcessFreeEnergyDr2ForSpeciesIndex:i in:D2dhMuExcessDr2 withT:t andWithP:p];

        if (i == mIndH2O) {
            for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=j; k<NR; k++) {
                d2gdr2[j][k] += nSpecies*dxSpeciesDr[mIndH2O][j]*DdhMuExcessDr[k]
                              + nSpecies*dxSpeciesDr[mIndH2O][k]*DdhMuExcessDr[j]
                              + nSpecies*xSpecies[mIndH2O]*D2dhMuExcessDr2[j][k];
                d2gdr2[k][j] = d2gdr2[j][k];
            }
        } else {
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            for (NSUInteger j=0; j<NR; j++) {
                double DmolalityDrj = dxSpeciesDr[i][j]*CapGam/xSpecies[mIndH2O]
                                    - xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0);
                for (NSUInteger k=j; k<NR; k++) {
                    double DmolalityDrk = dxSpeciesDr[i][k]*CapGam/xSpecies[mIndH2O]
                                        - xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0);
                    double D2molalityDrjk = - dxSpeciesDr[i][j]*CapGam*dxSpeciesDr[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0)
                                            - dxSpeciesDr[i][k]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0)
                                          + 2.0*xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]*dxSpeciesDr[mIndH2O][k]/pow(xSpecies[mIndH2O], 3.0);

                    d2gdr2[j][k] += nSpecies*dxSpeciesDr[i][j]*DmolalityDrk*R*t/molality
                                  + nSpecies*dxSpeciesDr[i][k]*R*t*DmolalityDrj/molality
                                  - nSpecies*xSpecies[i]*R*t*DmolalityDrj*DmolalityDrk/pow(molality, 2.0)
                                  + nSpecies*xSpecies[i]*R*t*D2molalityDrjk/molality;

                    d2gdr2[j][k] += nSpecies*dxSpeciesDr[i][j]*DdhMuExcessDr[k]
                                  + nSpecies*dxSpeciesDr[i][k]*DdhMuExcessDr[j]
                                  + nSpecies*xSpecies[i]*D2dhMuExcessDr2[j][k];

                    d2gdr2[k][j] = d2gdr2[j][k];
                }
            }
        }
    }
}

- (void)fillD2GDRDSin:(double [NR][NS])d2gdrds withT:(double)t andWithP:(double)p {
    memset(d2gdrds, 0, sizeof(d2gdrds[0][0])*NR*NS);
    initializeDHexcessTermsZerothOrder     = YES;
    initializeDHexcessTermsFirstOrderDr    = YES;
    initializeDHexcessTermsFirstOrderDs    = YES;
    initializeDHexcessTermsSecondOrderDrDs = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        double DdhMuExcessDs[NS], DdhMuExcessDr[NR], D2dhMuExcessDrDs[NR][NS];
        double gSS = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:DdhMuExcessDr withT:t andWithP:p];
        [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:DdhMuExcessDs withT:t andWithP:p];
        [self fillD2debyeHuckelExcessFreeEnergyDrDsForSpeciesIndex:i in:D2dhMuExcessDrDs withT:t andWithP:p];

        for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=0; k<NS; k++) {
            d2gdrds[j][k] += dnSpeciesds[k]*dxSpeciesDr[i][j]*gSS
                           + nSpecies*[self D2xSpeciesDrDsForIndex:i andDr:j andDs:k]*gSS;
        }

        if (i == mIndH2O) {
            for (NSUInteger j=0; j<NR; j++) for (NSUInteger k=0; k<NS; k++) {
                d2gdrds[j][k] += dnSpeciesds[k]*dxSpeciesDr[mIndH2O][j]*dhMuExcess
                               + nSpecies*[self D2xSpeciesDrDsForIndex:mIndH2O andDr:j andDs:k]*dhMuExcess
                               + nSpecies*dxSpeciesDr[mIndH2O][j]*DdhMuExcessDs[k]
                               + dnSpeciesds[k]*xSpecies[mIndH2O]*DdhMuExcessDr[j]
                               + nSpecies*dxSpecies[mIndH2O][k]*DdhMuExcessDr[j]
                               + nSpecies*xSpecies[mIndH2O]*D2dhMuExcessDrDs[j][k];
            }
        } else {
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            for (NSUInteger j=0; j<NR; j++) {
                double DmolalityDr = dxSpeciesDr[i][j]*CapGam/xSpecies[mIndH2O]
                                   - xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0);

                for (NSUInteger k=0; k<NS; k++) {
                    double DmolalityDs = dxSpecies[i][k]*CapGam/xSpecies[mIndH2O]
                                       - xSpecies[i]*CapGam*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0);
                    double d2xSpeciesDrDs = [self D2xSpeciesDrDsForIndex:i andDr:j andDs:k];
                    double d2xH2ODrDs     = [self D2xSpeciesDrDsForIndex:mIndH2O andDr:j andDs:k];
                    double D2molalityDrDs = d2xSpeciesDrDs*CapGam/xSpecies[mIndH2O]
                                            - dxSpeciesDr[i][j]*CapGam*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0)
                                            - dxSpecies[i][k]*CapGam*dxSpeciesDr[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0)
                                            - xSpecies[i]*CapGam*d2xH2ODrDs/pow(xSpecies[mIndH2O], 2.0)
                                          + 2.0*xSpecies[i]*CapGam*dxSpeciesDr[mIndH2O][j]*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 3.0);

                    d2gdrds[j][k] += dnSpeciesds[k]*dxSpeciesDr[i][j]*R*t*log(molality)
                                   + nSpecies*d2xSpeciesDrDs*R*t*log(molality)
                                   + nSpecies*dxSpeciesDr[i][j]*R*t*DmolalityDs/molality
                                   + dnSpeciesds[k]*xSpecies[i]*R*t*DmolalityDr/molality
                                   + nSpecies*dxSpecies[i][k]*R*t*DmolalityDr/molality
                                   - nSpecies*xSpecies[i]*R*t*DmolalityDr*DmolalityDs/pow(molality, 2.0)
                                   + nSpecies*xSpecies[i]*R*t*D2molalityDrDs/molality;
                    d2gdrds[j][k] += dnSpeciesds[k]*dxSpeciesDr[i][j]*dhMuExcess
                                   + nSpecies*d2xSpeciesDrDs*dhMuExcess
                                   + nSpecies*dxSpeciesDr[i][j]*DdhMuExcessDs[k]
                                   + dnSpeciesds[k]*xSpecies[i]*DdhMuExcessDr[j]
                                   + nSpecies*dxSpecies[i][k]*DdhMuExcessDr[j]
                                   + nSpecies*xSpecies[i]*D2dhMuExcessDrDs[j][k];
                }
            }
        }
    }

}

- (void)fillD2GDRDTin:(double [NR])d2gdrdt withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        d2gdrdt[k] = nSpecies*R*(log(xSpecies[k+1]) - log(xSpecies[0]));
    }
    NSLog(@"Function not yet implemented.");
}

- (void)fillD2GDRDPin:(double [NR])d2gdrdp withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        d2gdrdp[k] = 0.0;
    }
    NSLog(@"Function not yet implemented.");
}

- (void)fillD2GDS2in:(double [NS][NS])d2gds2 withT:(double)t andWithP:(double)p {
    memset(d2gds2, 0, sizeof(d2gds2[0][0])*NS*NS);
    initializeDHexcessTermsZerothOrder = YES;
    initializeDHexcessTermsFirstOrderDs  = YES;
    initializeDHexcessTermsSecondOrderDs2 = YES;

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        double DdhMuExcessDs[NS], D2dhMuExcessDs2[NS][NS];
        double gSS = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:DdhMuExcessDs withT:t andWithP:p];
        [self fillD2debyeHuckelExcessFreeEnergyDs2ForSpeciesIndex:i in:D2dhMuExcessDs2 withT:t andWithP:p];

        for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
            d2gds2[j][k] += (d2nSpeciesds2[j][k]*xSpecies[i] + dnSpeciesds[j]*dxSpecies[i][k]
                             + dnSpeciesds[k]*dxSpecies[i][j])*gSS;
            d2gds2[k][j] = d2gds2[j][k];
        }

        if (i == mIndH2O) {
            for (NSUInteger j=0; j<NS; j++) for (NSUInteger k=j; k<NS; k++) {
                d2gds2[j][k] += d2nSpeciesds2[j][k]*xSpecies[mIndH2O]*dhMuExcess
                              + dnSpeciesds[j]*dxSpecies[mIndH2O][k]*dhMuExcess
                              + dnSpeciesds[k]*dxSpecies[mIndH2O][j]*dhMuExcess
                              + dnSpeciesds[j]*xSpecies[mIndH2O]*DdhMuExcessDs[k]
                              + dnSpeciesds[k]*xSpecies[mIndH2O]*DdhMuExcessDs[j]
                              + nSpecies*dxSpecies[mIndH2O][j]*DdhMuExcessDs[k]
                              + nSpecies*dxSpecies[mIndH2O][k]*DdhMuExcessDs[j]
                              + nSpecies*xSpecies[mIndH2O]*D2dhMuExcessDs2[j][k];
                d2gds2[k][j] = d2gds2[j][k];
            }
        } else {
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            for (NSUInteger j=0; j<NS; j++) {
                double DmolalityDsj = dxSpecies[i][j]*CapGam/xSpecies[mIndH2O]
                                    - xSpecies[i]*CapGam*dxSpecies[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0);
                for (NSUInteger k=j; k<NS; k++) {
                    double DmolalityDsk = dxSpecies[i][k]*CapGam/xSpecies[mIndH2O]
                                        - xSpecies[i]*CapGam*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0);
                    double D2molalityDsjk = - dxSpecies[i][j]*CapGam*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 2.0)
                                            - dxSpecies[i][k]*CapGam*dxSpecies[mIndH2O][j]/pow(xSpecies[mIndH2O], 2.0)
                                          + 2.0*xSpecies[i]*CapGam*dxSpecies[mIndH2O][j]*dxSpecies[mIndH2O][k]/pow(xSpecies[mIndH2O], 3.0);

                    d2gds2[j][k] += d2nSpeciesds2[j][k]*xSpecies[i]*R*t*log(molality)
                                  + dnSpeciesds[j]*dxSpecies[i][k]*R*t*log(molality)
                                  + dnSpeciesds[k]*dxSpecies[i][j]*R*t*log(molality)
                                  + dnSpeciesds[j]*xSpecies[i]*DmolalityDsk*R*t/molality
                                  + dnSpeciesds[k]*xSpecies[i]*R*t*DmolalityDsj/molality
                                  + nSpecies*dxSpecies[i][j]*DmolalityDsk*R*t/molality
                                  + nSpecies*dxSpecies[i][k]*R*t*DmolalityDsj/molality
                                  - nSpecies*xSpecies[i]*R*t*DmolalityDsj*DmolalityDsk/pow(molality, 2.0)
                                  + nSpecies*xSpecies[i]*R*t*D2molalityDsjk/molality;

                    d2gds2[j][k] += d2nSpeciesds2[j][k]*xSpecies[i]*dhMuExcess
                                  + dnSpeciesds[j]*dxSpecies[i][k]*dhMuExcess
                                  + dnSpeciesds[k]*dxSpecies[i][j]*dhMuExcess
                                  + dnSpeciesds[j]*xSpecies[i]*DdhMuExcessDs[k]
                                  + dnSpeciesds[k]*xSpecies[i]*DdhMuExcessDs[j]
                                  + nSpecies*dxSpecies[i][j]*DdhMuExcessDs[k]
                                  + nSpecies*dxSpecies[i][k]*DdhMuExcessDs[j]
                                  + nSpecies*xSpecies[i]*D2dhMuExcessDs2[j][k];

                    d2gds2[k][j] = d2gds2[j][k];
                }
            }
        }
    }
}

- (void)fillD2GDSDTin:(double [NS])d2gdsdt withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) {
        d2gdsdt[i] = 0.0;
        for (NSUInteger j=0; j<NS; j++)
            d2gdsdt[i] += nSpecies*(/* - ss entropy i */ 0.0)*dxSpecies[NA+j][i]
                        + dnSpeciesds[i]*(/* - ss entropy j */ 0.0)*xSpecies[NA+j];
        for (NSUInteger j=0; j<NE; j++)
            d2gdsdt[i] += nSpecies*R*(1.0 + log(xSpecies[j]))*dxSpecies[j][i] + dnSpeciesds[i]*R*xSpecies[j]*log(xSpecies[j]);
    }
    NSLog(@"fillD2GDSDT not yet implemented.");
}

- (void)fillD2GDSDPin:(double [NS])d2gdsdp withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) {
        d2gdsdp[i] = 0.0;
        for (NSUInteger j=0; j<NS; j++)
            d2gdsdp[i] += nSpecies*(/* ss v of j */ 0.0)*dxSpecies[NA+j][i]
                        + dnSpeciesds[i]*(/* ss v of j */ 0.0)*xSpecies[NA+j];
    }
    NSLog(@"fillD2GDSDP not yet implemented.");
}

- (double)fillD2GDT2withT:(double)t andWithP:(double)p {
    NSLog(@"fillD2GDT2 not yet implemented.");
    return 0.0;
}

- (double)fillD2GDTDPwithT:(double)t andWithP:(double)p {
    NSLog(@"fillD2GDTDP not yet implemented.");
    return 0.0;
}

- (double)fillD2GDP2withT:(double)t andWithP:(double)p {
    NSLog(@"fillD2GDP2 not yet implemented.");
    return 0.0;
}

- (void)fillD3GDR2DSin:(double [NR][NR][NS])d3gdr2ds withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        for (NSUInteger l=0; l<NR; l++) {
            for (NSUInteger i=0; i<NS; i++) {
                d3gdr2ds[k][l][i] = -R*t*dxSpecies[0][i]/xSpecies[0]/xSpecies[0];
                if (k == l) d3gdr2ds[k][k][i] += -R*t*dxSpecies[k+1][i]/xSpecies[k+1]/xSpecies[k+1];
            }
        }
    }
    NSLog(@"fillD3GDR2DS not yet implemented.");
}

- (void)fillD3GDR2DTin:(double [NR][NR])d3gdr2dt withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        for (NSUInteger l=0; l<NR; l++) {
            d3gdr2dt[k][l] = R/xSpecies[0];
            if (k == l) d3gdr2dt[k][k] += R/xSpecies[k+1];
        }
    }
    NSLog(@"fillD3GDR2DT not yet implemented.");
}

- (void)fillD3GDR2DPin:(double [NR][NR])d3gdr2dp withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        for (NSUInteger l=0; l<NR; l++) {
            d3gdr2dp[k][l] = 0.0;
        }
    }
    NSLog(@"fillD3GDR2DP not yet implemented.");
}

- (void)fillD3GDRDS2in:(double [NR][NS][NS])d3gdrds2 withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        for (NSUInteger l=0; l<NS; l++) {
            for (NSUInteger m=0; m<NS; m++)
                d3gdrds2[k][l][m] = R*t*(-dxSpecies[k+1][l]*dxSpecies[k+1][m]/xSpecies[k+1]/xSpecies[k+1] + dxSpecies[0][l]*dxSpecies[0][m]/xSpecies[0]/xSpecies[0]);
        }
    }
    NSLog(@"fillD3GDRDS2 not yet implemented.");
}

- (void)fillD3GDRDSDTin:(double [NR][NS])d3gdrdsdt withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        for (NSUInteger l=0; l<NS; l++) {
            d3gdrdsdt[k][l] = R*(dxSpecies[k+1][l]/xSpecies[k+1] - dxSpecies[0][l]/xSpecies[0]);
        }
    }
    NSLog(@"fillD3GDRDSDT not yet implemented.");
}

- (void)fillD3GDRDSDPin:(double [NR][NS])d3gdrdsdp withT:(double)t andWithP:(double)p {
    for (NSUInteger k=0; k<NR; k++) {
        for (NSUInteger l=0; l<NS; l++) {
            d3gdrdsdp[k][l] = 0.0;
        }
    }
    NSLog(@"fillD3GDRDSDP not yet implemented.");
}

- (void)fillD3GDRDT2in:(double [NR])d3gdrdt2 withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NR; i++) d3gdrdt2[i] = 0.0;
    NSLog(@"fillD3GDRDT2 not yet implemented.");
}

- (void)fillD3GDRDTDPin:(double [NR])d3gdrdtdp withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NR; i++) d3gdrdtdp[i] = 0.0;
    NSLog(@"fillD3GDRDTDP not yet implemented.");
}

- (void)fillD3GDRDP2in:(double [NR])d3gdrdp2 withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NR; i++) d3gdrdp2[i] = 0.0;
    NSLog(@"fillD3GDRDP2 not yet implemented.");
}

// program this way to avoid declaring an array of dimension [NS][NS][NS]
- (double)fillD3GDS3atI:(NSUInteger)i andL:(NSUInteger)l andM:(NSUInteger)m withT:(double)t andWithP:(double)p {
//    - (double)fillD3nSpeciesDs3atI:(NSUInteger)i andL:(NSUInteger)l andM:(NSUInteger)m withS:(double [NS])s {
//        double coeff = [self speciationCoefficientWith:s];
//        if (coeff == 0.0) return 0.0;
//        else {
//            [self DspeciationCoefficientDsWith:s];
//            return -6.0*DcoeffDs[i]*DcoeffDs[l]*DcoeffDs[m]/coeff/coeff/coeff/coeff;
//        }
//    }
    double result = 0;
    for (NSUInteger j=0; j<NE; j++)
        result += -R*t*dxSpecies[j][l]*dxSpecies[j][i]*dxSpecies[j][m]/xSpecies[j]/xSpecies[j];
    if ((i == 0) && (l == 0) && (m == 0)) NSLog(@"fillD3GDS3 not yet implemented.");
    return result;
}

- (void)fillD3GDS2DTin:(double [NS][NS])d3gds2dt withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) {
        for (NSUInteger l=0; l<NS; l++) d3gds2dt[i][l] = 0.0;
        for (NSUInteger j=0; j<NE; j++) {
            for (NSUInteger l=0; l<NS; l++) d3gds2dt[i][l] += R*dxSpecies[j][l]*dxSpecies[j][i]/xSpecies[j];
        }
    }
    NSLog(@"fillD3GDS2DT not yet implemented.");
}

- (void)fillD3GDS2DPin:(double [NS][NS])d3gds2dp withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) {
        for (NSUInteger l=0; l<NS; l++) d3gds2dp[i][l] = 0.0;
        for (NSUInteger j=0; j<NE; j++) {
        }
    }
    NSLog(@"fillD3GDS2DP not yet implemented.");
}

- (void)fillD3GDSDT2in:(double [NS])d3gdsdt2 withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) d3gdsdt2[i] = 0.0;
    NSLog(@"fillD3GDSDT2 not yet implemented.");
}

- (void)fillD3GDSDTDPin:(double [NS])d3gdsdtdp withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) d3gdsdtdp[i] = 0.0;
    NSLog(@"fillD3GDSDTDP not yet implemented.");
}

- (void)fillD3GDSDP2in:(double [NS])d3gdsdp2 withT:(double)t andWithP:(double)p {
    for (NSUInteger i=0; i<NS; i++) d3gdsdp2[i] = 0.0;
    NSLog(@"fillD3GDSDP2 not yet implemented.");
}

- (double)fillD3GDT3withT:(double)t andWithP:(double)p {
    NSLog(@"fillD3GDT3 not yet implemented.");
    return 0.0;
}

- (double)fillD3GDT2DPwithT:(double)t andWithP:(double)p {
    NSLog(@"fillD3GDT2DP not yet implemented.");
    return 0.0;
}

- (double)fillD3GDTDP2withT:(double)t andWithP:(double)p {
    NSLog(@"fillD3GDTDP2 not yet implemented.");
    return 0.0;
}

- (double)fillD3GDP3withT:(double)t andWithP:(double)p {
    NSLog(@"fillD3GDP3 not yet implemented.");
    return 0.0;
}

#pragma mark -
#pragma mark Ordering (Speciation) routiunes

/* -------------------------------
 From SIMP1.C (Numerical Recipies)
 -------------------------------*/

- (void)simp1:(double [NA+4][NS+2])a mm:(NSInteger)mm ll:(NSInteger [NS+2])ll nll:(NSInteger)nll
         iabf:(NSInteger)iabf kp:(NSInteger *)kp bmax:(double *)bmax {
    NSInteger k;
    double test;

    *kp = ll[1];
    *bmax = a[mm+1][*kp+1];
    for (k=2; k<=nll; k++) {
        if (iabf == 0)
            test = a[mm+1][ll[k]+1] - (*bmax);
        else
            test = fabs(a[mm+1][ll[k]+1]) - fabs(*bmax);
        if (test > 0.0) {
            *bmax = a[mm+1][ll[k]+1];
            *kp=ll[k];
        }
    }
}

/* -------------------------------
 From SIMP2.C (Numerical Recipies)
  ------------------------------*/

#define EPS 1.0e-6
- (void)simp2:(double [NA+4][NS+2])a n:(NSInteger)n l2:(NSInteger [NA+4])l2 nl2:(NSInteger)nl2
           ip:(NSInteger *)ip kp:(NSInteger)kp q1:(double *)q1 {
    NSInteger k, ii, i;
    double qp = 0.0, q0 = 0.0, q;

    *ip = 0;
    for (i=1; i<=nl2; i++) {
        if (a[l2[i]+1][kp+1] < -EPS) {
            *q1 = -a[l2[i]+1][1]/a[l2[i]+1][kp+1];
            *ip = l2[i];
            for (i=i+1; i<=nl2; i++) {
                ii = l2[i];
                if (a[ii+1][kp+1] < -EPS) {
                    q = -a[ii+1][1]/a[ii+1][kp+1];
                    if (q < *q1) {
                        *ip = ii;
                        *q1 = q;
                    } else if (q == *q1) {
                        for (k=1; k<=n; k++) {
                            qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
                            q0 = -a[ii+1][k+1]/a[ii+1][kp+1];
                            if (q0 != qp) break;
                        }
                        if (q0 < qp) *ip=ii;
                    }
                }
            }
        }
    }
}
#undef EPS

/* -------------------------------
 From SIMP3.C (Numerical Recipies)
 -------------------------------*/

- (void)simp3:(double [NA+4][NS+2])a i1:(NSInteger)i1 k1:(NSInteger)k1 ip:(NSInteger)ip kp:(NSInteger)kp {
    NSInteger kk, ii;
    double piv;

    piv = 1.0/a[ip+1][kp+1];
    for (ii=1; ii<=i1+1; ii++)
        if (ii-1 != ip) {
            a[ii][kp+1] *= piv;
            for (kk=1; kk<=k1+1; kk++)
                if (kk-1 != kp)
                    a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
        }
    for (kk=1; kk<=k1+1; kk++)
        if (kk-1 != kp) a[ip+1][kk] *= -piv;
    a[ip+1][kp+1] = piv;
}

/* --------------------------------------------------
 From SIMPLX.C (Numerical Recipies)
 removed internal declarations of simp1, simp2, simp3
 --------------------------------------------------*/

#define EPS 1.0e-6
// m = NA+1, n = NS for this application
- (void)simplx:(double [NA+4][NS+2])a m:(NSInteger)m n:(NSInteger)n m1:(NSInteger)m1 m2:(NSInteger)m2
            m3:(NSInteger)m3 icasePt:(NSInteger *)icase izrov:(NSInteger [NS+2])izrov iposv:(NSInteger [NA+4])iposv
{
    NSInteger i, ip, ir, is, k, kh, kp, m12, nl1, nl2;
    NSInteger l1[NS+2],  l2[NA+4], l3[NA+4];
    double q1, bmax;

    NSAssert(m == m1+m2+m3, @"Bad input constraint counts in simplx");
    NSAssert(n >= 1, @"Bad input constraint counts in simplx");
    NSAssert(m >= 1, @"Bad input constraint counts in simplx");

    nl1 = n;
    for (k=1; k<=n; k++) l1[k] = izrov[k]=k;
    nl2 = m;
    for (i=1;i<=m;i++) {
        NSAssert(a[i+1][1] >= 0.0, @"Bad input tableau in simplx");
        l2[i] = i;
        iposv[i] = n+i;
    }
    for (i=1; i<=m2; i++) l3[i] = 1;
    if (m2 + m3) {
        ir = 1;
        for (k=1; k<=(n+1); k++) {
            q1 = 0.0;
            for (i=m1+1; i<=m; i++) q1 += a[i+1][k];
            a[m+2][k] = -q1;
        }
        do {
            [self simp1:a mm:m+1 ll:l1 nll:nl1 iabf:0 kp:&kp bmax:&bmax];
            if (bmax <= EPS && a[m+2][1] < -EPS) {
                *icase = -1;
                return;
            } else if (bmax <= EPS && a[m+2][1] <= EPS) {
                m12 = m1 + m2 + 1;
                if (m12 <= m) {
                    for (ip=m12; ip<=m; ip++) {
                        if (iposv[ip] == (ip+n)) {
                            [self simp1:a mm:ip ll:l1 nll:nl1 iabf:1 kp:&kp bmax:&bmax];
                            if (bmax > 0.0)
                                goto one;
                        }
                    }
                }
                --m12;
                if (m1+1 <= m12)
                    for (i=m1+1; i<=m12; i++)
                        if (l3[i-m1] == 1)
                            for (k=1; k<=n+1; k++)
                                a[i+1][k] = -a[i+1][k];
                break;
            }
            [self simp2:a n:n l2:l2 nl2:nl2 ip:&ip kp:kp q1:&q1];
            if (ip == 0) {
                *icase = -1;
                return;
            }
        one:
            [self simp3:a i1:m+1 k1:n ip:ip kp:kp];
            if (iposv[ip] >= (n+m1+m2+1)) {
                for (k=1; k<=nl1; k++)
                    if (l1[k] == kp) break;
                --nl1;
                for (is=k; is<=nl1; is++) l1[is] = l1[is+1];
                ++a[m+2][kp+1];
                for (i=1; i<=m+2; i++) a[i][kp+1] = -a[i][kp+1];
            } else {
                if (iposv[ip] >= (n+m1+1)) {
                    kh = iposv[ip] - m1 - n;
                    if (l3[kh]) {
                        l3[kh]=0;
                        ++a[m+2][kp+1];
                        for (i=1; i<=m+2; i++)
                            a[i][kp+1] = -a[i][kp+1];
                    }
                }
            }
            is = izrov[kp];
            izrov[kp] = iposv[ip];
            iposv[ip] = is;
        } while (ir);
    }
    for (;;) {
        [self simp1:a mm:0 ll:l1 nll:nl1 iabf:0 kp:&kp bmax:&bmax];
        if (bmax <= 0.0) {
            *icase = 0;
            return;
        }
        [self simp2:a n:n l2:l2 nl2:nl2 ip:&ip kp:kp q1:&q1];
        if (ip == 0) {
            *icase = 1;
            return;
        }
        [self simp3:a i1:m k1:n ip:ip kp:kp];
        is = izrov[kp];
        izrov[kp] = iposv[ip];
        iposv[ip] = is;
    }
}
#undef EPS

//    X1   CO2     r0
//    X2   O2      r1
//    X3   HF      r2
//    X4   NaOH    r3
//    X5   Mg(OH)2 r4
//    X6   HAlO2   r5
//    X7   SiO2    r6 - ignored, assumed neutral
//    X8   H3PO4   r7
//    X9   SO2     r8
//    X10  HCl     r9
//    X11  KOH     r10
//    X12  Ca(OH)2 r11
//    X13  H2CrO4  r12
//    X14  Mn(OH)2 r13
//    X15  Fe(OH)2 r14
//    X16  Co(OH)2 r15

- (void)guessOrderingParameters:(double [NS])s from:(double [NR])r atT:(double)t andP:(double)p {
    double rTotal = 0.0;
    for (NSUInteger i=0; i<NR; i++) rTotal += r[i];
    double XH2O = 1.0 - rTotal;
    double K;
    BOOL debug = self.debugV;

    double gH2O = [[endmembers objectAtIndex: 0] getGibbsFreeEnergyFromT:t andP:p];
    double gH   = [[endmembers objectAtIndex:17] getGibbsFreeEnergyFromT:t andP:p];
    double gOH  = [[endmembers objectAtIndex:18] getGibbsFreeEnergyFromT:t andP:p];
    double gO2  = [[endmembers objectAtIndex: 2] getGibbsFreeEnergyFromT:t andP:p];

                         // X0   H2O
    K = exp(-(gH+gOH-gH2O)/R/t);
    s[0] = sqrt(K*XH2O); // XH2O*1.0e-6;  // X17  H+  s0 (H2O = H+ + OH-)
    s[1] = sqrt(K*XH2O); // XH2O*1.0e-6;  // X18  OH- s1 (H2O = H+ + OH-)
    double XH  = s[0];
    double XOH = s[1];

    if (r[1] != 0.0) {   // X2 O2
    	double gH2 = [[endmembers objectAtIndex:19] getGibbsFreeEnergyFromT:t andP:p];
    	K = exp(-(gH2+gO2/2.0-gH2O)/R/t);
        s[2] = K*XH2O/sqrt(r[1]);  // X19  H2 s2 (H2O = H2 + O2/2)
        double total = s[2];
        total = (total >= XH2O*0.10) ? (1.0-100.0*DBL_EPSILON)*XH2O*0.10/total : 1.0;
        s[2] *= total;
        if (debug) NSLog(@"H2O %g H+ %g OH- %g O2 %g H2 %g", XH2O, XH, XOH, r[1], s[2]);
    }
    if (r[0] != 0.0) {   // X1   CO2      r0
    	double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
        double gCO3 = [[endmembers objectAtIndex:20] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCO3+2.0*gH-gH2O-gCO2)/R/t);
        s[3] = K*r[0]*XH2O/XH/XH; // X20  CO3=  s3 (CO2 + H2O = CO3= + 2H+)
        double gHCO3 = [[endmembers objectAtIndex:21] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHCO3+gH-gH2O-gCO2)/R/t);
        s[4] = K*r[0]/XH; // X21  HCO3- s4 (CO2 + H2O = HCO3- + H+)
        if (s[2] > 0.0) {
        	double gCO = [[endmembers objectAtIndex:22] getGibbsFreeEnergyFromT:t andP:p];
        	K = exp(-(gCO+gO2/2.0-gCO2)/R/t);
        	s[5] = K*r[0]/sqrt(s[2]); // X22  CO    s5 (CO2 = CO + O2/2)
        }
        double total = s[3] + s[4] + s[5];
        total = (total >= r[0]) ? (1.0-100.0*DBL_EPSILON)*r[0]/total : 1.0;
        s[3] *= total;
        s[4] *= total;
        s[5] *= total;
        if (debug) NSLog(@"CO2 %g CO3-2 %g HCO3- %g CO %g", r[0], s[3], s[4], s[5]);
    }
    if (r[2] > 0.0) {    // X3   HF r2
    	double gHF = [[endmembers objectAtIndex: 3] getGibbsFreeEnergyFromT:t andP:p];
    	double gF  = [[endmembers objectAtIndex:23] getGibbsFreeEnergyFromT:t andP:p];
    	K = exp(-(gF+gH-gHF)/R/t);
        s[6] = K*r[2]/XH; // X23  F- s6 (HF = H+ + F-)
        double total = s[6];
        total = (total >= r[2]) ? (1.0-100.0*DBL_EPSILON)*r[2]/total : 1.0;
        s[6] *= total;
        if (debug) NSLog(@"HF %g F- %g", r[2], s[6]);
    }
    if (r[3] > 0.0) {    // X4   NaOH     r3
        double gNaOH = [[endmembers objectAtIndex: 4] getGibbsFreeEnergyFromT:t andP:p];
        double gNa   = [[endmembers objectAtIndex:25] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gNa+gOH-gNaOH)/R/t);
        s[8] = K*r[3]/XOH; // X25  Na+ s8 (NaOH = Na+ + OH-)
        if (r[0] != 0.0) {
            double gCO2   = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
            double gNaCO3 = [[endmembers objectAtIndex:26] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gNaCO3+gH-gNaOH-gCO2)/R/t);
            s[ 9] = K*r[0]*r[3]/XH; // X26  NaCO3- s9  (NaOH + CO2 = NaCO3- + H+)
            double gNaHCO3 = [[endmembers objectAtIndex:27] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gNaHCO3-gNaOH-gCO2)/R/t);
            s[10] = K*r[0]*r[3]; // X27  NaHCO3 s10 (NaOH + CO2 = NaHCO3)
        }
        if (r[9] != 0.0) {
            double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
            double gNaCl = [[endmembers objectAtIndex:24] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gNaCl+gH2O-gNaOH-gHCl)/R/t);
            s[ 7] = K*r[9]*r[3]/XH2O; // X24  NaCl s7 (NaOH + HCl = NaCl + H2O)
        }
        if (r[6] != 0.0) {
            double gSiO2  = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
            double gHSiO3 = [[endmembers objectAtIndex:28] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gHSiO3-gNaOH-gSiO2)/R/t);
            s[11] = K*r[3]*r[6]; // X28  NaHSiO3 s11 (NaOH + SiO2 = NaHSiO3)
        }
        double total = s[8] + s[9] + s[10] + s[7] + s[11];
        total = (total >= r[3]) ? (1.0-100.0*DBL_EPSILON)*r[3]/total : 1.0;
        s[ 7] *= total;
        s[ 8] *= total;
        s[ 9] *= total;
        s[10] *= total;
        s[11] *= total;
        if (debug) NSLog(@"NaOH %g Na+ %g NaHSiO3 %g NaCl %g NaCO3- %g NaHCO3 %g", r[3], s[8], s[11], s[7], s[9], s[10]);
    }
    if (r[4] > 0.0) { // X5   Mg(OH)2 r4
    	double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
        double gMg     = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gMg+2.0*gOH-gMgO2H2)/R/t);
        s[15] = K*r[4]/XOH/XOH; // X32  Mg+2     s15 Mg(OH)2 = Mg++ + 2 OH-
        double gMgOH = [[endmembers objectAtIndex:17] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gMgOH+gOH-gMgO2H2)/R/t);
        s[17] = K*r[4]/XOH; // X34  MgOH+    s17 Mg(OH)2 = MgOH+ + OH-
        if (r[0] != 0.0) {
        	double gCO2   = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
            double gMgCO3 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMgCO3+gH2O-gMgO2H2-gCO2)/R/t);
            s[12] = K*r[4]*r[0]/XH2O; // X29  MgCO3    s12 Mg(OH)2 + CO2 = MgCO3 + H2O
            double gMgHCO3 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMgHCO3+gOH-gMgO2H2-gCO2)/R/t);
            s[14] = K*r[4]*r[0]/XOH; // X31  MgHCO3+  s14 Mg(OH)2 + CO2 = MgHCO3+ + OH-
        }
        if (r[9] != 0.0) {
            double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
            double gMgCl = [[endmembers objectAtIndex:33] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMgCl+2.0*gH2O-gMgO2H2-gHCl)/R/t);
            s[16] = K*r[4]*r[9]*XH/XH2O/XH2O; // X33  MgCl+    s16 Mg(OH)2 + HCl + H+ = MgCl+ + 2 H2O
        }
        if (r[6] != 0.0) {
            double gSiO2    = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
            double gMgHSiO3 = [[endmembers objectAtIndex:30] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMgHSiO3+gH2O-gMgO2H2-gSiO2)/R/t);
            s[13] = K*r[4]*r[6]*XH/XH2O; // X30  MgHSiO3+ s13 Mg(OH)2 + SiO2 + H+ = MgHSiO3+ + H2O
        }
        if (r[8] != 0.0) {
            double gSO2   = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
            double gMgSO4 = [[endmembers objectAtIndex:35] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMgSO4+gH2O-gMgO2H2-gSO2-gO2/2.0)/R/t);
            s[18] = K*r[4]*r[8]*sqrt(r[1])/XH2O; // X35  MgSO4    s18 Mg(OH)2 + SO2 + 1/2 O2 = MgSO4 + H2O
        }
        double total = s[12] + s[13] + s[14] + s[15] + s[16] + s[17] + s[18];
        total = (total >= r[4]) ? (1.0-100.0*DBL_EPSILON)*r[4]/total : 1.0;
        s[12] *= total;
        s[13] *= total;
        s[14] *= total;
        s[15] *= total;
        s[16] *= total;
        s[17] *= total;
        s[18] *= total;
        if (debug) NSLog(@"Mg(OH)2 %g MgHSiO3+ %g MgCO3 %g MgHCO3+ %g Mg++ %g MgCl+ %g MgOH+ %g MgSO4 %g",
            r[4], s[12], s[13], s[14], s[15], s[16], s[17], s[18]);
    }
    if (r[5] > 0.0) { // X6   HAlO2   r5
        double gHAlO2 = [[endmembers objectAtIndex: 6] getGibbsFreeEnergyFromT:t andP:p];
        double gAl    = [[endmembers objectAtIndex:36] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gAl+2.0*gH2O-gHAlO2-3.0*gH)/R/t);
        s[19] = K*r[5]*XH*XH*XH/XH2O/XH2O; // r[5]*0.1; // X36  Al+3  s19 (HAlO2 + 3H+ = Al+3 + 2H2O)
        double gAlO2  = [[endmembers objectAtIndex:37] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gAlO2+gH-gHAlO2)/R/t);
        s[20] = K*r[5]/XH; // r[5]*0.1; // X37  AlO2- s20 (HAlO2 = AlO2- + H+)
        double total = s[19] + s[20];
        total = (total >= r[5]) ? (1.0-100.0*DBL_EPSILON)*r[5]/total : 1.0;
        s[19] *= total;
        s[20] *= total;
        if (debug) NSLog(@"HAlO2 %g Al+3 %g AlO2- %g", r[4], s[19], s[20]);
    }
    if (r[6] > 0.0) { // X7   SiO2    r6
        double gSiO2  = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
        double gHSiO3 = [[endmembers objectAtIndex:38] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHSiO3+gH-gSiO2-gH2O)/R/t);
        s[21] = K*r[6]*XH2O/XH; // r[6]*0.2; // X38  HSiO3- s21 (SiO2 + H2O = HSiO3- + H+)
        double gSi2O4 = [[endmembers objectAtIndex:39] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gSi2O4-2.0*gSiO2)/R/t);
        s[22] = K*r[6]*r[6]; // r[6]*0.2; // X39  Si2O4  s22 (2SiO2 = Si2O4)
        double total = s[21] + 2.0*s[22];
        total = (total >= r[6]) ? (1.0-100.0*DBL_EPSILON)*r[6]/total : 1.0;
        s[21]  *= total;
        s[22] *= total;
        if (debug) NSLog(@"SiO2 %g HSiO3- %g Si2O4 %g", r[6], s[21], s[22]);
    }
    if (r[7] > 0.0) { // X8   H3PO4   r7
        double gH3PO4 = [[endmembers objectAtIndex: 8] getGibbsFreeEnergyFromT:t andP:p];
        double gH2PO4 = [[endmembers objectAtIndex:40] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gH2PO4+gH-gH3PO4)/R/t);
        s[23] = K*r[7]/XH; // X40  H2PO4-   s23 (H3PO4 = H2PO4- + H+)
        double gHPO4 = [[endmembers objectAtIndex:41] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHPO4+2.0*gH-gH3PO4)/R/t);
        s[24] = K*r[7]/XH/XH; // X41  HPO4-2   s24 (H3PO4 = HPO4-2 + 2H+)
        double gPO4 = [[endmembers objectAtIndex:42] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gPO4+3.0*gH-gH3PO4)/R/t);
        s[25] = K*r[7]/XH/XH/XH; // X42  PO4-3    s25 (H3PO4 = PO4-3 + 3H+)
        double gH3P2O7 = [[endmembers objectAtIndex:43] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gH3P2O7+gH+gH2O-2.0*gH3PO4)/R/t);
        s[26] = K*r[7]*r[7]/XH/XH2O; // X43  H3P2O7-  s26 (2H3PO4 = H3P2O7- + H+ + H2O)
        double gH2P2O7 = [[endmembers objectAtIndex:44] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gH2P2O7+2.0*gH+gH2O-2.0*gH3PO4)/R/t);
        s[27] = K*r[7]*r[7]/XH/XH/XH2O; // X44  H2P2O7-2 s27 (2H3PO4 = H2P2O7-2 + 2H+ + H2O)
        double total = s[23] + s[24] + s[25] + s[26] + s[27];
        total = (total >= r[7]) ? (1.0-100.0*DBL_EPSILON)*r[7]/total : 1.0;
        s[23] *= total;
        s[24] *= total;
        s[25] *= total;
        s[26] *= total;
        s[27] *= total;
        if (debug) NSLog(@"H3PO4 %g H2PO4- %g HPO4-2 %g PO4-3 %g H3P2O7- %g H2P2O7-2 %g",
            r[7], s[23], s[24], s[25], s[26], s[27]);
    }
    if (r[8] > 0.0) { // X9   SO2     r8
        double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
        double gH2S = [[endmembers objectAtIndex:45] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gH2S+3.0*gO2/2.0-gSO2-gH2O)/R/t);
        s[28] = K*r[8]*XH2O/pow(r[1], 3.0/2.0); // X45  H2S s28 (SO2 + H2O = H2S + 3/2 O2)
        double gHS = [[endmembers objectAtIndex:46] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHS+gH+3.0*gO2/2.0-gSO2-gH2O)/R/t);
        s[29] = K*r[8]*XH2O/XH/pow(r[1], 3.0/2.0); // X46  HS- s29 (SO2 + H2O = HS- + H+ + 3/2 O2)
        double gS2 = [[endmembers objectAtIndex:47] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS2+2.0*gH+5.0*gO2/2.0-2.0*gSO2-gH2O)/R/t);
        s[30] = K*r[8]*r[8]*XH2O/XH/XH/pow(r[1], 5.0/2.0); // X47  S2-2 s30 (2SO2 + H2O = S2-2 + 2H+ + 5/2 O2)
        double gS2O3 = [[endmembers objectAtIndex:48] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS2O3+2.0*gH+gO2-2.0*gSO2-gH2O)/R/t);
        s[31] = K*r[8]*r[8]*XH2O/XH/XH/r[1]; // X48  S2O3-2 s31 (2SO2 + H2O = S2O3-2 + 2H+ + O2)
        double gS2O4 = [[endmembers objectAtIndex:49] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS2O4+2.0*gH+gO2/2.0-2.0*gSO2-gH2O)/R/t);
        s[32] = K*r[8]*r[8]*XH2O/XH/XH/sqrt(r[1]); // X49  S2O4-2 s32 (2SO2 + H2O = S2O4-2 + 2H+ + 1/2 O2)
        double gS2O5 = [[endmembers objectAtIndex:50] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS2O5+2.0*gH-2.0*gSO2-gH2O)/R/t);
        s[33] = K*r[8]*r[8]*XH2O/XH/XH; // X50  S2O5-2 s33 (2SO2 + H2O = S2O5-2 + 2H+)
        double gS2O6 = [[endmembers objectAtIndex:51] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS2O6+2.0*gH-2.0*gSO2-gH2O-gO2/2.0)/R/t);
        s[34] = K*r[8]*r[8]*XH2O*sqrt(r[1])/XH/XH; // X51  S2O6-2 s34 (2SO2 + H2O + 1/2 O2 = S2O6-2 + 2H+)
        double gS2O8 = [[endmembers objectAtIndex:52] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS2O8+2.0*gH-2.0*gSO2-gH2O-3.0*gO2/2.0)/R/t);
        s[35] = K*r[8]*r[8]*XH2O*pow(r[1], 3.0/2.0)/XH/XH; // X52  S2O8-2 s35 (2SO2 + H2O + 3/2 O2 = S2O8-2 + 2H+)
        double gS3 = [[endmembers objectAtIndex:53] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS3+2.0*gH+7.0*gO2/2.0-3.0*gSO2-gH2O)/R/t);
        s[36] = K*r[8]*r[8]*r[8]*XH2O/XH/XH/pow(r[1], 7.0/2.0); // X53  S3-2 s36 (3SO2 + H2O = S3-2 + 2H+ + 7/2 O2)
        double gS3O6 = [[endmembers objectAtIndex:54] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS3O6+2.0*gH+gO2/2.0-3.0*gSO2-gH2O)/R/t);
        s[37] = K*r[8]*r[8]*r[8]*XH2O/XH/XH/sqrt(r[1]); // X54  S3O6-2   s37 (3SO2 + H2O = S3O6-2 + 2H+ + 1/2 O2)
        double gS4 = [[endmembers objectAtIndex:55] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS4+2.0*gH+9.0*gO2/2.0-4.0*gSO2-gH2O)/R/t);
        s[38] = K*r[8]*r[8]*r[8]*r[8]*XH2O/XH/XH/pow(r[1], 9.0/2.0); // X55  S4-2     s38 (4SO2 + H2O = S4-2 + 2H+ + 9/2 O2)
        double gS4O6 = [[endmembers objectAtIndex:56] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS4O6+2.0*gH+3.0*gO2/2.0-4.0*gSO2-gH2O)/R/t);
        s[39] = K*r[8]*r[8]*r[8]*r[8]*XH2O/XH/XH/pow(r[1], 3.0/2.0); // X56  S4O6-2   s39 (4SO2 + H2O = S4O6-2 + 2H+ + 3/2 O2)
        double gS5 = [[endmembers objectAtIndex:57] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS5+2.0*gH+11.0*gO2/2.0-5.0*gSO2-gH2O)/R/t);
        s[40] = K*r[8]*r[8]*r[8]*r[8]*r[8]*XH2O/XH/XH/pow(r[1], 11.0/2.0); // X57  S5-2 s40 (5SO2 + H2O = S5-2 + 2H+ + 11/2 O2)
        double gS5O6 = [[endmembers objectAtIndex:58] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gS5O6+2.0*gH+5.0*gO2/2.0-5.0*gSO2-gH2O)/R/t);
        s[41] = K*r[8]*r[8]*r[8]*r[8]*r[8]*XH2O/XH/XH/pow(r[1], 5.0/2.0); // X58  S5O6-2   s41 (5SO2 + H2O = S5O6-2 + 2H+ + 5/2 O2)
        double gSO3 = [[endmembers objectAtIndex:59] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gSO3+2.0*gH-gSO2-gH2O)/R/t);
        s[42] = K*r[8]*XH2O/XH/XH; // X59  SO3-2 s42 (SO2 + H2O = SO3-2 + 2H+)
        double gHSO3 = [[endmembers objectAtIndex:60] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHSO3+gH-gSO2-gH2O)/R/t);
        s[43] = K*r[8]*XH2O/XH; // X60  HSO3- s43 (SO2 + H2O = HSO3- + H+)
        double gSO4 = [[endmembers objectAtIndex:61] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gSO4+2.0*gH-gSO2-gH2O-gO2/2.0)/R/t);
        s[44] = K*r[8]*XH2O*sqrt(r[1])/XH/XH; // X61  SO4-2 s44 (SO2 + H2O + 1/2 O2 = SO4-2 + 2H+)
        double gHSO4 = [[endmembers objectAtIndex:62] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHSO4+gH-gSO2-gH2O-gO2/2.0)/R/t);
        s[45] = K*r[8]*XH2O*sqrt(r[1])/XH; // X62  HSO4-    s45 (SO2 + H2O + 1/2 O2 = HSO4- + H+)
        double gSO5 = [[endmembers objectAtIndex:63] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gSO5+gH-gSO2-gH2O-gO2)/R/t);
        s[46] = K*r[8]*XH2O*r[1]/XH; // X63  HSO5-    s46 (SO2 + H2O + O2 = HSO5- + H+)
        double total = 0.0;
        for (NSUInteger i=28; i<47; i++) total += s[i];
        total = (total >= r[8]) ? (1.0-100.0*DBL_EPSILON)*r[8]/total : 1.0;
        for (NSUInteger i=28; i<47; i++) s[i] *= total;
    }
    if (r[9] > 0.0) { // X10  HCl     r9
        double gHCl = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
        double gCl  = [[endmembers objectAtIndex:64] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCl+gH-gHCl)/R/t);
        s[47] = K*r[9]/XH; // X64  Cl- s47 HCl = Cl- + H+
        double total = s[47];
        total = (total >= r[9]) ? (1.0-100.0*DBL_EPSILON)*r[9]/total : 1.0;
        s[47] *= total;
        if (debug) NSLog(@"HCl %g Cl- %g", r[9], s[47]);
    }
    if (r[10] > 0.0) { // X11  KOH r10
        double gKOH = [[endmembers objectAtIndex:11] getGibbsFreeEnergyFromT:t andP:p];
        double gK   = [[endmembers objectAtIndex:65] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gK+gH2O-gKOH-gH)/R/t);
        s[48] = K*r[10]*XH/XH2O; // X65  K+ s48 (KOH + H+ = K+ + H2O)
        if (r[9] != 0.0) {
            double gHCl = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
            double gKCl = [[endmembers objectAtIndex:66] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gKCl+gH2O-gKOH-gHCl)/R/t);
            s[49] = K*r[10]*r[9]/XH2O; // X66  KCl   s49 (KOH + HCl = KCl + H2O)
        }
        if (r[8] != 0.0) {
            double gSO2  = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
            double gKSO4 = [[endmembers objectAtIndex:67] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gKSO4+gH-gKOH-gSO2-gO2/2.0)/R/t);
            s[50] = K*r[10]*r[8]*sqrt(r[1])/XH; // X67  KSO4- s50 (KOH + SO2 + 1/2 O2 = KSO4- + H+)
        }
        double total = s[48] + s[49] + s[50];
        total = (total >= r[10]) ? (1.0-100.0*DBL_EPSILON)*r[10]/total : 1.0;
        s[48] *= total;
        s[49] *= total;
        s[50] *= total;
        if (debug) NSLog(@"KOH %g K+ %g KCl %g KSO4- %g", r[10], s[48], s[49], s[50]);
    }
    if (r[11] > 0.0) { // X12  Ca(OH)2 r11
        double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
        double gCaOH   = [[endmembers objectAtIndex:70] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCaOH+gH2O-gCaO2H2-gH)/R/t);
        s[53] = K*r[11]*XH/XH2O; // X70  CaOH+ s53 (Ca(OH)2 + H+ = CaOH+ + H2O)
        double gCa = [[endmembers objectAtIndex:71] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCa+2.0*gH2O-gCaO2H2-2.0*gH)/R/t);
        s[54] = K*r[11]*XH*XH/XH2O/XH2O; // X71  Ca+2  s54 (Ca(OH)2 + 2H+ = Ca+2 + 2 H2O)
        if (r[0] != 0.0) {
            double gCO2   = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
            double gCaCO3 = [[endmembers objectAtIndex:68] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gCaCO3+gH2O-gCaO2H2-gCO2)/R/t);
            s[51] = K*r[0]*r[11]/XH2O; // X68  CaCO3   s51 (Ca(OH)2 + CO2 = CaCO3 + H2O)
            double gCaHCO3 = [[endmembers objectAtIndex:69] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gCaHCO3+gOH-gCaO2H2-gCO2)/R/t);
            s[52] = K*r[0]*r[11]/XOH; // X69  CaHCO3+ s52 (Ca(OH)2 + CO2 = CaHCO3+ + OH-)
        }
        if (r[9] != 0.0) {
            double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
            double gCaCl = [[endmembers objectAtIndex:72] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gCaCl+2.0*gH2O-gCaO2H2-gHCl-gH)/R/t);
            s[55] = K*r[9]*r[11]*XH/XH2O/XH2O; // X72  CaCl+ s55 (Ca(OH)2 + HCl + H+ = CaCl+ + 2 H2O)
            double gCaCl2 = [[endmembers objectAtIndex:73] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gCaCl2+2.0*gH2O-gCaO2H2-2.0*gHCl)/R/t);
            s[56] = K*r[9]*r[9]*r[11]/XH2O/XH2O; // X73  CaCl2 s56 (Ca(OH)2 + 2 HCl = CaCl2 + 2 H2O)
        }
        if (r[8] != 0.0) {
            double gSO2   = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
            double gCaSO4 = [[endmembers objectAtIndex:74] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gCaSO4+gH2O-gCaO2H2-gSO2-gO2/2.0)/R/t);
            s[57] = K*r[8]*r[11]*sqrt(r[1])/XH2O; // X74  CaSO4 s57 (Ca(OH)2 + SO2 + 1/2 O2 = CaSO4 + H2O)
        }
        double total = s[51] + s[52] + s[53] + s[54] + s[55] + s[56] + s[57];
        total = (total >= r[11]) ? (1.0-100.0*DBL_EPSILON)*r[11]/total : 1.0;
        s[51] *= total;
        s[52] *= total;
        s[53] *= total;
        s[54] *= total;
        s[55] *= total;
        s[56] *= total;
        s[57] *= total;
        if (debug) NSLog(@"Ca(OH)2 %g CaCO3 %g CaHCO3+ %g CaOH+ %g Ca++ %g CaCl+ %g CaCl2 %g CaSO4 %g",
            r[11], s[51], s[52], s[53], s[54], s[55], s[56], s[57]);
    }
    if (r[12] > 0.0) { // X13  H2CrO4  r12
        double gH2CrO4 = [[endmembers objectAtIndex:13] getGibbsFreeEnergyFromT:t andP:p];
        double gCr     = [[endmembers objectAtIndex:75] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCr+2.0*gH2O+gO2-gH2CrO4-2.0*gH)/R/t);
        s[58] = K*r[12]*XH*XH/XH2O/XH2O/r[1]; // X75  Cr+2 s58 (H2CrO4 + 2 H+ = Cr+2 + 2 H2O + O2)
        gCr = [[endmembers objectAtIndex:76] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCr+5.0*gH2O/2.0+3.0*gO2/2.0-gH2CrO4-3.0*gH)/R/t);
        s[59] = K*r[12]*XH*XH*XH/pow(XH2O,5.0/2.0)/pow(r[1],3.0/2.0); // X76  Cr+3 s59 (H2CrO4 + 3 H+ = Cr+3 + 5/2 H2O + 3/2 O2)
        double gCr2O7 = [[endmembers objectAtIndex:77] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCr2O7+2.0*gH+gH2O-2.0*gH2CrO4)/R/t);
        s[60] = K*r[12]*r[12]/XH/XH/XH2O; // X77  Cr2O7-2  s60 (2H2CrO4 = Cr2O7-2 + 2H+ + H2O)
        double gCrO4 = [[endmembers objectAtIndex:78] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCrO4+2.0*gH-gH2CrO4)/R/t);
        s[61] = K*r[12]/XH/XH; // X78  CrO4-2 s61 (H2CrO4 = CrO4-2 + 2H+)
        double gHCrO4 = [[endmembers objectAtIndex:79] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gHCrO4+gH-gH2CrO4)/R/t);
        s[62] = K*r[12]/XH; // X79  HCrO4- s62 (H2CrO4 = HCrO4- + H+)
        double total = s[58] + s[59] + s[60] + s[61] + s[62];
        total = (total >= r[12]) ? (1.0-100.0*DBL_EPSILON)*r[12]/total : 1.0;
        s[58] *= total;
        s[59] *= total;
        s[60] *= total;
        s[61] *= total;
        s[62] *= total;
        if (debug) NSLog(@"H2CrO4 %g Cr2+ %g Cr3+ %g Cr2O7-2 %g CrO4-2 %g HCrO4- %g",
            r[12], s[58], s[59], s[60], s[61], s[62]);
    }
    if (r[13] > 0.0) { // X14  Mn(OH)2 r13
        double gMnO2H2 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
        double gMn     = [[endmembers objectAtIndex:80] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gMn+2.0*gOH-gMnO2H2)/R/t);
        s[63] = K*r[13]/XOH/XOH; // X80  Mn+2 s63 (Mn(OH)2 = Mn+2 + 2 OH-)
        double gMnO4 = [[endmembers objectAtIndex:82] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gMnO4+gH+gH2O/2.0-gMnO2H2-5.0*gO2/4.0)/R/t);
        s[65] = K*r[13]*pow(r[1],5.0/4.0)/XH/sqrt(XH2O); // X82  MnO4- s65 (Mn(OH)2 + 5/4 O2 = MnO4- + H+ + 1/2 H2O)
        gMnO4 = [[endmembers objectAtIndex:83] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gMnO4+2.0*gH-gMnO2H2-gO2)/R/t);
        s[66] = K*r[13]*r[1]/XH/XH; // X83  MnO4-2   s66 (Mn(OH)2 + O2 = MnO4-2 + 2 H+)
        if (r[9] != 0.0) {
            double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
            double gMnCl = [[endmembers objectAtIndex:81] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMnCl+2.0*gH2O-gMnO2H2-gHCl-gH)/R/t);
            s[64] = K*r[9]*r[13]*XH/XH2O/XH2O; // X81  MnCl+ s64 (Mn(OH)2 + HCl + H+ = MnCl+ + 2 H2O)
        }
        if (r[8] != 0.0) {
            double gSO2   = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
            double gMnSO4 = [[endmembers objectAtIndex:84] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gMnSO4+gH2O-gMnO2H2-gSO2-gO2/2.0)/R/t);
            s[67] = K*r[8]*r[13]*sqrt(r[1])/XH2O; // X84  MnSO4 s67 (Mn(OH)2 + SO2 + 1/2 O2 = MnSO4 + H2O)
        }
        double total = s[63] + s[64] + s[65] + s[66] + s[67];
        total = (total >= r[13]) ? (1.0-100.0*DBL_EPSILON)*r[13]/total : 1.0;
        s[63] *= total;
        s[64] *= total;
        s[65] *= total;
        s[66] *= total;
        s[67] *= total;
        if (debug) NSLog(@"Mn(OH)2 %g Mn++ %g MnCl+ %g MnO4- %g MnO4-2 %g MnSO4 %g",
            r[13], s[63], s[64], s[65], s[66], s[67]);
    }
    if (r[14] > 0.0) { // X15  Fe(OH)2 r14
        double gFeO2H2 = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
        double gFe     = [[endmembers objectAtIndex:85] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gFe+gH2O-gFeO2H2-2.0*gH)/R/t);
        s[68] = K*r[14]*XH*XH/XH2O/XH2O; // X85  Fe2+  s68 (Fe(OH)2 + 2 H+ = Fe2+ + 2 H2O)
        gFe     = [[endmembers objectAtIndex:86] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gFe+5.0*gH2O/2.0-gFeO2H2-3.0*gH-gO2/2.0)/R/t);
        s[69] = K*r[14]*XH*XH*XH*pow(r[1],1.0/4.0)/pow(XH2O,5.0/2.0); // X86  Fe3+  s69 (Fe(OH)2 + 3 H+ + 1/4 O2 = Fe3+ + 5/2 H2O)
        if (r[9] != 0.0) {
            double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
            double gFeCl = [[endmembers objectAtIndex:87] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gFeCl+2.0*gH2O-gFeO2H2-gHCl-gH)/R/t);
            s[70] = K*r[9]*r[14]*XH/XH2O/XH2O; // X87  FeCl+    s70 (Fe(OH)2 + HCl + H+ = FeCl+ + 2 H2O)
            gFeCl = [[endmembers objectAtIndex:88] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gFeCl+5.0*gH2O/2.0-gFeO2H2-gHCl-2.0*gH-gO2/4.0)/R/t);
            s[71] = K*r[9]*r[14]*XH*XH/pow(r[1],1.0/4.0)/pow(XH2O,5.0/2.0); // X88  FeCl+2 s71 (Fe(OH)2 + HCl + 2 H+ + 1/4 O2 = FeCl+2 + 5/2 H2O)
            double gFeCl2 = [[endmembers objectAtIndex:89] getGibbsFreeEnergyFromT:t andP:p];
            K = exp(-(gFeCl2+2.0*gH2O-gFeO2H2-2.0*gHCl)/R/t);
            s[72] = K*r[14]*r[9]*r[9]/XH2O/XH2O; // X89  FeCl2 s72 (Fe(OH)2 + 2 HCl = FeCl2 + 2 H2O)
        }
        double total = s[68] + s[69] + s[70] + s[71] + s[72];
        total = (total >= r[14]) ? (1.0-100.0*DBL_EPSILON)*r[14]/total : 1.0;
        s[68] *= total;
        s[69] *= total;
        s[70] *= total;
        s[71] *= total;
        s[72] *= total;
        if (debug) NSLog(@"Fe(OH)2 %g Fe+2 %g Fe+3 %g FeCl+ %g FeCl+2 %g FeCl2 %g",
            r[14], s[68], s[69], s[70], s[71], s[72]);
    }
    if (r[15] > 0.0) { // X16  Co(OH)2 r15
        double gCoO2H2 = [[endmembers objectAtIndex:16] getGibbsFreeEnergyFromT:t andP:p];
        double gCo     = [[endmembers objectAtIndex:90] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCo+2.0*gH2O-gCoO2H2-2.0*gH)/R/t);
        s[73] = K*r[15]*XH*XH/XH2O/XH2O; // X90  Co+2 s73 (Co(OH)2 + 2 H+ = Co2+ + 2 H2O)
        gCo = [[endmembers objectAtIndex:91] getGibbsFreeEnergyFromT:t andP:p];
        K = exp(-(gCo+5.0*gH2O/2.0-gCoO2H2-3.0*gH-gO2/4.0)/R/t);
        s[74] = K*r[15]*XH*XH*XH*pow(r[1],1.0/4.0)/pow(XH2O,5.0/2.0); // X91  Co+3 s74 (Co(OH)2 + 3 H+ + 1/4 O2 = Co3+ + 5/2 H2O)
        double total = s[73] + s[74];
        total = (total >= r[15]) ? (1.0-100.0*DBL_EPSILON)*r[15]/total : 1.0;
        s[73] *= total;
        s[74] *= total;
        if (debug) NSLog(@"Co(OH)2 %g Co+2 %g Co+3 %g", r[15], s[73], s[74]);
    }

    for (NSUInteger i=0; i<NS; i++)
    	if ((s[i] > 0.0) && (s[i] <= DBL_EPSILON)) {
    		lowConcFlagForS[i] = YES;
    		s[i] = 0.0;
    	}

    double cations = s[0] + s[8] + s[13] + s[14] + 2.0*s[15] + s[16] + s[17] + 3.0*s[19] + s[48] + s[52] + s[53] + 2.0*s[54]
    + s[55] + 2.0*s[58] + 3.0*s[59] + 2.0*s[63] + s[64] + 2.0*s[68] + 3.0*s[69] + s[70] + 2.0*s[71]
    + 2.0*s[73] + 3.0*s[74];
    double anions = s[1] + 2.0*s[3] + s[4] + s[6] + s[9] + s[20] + s[21] + s[23] + 2.0*s[24] + 3.0*s[25] + s[26] + 2.0*s[27]
    + s[29] + 2.0*s[30] + 2.0*s[31] + 2.0*s[32] + 2.0*s[33] + 2.0*s[34] + 2.0*s[35] + 2.0*s[36] + 2.0*s[37]
    + 2.0*s[38] + 2.0*s[39] + 2.0*s[40] + 2.0*s[41] + 2.0*s[42] + s[43] + 2.0*s[44] + s[45] + s[46] + s[47]
    + s[50] + 2.0*s[60] + 2.0*s[61] + s[62] + s[65] + 2.0*s[66];

    if (cations > anions) {
        // excess cations
        s[1] += cations-anions;  // X18  OH-      s1
    } else if (anions > cations) {
        // excess anions
        s[0] += anions-cations;  // X17  H+       s0
    }

    BOOL feasible = [self fillxSpeciesWith:r andWith:s];
    NSUInteger loopNumber = 0;
retry:
    if (!feasible) {
        if (debug) NSLog(@"Estimate from initialGuessOrdering is not feasible. Correcting...");
        double coeff = [self speciationCoefficientWith:s];
        if (debug) NSLog(@"... mole correction coefficient is %g", coeff);
        for (NSUInteger i=0; i<NA; i++) {
            if (xSpecies[i] < 0.0) {
                switch (i) {
                    case 1:
                    {
                        double total = s[3] + s[4] + s[5] + s[9] + s[10] + s[12] + s[14] + s[51] + s[52];
                        total = (total >= r[0]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[0]*coeff/total : 1.0;
                        s[ 3] *= total;
                        s[ 4] *= total;
                        s[ 5] *= total;
                        s[ 9] *= total;
                        s[10] *= total;
                        s[12] *= total;
                        s[14] *= total;
                        s[51] *= total;
                        s[52] *= total;
                        if (debug) NSLog(@"CO2 %g CO3-2 %g HCO3- %g CO %g", r[0], s[3], s[4], s[5]);
                        if (debug) NSLog(@"NaCO3- %g NaHCO3 %g MgCO3 %g MgHCO3+ %g CaCO3 %g CaHCO3+ %g", s[9], s[10], s[12], s[14], s[51], s[52]);
                    }
                    case 4:
                    {
                        double total = s[8] + s[9] + s[10] + s[7] + s[11];
                        total = (total >= r[3]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[3]*coeff/total : 1.0;
                        s[ 7] *= total;
                        s[ 8] *= total;
                        s[ 9] *= total;
                        s[10] *= total;
                        s[11] *= total;
                        if (debug) NSLog(@"NaOH %g Na+ %g NaHSiO3 %g NaCl %g NaCO3- %g NaHCO3 %g", r[3], s[8], s[11], s[7], s[9], s[10]);
                        break;
                    }
                    case 5:
                    {
                        double total = s[13] + s[14] + s[15] + s[16] + s[17] + s[18];
                        total = (total >= r[4]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[4]*coeff/total : 1.0;
                        s[13] *= total;
                        s[14] *= total;
                        s[15] *= total;
                        s[16] *= total;
                        s[17] *= total;
                        s[18] *= total;
                        if (debug) NSLog(@"Mg(OH)2 %g MgHSiO3+ %g MgHCO3+ %g Mg++ %g MgCl+ %g MgOH+ %g MgSO4 %g",
                            r[4], s[13], s[14], s[15], s[16], s[17], s[18]);
                        break;
                    }
                    case 7:
                    {
                        double total = s[21] + 2.0*s[22] + s[11] +s[13];
                        total = (total >= r[6]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[6]*coeff/total : 1.0;
                        s[21] *= total;
                        s[22] *= total;
                        s[11] *= total;
                        s[13] *= total;
                        if (debug) NSLog(@"SiO2 %g HSiO3- %g Si2O4 %g NaHSiO3 %g Mg(HSiO3)+ %g, diff = %g",
                            r[6], s[21], s[22], s[11], s[13], r[6]*coeff-s[21]-2.0*s[22]-s[11]-s[13]);
                        break;
                    }
                    case 8: {
                        double total = s[23] + s[24] + s[25] + s[26] + s[27];
                        total = (total >= r[7]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[7]*coeff/total : 1.0;
                        s[23] *= total;
                        s[24] *= total;
                        s[25] *= total;
                        s[26] *= total;
                        s[27] *= total;
                        if (debug) NSLog(@"H3PO4 %g H2PO4- %g HPO4-2 %g PO4-3 %g H3P2O7- %g H2P2O7-2 %g",
                            r[7], s[23], s[24], s[25], s[26], s[27]);
                        break;
                    }
                    case 11:
                    {
                        double total = s[48] + s[49] + s[50];
                        total = (total >= r[10]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[10]*coeff/total : 1.0;
                        s[48] *= total;
                        s[49] *= total;
                        s[50] *= total;
                        if (debug) NSLog(@"KOH %g K+ %g KCl %g KSO4- %g", r[10], s[48], s[49], s[50]);
                        break;
                    }
                    case 12:
                    {
                        double total = s[51] + s[52] + s[53] + s[54] + s[55] + s[56] + s[57];
                        total = (total >= r[11]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[11]*coeff/total : 1.0;
                        s[51] *= total;
                        s[52] *= total;
                        s[53] *= total;
                        s[54] *= total;
                        s[55] *= total;
                        s[56] *= total;
                        s[57] *= total;
                        if (debug) NSLog(@"Ca(OH)2 %g CaCO3 %g CaHCO3+ %g CaOH+ %g Ca++ %g CaCl+ %g CaCl2 %g CaSO4 %g",
                            r[11], s[51], s[52], s[53], s[54], s[55], s[56], s[57]);
                        break;
                    }
                    case 13:
                    {
                        double total = s[58] + s[59] + s[60] + s[61] + s[62];
                        total = (total >= r[12]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[12]*coeff/total : 1.0;
                        s[58] *= total;
                        s[59] *= total;
                        s[60] *= total;
                        s[61] *= total;
                        s[62] *= total;
                        if (debug) NSLog(@"H2CrO4 %g Cr2+ %g Cr3+ %g Cr2O7-2 %g CrO4-2 %g HCrO4- %g",
                            r[12], s[58], s[59], s[60], s[61], s[62]);
                        break;
                    }
                    case 15:
                    {
                        double total = s[68] + s[69] + s[70] + s[71] + s[72];
                        total = (total >= r[14]*coeff) ? (1.0-100.0*DBL_EPSILON)*r[14]*coeff/total : 1.0;
                        s[68] *= total;
                        s[69] *= total;
                        s[70] *= total;
                        s[71] *= total;
                        s[72] *= total;
                        if (debug) NSLog(@"Fe(OH)2 %g Fe+2 %g Fe+3 %g FeCl+ %g FeCl+2 %g FeCl2 %g",
                            r[14], s[68], s[69], s[70], s[71], s[72]);
                        break;
                    }
                    default:
                        if (debug) NSLog(@"<><><> Case %lu not properly accounted for.", i);
                        break;
                }
            }
        }
        double cations = s[0] + s[8] + s[13] + s[14] + 2.0*s[15] + s[16] + s[17] + 3.0*s[19] + s[48] + s[52] + s[53] + 2.0*s[54]
            + s[55] + 2.0*s[58] + 3.0*s[59] + 2.0*s[63] + s[64] + 2.0*s[68] + 3.0*s[69] + s[70] + 2.0*s[71]
            + 2.0*s[73] + 3.0*s[74];
        double anions = s[1] + 2.0*s[3] + s[4] + s[6] + s[9] + s[20] + s[21] + s[23] + 2.0*s[24] + 3.0*s[25] + s[26] + 2.0*s[27]
            + s[29] + 2.0*s[30] + 2.0*s[31] + 2.0*s[32] + 2.0*s[33] + 2.0*s[34] + 2.0*s[35] + 2.0*s[36] + 2.0*s[37]
            + 2.0*s[38] + 2.0*s[39] + 2.0*s[40] + 2.0*s[41] + 2.0*s[42] + s[43] + 2.0*s[44] + s[45] + s[46] + s[47]
            + s[50] + 2.0*s[60] + 2.0*s[61] + s[62] + s[65] + 2.0*s[66];

        if (cations > anions) {
            // excess cations
            s[1] += cations-anions;  // X18  OH-      s1
        } else if (anions > cations) {
            // excess anions
            s[0] += anions-cations;  // X17  H+       s0
        }
    }
    feasible = [self fillxSpeciesWith:r andWith:s];
    if (!feasible) {
        if (debug) {
            NSLog(@"Estimate from initialGuessOrdering is not feasible even after correction.");
            for (NSUInteger i=0; i<NE; i++)
                NSLog(@"<><><> Species %@ [%2.2lu] = %13.6g %@",
                    [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0], i, xSpecies[i],
                    (xSpecies[i] < 0.0) ? @"has negative concentration" : @"");
        }
        loopNumber++;
        if (loopNumber < 10) goto retry;
        NSAssert(NO, @"ERROR");
    }

}

- (void)adjustOrderingParameters:(double [NS])s from:(double [NR])r atT:(double)t andP:(double)p {
    double rTotal = 0.0;
    for (NSUInteger i=0; i<NR; i++) rTotal += r[i];
    double XH2O = 1.0 - rTotal;
    double K;

    double gH2O = [[endmembers objectAtIndex: 0] getGibbsFreeEnergyFromT:t andP:p];
    double gH   = [[endmembers objectAtIndex:17] getGibbsFreeEnergyFromT:t andP:p];
    double gOH  = [[endmembers objectAtIndex:18] getGibbsFreeEnergyFromT:t andP:p];
    double gO2  = [[endmembers objectAtIndex: 2] getGibbsFreeEnergyFromT:t andP:p];

    for (NSUInteger i=0; i<NS; i++) if (lowConcFlagForS[i]) {
    	switch (i) {
    		case 0:
            {
                K = exp(-(gH+gOH-gH2O)/R/t);
                s[0] = sqrt(K*XH2O); // XH2O*1.0e-6;  // X17  H+  s0 (H2O = H+ + OH-)
                break;
            }
    		case 1:
    		{
   				K = exp(-(gH+gOH-gH2O)/R/t);
    			s[1] = sqrt(K*XH2O); // XH2O*1.0e-6;  // X18  OH- s1 (H2O = H+ + OH-)
    			break;
    		}
    		case 3:
            {
                double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                double gCO3 = [[endmembers objectAtIndex:20] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCO3+2.0*gH-gH2O-gCO2)/R/t);
                s[3] = K*r[0]*XH2O/s[0]/s[0]; // X20  CO3=  s3 (CO2 + H2O = CO3= + 2H+)
                break;
            }
    		case 4:
            {
                double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                double gHCO3 = [[endmembers objectAtIndex:21] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gHCO3+gH-gH2O-gCO2)/R/t);
                s[4] = K*r[0]/s[0]; // X21  HCO3- s4 (CO2 + H2O = HCO3- + H+)
                break;
            }
    		case 5:
    		{
        		double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                double gCO = [[endmembers objectAtIndex:22] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCO+gO2/2.0-gCO2)/R/t);
                s[5] = K*r[0]/sqrt(s[2]); // X22  CO    s5 (CO2 = CO + O2/2)
        		break;
    		}
    		case 2:
    		{
        		double gH2 = [[endmembers objectAtIndex:19] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gH2+gO2/2.0-gH2O)/R/t);
                s[2] = K*XH2O/sqrt(r[1]);  // X19  H2 s2 (H2O = H2 + O2/2)
        		break;
    		}
    		case 6:
			{
        		double gHF = [[endmembers objectAtIndex: 3] getGibbsFreeEnergyFromT:t andP:p];
                double gF  = [[endmembers objectAtIndex:23] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gF+gH-gHF)/R/t);
                s[6] = K*r[2]/s[0]; // X23  F- s6 (HF = H+ + F-)
        		break;
    		}
    		case 8:
    		{
        		double gNaOH = [[endmembers objectAtIndex: 4] getGibbsFreeEnergyFromT:t andP:p];
        		double gNa   = [[endmembers objectAtIndex:25] getGibbsFreeEnergyFromT:t andP:p];
        		K = exp(-(gNa+gOH-gNaOH)/R/t);
        		s[8] = K*r[3]/s[1]; // r[3]*0.1; // X25  Na+ s8 (NaOH = Na+ + OH-)
        		break;
        	}
        	case 9:
        	{
                double gNaOH = [[endmembers objectAtIndex: 4] getGibbsFreeEnergyFromT:t andP:p];
            	double gNaCO3 = [[endmembers objectAtIndex:26] getGibbsFreeEnergyFromT:t andP:p];
                double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gNaCO3+gH-gNaOH-gCO2)/R/t);
                s[ 9] = K*r[0]*r[3]/s[0]; // X26  NaCO3- s9  (NaOH + CO2 = NaCO3- + H+)
            	break;
            }
            case 10:
            {
                double gNaOH = [[endmembers objectAtIndex: 4] getGibbsFreeEnergyFromT:t andP:p];
            	double gNaHCO3 = [[endmembers objectAtIndex:27] getGibbsFreeEnergyFromT:t andP:p];
                double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gNaHCO3-gNaOH-gCO2)/R/t);
                s[10] = K*r[0]*r[3]; // X27  NaHCO3 s10 (NaOH + CO2 = NaHCO3)
            	break;
            }
        	case 7:
        	{
                double gNaOH = [[endmembers objectAtIndex: 4] getGibbsFreeEnergyFromT:t andP:p];
            	double gNaCl = [[endmembers objectAtIndex:24] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gNaCl+gH2O-gNaOH-gHCl)/R/t);
                s[ 7] = K*r[9]*r[3]/XH2O; // X24  NaCl s7 (NaOH + HCl = NaCl + H2O)
            	break;
        	}
        	case 11:
        	{
            	double gSiO2  = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
            	double gHSiO3 = [[endmembers objectAtIndex:28] getGibbsFreeEnergyFromT:t andP:p];
            	double gNaOH = [[endmembers objectAtIndex: 4] getGibbsFreeEnergyFromT:t andP:p];
            	K = exp(-(gHSiO3-gNaOH-gSiO2)/R/t);
            	s[11] = K*r[3]*r[6]; // 0.1*MIN(r[3], r[6]); // X28  NaHSiO3 s11 (NaOH + SiO2 = NaHSiO3)
            	break;
        	}
        	case 12:
            {
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
                double gMgCO3 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMgCO3+gH2O-gMgO2H2-gCO2)/R/t);
                s[12] = K*r[4]*r[0]/XH2O; // X29  MgCO3    s12 Mg(OH)2 + CO2 = MgCO3 + H2O
                break;
            }
        	case 13:
            {
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
                double gMgHSiO3 = [[endmembers objectAtIndex:30] getGibbsFreeEnergyFromT:t andP:p];
                double gSiO2  = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMgHSiO3+gH2O-gMgO2H2-gSiO2)/R/t);
                s[13] = K*r[4]*r[6]*s[0]/XH2O; // X30  MgHSiO3+ s13 Mg(OH)2 + SiO2 + H+ = MgHSiO3+ + H2O
                break;
            }
        	case 14:
            {
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
                double gMgHCO3 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
                double gCO2 = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMgHCO3+gOH-gMgO2H2-gCO2)/R/t);
                s[14] = K*r[4]*r[0]/s[1]; // X31  MgHCO3+  s14 Mg(OH)2 + CO2 = MgHCO3+ + OH-
                break;
            }
        	case 15:
            {
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
                double gMg     = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMg+2.0*gOH-gMgO2H2)/R/t);
                s[15] = K*r[4]/s[1]/s[1]; // X32  Mg+2     s15 Mg(OH)2 = Mg++ + 2 OH-
                break;
            }
        	case 16:
            {
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
                double gMgCl = [[endmembers objectAtIndex:33] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMgCl+2.0*gH2O-gMgO2H2-gHCl)/R/t);
                s[16] = K*r[4]*r[9]*s[0]/XH2O/XH2O; // X33  MgCl+    s16 Mg(OH)2 + HCl + H+ = MgCl+ + 2 H2O
                break;
            }
        	case 17:
            {
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
                double gMgOH = [[endmembers objectAtIndex:17] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMgOH+gOH-gMgO2H2)/R/t);
                s[17] = K*r[4]/s[1]; // X34  MgOH+    s17 Mg(OH)2 = MgOH+ + OH-
                break;
            }
        	case 18:
        	{
                double gMgO2H2 = [[endmembers objectAtIndex: 5] getGibbsFreeEnergyFromT:t andP:p];
        		double gMgSO4 = [[endmembers objectAtIndex:35] getGibbsFreeEnergyFromT:t andP:p];
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMgSO4+gH2O-gMgO2H2-gSO2-gO2/2.0)/R/t);
                s[18] = K*r[4]*r[8]*sqrt(r[1])/XH2O; // X35  MgSO4    s18 Mg(OH)2 + SO2 + 1/2 O2 = MgSO4 + H2O
        		break;
    		}
    		case 19:
    		{
        		double gHAlO2 = [[endmembers objectAtIndex: 6] getGibbsFreeEnergyFromT:t andP:p];
        		double gAl    = [[endmembers objectAtIndex:36] getGibbsFreeEnergyFromT:t andP:p];
        		K = exp(-(gAl+2.0*gH2O-gHAlO2-3.0*gH)/R/t);
        		s[19] = K*r[5]*s[0]*s[0]*s[0]/XH2O/XH2O; // r[5]*0.1; // X36  Al+3  s19 (HAlO2 + 3H+ = Al+3 + 2H2O)
        		break;
        	}
        	case 20:
        	{
        		double gHAlO2 = [[endmembers objectAtIndex: 6] getGibbsFreeEnergyFromT:t andP:p];
        		double gAlO2  = [[endmembers objectAtIndex:37] getGibbsFreeEnergyFromT:t andP:p];
        		K = exp(-(gAlO2+gH-gHAlO2)/R/t);
        		s[20] = K*r[5]/s[0]; // r[5]*0.1; // X37  AlO2- s20 (HAlO2 = AlO2- + H+)
        		break;
        	}
    		case 21:
    		{
        		double gSiO2  = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
        		double gHSiO3 = [[endmembers objectAtIndex:38] getGibbsFreeEnergyFromT:t andP:p];
        		K = exp(-(gHSiO3+gH-gSiO2-gH2O)/R/t);
        		s[21] = K*r[6]*XH2O/s[0]; // r[6]*0.2; // X38  HSiO3- s21 (SiO2 + H2O = HSiO3- + H+)
        		break;
        	}
        	case 22:
        	{
        		double gSiO2  = [[endmembers objectAtIndex: 7] getGibbsFreeEnergyFromT:t andP:p];
        		double gSi2O4 = [[endmembers objectAtIndex:39] getGibbsFreeEnergyFromT:t andP:p];
        		K = exp(-(gSi2O4-2.0*gSiO2)/R/t);
        		s[22] = K*r[6]*r[6]; // r[6]*0.2; // X39  Si2O4  s22 (2SiO2 = Si2O4)
        		break;
    		}
    		case 23:
            {
                double gH3PO4 = [[endmembers objectAtIndex: 8] getGibbsFreeEnergyFromT:t andP:p];
                double gH2PO4 = [[endmembers objectAtIndex:40] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gH2PO4+gH-gH3PO4)/R/t);
                s[23] = K*r[7]/s[0]; // X40  H2PO4-   s23 (H3PO4 = H2PO4- + H+)
                break;
            }
    		case 24:
            {
                double gH3PO4 = [[endmembers objectAtIndex: 8] getGibbsFreeEnergyFromT:t andP:p];
                double gHPO4 = [[endmembers objectAtIndex:41] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gHPO4+2.0*gH-gH3PO4)/R/t);
                s[24] = K*r[7]/s[0]/s[0]; // X41  HPO4-2   s24 (H3PO4 = HPO4-2 + 2H+)
                break;
            }
    		case 25:
            {
                double gH3PO4 = [[endmembers objectAtIndex: 8] getGibbsFreeEnergyFromT:t andP:p];
                double gPO4 = [[endmembers objectAtIndex:42] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gPO4+3.0*gH-gH3PO4)/R/t);
                s[25] = K*r[7]/s[0]/s[0]/s[0]; // X42  PO4-3    s25 (H3PO4 = PO4-3 + 3H+)
                break;
            }
    		case 26:
            {
                double gH3PO4 = [[endmembers objectAtIndex: 8] getGibbsFreeEnergyFromT:t andP:p];
                double gH3P2O7 = [[endmembers objectAtIndex:43] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gH3P2O7+gH+gH2O-2.0*gH3PO4)/R/t);
                s[26] = K*r[7]*r[7]/s[0]/XH2O; // X43  H3P2O7-  s26 (2H3PO4 = H3P2O7- + H+ + H2O)
                break;
            }
    		case 27:
    		{
                double gH3PO4 = [[endmembers objectAtIndex: 8] getGibbsFreeEnergyFromT:t andP:p];
        		double gH2P2O7 = [[endmembers objectAtIndex:44] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gH2P2O7+2.0*gH+gH2O-2.0*gH3PO4)/R/t);
                s[27] = K*r[7]*r[7]/s[0]/s[0]/XH2O; // X44  H2P2O7-2 s27 (2H3PO4 = H2P2O7-2 + 2H+ + H2O)
        		break;
    		}
    		case 28:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gH2S = [[endmembers objectAtIndex:45] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gH2S+3.0*gO2/2.0-gSO2-gH2O)/R/t);
                s[28] = K*r[8]*XH2O/pow(r[1], 3.0/2.0); // X45  H2S s28 (SO2 + H2O = H2S + 3/2 O2)
                break;
            }
    		case 29:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gHS = [[endmembers objectAtIndex:46] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gHS+gH+3.0*gO2/2.0-gSO2-gH2O)/R/t);
                s[29] = K*r[8]*XH2O/s[0]/pow(r[1], 3.0/2.0); // X46  HS- s29 (SO2 + H2O = HS- + H+ + 3/2 O2)
                break;
            }
    		case 30:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS2 = [[endmembers objectAtIndex:47] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS2+2.0*gH+5.0*gO2/2.0-2.0*gSO2-gH2O)/R/t);
                s[30] = K*r[8]*r[8]*XH2O/s[0]/s[0]/pow(r[1], 5.0/2.0); // X47  S2-2 s30 (2SO2 + H2O = S2-2 + 2H+ + 5/2 O2)
                break;
            }
    		case 31:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS2O3 = [[endmembers objectAtIndex:48] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS2O3+2.0*gH+gO2-2.0*gSO2-gH2O)/R/t);
                s[31] = K*r[8]*r[8]*XH2O/s[0]/s[0]/r[1]; // X48  S2O3-2 s31 (2SO2 + H2O = S2O3-2 + 2H+ + O2)
                break;
            }
    		case 32:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS2O4 = [[endmembers objectAtIndex:49] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS2O4+2.0*gH+gO2/2.0-2.0*gSO2-gH2O)/R/t);
                s[32] = K*r[8]*r[8]*XH2O/s[0]/s[0]/sqrt(r[1]); // X49  S2O4-2 s32 (2SO2 + H2O = S2O4-2 + 2H+ + 1/2 O2)
                break;
            }
    		case 33:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS2O5 = [[endmembers objectAtIndex:50] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS2O5+2.0*gH-2.0*gSO2-gH2O)/R/t);
                s[33] = K*r[8]*r[8]*XH2O/s[0]/s[0]; // X50  S2O5-2 s33 (2SO2 + H2O = S2O5-2 + 2H+)
                break;
            }
    		case 34:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS2O6 = [[endmembers objectAtIndex:51] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS2O6+2.0*gH-2.0*gSO2-gH2O-gO2/2.0)/R/t);
                s[34] = K*r[8]*r[8]*XH2O*sqrt(r[1])/s[0]/s[0]; // X51  S2O6-2 s34 (2SO2 + H2O + 1/2 O2 = S2O6-2 + 2H+)
                break;
            }
    		case 35:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS2O8 = [[endmembers objectAtIndex:52] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS2O8+2.0*gH-2.0*gSO2-gH2O-3.0*gO2/2.0)/R/t);
                s[35] = K*r[8]*r[8]*XH2O*pow(r[1], 3.0/2.0)/s[0]/s[0]; // X52  S2O8-2 s35 (2SO2 + H2O + 3/2 O2 = S2O8-2 + 2H+)
                break;
            }
    		case 36:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS3 = [[endmembers objectAtIndex:53] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS3+2.0*gH+7.0*gO2/2.0-3.0*gSO2-gH2O)/R/t);
                s[36] = K*r[8]*r[8]*r[8]*XH2O/s[0]/s[0]/pow(r[1], 7.0/2.0); // X53  S3-2 s36 (3SO2 + H2O = S3-2 + 2H+ + 7/2 O2)
                break;
            }
    		case 37:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS3O6 = [[endmembers objectAtIndex:54] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS3O6+2.0*gH+gO2/2.0-3.0*gSO2-gH2O)/R/t);
                s[37] = K*r[8]*r[8]*r[8]*XH2O/s[0]/s[0]/sqrt(r[1]); // X54  S3O6-2   s37 (3SO2 + H2O = S3O6-2 + 2H+ + 1/2 O2)
                break;
            }
    		case 38:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS4 = [[endmembers objectAtIndex:55] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS4+2.0*gH+9.0*gO2/2.0-4.0*gSO2-gH2O)/R/t);
                s[38] = K*r[8]*r[8]*r[8]*r[8]*XH2O/s[0]/s[0]/pow(r[1], 9.0/2.0); // X55  S4-2     s38 (4SO2 + H2O = S4-2 + 2H+ + 9/2 O2)
                break;
            }
    		case 39:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS4O6 = [[endmembers objectAtIndex:56] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS4O6+2.0*gH+3.0*gO2/2.0-4.0*gSO2-gH2O)/R/t);
                s[39] = K*r[8]*r[8]*r[8]*r[8]*XH2O/s[0]/s[0]/pow(r[1], 3.0/2.0); // X56  S4O6-2   s39 (4SO2 + H2O = S4O6-2 + 2H+ + 3/2 O2)
                break;
            }
    		case 40:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS5 = [[endmembers objectAtIndex:57] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS5+2.0*gH+11.0*gO2/2.0-5.0*gSO2-gH2O)/R/t);
                s[40] = K*r[8]*r[8]*r[8]*r[8]*r[8]*XH2O/s[0]/s[0]/pow(r[1], 11.0/2.0); // X57  S5-2 s40 (5SO2 + H2O = S5-2 + 2H+ + 11/2 O2)
                break;
            }
    		case 41:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gS5O6 = [[endmembers objectAtIndex:58] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gS5O6+2.0*gH+5.0*gO2/2.0-5.0*gSO2-gH2O)/R/t);
                s[41] = K*r[8]*r[8]*r[8]*r[8]*r[8]*XH2O/s[0]/s[0]/pow(r[1], 5.0/2.0); // X58  S5O6-2   s41 (5SO2 + H2O = S5O6-2 + 2H+ + 5/2 O2)
                break;
            }
    		case 42:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gSO3 = [[endmembers objectAtIndex:59] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gSO3+2.0*gH-gSO2-gH2O)/R/t);
                s[42] = K*r[8]*XH2O/s[0]/s[0]; // X59  SO3-2 s42 (SO2 + H2O = SO3-2 + 2H+)
                break;
            }
    		case 43:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gHSO3 = [[endmembers objectAtIndex:60] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gHSO3+gH-gSO2-gH2O)/R/t);
                s[43] = K*r[8]*XH2O/s[0]; // X60  HSO3- s43 (SO2 + H2O = HSO3- + H+)
                break;
            }
    		case 44:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gSO4 = [[endmembers objectAtIndex:61] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gSO4+2.0*gH-gSO2-gH2O-gO2/2.0)/R/t);
                s[44] = K*r[8]*XH2O*sqrt(r[1])/s[0]/s[0]; // X61  SO4-2 s44 (SO2 + H2O + 1/2 O2 = SO4-2 + 2H+)
                break;
            }
    		case 45:
            {
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gHSO4 = [[endmembers objectAtIndex:62] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gHSO4+gH-gSO2-gH2O-gO2/2.0)/R/t);
                s[45] = K*r[8]*XH2O*sqrt(r[1])/s[0]; // X62  HSO4-    s45 (SO2 + H2O + 1/2 O2 = HSO4- + H+)
                break;
            }
    		case 46:
    		{
                double gSO2 = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
        		double gSO5 = [[endmembers objectAtIndex:63] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gSO5+gH-gSO2-gH2O-gO2)/R/t);
                s[46] = K*r[8]*XH2O*r[1]/s[0]; // X63  HSO5-    s46 (SO2 + H2O + O2 = HSO5- + H+)
        		break;
    		}
    		case 47:
    		{
        		double gHCl = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gCl  = [[endmembers objectAtIndex:64] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCl+gH-gHCl)/R/t);
                s[47] = K*r[9]/s[0]; // X64  Cl- s47 HCl = Cl- + H+
        		break;
    		}
    		case 48:
            {
                double gKOH = [[endmembers objectAtIndex:11] getGibbsFreeEnergyFromT:t andP:p];
                double gK   = [[endmembers objectAtIndex:65] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gK+gH2O-gKOH-gH)/R/t);
                s[48] = K*r[10]*s[0]/XH2O; // X65  K+ s48 (KOH + H+ = K+ + H2O)
                break;
            }
    		case 49:
            {
                double gKOH = [[endmembers objectAtIndex:11] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gKCl = [[endmembers objectAtIndex:66] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gKCl+gH2O-gKOH-gHCl)/R/t);
                s[49] = K*r[10]*r[9]/XH2O; // X66  KCl   s49 (KOH + HCl = KCl + H2O)
                break;
            }
    		case 50:
    		{
                double gKOH = [[endmembers objectAtIndex:11] getGibbsFreeEnergyFromT:t andP:p];
        		double gSO2  = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gKSO4 = [[endmembers objectAtIndex:67] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gKSO4+gH-gKOH-gSO2-gO2/2.0)/R/t);
                s[50] = K*r[10]*r[8]*sqrt(r[1])/s[0]; // X67  KSO4- s50 (KOH + SO2 + 1/2 O2 = KSO4- + H+)
        		break;
    		}
    		case 51:
            {
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gCO2   = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                double gCaCO3 = [[endmembers objectAtIndex:68] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCaCO3+gH2O-gCaO2H2-gCO2)/R/t);
                s[51] = K*r[0]*r[11]/XH2O; // X68  CaCO3   s51 (Ca(OH)2 + CO2 = CaCO3 + H2O)
                break;
            }
    		case 52:
            {
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gCO2   = [[endmembers objectAtIndex: 1] getGibbsFreeEnergyFromT:t andP:p];
                double gCaHCO3 = [[endmembers objectAtIndex:69] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCaHCO3+gOH-gCaO2H2-gCO2)/R/t);
                s[52] = K*r[0]*r[11]/s[1]; // X69  CaHCO3+ s52 (Ca(OH)2 + CO2 = CaHCO3+ + OH-)
                break;
            }
    		case 53:
            {
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gCaOH   = [[endmembers objectAtIndex:70] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCaOH+gH2O-gCaO2H2-gH)/R/t);
                s[53] = K*r[11]*s[0]/XH2O; // X70  CaOH+ s53 (Ca(OH)2 + H+ = CaOH+ + H2O)
                break;
            }
    		case 54:
            {
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gCa = [[endmembers objectAtIndex:71] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCa+2.0*gH2O-gCaO2H2-2.0*gH)/R/t);
                s[54] = K*r[11]*s[0]*s[0]/XH2O/XH2O; // X71  Ca+2  s54 (Ca(OH)2 + 2H+ = Ca+2 + 2 H2O)
                break;
            }
    		case 55:
            {
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gCaCl = [[endmembers objectAtIndex:72] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCaCl+2.0*gH2O-gCaO2H2-gHCl-gH)/R/t);
                s[55] = K*r[9]*r[11]*s[0]/XH2O/XH2O; // X72  CaCl+ s55 (Ca(OH)2 + HCl + H+ = CaCl+ + 2 H2O)
                break;
            }
    		case 56:
            {
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gCaCl2 = [[endmembers objectAtIndex:73] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCaCl2+2.0*gH2O-gCaO2H2-2.0*gHCl)/R/t);
                s[56] = K*r[9]*r[9]*r[11]/XH2O/XH2O; // X73  CaCl2 s56 (Ca(OH)2 + 2 HCl = CaCl2 + 2 H2O)
                break;
            }
    		case 57:
    		{
                double gCaO2H2 = [[endmembers objectAtIndex:12] getGibbsFreeEnergyFromT:t andP:p];
        		double gSO2   = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gCaSO4 = [[endmembers objectAtIndex:74] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCaSO4+gH2O-gCaO2H2-gSO2-gO2/2.0)/R/t);
                s[57] = K*r[8]*r[11]*sqrt(r[1])/XH2O; // X74  CaSO4 s57 (Ca(OH)2 + SO2 + 1/2 O2 = CaSO4 + H2O)
            	break;
        	}
        	case 58:
            {
                double gH2CrO4 = [[endmembers objectAtIndex:13] getGibbsFreeEnergyFromT:t andP:p];
                double gCr     = [[endmembers objectAtIndex:75] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCr+2.0*gH2O+gO2-gH2CrO4-2.0*gH)/R/t);
                s[58] = K*r[12]*s[0]*s[0]/XH2O/XH2O/r[1]; // X75  Cr+2 s58 (H2CrO4 + 2 H+ = Cr+2 + 2 H2O + O2)
                break;
            }
        	case 59:
            {
                double gH2CrO4 = [[endmembers objectAtIndex:13] getGibbsFreeEnergyFromT:t andP:p];
                double gCr = [[endmembers objectAtIndex:76] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCr+5.0*gH2O/2.0+3.0*gO2/2.0-gH2CrO4-3.0*gH)/R/t);
                s[59] = K*r[12]*s[0]*s[0]*s[0]/pow(XH2O,5.0/2.0)/pow(r[1],3.0/2.0); // X76  Cr+3 s59 (H2CrO4 + 3 H+ = Cr+3 + 5/2 H2O + 3/2 O2)
                break;
            }
        	case 60:
            {
                double gH2CrO4 = [[endmembers objectAtIndex:13] getGibbsFreeEnergyFromT:t andP:p];
                double gCr2O7 = [[endmembers objectAtIndex:77] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCr2O7+2.0*gH+gH2O-2.0*gH2CrO4)/R/t);
                s[60] = K*r[12]*r[12]/s[0]/s[0]/XH2O; // X77  Cr2O7-2  s60 (2H2CrO4 = Cr2O7-2 + 2H+ + H2O)
                break;
            }
        	case 61:
            {
                double gH2CrO4 = [[endmembers objectAtIndex:13] getGibbsFreeEnergyFromT:t andP:p];
                double gCrO4 = [[endmembers objectAtIndex:78] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCrO4+2.0*gH-gH2CrO4)/R/t);
                s[61] = K*r[12]/s[0]/s[0]; // X78  CrO4-2 s61 (H2CrO4 = CrO4-2 + 2H+)
                break;
            }
        	case 62:
        	{
                double gH2CrO4 = [[endmembers objectAtIndex:13] getGibbsFreeEnergyFromT:t andP:p];
        		double gHCrO4 = [[endmembers objectAtIndex:79] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gHCrO4+gH-gH2CrO4)/R/t);
                s[62] = K*r[12]/s[0]; // X79  HCrO4- s62 (H2CrO4 = HCrO4- + H+)
        		break;
    		}
    		case 63:
            {
                double gMnO2H2 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
                double gMn     = [[endmembers objectAtIndex:80] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMn+2.0*gOH-gMnO2H2)/R/t);
                s[63] = K*r[13]/s[1]/s[1]; // X80  Mn+2 s63 (Mn(OH)2 = Mn+2 + 2 OH-)
                break;
            }
    		case 64:
            {
                double gMnO2H2 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gMnCl = [[endmembers objectAtIndex:81] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMnCl+2.0*gH2O-gMnO2H2-gHCl-gH)/R/t);
                s[64] = K*r[9]*r[13]*s[0]/XH2O/XH2O; // X81  MnCl+ s64 (Mn(OH)2 + HCl + H+ = MnCl+ + 2 H2O)
                break;
            }
    		case 65:
            {
                double gMnO2H2 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
                double gMnO4 = [[endmembers objectAtIndex:82] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMnO4+gH+gH2O/2.0-gMnO2H2-5.0*gO2/4.0)/R/t);
                s[65] = K*r[13]*pow(r[1],5.0/4.0)/s[0]/sqrt(XH2O); // X82  MnO4- s65 (Mn(OH)2 + 5/4 O2 = MnO4- + H+ + 1/2 H2O)
                break;
            }
    		case 66:
            {
                double gMnO2H2 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
                double gMnO4 = [[endmembers objectAtIndex:83] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMnO4+2.0*gH-gMnO2H2-gO2)/R/t);
                s[66] = K*r[13]*r[1]/s[0]/s[0]; // X83  MnO4-2   s66 (Mn(OH)2 + O2 = MnO4-2 + 2 H+)
                break;
            }
    		case 67:
    		{
                double gMnO2H2 = [[endmembers objectAtIndex:14] getGibbsFreeEnergyFromT:t andP:p];
    			double gSO2   = [[endmembers objectAtIndex: 9] getGibbsFreeEnergyFromT:t andP:p];
                double gMnSO4 = [[endmembers objectAtIndex:84] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gMnSO4+gH2O-gMnO2H2-gSO2-gO2/2.0)/R/t);
                s[67] = K*r[8]*r[13]*sqrt(r[1])/XH2O; // X84  MnSO4 s67 (Mn(OH)2 + SO2 + 1/2 O2 = MnSO4 + H2O)
        		break;
    		}
    		case 68:
            {
                double gFeO2H2 = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
                double gFe     = [[endmembers objectAtIndex:85] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gFe+gH2O-gFeO2H2-2.0*gH)/R/t);
                s[68] = K*r[14]*s[0]*s[0]/XH2O/XH2O; // X85  Fe2+  s68 (Fe(OH)2 + 2 H+ = Fe2+ + 2 H2O)
                break;
            }
    		case 69:
            {
                double gFeO2H2 = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
                double gFe     = [[endmembers objectAtIndex:86] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gFe+5.0*gH2O/2.0-gFeO2H2-3.0*gH-gO2/2.0)/R/t);
                s[69] = K*r[14]*s[0]*s[0]*s[0]*pow(r[1],1.0/4.0)/pow(XH2O,5.0/2.0); // X86  Fe3+  s69 (Fe(OH)2 + 3 H+ + 1/4 O2 = Fe3+ + 5/2 H2O)
                break;
            }
    		case 70:
            {
                double gFeO2H2 = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gFeCl = [[endmembers objectAtIndex:87] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gFeCl+2.0*gH2O-gFeO2H2-gHCl-gH)/R/t);
                s[70] = K*r[9]*r[14]*s[0]/XH2O/XH2O; // X87  FeCl+    s70 (Fe(OH)2 + HCl + H+ = FeCl+ + 2 H2O)
                break;
            }
    		case 71:
            {
                double gFeO2H2 = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
                double gFeCl = [[endmembers objectAtIndex:88] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gFeCl+5.0*gH2O/2.0-gFeO2H2-gHCl-2.0*gH-gO2/4.0)/R/t);
                s[71] = K*r[9]*r[14]*s[0]*s[0]/pow(r[1],1.0/4.0)/pow(XH2O,5.0/2.0); // X88  FeCl+2 s71 (Fe(OH)2 + HCl + 2 H+ + 1/4 O2 = FeCl+2 + 5/2 H2O)
                break;
            }
    		case 72:
    		{
                double gFeO2H2 = [[endmembers objectAtIndex:15] getGibbsFreeEnergyFromT:t andP:p];
                double gHCl  = [[endmembers objectAtIndex:10] getGibbsFreeEnergyFromT:t andP:p];
    			double gFeCl2 = [[endmembers objectAtIndex:89] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gFeCl2+2.0*gH2O-gFeO2H2-2.0*gHCl)/R/t);
                s[72] = K*r[14]*r[9]*r[9]/XH2O/XH2O; // X89  FeCl2 s72 (Fe(OH)2 + 2 HCl = FeCl2 + 2 H2O)
            	break;
        	}
        	case 73:
            {
                double gCoO2H2 = [[endmembers objectAtIndex:16] getGibbsFreeEnergyFromT:t andP:p];
                double gCo     = [[endmembers objectAtIndex:90] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCo+2.0*gH2O-gCoO2H2-2.0*gH)/R/t);
                s[73] = K*r[15]*s[0]*s[0]/XH2O/XH2O; // X90  Co+2 s73 (Co(OH)2 + 2 H+ = Co2+ + 2 H2O)
                break;
            }
        	case 74:
        	{
                double gCoO2H2 = [[endmembers objectAtIndex:16] getGibbsFreeEnergyFromT:t andP:p];
    			double gCo = [[endmembers objectAtIndex:91] getGibbsFreeEnergyFromT:t andP:p];
                K = exp(-(gCo+5.0*gH2O/2.0-gCoO2H2-3.0*gH-gO2/4.0)/R/t);
                s[74] = K*r[15]*s[0]*s[0]*s[0]*pow(r[1],1.0/4.0)/pow(XH2O,5.0/2.0); // X91  Co+3 s74 (Co(OH)2 + 3 H+ + 1/4 O2 = Co3+ + 5/2 H2O)
        		break;
        	}
        	default:
        		break;
        }
    }
}

- (void)initialGuessOrderingGiven:(double [NR])r and:(double [NS])s atT:(double)t andP:(double)p {
    BOOL finiteOrderParameters[NS];
    BOOL debug = self.debugV;
    BOOL alwaysUseSimplex = NO;
    if (NS == 0) return;
    memset(s, 0, sizeof(s[0])*NS);

    [self guessOrderingParameters:s from:r atT:t andP:p];
    if (debug) NSLog(@"   Charge balance condition: %g", [self chargeBalanceWith:s]);
    alwaysUseSimplex |= ![self fillxSpeciesWith:r andWith:s];
    [self DxSpeciesDsWith:r andWith:s];
    [self fillnSpeciesWith:s];
    [self DnSpeciesDsWith:s];
    [self D2nSpeciesDs2With:s];
    for (NSUInteger i=0; i<NS; i++) finiteOrderParameters[i] = (xSpecies[NA+i] != 0.0);

    if (alwaysUseSimplex) {
        double a[NA+4][NS+2], sMin[NS], sMax[NS], tolerance = pow(DBL_EPSILON, (double) (2.0/3.0));
        NSInteger iposv[NA+4], izrov[NS+2];
        NSInteger m1 = 0, m2 = 0, m3 = 0, m, n, icase;

        if(debug) {
            NSLog(@"Setting up for a call to the Simplex method:");
            for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] != 0.0) {
                NSLog(@"[%3.3ld] %@ %g ", i,
                      [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0], xSpecies[i]);
            }
        }

        /***********************************/
        /* Find minimum length Soln vector */
        /***********************************/

        /* objective fucntion row for finding minimum length solution vector that satisfies constraints */
        a[1][1] = 0.0; for (NSUInteger i=0; i<NS; i++) a[1][2+i] = -1.0;
        /* remainder of the tableau */
        for (NSUInteger i=0; i<NA; i++) if (-xSpecies[i]  < 0.0) {
            m1++; a[1+m1][1] = xSpecies[i];  for (NSUInteger j=0; j<NS; j++) a[1+m1][2+j] = dxSpecies[i][j];
            if (debug) {
                id object = [endmembers objectAtIndex:i];
                NSString *name;
                if ([[object class] instancesRespondToSelector:@selector(phaseName)]) name = [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                else name = [[@"Species-" stringByAppendingFormat:@"%ld", i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                NSLog(@"m1: xSpecies[%@] = %e", name, a[1+m1][1]);
            }
        }
        for (NSUInteger i=0; i<NA; i++) if (-xSpecies[i]  > 0.0) {
            m2++; a[1+m1+m2][1] = -xSpecies[i]; for (NSUInteger j=0; j<NS; j++) a[1+m1+m2][2+j] = -dxSpecies[i][j];
            if (debug) {
                id object = [endmembers objectAtIndex:i];
                NSString *name;
                if ([[object class] instancesRespondToSelector:@selector(phaseName)]) name = [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                else name = [[@"Species-" stringByAppendingFormat:@"%ld", i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                NSLog(@"m2: xSpecies[%@] = %e", name, a[1+m1+m2][1]);
            }
        }
        for (NSUInteger i=0; i<NA; i++) if (-xSpecies[i] == 0.0) {
            m3++; a[1+m1+m2+m3][1] = 0.0; for (NSUInteger j=0; j<NS; j++) a[1+m1+m2+m3][2+j] = -dxSpecies[i][j];
            if (debug) {
                id object = [endmembers objectAtIndex:i];
                NSString *name;
                if ([[object class] instancesRespondToSelector:@selector(phaseName)]) name = [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                else name = [[@"Species-" stringByAppendingFormat:@"%ld", i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                NSLog(@"m3: xSpecies[%@] = %e", name, a[1+m1+m2+m3][1]);
            }
        }
        [self DchargeBalanceDsWith:DcbDs];
        m3++; a[1+m1+m2+m3][1] = 0.0; for (NSUInteger j=0; j<NS; j++) a[1+m1+m2+m3][2+j] = DcbDs[j];

        m = NA+1; n = NS;

        if (debug) {
//            for (NSInteger i=1; i<=(m+1); i++) {
//                for (NSInteger j=1; j<=(n+1); j++) NSLog(@"%6.3lf ", a[i][j]);
//            }
            NSLog(@"Making call to the Simplex method: m1 = %ld, m2 = %ld, m3 = %ld, m = %ld, n = %ld", m1, m2, m3, m, n);
        }

        [self simplx:a m:m n:n m1:m1 m2:m2 m3:m3 icasePt:&icase izrov:izrov iposv:iposv];

        if (debug) {
            NSLog(@"Tableau output from the Simplex method - minimum solution:");
//            for (NSUInteger i=1; i<=(m+1); i++) {
//                for (NSUInteger j=1; j<=(n+1); j++) NSLog(@"%6.3lf ", a[i][j]);
//            }
            NSLog(@"icase = %ld", icase);

            NSLog(@"Solution:");
            for (NSUInteger j=1,i; j<=m; j++) if ((i = iposv[j]) <= n)
                NSLog(@"s[%@][iposv %3.3ld] = %g",
                      [[[endmembers objectAtIndex:NA+i-1] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                      i, a[j+1][1]);
            for (NSUInteger j=1,i; j<=n; j++) if ((i = izrov[j]) <= n)
                NSLog(@"s[%@][izrov %3.3ld] = %g",
                      [[[endmembers objectAtIndex:NA+i-1] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                      i, 0.0);
        }

        /* if (icase != 0) return; */

        for (NSUInteger j=1,i; j<=m; j++) if ((i = iposv[j]) <= n) sMin[i-1] = a[j+1][1];
        for (NSUInteger j=1,i; j<=n; j++) if ((i = izrov[j]) <= n) sMin[i-1] = 0.0;

        /***********************************/
        /* Find maximum length Soln vector */
        /***********************************/

        /* objective fucntion row for finding maximum length solution vector that satisfies constraints */
        a[1][1] = 0.0; for (NSUInteger i=0; i<NS; i++) a[1][2+i] = 1.0;
        /* remainder of the tableau */
        m1 = 0; m2 = 0; m3 = 0;
        for (NSUInteger i=0; i<NA; i++) if (-xSpecies[i]  < 0.0) {
            m1++; a[1+m1][1] = xSpecies[i];  for (NSUInteger j=0; j<NS; j++) a[1+m1][2+j] = dxSpecies[i][j];
        }
        for (NSUInteger i=0; i<NA; i++) if (-xSpecies[i]  > 0.0) {
            m2++; a[1+m1+m2][1] = -xSpecies[i]; for (NSUInteger j=0; j<NS; j++) a[1+m1+m2][2+j] = -dxSpecies[i][j];
        }
        for (NSUInteger i=0; i<NA; i++) if (-xSpecies[i] == 0.0) {
            m3++; a[1+m1+m2+m3][1] = 0.0; for (NSUInteger j=0; j<NS; j++) a[1+m1+m2+m3][2+j] = -dxSpecies[i][j];
        }
        [self DchargeBalanceDsWith:DcbDs];
        m3++; a[1+m1+m2+m3][1] = 0.0; for (NSUInteger j=0; j<NS; j++) a[1+m1+m2+m3][2+j] = DcbDs[j];

        m = NA+1; n = NS;

        if (debug) {
//            for (NSUInteger i=1; i<=(m+1); i++) {
//                for (NSUInteger j=1; j<=(n+1); j++) NSLog(@"%6.3lf ", a[i][j]);
//            }
            NSLog(@"Making call to the Simplex method: m1 = %ld, m2 = %ld, m3 = %ld, m = %ld, n = %ld\n", m1, m2, m3, m, n);
        }

        [self simplx:a m:m n:n m1:m1 m2:m2 m3:m3 icasePt:&icase izrov:izrov iposv:iposv];

        if (debug) {
            NSLog(@"Tableau output from the Simplex method - maximum method:");
//            for (NSUInteger i=1; i<=(m+1); i++) {
//                for (NSUInteger j=1; j<=(n+1); j++) NSLog(@"%6.3lf ", a[i][j]);
//            }
            NSLog(@"icase = %ld", icase);

            NSLog(@"Solution:");
            for (NSUInteger j=1,i; j<=m; j++) if ((i = iposv[j]) <= n)
                NSLog(@"s[%@][iposv %3.3ld] = %g",
                      [[[endmembers objectAtIndex:NA+i-1] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                      i, a[j+1][1]);
            for (NSUInteger j=1,i; j<=n; j++) if ((i = izrov[j]) <= n)
                NSLog(@"s[%@][izrov %3.3ld] = %g",
                      [[[endmembers objectAtIndex:NA+i-1] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                      i, 0.0);
        }

        /* if (icase != 0) return; */

        for (NSUInteger j=1,i; j<=m; j++) if ((i = iposv[j]) <= n) sMax[i-1] = a[j+1][1];
        for (NSUInteger j=1,i; j<=n; j++) if ((i = izrov[j]) <= n) sMax[i-1] = 0.0;

        /***************************************************************/
        /* Return average of two solutions - should always be feasible */
        /***************************************************************/

        for (NSUInteger i=0; i<NS; i++) {
            s[i] = (sMin[i] + sMax[i])/2.0;
            if (finiteOrderParameters[i] && (fabs(s[i]) < tolerance)) s[i] = sqrt(DBL_EPSILON);
            else if (fabs(s[i]) < tolerance) s[i] = 0.0;
        }
        double cations = s[0] + s[8] + s[13] + s[14] + 2.0*s[15] + s[16] + s[17] + 3.0*s[19] + s[48] + s[52] + s[53] + 2.0*s[54]
        + s[55] + 2.0*s[58] + 3.0*s[59] + 2.0*s[63] + s[64] + 2.0*s[68] + 3.0*s[69] + s[70] + 2.0*s[71]
        + 2.0*s[73] + 3.0*s[74];
        double anions = s[1] + 2.0*s[3] + s[4] + s[6] + s[9] + s[20] + s[21] + s[23] + 2.0*s[24] + 3.0*s[25] + s[26] + 2.0*s[27]
        + s[29] + 2.0*s[30] + 2.0*s[31] + 2.0*s[32] + 2.0*s[33] + 2.0*s[34] + 2.0*s[35] + 2.0*s[36] + 2.0*s[37]
        + 2.0*s[38] + 2.0*s[39] + 2.0*s[40] + 2.0*s[41] + 2.0*s[42] + s[43] + 2.0*s[44] + s[45] + s[46] + s[47]
        + s[50] + 2.0*s[60] + 2.0*s[61] + s[62] + s[65] + 2.0*s[66];

        if (cations > anions) {
            // excess cations
            s[1] += cations-anions;  // X18  OH-      s1
        } else if (anions > cations) {
            // excess anions
            s[0] += anions-cations;  // X17  H+       s0
        }

        if(![self fillxSpeciesWith:r andWith:s])
            if (debug) NSLog(@"Simplex method suceeded but failed to find feasible solution in initialGuessOrdering.");

        if (debug) {
            NSLog(@"Results of call to initialGuessOrdering (after simplex):");
            NSLog(@"   %@ %@ %@ %@",
                  [@"Species" stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                  [@"Mole frac" stringByPaddingToLength:13 withString:@" " startingAtIndex:0],
                  [@"r" stringByPaddingToLength:13 withString:@" " startingAtIndex:0],
                  [@"s" stringByPaddingToLength:13 withString:@" " startingAtIndex:0]);
            NSLog(@"   %@ %13.6g",
                  [[[endmembers objectAtIndex:0] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                  xSpecies[0]);
            for (NSUInteger i=0;  i<NR; i++) {
                id object = [endmembers objectAtIndex:i+1];
                if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                    NSLog(@"   %@ %13.6g %13.6g",
                          [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                          xSpecies[i+1], r[i]);
                } else {
                    NSLog(@"   %@ %13.6g %13.6g",
                          [[@"Species-" stringByAppendingFormat:@"%ld", i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                          xSpecies[i+1], r[i]);
                }

            }
            for (NSUInteger i=0;  i<NS; i++) {
                id object = [endmembers objectAtIndex:i+NA];
                if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                    NSLog(@"   %@ %13.6g %13.13s %13.6g",
                          [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                          xSpecies[i+NA], "", s[i]);
                } else {
                    NSLog(@"   %@ %13.6g %13.13s %13.6g",
                          [[@"Species-" stringByAppendingFormat:@"%ld", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                          xSpecies[i+NA], "", s[i]);
                }
            }
            NSLog(@"   sMin: Charge balance condition: %g", [self chargeBalanceWith:sMin]);
            NSLog(@"   sMax: Charge balance condition: %g", [self chargeBalanceWith:sMax]);
            NSLog(@"   sAve: Charge balance condition: %g", [self chargeBalanceWith:s]);
        }

        return;
    } /* end block on simplex method */

    /*
    double sCorr[NS];
    memset(sCorr, 0, sizeof(sCorr));
    for (NSUInteger i=0; i<NS; i++) {
        s[i] = sqrt(DBL_EPSILON);
        if (![self fillxSpeciesWith:r andWith:s]) s[i] = 0.0;
        else {
            sCorr[i] = 0.5; s[i] = 1.0;
            while (sCorr[i] > sqrt(DBL_EPSILON)) {
                if(![self fillxSpeciesWith:r andWith:s]) s[i] -= sCorr[i]; else s[i] += sCorr[i];
                sCorr[i] /= 2.0;
            }
        }
        sCorr[i] = s[i]/2.0;
        s[i]     = 0.0;
    }

    for (NSUInteger i=0; i<NS; i++) s[i] = sCorr[i];
    double factor = 1.0;
    while (factor > sqrt(DBL_EPSILON) && ![self fillxSpeciesWith:r andWith:s]) {
        factor /= 2.0;
        for (NSUInteger i=0; i<NS; i++) s[i] = sCorr[i]*factor;
    }
     */

    if (debug) {
        NSLog(@"Results of call to initialGuessOrdering (bypass simplex):");
        NSLog(@"   %@ %@ %@ %@",
              [@"Species" stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
              [@"Mole frac" stringByPaddingToLength:13 withString:@" " startingAtIndex:0],
              [@"r" stringByPaddingToLength:13 withString:@" " startingAtIndex:0],
              [@"s" stringByPaddingToLength:13 withString:@" " startingAtIndex:0]);
        NSLog(@"   %@ %13.6g",
              [[[endmembers objectAtIndex:0] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
              xSpecies[0]);
        for (NSUInteger i=0;  i<NR; i++) {
            id object = [endmembers objectAtIndex:i+1];
            if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                NSLog(@"   %@ %13.6g %13.6g",
                      [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                      xSpecies[i+1], r[i]);
            } else {
                NSLog(@"   %@ %13.6g %13.6g",
                      [[@"Species-" stringByAppendingFormat:@"%ld", i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                      xSpecies[i+1], r[i]);
            }

        }
        for (NSUInteger i=0;  i<NS; i++) {
            id object = [endmembers objectAtIndex:i+NA];
            if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                NSLog(@"   %@ %13.6g %13.13s %13.6g",
                      [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                      xSpecies[i+NA], "", s[i]);
            } else {
                NSLog(@"   %@ %13.6g %13.13s %13.6g",
                      [[@"Species-" stringByAppendingFormat:@"%ld", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                      xSpecies[i+NA], "", s[i]);
            }
        }
        NSLog(@"   Charge balance condition: %g", [self chargeBalanceWith:s]);
    }
}

- (void)speciation:(int)mask
                 t:(double)t
                 p:(double)p
                 r:(double [NR])r
                 s:(double [NS])s           // s[NS]                BINARY MASK: 0000000001
                dr:(double [NS][NR])dr      // ds[NS]/dr[NR]        BINARY MASK: 0000000010
                dt:(double [NS])dt          // ds[NS]/dt            BINARY MASK: 0000000100
                dp:(double [NS])dp          // ds[NS]/dp            BINARY MASK: 0000001000
               dr2:(double [NS][NR][NR])dr2 // d2s[NS]/dr[NR]dr[NR] BINARY MASK: 0000010000
               drt:(double [NS][NR])drt     // d2s[NS]/dr[NR]dt     BINARY MASK: 0000100000
               drp:(double [NS][NR])drp     // d2s[NS]/dr[NR]dp     BINARY MASK: 0001000000
               dt2:(double [NS])dt2         // d2s[NS]/dt2          BINARY MASK: 0010000000
               dtp:(double [NS])dtp         // d2s[NS]/dtp          BINARY MASK: 0100000000
               dp2:(double [NS])dp2         // d2s[NS]/dp2          BINARY MASK: 1000000000
{
    int i, j, k, l, kk, ll, iter = 0;
    BOOL debugV = self.debugV;
    BOOL debugS = self.debugS;

    /* look-up or compute the current ordering state */
    if ((t != tOld) || (p != pOld) ||
        (r[ 0] != rOld[ 0]) || (r[ 1] != rOld[ 1]) || (r[ 2] != rOld[ 2]) || (r[ 3] != rOld[ 3]) || (r[ 4] != rOld[ 4]) ||
        (r[ 5] != rOld[ 5]) || (r[ 6] != rOld[ 6]) || (r[ 7] != rOld[ 7]) || (r[ 8] != rOld[ 8]) || (r[ 9] != rOld[ 9]) ||
        (r[10] != rOld[10]) || (r[11] != rOld[11]) || (r[12] != rOld[12]) || (r[13] != rOld[13]) || (r[14] != rOld[14]) ||
        (r[15] != rOld[15]) ) {

        double dgds[NS], sNew[NS], charge[NS], rTotal = 0.0;
        BOOL doLoop = YES;

        for (i=0; i<NS; i++) sOld[i] = 2.0;
        for (i=0; i<NR; i++) rTotal += r[i];

        [self initialGuessOrderingGiven:r and:sNew atT:t andP:p];
        if ((r[ 0] != rOld[ 0]) || (r[ 1] != rOld[ 1]) || (r[ 2] != rOld[ 2]) || (r[ 3] != rOld[ 3]) || (r[ 4] != rOld[ 4]) ||
            (r[ 5] != rOld[ 5]) || (r[ 6] != rOld[ 6]) || (r[ 7] != rOld[ 7]) || (r[ 8] != rOld[ 8]) || (r[ 9] != rOld[ 9]) ||
            (r[10] != rOld[10]) || (r[11] != rOld[11]) || (r[12] != rOld[12]) || (r[13] != rOld[13]) || (r[14] != rOld[14]) ||
            (r[15] != rOld[15])) {
            // azero is a function ONLY of bulk composition, NOT speciation and NOT t and p
            [self azeroWithR:r];
        }

        while (doLoop && (iter < MAX_ITER)) {
            double s[NS], deltaS[NS], lambda, d2gds2[NS][NS];

            for (i=0; i<NS; i++) s[i] = sNew[i];

            [self fillxSpeciesWith:r andWith:s];
            [self DxSpeciesDsWith:r andWith:s];
            [self DxSpeciesDrWith:r andWith:s];
            [self fillnSpeciesWith:s];
            [self DnSpeciesDsWith:s];
            [self D2nSpeciesDs2With:s];

            [self fillDGDSin:dgds withT:t andWithP:p];
            if (debugV) {
                for (NSUInteger i=0; i<NS; i++) {
                    if (s[i] > 0.0) {
                        id object = [endmembers objectAtIndex:i+NA];
                        if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                            NSLog(@"dgds[%@, s:%13.6e] = %g",
                                  [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                                  s[i], dgds[i]);
                        } else {
                            NSLog(@"dgds[%@, s:%13.6e] = %g",
                                  [[@"Species-" stringByAppendingFormat:@"%ld", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                                  s[i], dgds[i]);
                        }
                    }
                }
            }

            [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
            if (debugV) {
                for (NSUInteger i=0; i<NS; i++) {
                    if (s[i] > 0.0) {
                        for (NSUInteger j=i; j<NS; j++) {
                            if (s[j] > 0.0) {
                                id objectForI = [endmembers objectAtIndex:i+NA];
                                NSString *nameForI = @"";
                                if ([[objectForI class] instancesRespondToSelector:@selector(phaseName)])
                                    nameForI = [[objectForI phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                                else
                                    nameForI = [[@"Species-" stringByAppendingFormat:@"%ld", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];

                                id objectForJ = [endmembers objectAtIndex:j+NA];
                                NSString *nameForJ = @"";
                                if ([[objectForJ class] instancesRespondToSelector:@selector(phaseName)])
                                    nameForJ = [[objectForJ phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];
                                else
                                    nameForJ = [[@"Species-" stringByAppendingFormat:@"%ld", j+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0];

                                NSLog(@"d2gds2[%@, %@] = %g", nameForI, nameForJ, d2gds2[i][j]);
                            }
                        }
                    }
                }
            }

            for (j=0, k=0; j<NS; j++) {
                charge[j] = 0.0;
                sOld[j] = s[j];
                if (xSpecies[j+NA] > 0.0) {
                    incOfxSpecies[j] = YES;
                    id object = [endmembers objectAtIndex:j+NA];
                    if ([object isKindOfClass:[HKFspeciesProperties class]]) charge[j] = [(HKFspeciesProperties *)object charge];
                    k++;
                } else incOfxSpecies[j] = NO;
            }

            dgglse_m     = k-1;    // matrix order - rows, there is no row for the hydrogen ion
            dgglse_n     = k;      // matrix order - columns, but there is a column for the hydrogen ion
            dgglse_p     = 1;      // number of constraints for charge balance
            dgglse_lda   = k-1;    // row-dimension of matrix "A"
            dgglse_ldb   = 1;      // row-dimension of matrix "B"
            dgglse_lwork = 2*NS+1; // dimension of work array

            // LAPACK normal matrx storage for dgglse
            dgglse_d[0] = 0.0;
            [self DchargeBalanceDsWith:s];
            for (i=0, k=0, kk=0; i<NS; i++) {
                if (i != indexOfHydrogenIon) {
                    if (xSpecies[i+NA] > 0.0) {
                        for (j=0, l=0, ll=0; j<NS; j++) if (xSpecies[j+NA] > 0.0) {
                            if (charge[i] != 0.0) {
                                dgglse_a[k+dgglse_lda*l] = d2gds2[i][j] - charge[i]*d2gds2[indexOfHydrogenIon][j];
                            } else {
                                dgglse_a[k+dgglse_lda*l] = d2gds2[i][j];
                            }
                            l++;
                            if (j != indexOfHydrogenIon) {
                                if (i <= j) { // indexing of ap(i,j) = i + (j-1)*j/2 1 <=i<=j
                                    dspsvx_ap[k+(ll+1)*ll/2]  = d2gds2[i][j];  // lower right symmetric matrix
                                    dspsvx_ap[k+(ll+1)*ll/2] += (charge[i] != 0.0) ? -charge[i]*d2gds2[indexOfHydrogenIon][j] :0.0;
                                    dspsvx_ap[k+(ll+1)*ll/2] += (charge[j] != 0.0) ? -charge[j]*d2gds2[indexOfHydrogenIon][i] :0.0;
                                    dspsvx_ap[k+(ll+1)*ll/2] += ((charge[i] != 0.0) && (charge[j] != 0.0)) ?
                                                                charge[j]*charge[i]*d2gds2[indexOfHydrogenIon][indexOfHydrogenIon] : 0.0;
                                }
                                ll++;
                            }
                        }
                        // if anion,  then dgds[i] + |charge|*dgds[H+]
                        // if cation, then dgds[i] - |charge|*dgds[H+]
                        if      (charge[i] != 0.0) dgglse_c[k] = -(dgds[i] - charge[i]*dgds[indexOfHydrogenIon]);
                        else                       dgglse_c[k] = -dgds[i];
                        dspsvx_b[k]  = 0.0;  // no rhs vector is referenced
                        dgglse_b[kk] = DcbDs[i];
                        k++;
                        kk++;
                    }
                } else {
                    // Hydrogen ion is in the charge balance constraint, but NOT in the Jacobian rows
                    dgglse_b[kk] = DcbDs[indexOfHydrogenIon];
                    kk++;
                }
            }

            if (debugV) {
                NSLog(@"On entry to dgglse (d = %e):", dgglse_d[0]);
                for (NSUInteger i=0; i<dgglse_n; i++) NSLog(@"b[%lu] = %e", i, dgglse_b[i]);
            }

            dgglse_(&dgglse_m, &dgglse_n, &dgglse_p, dgglse_a, &dgglse_lda, dgglse_b, &dgglse_ldb, dgglse_c, dgglse_d,
                    dgglse_x, dgglse_work, &dgglse_lwork, &dgglse_info);

            if ((dgglse_info != 0) && (debugV || debugS)) {
                NSLog(@"DEW:order() %d argument to LAPACK dgglse is incorrect.", -dgglse_info);
                NSAssert(NO, @"Error returned by LAPACK dggglse in DEW:order().");
            }
            if (debugV) {
                double SSresiduals = 0.0;
                for (i=0; i<dgglse_n-dgglse_p; i++) SSresiduals += pow(dgglse_c[i], 2.0);
                NSLog(@"DEW:order() Iter:%3.3d Resiudal SS: %g Remainder SS: %g", iter, SSresiduals, pow(dgglse_c[dgglse_n-dgglse_p], 2.0));
            }

            if (debugV) NSLog(@"<>Ordering parameters before lambda correction ...");
            double cbOld=0.0, cbNew=0.0, cbDelta=0.0;
            for (i=0, k=0; i<NS; i++) {
                if (incOfxSpecies[i]) {
                    s[i] += dgglse_x[k];
                    if (debugV) {
                        id object = [endmembers objectAtIndex:i+NA];
                        if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                            NSLog(@"...s[%@] = %13.6e sOld[] = %13.6e delta = %13.6e",
                                  [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                                  s[i], sOld[i], dgglse_x[k]);
                        } else {
                            NSLog(@"...s[%@] = %13.6e sOld[] = %13.6e delta = %13.6e",
                                  [[@"Species-" stringByAppendingFormat:@"%d", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                                  s[i], sOld[i], dgglse_x[k]);
                        }
                    }
                    k++;
                } else s[i] = 0.0;
                deltaS[i] = s[i] - sOld[i];
                if (debugV) {
                    cbOld   += sOld[i]*charge[i];
                    cbNew   += s[i]*charge[i];
                    cbDelta += deltaS[i]*charge[i];
                }
            }
            if (debugV) NSLog(@"...charge balnace old: %e, new: %e, delta %e ", cbOld, cbNew, cbDelta);

            lambda = 1.0;
            BOOL guessIsOkay = NO;
            while (!guessIsOkay) {
                for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*deltaS[i];
                guessIsOkay = [self fillxSpeciesWith:r andWith:s];
                if (!guessIsOkay) lambda = lambda/2.0;
                if (lambda <= DBL_EPSILON) {
                    for (i=0; i<NS; i++) s[i] = sOld[i];
                    guessIsOkay = YES;
                }
            }
            if (debugV) NSLog(@"...lambda = %13.6e", lambda);

            if (debugV) NSLog(@"<>Ordering parameters after lambda correction ...");
            for (i=0; i<NS; i++) {
                sNew[i] = s[i];
                if (debugV && incOfxSpecies[i]) {
                    id object = [endmembers objectAtIndex:i+NA];
                    if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                        NSLog(@"...s[%@] = %13.6e sOld[] = %13.6e delta = %13.6e",
                            [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                            sNew[i], sOld[i], sNew[i]-sOld[i]);
                    } else {
                        NSLog(@"...s[%@] = %13.6e sOld[] = %13.6e delta = %13.6e",
                            [[@"Species-" stringByAppendingFormat:@"%d", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                            sNew[i], sOld[i], sNew[i]-sOld[i]);
                    }
                }
            }

            double cb = [self chargeBalanceWith:s];
            if (debugV) {
                if (cb < -DBL_EPSILON*100.0)
                    NSLog(@"... anion excess of %g with H+ concentration at %g", -cb, xSpecies[mIndHplus]);
                else if (cb > DBL_EPSILON*100.0)
                    NSLog(@"... anion excess of %g with OH- concentration at %g",  cb, xSpecies[mIndOHminus]);
                else
                    NSLog(@"... system is in charge balance (%g)", cb);
            }

            doLoop = NO;
            for (i=0; i<NS; i++) doLoop |= (fabs(sNew[i]-sOld[i]) > 10.0*DBL_EPSILON);

            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

        if (debugS || (iter >= 200)) {
            if (iter >= MAX_ITER) NSLog(@"ERROR in DEWFluid Class (function ORDER). Failed to converge!\n  Iteration limit (%4d) exceeded.", iter);
            else NSLog(@"DEBUG: DEWFluid Class (function ORDER) Iter = %4d.", iter);
            for (i=0; i<NS; i++) {
                if (fabs(dgds[i]) > sqrt(DBL_EPSILON) && fabs(sOld[i]) > DBL_EPSILON) {
                    id object = [endmembers objectAtIndex:i+NA];
                    if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                        NSLog(@"  s[%@] = %13.6g  dgds[] = %13.6g  dcbds[] = %13.6g",
                              [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                              sOld[i], dgds[i], DcbDs[i]);
                    } else {
                        NSLog(@"  s[%@] = %13.6g  dgds[] = %13.6g  dcbds[] = %13.6g",
                              [[@"Species-" stringByAppendingFormat:@"%d", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                              sOld[i], dgds[i], DcbDs[i]);
                    }
                }
            }
        }

        [self adjustOrderingParameters:sOld from:r atT:t andP:p];
        if (debugV) {
        	NSLog(@"DEWFluid Class (after adjustment):");
        	for (i=0; i<NS; i++) {
                if (fabs(sOld[i]) > 0.0) {
                    id object = [endmembers objectAtIndex:i+NA];
                    if ([[object class] instancesRespondToSelector:@selector(phaseName)]) {
                        NSLog(@"  s[%@] = %13.6g",
                              [[object phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0], sOld[i]);
                    } else {
                        NSLog(@"  s[%@] = %13.6g",
                              [[@"Species-" stringByAppendingFormat:@"%d", i+NA] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                              sOld[i]);
                    }
                }
            }
        }

        //
        // Compute factorization at the solution for use in order parameter derivatives
        //
        dspsvx_fact = 'N';      // factor the system; 'F' denotes matrix is factored already and stored in AFP
        dspsvx_uplo = 'U';      // stored in upper triangle of A; 'L' would be lower triangle
        dspsvx_ldb  = NS;       // leading dimension of B
        dspsvx_ldx  = NS;       // leading dimension of X
        dspsvx_n    = dgglse_m; // matrix order define above
        dspsvx_nrhs =  0;       // number of right-hand-sides (just factor the matrix)

        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs,
                dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb,
                dspsvx_x, &dspsvx_ldx,
                &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);

        if ((dspsvx_info != 0) && debugS) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            doLoop = NO;
            // NSAssert(NO, @"Error returned by LAPACK dppsvx in DEW:order().");
        }

        if (debugS) NSLog(@"... Condition number = %g", dspsvx_rcond);

        // set up flags for computation with factorized matrixes
        dspsvx_fact = 'F';
        dspsvx_nrhs = 1;
        dspsvx_ldb  = dspsvx_n;
        dspsvx_ldx  = dspsvx_n;

    }

    // NSAssert(NO, @"ERROR");

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
        double d2gdrds[NR][NS], rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self DxSpeciesDrWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];

        dspsvx_nrhs = 0;
        memset(dr, 0, sizeof(dr[0][0])*NR*NS);
        for (NSUInteger j=0, k=0; j<NR; j++) if (r[j] != 0.0) {
            dspsvx_nrhs++;
            NSString *nameOfJ = (debugV) ? [[[endmembers objectAtIndex:j+1] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0] : @"";
            for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
                dspsvx_b[k] = - d2gdrds[j][i];
                id object = [endmembers objectAtIndex:i+NA];
                if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                    double z = [(HKFspeciesProperties *)object charge];
                    dspsvx_b[k] += z*d2gdrds[j][indexOfHydrogenIon]; // this is -(-z)
                }
                NSString *nameOfI = (debugV) ? [[object phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0] : @"";
                if (debugV) NSLog(@"r[%@] s[%@] dspsvx_b[%3.3lu] = %g", nameOfJ, nameOfI, k, dspsvx_b[k]);
                k++;
            }
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger j=0, k=0; j<NR; j++) if (r[j] != 0.0) {
            NSString *nameOfJ = (debugV) ? [[[endmembers objectAtIndex:j+1] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0] : @"";
            for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
                dr[i][j] = dspsvx_x[k++];
                id object = [endmembers objectAtIndex:i+NA];
                if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                    double z = [(HKFspeciesProperties *)object charge];
                    dr[indexOfHydrogenIon][j] += -z*dr[i][j];
                }
                NSString *nameOfI = (debugV) ? [[object phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0] : @"";
                if (debugV) NSLog(@"r[%@] s[%@] dspsvx_x[%3.3lu] = %13.6g H+ = %g", nameOfJ, nameOfI, k, dr[i][j], dr[indexOfHydrogenIon][j]);
            }
        }

//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                dr[i][j] = 0.0;
//                for (k=0; k<NS; k++) dr[i][j] += (s[i] > 0.0) ? - invd2gds2[i][k]*d2gdrds[j][k] : 0.0;
//            }
//        }
    }
    if (mask & THIRD) {   /* compute ds/dt:  */
        double d2gdsdt[NS], rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];

        dspsvx_nrhs =  1;
        memset(dt, 0, sizeof(dt[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) {
            if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
                dspsvx_b[k] = -d2gdsdt[i];
                id object = [endmembers objectAtIndex:i+NA];
                if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                    double z = [(HKFspeciesProperties *)object charge];
                    dspsvx_b[k] += z*d2gdsdt[indexOfHydrogenIon]; // this is -(-z)
                }
                k++;
            }
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dt[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dt[indexOfHydrogenIon] += -z*dt[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dt[i] = 0.0;
//            for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
//        }
    }
    if (mask & FOURTH) {   /* compute ds/dp:  */
        double d2gdsdp[NS], rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];

        dspsvx_nrhs =  1;
        memset(dp, 0, sizeof(dp[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdp[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdp[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        //NSLog(@"%@", [NSThread callStackSymbols]);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dp[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dp[indexOfHydrogenIon] += -z*dp[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dp[i] = 0.0;
//            for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
//        }
    }
    if (mask & FIFTH) {   /* compute d2s/dr2 */
        double d2gdrds[NR][NS], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS], dsdr[NS][NR], temp[NS];
        int m, n, ii, ll;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self DxSpeciesDrWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD3GDR2DSin:d3gdr2ds withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];

        /* compute dsdr matrix */
        dspsvx_nrhs =  NR;
        memset(dsdr, 0, sizeof(dsdr[0][0])*NR*NS);
        for (NSUInteger j=0, k=0; j<NR; j++) for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdrds[j][i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdrds[j][indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger j=0, k=0; j<NR; j++) for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdr[i][j] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dt[indexOfHydrogenIon] += -z*dsdr[i][j];
            }
        }

//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                dsdr[i][j] = 0.0;
//                for (k=0; k<NS; k++) dsdr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
//            }
//        }

        /* compute dsdr2 cube */
        dspsvx_nrhs =  1;
        memset(dr2, 0, sizeof(dr2[0][0][0])*NS*NR*NR);
        for (j=0; j<NR; j++) {
            for (k=0; k<NR; k++) {
                for (i=0, ii=0; i<NS; i++) if(incOfxSpecies[i]) {
                    for (l=0, ll=0; l<NS; l++) if (incOfxSpecies[l]) {
                        temp[ll] = d3gdr2ds[j][k][l];
                        for (m=0; m<NS; m++) if (incOfxSpecies[m]) {
                            temp[ll] += d3gdrds2[j][l][m]*dsdr[m][k] + d3gdrds2[k][l][m]*dsdr[m][j];
                            for (n=0; n<NS; n++) if (incOfxSpecies[n])
                                temp[ll] += [self fillD3GDS3atI:l andL:m andM:n withT:t andWithP:p]*dsdr[m][j]*dsdr[n][k];
                        }
                        ll++;
                    }
                    for (ll=0; ll<dspsvx_n; ll++) dspsvx_b[ll] = -temp[ll];
                    NSLog(@"Call to DSPSVX for d2s/dr2 has not been properly configured.");
                    if (NO) {
                        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                            dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                            dspsvx_work, dspsvx_iwork, &dspsvx_info);
                        if (dspsvx_info != 0) {
                            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
                            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
                            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
                            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
                        }
                        dr2[i][j][k] = dspsvx_x[ii];
                    } else {
                        dr2[i][j][k] = 0.0;
                    }
                    ii++;
                }
            }
        }

//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                for (k=0; k<NR; k++) {
//                    for (l=0; l<NS; l++) {
//                        temp[l] = d3gdr2ds[j][k][l];
//                        for (m=0; m<NS; m++) {
//                            temp[l] += d3gdrds2[j][l][m]*dsdr[m][k] + d3gdrds2[k][l][m]*dsdr[m][j];
//                            for (n=0; n<NS; n++)
//                                temp[l] += d3gds3[l][m][n]*dsdr[m][j]*dsdr[n][k];
//                        }
//                    }
//                    dr2[i][j][k] = 0.0;
//                    for (l=0; l<NS; l++) dr2[i][j][k] += - invd2gds2[i][l]*temp[l];
//                }
//            }
//        }

    }
    if (mask & SIXTH) {   /* compute d2s/drt */
        double d2gdrds[NR][NS], d2gdsdt[NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS], d3gds2dt[NS][NS];
        double dsdr[NS][NR], dsdt[NS], temp[NS];
        int k, l, m, ii, kk;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self DxSpeciesDrWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDTin:d3gdrdsdt withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];

        /* compute dsdr matrix */
        dspsvx_nrhs =  NR;
        memset(dsdr, 0, sizeof(dsdr[0][0])*NR*NS);
        for (NSUInteger j=0, k=0; j<NR; j++) for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdrds[j][i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdrds[j][indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger j=0, k=0; j<NR; j++) for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdr[i][j] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdr[indexOfHydrogenIon][j] += -z*dsdr[i][j];
            }
        }

//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                dsdr[i][j] = 0.0;
//                for (k=0; k<NS; k++) dsdr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
//            }
//        }

        /* compute dsdt vector */
        dspsvx_nrhs =  1;
        memset(dsdt, 0, sizeof(dsdt[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdt[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdt[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdt[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdt[indexOfHydrogenIon] += -z*dsdt[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dsdt[i] = 0.0;
//            for (j=0; j<NS; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
//        }

        /* compute dsdrdt matrix */
        dspsvx_nrhs =  1;
        memset(drt, 0, sizeof(drt[0][0])*NR*NS);
        for (j=0; j<NR; j++) {
            for (i=0, ii=0; i<NS; i++) if (incOfxSpecies[i]) {
                for (k=0, kk=0; k<NS; k++) if (incOfxSpecies[k]) {
                    temp[kk] = d3gdrdsdt[j][k];
                    for (l=0; l<NS; l++) {
                        temp[kk] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
                        for (m=0; m<NS; m++) temp[kk] += [self fillD3GDS3atI:k andL:l andM:m withT:t andWithP:p]*dsdr[l][j]*dsdt[m];
                    }
                    kk++;
                }
                for (kk=0; kk<dspsvx_n; kk++) dspsvx_b[kk] = -temp[kk];
                if (j == 0) NSLog(@"Call to DSPSVX for d2s/drt has not been properly configured.");
                if (NO) {
                    dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                        dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                        dspsvx_work, dspsvx_iwork, &dspsvx_info);
                    if (dspsvx_info != 0) {
                        if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
                        else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
                        else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
                        // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
                    }
                    drt[i][j] = dspsvx_x[ii];
                } else {
                   drt[i][j] = 0.0;
                }
                ii++;
            }
        }
//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                for (k=0; k<NS; k++) {
//                    temp[k] = d3gdrdsdt[j][k];
//                    for (l=0; l<NS; l++) {
//                        temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
//                        for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
//                    }
//                }
//                drt[i][j] = 0.0;
//                for (k=0; k<NS; k++) drt[i][j] += - invd2gds2[i][k]*temp[k];
//            }
//        }

    }
    if (mask & SEVENTH) {   /* compute d2s/drp */
        double d2gdrds[NR][NS], d2gdsdp[NS], d3gdrds2[NR][NS][NS], d3gdrdsdp[NR][NS], d3gds2dp[NS][NS];
        double dsdr[NS][NR], dsdp[NS], temp[NS];
        int k, l, m, ii, kk;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self DxSpeciesDrWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDPin:d3gdrdsdp withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];

        /* compute dsdr matrix */
        dspsvx_nrhs =  NR;
        memset(dsdr, 0, sizeof(dsdr[0][0])*NR*NS);
        for (NSUInteger j=0, k=0; j<NR; j++) for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdrds[j][i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdrds[j][indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger j=0, k=0; j<NR; j++) for (NSUInteger i=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdr[i][j] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdr[indexOfHydrogenIon][j] += -z*dsdr[i][j];
            }
        }

//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                dsdr[i][j] = 0.0;
//                for (k=0; k<NS; k++) dsdr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
//            }
//        }

        /* compute dsdp vector */
        dspsvx_nrhs =  1;
        memset(dsdp, 0, sizeof(dsdp[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdp[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdp[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdp[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdp[indexOfHydrogenIon] += -z*dsdp[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dsdp[i] = 0.0;
//            for (j=0; j<NS; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
//        }

        /* compute dsdrdp matrix */
        dspsvx_nrhs =  1;
        memset(drt, 0, sizeof(drp[0][0])*NR*NS);
        for (j=0; j<NR; j++) {
            for (i=0, ii=0; i<NS; i++) if (incOfxSpecies[i]) {
                for (k=0, kk=0; k<NS; k++) if (incOfxSpecies[k]) {
                    temp[kk] = d3gdrdsdp[j][k];
                    for (l=0; l<NS; l++) {
                        temp[kk] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
                        for (m=0; m<NS; m++) temp[kk] += [self fillD3GDS3atI:k andL:l andM:m withT:t andWithP:p]*dsdr[l][j]*dsdp[m];
                    }
                    kk++;
                }
                for (kk=0; kk<dspsvx_n; kk++) dspsvx_b[kk] = -temp[kk];
                if (j == 0) NSLog(@"Call to DSPSVX for d2s/drp has not been properly configured.");
                if (NO) {
                    dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                        dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                        dspsvx_work, dspsvx_iwork, &dspsvx_info);
                    if (dspsvx_info != 0) {
                        if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
                        else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
                        else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
                        // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
                    }
                    drp[i][j] = dspsvx_x[ii];
                } else {
                    drp[i][j] = 0.0;
                }
                ii++;
            }
        }
//        for (i=0; i<NS; i++) {
//            for (j=0; j<NR; j++) {
//                for (k=0; k<NS; k++) {
//                    temp[k] = d3gdrdsdp[j][k];
//                    for (l=0; l<NS; l++) {
//                        temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
//                        for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
//                    }
//                }
//                drp[i][j] = 0.0;
//                for (k=0; k<NS; k++) drp[i][j] += - invd2gds2[i][k]*temp[k];
//            }
//        }

    }
    if (mask & EIGHTH ) {   /* compute d2s/dt2 */
        double d2gdsdt[NS], d3gds2dt[NS][NS], d3gdsdt2[NS], dsdt[NS], temp[NS];
        int k, l, ii, jj;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDSDT2in:d3gdsdt2 withT:t andWithP:p];

        /* compute dsdt vector */
        dspsvx_nrhs =  1;
        memset(dsdt, 0, sizeof(dsdt[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdt[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdt[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdt[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdt[indexOfHydrogenIon] += -z*dsdt[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dsdt[i] = 0.0;
//            for (j=0; j<NS; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
//        }

        /* compute dsdt2 vector */
        dspsvx_nrhs =  1;
        memset(dt2, 0, sizeof(dt2[0])*NS);
        for (i=0, ii=0; i<NS; i++) if (incOfxSpecies[i]) {
            for (j=0, jj=0; j<NS; j++) if (incOfxSpecies[j]) {
                temp[jj] = d3gdsdt2[j];
                for (k=0; k<NS; k++) {
                    temp[jj] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                    for (l=0; l<NS; l++) temp[jj] += [self fillD3GDS3atI:j andL:k andM:l withT:t andWithP:p]*dsdt[k]*dsdt[l];
                }
                jj++;
            }
            for (jj=0; jj<dspsvx_n; jj++) dspsvx_b[jj] = -temp[jj];
            if (i == 0) NSLog(@"Call to DSPSVX for d2s/dt2 has not been properly configured.");
            if (NO) {
                dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                    dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                    dspsvx_work, dspsvx_iwork, &dspsvx_info);
                if (dspsvx_info != 0) {
                    if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
                    else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
                    else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
                    // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
                }
                dt2[i] = dspsvx_x[ii];
            } else {
                dt2[i] = 0.0;
            }
            ii++;
        }
//        for (i=0; i<NS; i++) {
//            for (j=0; j<NS; j++) {
//                temp[j] = d3gdsdt2[j];
//                for (k=0; k<NS; k++) {
//                    temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
//                    for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
//                }
//            }
//            dt2[i] = 0.0;
//            for (j=0; j<NS; j++) dt2[i] += - invd2gds2[i][j]*temp[j];
//        }

    }
    if (mask & NINTH  ) {   /* compute d2s/dtp */
        double d2gdsdt[NS], d2gdsdp[NS], d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS];
        double dsdt[NS], dsdp[NS], temp[NS];
        int k, l, ii, jj;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];
        [self fillD3GDSDTDPin:d3gdsdtdp withT:t andWithP:p];

        /* compute dsdt vector */
        dspsvx_nrhs =  1;
        memset(dsdt, 0, sizeof(dsdt[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdt[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdt[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdt[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdt[indexOfHydrogenIon] += -z*dsdt[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dsdt[i] = 0.0;
//            for (j=0; j<NS; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
//        }

        /* compute dsdp vector */
        dspsvx_nrhs =  1;
        memset(dsdp, 0, sizeof(dsdp[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdp[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdp[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdp[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdp[indexOfHydrogenIon] += -z*dsdp[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dsdp[i] = 0.0;
//            for (j=0; j<NS; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
//        }

        /* compute dsdtp vector */
        dspsvx_nrhs =  1;
        memset(dtp, 0, sizeof(dtp[0])*NS);
        for (i=0, ii=0; i<NS; i++) if (incOfxSpecies[i]) {
            for (j=0, jj=0; j<NS; j++) if (incOfxSpecies[j]) {
                temp[jj] = d3gdsdtdp[j];
                for (k=0; k<NS; k++) {
                    temp[jj] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                    for (l=0; l<NS; l++) temp[jj] += [self fillD3GDS3atI:j andL:k andM:l withT:t andWithP:p]*dsdt[k]*dsdp[l];
                }
                jj++;
            }
            for (jj=0; jj<dspsvx_n; jj++) dspsvx_b[jj] = -temp[jj];
            if (i == 0) NSLog(@"Call to DSPSVX for d2s/dtp has not been properly configured.");
            if (NO) {
                dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                    dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                    dspsvx_work, dspsvx_iwork, &dspsvx_info);
                if (dspsvx_info != 0) {
                    if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
                    else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
                    else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
                    // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
                }
                dtp[i] = dspsvx_x[ii];
            } else {
                dtp[i] = 0.0;
            }
            ii++;
        }
//        for (i=0; i<NS; i++) {
//            for (j=0; j<NS; j++) {
//                temp[j] = d3gdsdtdp[j];
//                for (k=0; k<NS; k++) {
//                    temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
//                    for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
//                }
//            }
//            dtp[i] = 0.0;
//            for (j=0; j<NS; j++) dtp[i] += - invd2gds2[i][j]*temp[j];
//        }

    }
    if (mask & TENTH  ) {   /* compute d2s/dp2 */
        double d2gdsdp[NS], d3gds2dp[NS][NS], d3gdsdp2[NS], dsdp[NS], temp[NS];
        int k, l, ii, jj;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        [self fillxSpeciesWith:r andWith:s];
        [self DxSpeciesDsWith:r andWith:s];
        [self fillnSpeciesWith:s];
        [self DnSpeciesDsWith:s];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];
        [self fillD3GDSDP2in:d3gdsdp2 withT:t andWithP:p];

        /* compute dsdp vector */
        dspsvx_nrhs =  1;
        memset(dsdp, 0, sizeof(dsdp[0])*NS);
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dspsvx_b[k] = -d2gdsdp[i];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dspsvx_b[k] += z*d2gdsdp[indexOfHydrogenIon]; // this is -(-z)
            }
            k++;
        }
        dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                dspsvx_work, dspsvx_iwork, &dspsvx_info);
        if (dspsvx_info != 0) {
            if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
            else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
            else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
            // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
        }
        for (NSUInteger i=0, k=0; i<NS; i++) if ((i != indexOfHydrogenIon) && incOfxSpecies[i]) {
            dsdp[i] = dspsvx_x[k++];
            id object = [endmembers objectAtIndex:i+NA];
            if ([object isKindOfClass:[HKFspeciesProperties class]]) {
                double z = [(HKFspeciesProperties *)object charge];
                dsdp[indexOfHydrogenIon] += -z*dsdp[i];
            }
        }

//        for (i=0; i<NS; i++) {
//            dsdp[i] = 0.0;
//            for (j=0; j<NS; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
//        }

        /* compute dsdp2 vector */
        dspsvx_nrhs =  1;
        memset(dp2, 0, sizeof(dp2[0])*NS);
        for (i=0, ii=0; i<NS; i++) if (incOfxSpecies[i]) {
            for (j=0, jj=0; j<NS; j++) if (incOfxSpecies[j]) {
                temp[jj] = d3gdsdp2[j];
                for (k=0; k<NS; k++) {
                    temp[jj] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                    for (l=0; l<NS; l++) temp[jj] += [self fillD3GDS3atI:j andL:k andM:l withT:t andWithP:p]*dsdp[k]*dsdp[l];
                }
                jj++;
            }
            for (jj=0; jj<dspsvx_n; jj++) dspsvx_b[jj] = -temp[jj];
            if (i == 0) NSLog(@"Call to DSPSVX for d2s/dp2 has not been properly configured.");
            if (NO) {
                dspsvx_(&dspsvx_fact, &dspsvx_uplo, &dspsvx_n, &dspsvx_nrhs, dspsvx_ap, dspsvx_afp, dspsvx_ipiv,
                    dspsvx_b, &dspsvx_ldb, dspsvx_x, &dspsvx_ldx, &dspsvx_rcond, dspsvx_ferr, dspsvx_berr,
                    dspsvx_work, dspsvx_iwork, &dspsvx_info);
                if (dspsvx_info != 0) {
                    if (dspsvx_info < 0)              NSLog(@"DEW:order() %d argument to LAPACK dppsvx is incorrect.", -dspsvx_info);
                    else if (dspsvx_info <= dspsvx_n) NSLog(@"DEW:order() %d D(i,i) is zero.", dspsvx_info);
                    else NSLog(@"DEW:order() Condition number (%g) of matrix is below machine precision.", dspsvx_rcond);
                    // NSAssert(NO, @"ERROR: DEWspecies(order)(%d). Internal error in dppsvx (%d)", mask, dspsvx_info);
                }
                dp2[i] = dspsvx_x[ii];
            } else {
                dp2[i] = 0.0;
            }
            ii++;
        }
//        for (i=0; i<NS; i++) {
//            for (j=0; j<NS; j++) {
//                temp[j] = d3gdsdp2[j];
//                for (k=0; k<NS; k++) {
//                    temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
//                    for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
//                }
//            }
//            dp2[i] = 0.0;
//            for (j=0; j<NS; j++) dp2[i] += - invd2gds2[i][j]*temp[j];
//        }

    }
}

/*
 *=============================================================================
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    r[NR] -  Array of independent compositional variables
 */

-(BOOL)test:(int)mask
          t:(double)t
          p:(double)p
         na:(int)na           // Expected number of endmember components
         nr:(int)nr           // Expected number of independent compositional variables
      names:(char **)names    // array of strings of names of endmember components
   formulas:(char **)formulas // array of strings of formulas of endmember components
          r:(double [NR])r    // array of indepependent compos variables, check bounds
          m:(double [NA])m    // array of moles of endmember components, check bounds
{
    const char *NAMES[NA] = { "H2O", "CO2,aq", "O2,aq", "HF,aq", "NaOH,aq", "Mg(OH)2,aq", "HAlO2,aq", "SiO2,aq",
        "H3PO4,aq", "SO2,aq", "HCl,aq", "KOH,aq", "Ca(OH)2,aq", "H2CrO4,aq", "Mn(OH)2,aq", "Fe(OH)2,aq", "Co(OH)2,aq" };
    const char *FORMULAS[NA] = { "H2O", "CO2,aq", "O2,aq", "HF,aq", "NaOH,aq", "Mg(OH)2,aq", "HAlO2,aq", "SiO2,aq",
        "H3PO4,aq", "SO2,aq", "HCl,aq", "KOH,aq", "Ca(OH)2,aq", "H2CrO4,aq", "Mn(OH)2,aq", "Fe(OH)2,aq", "Co(OH)2,aq" };
    BOOL result = YES;
    int i;
    double sum;

    if (mask & FIRST) {
        result = result && (na == NA);
    }
    if (mask & SECOND) {
        result = result && (nr == NR);
    }
    if (mask & THIRD) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(names[i],NAMES[i]) == 0);
        }
    }
    if (mask & FOURTH) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(formulas[i],FORMULAS[i]) == 0);
        }
    }
    /* Check bounds on the independent compositional variables */
    if (mask & FIFTH) {
        for (i=0, sum=0.0; i<NR; i++) {
            result = result && (r[i] >= 0.0) && (r[i] <= 1.0);
            sum += r[i];
        }
        result = result && (sum <= 1.0);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0; i<NA; i++) result = result && (m[i] >= 0.0);
    }

    return result;
}

-(void)convert:(int)inpMask
       outMask:(int)outMask
             t:(double)t
             p:(double)p
             e:(double [107])e              // comp of olivine in moles of elements
             m:(double [NA])m               // comp of olivine in moles of endmember components
             r:(double [NR])r               // comp of olivine in terms of the independent comp var
             x:(double [NA])x               // comp of olivine in mole fractions of endmember comp
            dm:(double [NR][NA])dm          // Jacobian matrix: dm[i][j] = dr[i]/dm[j]
           d2m:(double [NR][NA][NA])d2m     // vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]
            dr:(double [NA][NR])dr          // Jacobian matrix: dr[i][j] = dx[i]/dr[j]
           d3m:(double [NR][NA][NA][NA])d3m // 3rd deriv matrix: d3m[i][j][k][l]=d3r[i]/dm[j]dm[k]dm[l]
{
    /*---------------------------------------------------------------------------
     Not all combinations of inpMask and outMask are feasible. Valid
     combinations are:

     inpMask          outMask
     (1)  FIRST            SECOND
     (2)  SECOND           THIRD | FOURTH | FIFTH | SIXTH | EIGHTH
     (3)  THIRD            FOURTH | SEVENTH

     (1) converts a vector of moles of elements into a vector of moles of
     endmember olivine components.
     (2) calculates from a vector of moles of endmember components, one or
     all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
     (3) calculates from a vector of independent compositional variables
     mole fractions of endmember components and/or the Jacobian matrix
     dx[]/dr[]
     ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {

        /* Converts a vector of moles of elements into a vector of moles of
         end-member components. */
        static const int  H =  1;
        static const int  C =  6;
        static const int  O =  8;
        static const int  F =  9;
        static const int Na = 11;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Si = 14;
        static const int  P = 15;
        static const int  S = 16;
        static const int Cl = 17;
        static const int  K = 19;
        static const int Ca = 20;
        static const int Cr = 24;
        static const int Mn = 25;
        static const int Fe = 26;
        static const int Co = 27;

                        // H2O
        m[ 0] = (e[H]-e[F]-e[Na]-2.0*e[Mg]-e[Al]-3.0*e[P]-e[Cl]-e[K]-2.0*e[Ca]-2.0*e[Cr]-2.0*e[Mn]-2.0*e[Fe]-2.0*e[Co])/2.0;
        m[ 1] = e[C] ;  // CO2,aq
                        // O2,aq
        m[ 2] = (e[O]-e[H]/2.0-2.0*e[C]-e[Na]-2.0*e[Mg]-2.0*e[Al]-2.0*e[Si]-4.0*e[P]-2.0*e[S]-e[K]-2.0*e[Ca]-4.0*e[Cr]-2.0*e[Mn]-2.0*e[Fe]-2.0*e[Co])/2.0;
        m[ 3] = e[F];   // HF,aq
        m[ 4] = e[Na];  // NaOH,aq
        m[ 5] = e[Mg];  // Mg(OH)2,aq
        m[ 6] = e[Al];  // HAlO2,aq
        m[ 7] = e[Si];  // SiO2,aq
        m[ 8] = e[P];   // H3PO4,aq
        m[ 9] = e[S];   // SO2,aq
        m[10] = e[Cl];  // HCl,aq
        m[11] = e[K];   // KOH,aq
        m[12] = e[Ca];  // Ca(OH)2,aq
        m[13] = e[Cr];  // H2CrO4,aq
        m[14] = e[Mn];  // Mn(OH)2,aq
        m[15] = e[Fe];  // Fe(OH)2,aq
        m[16] = e[Co];  // Co(OH)2,aq

        if (m[2] < 10*DBL_EPSILON) m[2] = 0.0;

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            NSLog(@"Illegal call to convert with inpMask = %o and outMask = %o", inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
             of independent compositional variables (r) required as input for the
             remaining public functions.
             The dependent variable is taken to be SiO2 (1st component), as this
             component will never have a mole fraction of zero.                   */

            for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
        }

        if (outMask & FOURTH) {
            /* Converts a vector of moles of end-member components (m) into a vector
             of mole fractions of endmember components                            */

            for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
        }

        if (outMask & FIFTH) {
            /* Calculates the matrix dr[i]/dm[j] using m[] as input                 */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) { for (j=0; j<NA; j++) dm[i][j] = 0.0; }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)
                        dm[i][j] = (i+1 == j) ? (1.0-m[i+1]/sum)/sum : - m[i+1]/SQUARE(sum);
                }
            }
        }

        if (outMask & SIXTH) {
            /* Calculates the matrix d2r[i]/dm[j]dm[k] using m[] as input           */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++) d2m[i][j][k] = 0.0;
                    }
                }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++) {
                            d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
                            d2m[i][j][k] -= (i+1 == j) ? 1.0/SQUARE(sum) : 0.0;
                            d2m[i][j][k] -= (i+1 == k) ? 1.0/SQUARE(sum) : 0.0;
                        }
                    }
                }
            }

        }

        if (outMask & EIGHTH) {
            /* Calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input      */
            int l;

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++)  {
                            for (l=0; l<NA; l++) d3m[i][j][k][l] = 0.0;
                        }
                    }
                }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++) {
                        for (k=0; k<NA; k++) {
                            for (l=0; l<NA; l++) {
                                d3m[i][j][k][l] = -12.0*m[i]/QUARTIC(sum);
                                d3m[i][j][k][l] += (i == j) ? 4.0/CUBE(sum) : 0.0;
                                d3m[i][j][k][l] += (i == k) ? 4.0/CUBE(sum) : 0.0;
                                d3m[i][j][k][l] += (i == l) ? 4.0/CUBE(sum) : 0.0;
                                if (i == 4) d3m[i][j][k][l] /= 2.0;
                            }
                        }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            NSLog(@"Illegal call to convert with inpMask = %o and outMask = %o", inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r)
             into a vector of mole fractions of end-member components (x)            */

            for (i=0, x[0] = 1.0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (j=0; j<NR; j++) dr[0][j] = -1.0;
            for (i=1; i<NA; i++) {
                for (j=0; j<NR; j++) dr[i][j] = 0.0;
                dr[i][i-1] = 1.0;
            }
        }

    } else {
        NSLog(@"Illegal call to DEWFluid with inpMask = %o and outMask = %o", inpMask, outMask);
    }

}

-(NSString *)displayFormula:(double)t
                          p:(double)p
                          r:(double [NA])r
{
    return [NSString stringWithFormat:@"DEWFluid - Formula NYI"];
}

-(void)activity:(int)mask
              t:(double)t
              p:(double)p
              r:(double [NA])r
              a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
             mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
             dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
    double s[NS], g, dgdr[NR];
    double fr[NA][NR];
    int i, j;
    double rTotal = 0.0;
    for (i=0; i<NR; i++) rTotal += r[i];

    for (i=0; i<NR; i++) fr[0][i] = -r[i];
    for (j=1; j<NA; j++) for (i=0; i<NR; i++) fr[j][i] = (i+1 == j) ? 1.0 - r[i] : - r[i];

//    NSLog(@"Just before call to speciation.");
//    SEL spPt = @selector(speciation:t:p:r:s:dr:dt:dp:dr2:drt:drp:dt2:dtp:dp2:);
//    if ([self respondsToSelector:spPt]) NSLog(@"speciation selector exists");
//    else NSLog(@"speciation selector does not exists");
//    IMP spAd = [self methodForSelector:spPt];
//    NSLog(@"address of speciation is %p", spAd);
    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

    [self fillxSpeciesWith:r andWith:s];
    [self DxSpeciesDsWith:r andWith:s];
    [self DxSpeciesDrWith:r andWith:s];
    [self fillnSpeciesWith:s];
    [self DnSpeciesDsWith:s];
    g = [self fillGwithT:t andwithP:p];
    [self fillDGDRin:dgdr withT:t andWithP:p];

    if (mask & FIRST) {
        for (i=0; i<NA; i++) {
            for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
            id component = [endmembers objectAtIndex:i];
            a[i] -= [component getGibbsFreeEnergyFromT:t andP:p];
            a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & SECOND) {
        for (i=0; i<NA; i++) {
            for (mu[i]=g, j=0; j<NR; j++) mu[i] += (((i > 0) && (i == j+1)) ? 1.0-r[j] : - r[j])*dgdr[j];
        }
    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
        dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
        int k, l;

        [self D2nSpeciesDs2With:s];
        [self fillD2GDR2in:d2gdr2 withT:t andWithP:p];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];

        for (j=0; j<NA; j++) for (i=0; i<NR; i++) dfrdr[j][i] = -1.0;
        for (j=0; j<NA; j++) for (i=0; i<NS; i++) gs[j][i]    = -s[i];
        for (j=0; j<NA; j++) for (i=0; i<NS; i++) dgsds[j][i] = -1.0;

        [self speciation:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                /* compute activity of the i-th component */
                for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
                dx[i][k] = exp(dx[i][k]/(R*t));

                /* compute derivative of i-th activity with respect to r(k) */
                sum = (1.0+dfrdr[i][k])*dgdr[k];
                for (j=0; j<NR; j++) {
                    sum += fr[i][j]*d2gdr2[j][k];
                    for (l=0; l<NS; l++) if (incOfxSpecies[l]) sum += fr[i][j]*d2gdrds[j][l]*dsdr[l][k];
                }
                for (j=0; j<NS; j++) if (incOfxSpecies[j]) {
                    sum += gs[i][j]*d2gdrds[k][j];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l]) sum += gs[i][j]*d2gds2[j][l]*dsdr[l][k];
                }
                dx[i][k] *= sum/(R*t);
            }
        }

    }

}

-(void)gmix:(int)mask
          t:(double)t
          p:(double)p
          r:(double [NR])r
       gmix:(double *)gmix           // Gibbs energy of mixing             BINARY MASK: 0001
         dx:(double [NR])dr          // (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010
        dx2:(double [NR][NR])dr2     // (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100
        dx3:(double [NR][NR][NR])dx3 // (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000
{
    double s[NS];
    double rTotal = 0.0;
    int i;
    for (i=0; i<NR; i++) rTotal += r[i];

    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillxSpeciesWith:r andWith:s];
    [self DxSpeciesDsWith:r andWith:s];
    [self DxSpeciesDrWith:r andWith:s];
    [self fillnSpeciesWith:s];

    if (mask & FIRST) {
        double g;

        g = [self fillGwithT:t andwithP:p];
        *gmix = g;
    }

    if(mask & SECOND) {
        double dgdr[NR];

        [self fillDGDRin:dgdr withT:t andWithP:p];
        memset(dr, 0, sizeof(dr[0])*NR);
        for (i=0; i<NR; i++) if (r[i] != 0.0) dr[i] = dgdr[i];
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR];
        int j, k, l;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDR2in:d2gdr2 withT:t andWithP:p];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];

        [self speciation:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        memset(dr2, 0, sizeof(dr2[0][0])*NR*NR);
        for (i=0; i<NR; i++) if (r[i] != 0.0) {
            for (j=0; j<NR; j++) if (r[j] != 0.0) {
                dr2[i][j] = d2gdr2[i][j];
                for (k=0; k<NS; k++) if (incOfxSpecies[k]) {
                    dr2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l]) dr2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
                }
            }
        }
    }

    // To be completed
    if(mask & FOURTH) {
        double d3gdr3[NR][NR][NR], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS];
        double dsdr[NS][NR];
        int j, k, l, m, n;

        for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) d3gdr3[i][j][k] = 0.0;
        [self fillD3GDR2DSin:d3gdr2ds withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];

        [self speciation:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                for (k=0; k<NR; k++) {
                    dx3[i][j][k] = d3gdr3[i][j][k];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l]) {
                        dx3[i][j][k] += d3gdr2ds[i][j][l]*dsdr[l][k] +
                        d3gdr2ds[j][k][l]*dsdr[l][i] + d3gdr2ds[k][i][l]*dsdr[l][j];
                        for (m=0; m<NS; m++) if(incOfxSpecies[m]) {
                            dx3[i][j][k] +=
                            d3gdrds2[i][l][m]*dsdr[l][j]*dsdr[m][k] +
                            d3gdrds2[j][l][m]*dsdr[l][k]*dsdr[m][i] +
                            d3gdrds2[k][l][m]*dsdr[l][i]*dsdr[m][j];
                            for (n=0; n<NS; n++) if (incOfxSpecies[n])
                                dx3[i][j][k] +=
                                [self fillD3GDS3atI:l andL:m andM:n withT:t andWithP:p]*dsdr[l][i]*dsdr[m][j]*dsdr[n][k];
                        }
                    }
                }
            }
        }
    }
}

#pragma mark -
#pragma mark Testing routines for derivatives

-(void)testD2GDR2withT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR];
    double EPS = sqrt(DBL_EPSILON);
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];
    double g, dgdr[NR], d2gdr2[NR][NR];

    NSLog(@"<><><> Test of dgdr and d2gdr2 derivatives calling gmix <><><>");
    [self gmix:(FIRST | SECOND | THIRD) t:t p:p r:r gmix:&g dx:dgdr dx2:d2gdr2 dx3:NULL];
    for (NSUInteger i=0; i<NR; i++) if (r[i] != 0.0) {
        double gForward, dgdrForward[NR], g2Forward, dgdr2Forward[NR];
        double gBackward, dgdrBackward[NR], g2Backward, dgdr2Backward[NR];
        rStep[i] = r[i]*(1.0+EPS);
        [self gmix:(FIRST | SECOND) t:t p:p r:rStep gmix:&gForward dx:dgdrForward dx2:NULL dx3:NULL];
        rStep[i] = r[i]*(1.0-EPS);
        [self gmix:(FIRST | SECOND) t:t p:p r:rStep gmix:&gBackward dx:dgdrBackward dx2:NULL dx3:NULL];
        rStep[i] = r[i]*(1.0+2.0*EPS);
        [self gmix:(FIRST | SECOND) t:t p:p r:rStep gmix:&g2Forward dx:dgdr2Forward dx2:NULL dx3:NULL];
        rStep[i] = r[i]*(1.0-2.0*EPS);
        [self gmix:(FIRST | SECOND) t:t p:p r:rStep gmix:&g2Backward dx:dgdr2Backward dx2:NULL dx3:NULL];
        rStep[i] = r[i];
        double gEst = (-g2Forward + 8.0*gForward - 8.0*gBackward + g2Backward)/12.0/r[i]/EPS;
        NSLog(@"dgdr[%lu]      est %23.15e act %23.15e dif %10.3f %%", i, gEst, dgdr[i], 100.0*(gEst - dgdr[i])/dgdr[i]);
        for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) {
            double d2gdr2Est = (-dgdr2Forward[j] + 8.0*dgdrForward[j] - 8.0*dgdrBackward[j] + dgdr2Backward[j])/r[i]/EPS/12.0;
            NSLog(@"d2gdr2[%lu][%lu] est %23.15e act %23.15e dif %10.3f %%", i, j, d2gdr2Est, d2gdr2[i][j],
                  100.0*(d2gdr2Est - d2gdr2[i][j])/d2gdr2[i][j]);
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testPrimD2GDR2withT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR];
    double EPS = sqrt(DBL_EPSILON);
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];
    double dgdr[NR], d2gdr2[NR][NR], sStep[NS];

    NSLog(@"<><><> Test of dgdr and d2gdr2 derivatives calling primitive functions <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self azeroWithR:rStep];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];
    [self DxSpeciesDrWith:rStep andWith:sStep];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];
    [self fillDGDRin:dgdr withT:t andWithP:p];
    [self fillD2GDR2in:d2gdr2 withT:t andWithP:p];

    for (NSUInteger i=0; i<NR; i++) if (r[i] != 0.0) {
        double gForward, dgdrForward[NR], g2Forward, dgdr2Forward[NR];
        double gBackward, dgdrBackward[NR], g2Backward, dgdr2Backward[NR];
        rStep[i] = r[i]*(1.0+EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self azeroWithR:rStep];
        gForward = [self fillGwithT:t andwithP:p];
        [self fillDGDRin:dgdrForward withT:t andWithP:p];
        rStep[i] = r[i]*(1.0-EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self azeroWithR:rStep];
        gBackward = [self fillGwithT:t andwithP:p];
        [self fillDGDRin:dgdrBackward withT:t andWithP:p];
        rStep[i] = r[i]*(1.0+2.0*EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self azeroWithR:rStep];
        g2Forward = [self fillGwithT:t andwithP:p];
        [self fillDGDRin:dgdr2Forward withT:t andWithP:p];
        rStep[i] = r[i]*(1.0-2.0*EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self azeroWithR:rStep];
        g2Backward = [self fillGwithT:t andwithP:p];
        [self fillDGDRin:dgdr2Backward withT:t andWithP:p];
        rStep[i] = r[i];
        double gEst = (-g2Forward + 8.0*gForward - 8.0*gBackward + g2Backward)/12.0/r[i]/EPS;
        NSLog(@"dgdr[%lu]      est %23.15e act %23.15e dif %10.3f %%", i, gEst, dgdr[i], 100.0*(gEst - dgdr[i])/dgdr[i]);
        for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) {
            double d2gdr2Est = (-dgdr2Forward[j] + 8.0*dgdrForward[j] - 8.0*dgdrBackward[j] + dgdr2Backward[j])/r[i]/EPS/12.0;
            NSLog(@"d2gdr2[%lu][%lu] est %23.15e act %23.15e dif %10.3f %%", i, j, d2gdr2Est, d2gdr2[i][j],
                  100.0*(d2gdr2Est - d2gdr2[i][j])/d2gdr2[i][j]);
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testPrimD2GDS2withT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS], dgds[NS], d2gds2[NS][NS];
    double EPS = 0.001;
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of dgds and d2gds2 derivatives calling primitive functions <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self azeroWithR:rStep];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];
    [self DxSpeciesDrWith:rStep andWith:sStep];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];

    [self fillDGDSin:dgds withT:t andWithP:p];
    [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
    for (NSUInteger i=0; i<NS; i++) if (sStep[i] != 0.0) {
        double gForward, dgdsForward[NS], g2Forward, dgds2Forward[NS];
        double gBackward, dgdsBackward[NS], g2Backward, dgds2Backward[NS];
        NSString *iName = [[[endmembers objectAtIndex:i+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        double sCur = sStep[i];

        sStep[i] = sCur*(1.0+EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        gForward = [self fillGwithT:t andwithP:p];
        [self fillDGDSin:dgdsForward withT:t andWithP:p];

        sStep[i] = sCur*(1.0-EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        gBackward = [self fillGwithT:t andwithP:p];
        [self fillDGDSin:dgdsBackward withT:t andWithP:p];

        sStep[i] = sCur*(1.0+2.0*EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        g2Forward = [self fillGwithT:t andwithP:p];
        [self fillDGDSin:dgds2Forward withT:t andWithP:p];

        sStep[i] = sCur*(1.0-2.0*EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        g2Backward = [self fillGwithT:t andwithP:p];
        [self fillDGDSin:dgds2Backward withT:t andWithP:p];

        sStep[i] = sCur;
        double gEst = (-g2Forward + 8.0*gForward - 8.0*gBackward + g2Backward)/12.0/sStep[i]/EPS;
        NSLog(@"dgds[%@]               est %23.15e act %23.15e dif %10.3f %%", iName, gEst, dgds[i], 100.0*(gEst - dgds[i])/dgds[i]);
        for (NSUInteger j=0; j<NS; j++) if (sStep[j] != 0.0) {
            NSString *jName = [[[endmembers objectAtIndex:j+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
            double d2gds2Est = (-dgds2Forward[j] + 8.0*dgdsForward[j] - 8.0*dgdsBackward[j] + dgds2Backward[j])/sStep[i]/EPS/12.0;
            NSLog(@"d2gds2[%@][%@] est %23.15e act %23.15e dif %10.3f %%", iName, jName, d2gds2Est, d2gds2[i][j],
                  100.0*(d2gds2Est - d2gds2[i][j])/d2gds2[i][j]);
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testPrimD2GDRDSwithT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS], d2gdrds[NR][NS];
    double EPS = 0.001;
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of d2gdrds derivative calling primitive functions <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self azeroWithR:rStep];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];
    [self DxSpeciesDrWith:rStep andWith:sStep];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];
    [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
    for (NSUInteger i=0; i<NS; i++) if (sStep[i] != 0.0) {
        double dgdrForward[NR], dgdr2Forward[NR];
        double dgdrBackward[NR], dgdr2Backward[NR];
        NSString *iName = [[[endmembers objectAtIndex:i+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        double sCur = sStep[i];

        sStep[i] = sCur*(1.0+EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        [self fillDGDRin:dgdrForward withT:t andWithP:p];

        sStep[i] = sCur*(1.0-EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        [self fillDGDRin:dgdrBackward withT:t andWithP:p];

        sStep[i] = sCur*(1.0+2.0*EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        [self fillDGDRin:dgdr2Forward withT:t andWithP:p];

        sStep[i] = sCur*(1.0-2.0*EPS);
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        [self fillDGDRin:dgdr2Backward withT:t andWithP:p];

        sStep[i] = sCur;
        for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) {
            double d2gdrdsEst = (-dgdr2Forward[j] + 8.0*dgdrForward[j] - 8.0*dgdrBackward[j] + dgdr2Backward[j])/sStep[i]/EPS/12.0;
            NSLog(@"d2gdrds[%lu][%@] est %23.15e act %23.15e dif %10.3f %%", j, iName, d2gdrdsEst, d2gdrds[j][i],
                  100.0*(d2gdrdsEst - d2gdrds[j][i])/d2gdrds[j][i]);
        }

    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testAzerowithT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR];
    double EPS = sqrt(DBL_EPSILON);
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of dazerodr and d2azerodr2 derivatives calling azero <><><>");
    [self azeroWithR:rStep];
    double azeroRef = aZero;
    NSLog(@"azero = %23.15e", azeroRef);
    double dazerodrRef[NR];
    for (NSUInteger i=0; i<NR; i++) if (r[i] != 0.0) dazerodrRef[i] = [self DazeroDrForIndex:i];

    for (NSUInteger i=0; i<NR; i++) if (rStep[i] != 0.0) {
        double dazerodrForward[NR];
        rStep[i] = r[i]*(1.0+EPS);
        [self azeroWithR:rStep];
        for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) dazerodrForward[j] = [self DazeroDrForIndex:j];
        NSLog(@"dazerodr[%lu]      est %23.15e act %23.15e dif %10.3f %%", i, (aZero-azeroRef)/EPS/r[i], dazerodrRef[i],
              100.0*((aZero-azeroRef)/EPS/r[i] - dazerodrRef[i])/dazerodrRef[i]);
        rStep[i] = r[i];
        [self azeroWithR:rStep]; // must reset aZeroDifference
        for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) {
            double d2azerodr2 = [self D2azeroDr2ForIndex:i andIndex:j];
            double d2azerodr2Est = (dazerodrForward[j]-dazerodrRef[j])/r[i]/EPS;
            NSLog(@"d2azerodr2[%lu][%lu] est %23.15e act %23.15e dif %10.3f %%", i, j, d2azerodr2Est, d2azerodr2,
                  100.0*(d2azerodr2Est-d2azerodr2)/d2azerodr2);
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testD2DHDR2withT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS];
    double EPS = sqrt(DBL_EPSILON);
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of fillDebyeHuckelExcessFreeEnergy and derivatives <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self azeroWithR:rStep];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];
    [self DxSpeciesDrWith:rStep andWith:sStep];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];
    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] != 0.0) {
        NSString *spName = [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        double dgDHdr[NR], d2gDHdr2[NR][NR];
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        initializeDHexcessTermsZerothOrder = YES;
        initializeDHexcessTermsFirstOrderDr = YES;
        [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdr withT:t andWithP:p];
        [self fillD2debyeHuckelExcessFreeEnergyDr2ForSpeciesIndex:i in:d2gDHdr2 withT:t andWithP:p];
        for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) {
            double dgDHdrForward[NR];
            rStep[j] = r[j]*(1.0+EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            double gDHForward = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrForward withT:t andWithP:p];

            double dgDHdrBackward[NR];
            rStep[j] = r[j]*(1.0-EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            double gDHBackward = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrBackward withT:t andWithP:p];

            double dgDHdrForward2[NR];
            rStep[j] = r[j]*(1.0+2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            double gDHForward2 = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrForward2 withT:t andWithP:p];

            double dgDHdrBackward2[NR];
            rStep[j] = r[j]*(1.0-2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            double gDHBackward2 = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrBackward2 withT:t andWithP:p];

            rStep[j] = r[j];
            double dgDHdrEst = (-gDHForward2 + 8.0*gDHForward - 8.0*gDHBackward + gDHBackward2)/r[j]/EPS/12.0;
            NSLog(@"dDHdr  [%@][%lu]    est %23.15e act %23.15e dif %10.3f %%", spName, j, dgDHdrEst, dgDHdr[j], 100.0*(dgDHdrEst - dgDHdr[j])/dgDHdr[j]);
            for (NSUInteger k=0; k<NR; k++) if (r[k] != 0.0) {
                double d2DHdr2Est = (-dgDHdrForward2[k] + 8.0*dgDHdrForward[k] - 8.0*dgDHdrBackward[k] + dgDHdrBackward2[k])/r[j]/EPS/12.0;
                NSLog(@"d2DHdr2[%@][%lu][%lu] est %23.15e act %23.15e dif %10.3f %%", spName, j, k, d2DHdr2Est, d2gDHdr2[j][k],
                      100.0*(d2DHdr2Est - d2gDHdr2[j][k])/d2gDHdr2[j][k]);
            }
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testD2DHDRDSwithT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS];
    double EPS = 0.001;
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of fillDebyeHuckelExcessFreeEnergy and drds derivatives <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self azeroWithR:rStep];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];
    [self DxSpeciesDrWith:rStep andWith:sStep];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];
    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] != 0.0) {
        NSString *spName = [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        double dgDHdr[NR], d2gDHdrds[NR][NS];
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        initializeDHexcessTermsZerothOrder = YES;
        initializeDHexcessTermsFirstOrderDr = YES;
        initializeDHexcessTermsFirstOrderDs = YES;
        initializeDHexcessTermsSecondOrderDrDs = YES;
        [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdr withT:t andWithP:p];
        [self fillD2debyeHuckelExcessFreeEnergyDrDsForSpeciesIndex:i in:d2gDHdrds withT:t andWithP:p];
        for (NSUInteger j=0; j<NS; j++) if (sStep[j] != 0.0) {
            double sCur = sStep[j];
            NSString *jName = [[[endmembers objectAtIndex:j+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];

            double dgDHdrForward[NR];
            sStep[j] = sCur*(1.0+EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrForward withT:t andWithP:p];

            double dgDHdrBackward[NR];
            sStep[j] = sCur*(1.0-EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrBackward withT:t andWithP:p];

            double dgDHdrForward2[NR];
            sStep[j] = sCur*(1.0+2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrForward2 withT:t andWithP:p];

            double dgDHdrBackward2[NR];
            sStep[j] = sCur*(1.0-2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            initializeDHexcessTermsFirstOrderDr = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDrForSpeciesIndex:i in:dgDHdrBackward2 withT:t andWithP:p];

            sStep[j] = sCur;
            for (NSUInteger k=0; k<NR; k++) if (r[k] != 0.0) {
                double d2DHdrdsEst = (-dgDHdrForward2[k] + 8.0*dgDHdrForward[k] - 8.0*dgDHdrBackward[k] + dgDHdrBackward2[k])/sStep[j]/EPS/12.0;
                NSLog(@"d2DHdrds[%@][%lu][%@] est %23.15e act %23.15e dif %10.3f %%", spName, k, jName, d2DHdrdsEst, d2gDHdrds[k][j],
                      (fabs(d2gDHdrds[k][j]) >= 10.0*DBL_EPSILON) ? 100.0*(d2DHdrdsEst - d2gDHdrds[k][j])/d2gDHdrds[k][j] : 0.0);
            }
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testD2DHDS2withT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS];
    double EPS = 0.001;
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of fillDebyeHuckelExcessFreeEnergy and ordering derivatives <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self azeroWithR:rStep];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];
    [self DxSpeciesDrWith:rStep andWith:sStep];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];
    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] != 0.0) {
        NSString *spName = [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        double dgDHds[NS], d2gDHds2[NS][NS];
        [self fillxSpeciesWith:rStep andWith:sStep];
        [self DxSpeciesDsWith:rStep andWith:sStep];
        [self DxSpeciesDrWith:rStep andWith:sStep];
        [self fillnSpeciesWith:sStep];
        [self DnSpeciesDsWith:sStep];
        [self D2nSpeciesDs2With:sStep];
        [self azeroWithR:rStep];
        initializeDHexcessTermsZerothOrder    = YES;
        initializeDHexcessTermsFirstOrderDs   = YES;
        [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:dgDHds withT:t andWithP:p];
        initializeDHexcessTermsSecondOrderDs2 = YES;
        [self fillD2debyeHuckelExcessFreeEnergyDs2ForSpeciesIndex:i in:d2gDHds2 withT:t andWithP:p];
        for (NSUInteger j=0; j<NS; j++) if (sStep[j] != 0.0) {
            double sCur = sStep[j];
            NSString *jName = [[[endmembers objectAtIndex:j+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];

            double dgDHdsForward[NS];
            sStep[j] = sCur*(1.0+EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder  = YES;
            double gDHForward = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            initializeDHexcessTermsFirstOrderDs = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:dgDHdsForward withT:t andWithP:p];

            double dgDHdsBackward[NS];
            sStep[j] = sCur*(1.0-EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            double gDHBackward = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            initializeDHexcessTermsFirstOrderDs = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:dgDHdsBackward withT:t andWithP:p];

            double dgDHdsForward2[NS];
            sStep[j] = sCur*(1.0+2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            double gDHForward2 = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            initializeDHexcessTermsFirstOrderDs = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:dgDHdsForward2 withT:t andWithP:p];

            double dgDHdsBackward2[NS];
            sStep[j] = sCur*(1.0-2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            [self DxSpeciesDsWith:rStep andWith:sStep];
            [self DxSpeciesDrWith:rStep andWith:sStep];
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            [self D2nSpeciesDs2With:sStep];
            [self azeroWithR:rStep];
            initializeDHexcessTermsZerothOrder = YES;
            double gDHBackward2 = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
            initializeDHexcessTermsFirstOrderDs = YES;
            [self fillDdebyeHuckelExcessFreeEnergyDsForSpeciesIndex:i in:dgDHdsBackward2 withT:t andWithP:p];

            sStep[j] = sCur;
            double dgDHdsEst = (-gDHForward2 + 8.0*gDHForward - 8.0*gDHBackward + gDHBackward2)/sStep[j]/EPS/12.0;
            NSLog(@"dDHds  [%@][%@]             est %23.15e act %23.15e dif %10.3f %%", spName, jName, dgDHdsEst, dgDHds[j],
                  100.0*(dgDHdsEst - dgDHds[j])/dgDHds[j]);
            for (NSUInteger k=0; k<NS; k++) if (sStep[k] != 0.0) {
                NSString *kName = [[[endmembers objectAtIndex:k+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
                double d2DHds2Est = (-dgDHdsForward2[k] + 8.0*dgDHdsForward[k] - 8.0*dgDHdsBackward[k] + dgDHdsBackward2[k])/sStep[j]/EPS/12.0;
                NSLog(@"d2DHds2[%@][%@][%@] est %23.15e act %23.15e dif %10.3f %%", spName, jName, kName, d2DHds2Est, d2gDHds2[j][k],
                      100.0*(d2DHds2Est - d2gDHds2[j][k])/d2gDHds2[j][k]);
            }
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testFillnSpecieswithT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS];
    double EPS = 0.001;
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of nSpecies and derivatives <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillnSpeciesWith:sStep];
    [self DnSpeciesDsWith:sStep];
    [self D2nSpeciesDs2With:sStep];
    double DnSpeciesRef[NS];
    for (NSUInteger k=0; k<NS; k++) DnSpeciesRef[k] = dnSpeciesds[k];

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] != 0.0) {
        NSString *spName = [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        for (NSUInteger j=0; j<NS; j++) if (sStep[j] != 0.0) {
            NSString *jName = [[[endmembers objectAtIndex:j+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
            double sCur = sStep[j];

            double DnSpeciesForward[NS];
            sStep[j] = sCur*(1.0+EPS);
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            double nSpeciesForward = nSpecies;
            for (NSUInteger k=0; k<NS; k++) DnSpeciesForward[k] = dnSpeciesds[k];

            double DnSpeciesBackward[NS];
            sStep[j] = sCur*(1.0-EPS);
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            double nSpeciesBackward = nSpecies;
            for (NSUInteger k=0; k<NS; k++) DnSpeciesBackward[k] = dnSpeciesds[k];

            double DnSpeciesForward2[NS];
            sStep[j] = sCur*(1.0+2.0*EPS);
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            double nSpeciesForward2 = nSpecies;
            for (NSUInteger k=0; k<NS; k++) DnSpeciesForward2[k] = dnSpeciesds[k];

            double DnSpeciesBackward2[NS];
            sStep[j] = sCur*(1.0-2.0*EPS);
            [self fillnSpeciesWith:sStep];
            [self DnSpeciesDsWith:sStep];
            double nSpeciesBackward2 = nSpecies;
            for (NSUInteger k=0; k<NS; k++) DnSpeciesBackward2[k] = dnSpeciesds[k];

            sStep[j] = sCur;
            double DnSpeciesEst = (-nSpeciesForward2 + 8.0*nSpeciesForward - 8.0*nSpeciesBackward + nSpeciesBackward2)/sStep[j]/EPS/12.0;
            NSLog(@"dnSpecies [%@][%@]             est %23.15e act %23.15e dif %10.3f %%", spName, jName, DnSpeciesEst, DnSpeciesRef[j],
                  (fabs(DnSpeciesRef[j]) >= 10.0*DBL_EPSILON) ? 100.0*(DnSpeciesEst - DnSpeciesRef[j])/DnSpeciesRef[j] : 0.0);
            for (NSUInteger k=0; k<NS; k++) if (sStep[k] != 0.0) {
                NSString *kName = [[[endmembers objectAtIndex:k+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
                double D2nSpeciesEst = (-DnSpeciesForward2[k] + 8.0*DnSpeciesForward[k] - 8.0*DnSpeciesBackward[k]
                                        + DnSpeciesBackward2[k])/sStep[j]/EPS/12.0;
                NSLog(@"d2nSpecies[%@][%@][%@] est %23.15e act %23.15e dif %10.3f %%", spName, jName, kName, D2nSpeciesEst, d2nSpeciesds2[j][k],
                      (fabs(d2nSpeciesds2[j][k]) >= 10.0*DBL_EPSILON) ? 100.0*(D2nSpeciesEst - d2nSpeciesds2[j][k])/d2nSpeciesds2[j][k] : 0.0);
            }
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testFillxSpecieswithT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], sStep[NS];
    double EPS = 0.001; //sqrt(DBL_EPSILON);
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of xSpecies and derivatives <><><>");
    [self speciation:FIRST t:t p:p r:rStep s:sStep dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillxSpeciesWith:rStep andWith:sStep];
    [self DxSpeciesDsWith:rStep andWith:sStep];

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] != 0.0) {
        NSString *spName = [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        for (NSUInteger j=0; j<NS; j++) if (sStep[j] != 0.0) {
            NSString *jName = [[[endmembers objectAtIndex:j+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
            double sCur = sStep[j];

            sStep[j] = sCur*(1.0+EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            double xSpeciesForward = xSpecies[i];

            sStep[j] = sCur*(1.0-EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            double xSpeciesBackward = xSpecies[i];

            sStep[j] = sCur*(1.0+2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            double xSpeciesForward2 = xSpecies[i];

            sStep[j] = sCur*(1.0-2.0*EPS);
            [self fillxSpeciesWith:rStep andWith:sStep];
            double xSpeciesBackward2 = xSpecies[i];

            sStep[j] = sCur;
            double DxSpeciesEst = (-xSpeciesForward2 + 8.0*xSpeciesForward - 8.0*xSpeciesBackward + xSpeciesBackward2)/sStep[j]/EPS/12.0;
            NSLog(@"dxSpecies [%@][%@] est %23.15e act %23.15e dif %10.3f %%", spName, jName, DxSpeciesEst, dxSpecies[i][j],
                  (fabs(dxSpecies[i][j]) >= 10.0*DBL_EPSILON) ? 100.0*(DxSpeciesEst - dxSpecies[i][j])/dxSpecies[i][j] : 0.0);
        }
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)testDRDSwithT:(double)t p:(double)p r:(double [16])r {
    double rStep[NR], dsdrRef[NS][NR], sRef[NS];
    double EPS = 0.000001; //sqrt(DBL_EPSILON);
    for (NSUInteger i=0; i<NR; i++) rStep[i] = r[i];

    NSLog(@"<><><> Test of ordering function dsdr <><><>");
    [self speciation:FIRST | SECOND t:t p:p r:r s:sRef dr:dsdrRef dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    for (NSUInteger i=0; i<NS; i++) if (sRef[i] != 0.0) {
        NSString *spName = [[[endmembers objectAtIndex:i+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
        NSLog(@"   s[%@]    act %23.15e", spName, sRef[i]);
    }

    for (NSUInteger j=0; j<NR; j++) if (r[j] != 0.0) {
        double sForward[NS];
        rStep[j] = r[j]*(1.0+EPS);
        [self speciation:FIRST t:t p:p r:rStep s:sForward dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        double sBackward[NS];
        rStep[j] = r[j]*(1.0-EPS);
        [self speciation:FIRST t:t p:p r:rStep s:sBackward dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        double sForward2[NS];
        rStep[j] = r[j]*(1.0+2.0*EPS);
        [self speciation:FIRST t:t p:p r:rStep s:sForward2 dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        double sBackward2[NS];
        rStep[j] = r[j]*(1.0-2.0*EPS);
        [self speciation:FIRST t:t p:p r:rStep s:sBackward2 dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (NSUInteger i=0; i<NS; i++) if (sRef[i] != 0.0) {
            NSString *spName = [[[endmembers objectAtIndex:i+NA] phaseName] stringByPaddingToLength:10 withString:@" " startingAtIndex:0];
            double dsdrEst = (-sForward2[i] + 8.0*sForward[i] - 8.0*sBackward[i] + sBackward2[i])/r[j]/EPS/12.0;
            NSLog(@"dsdr[%@][%lu] est %23.15e act %23.15e dif %10.3f %%", spName, j, dsdrEst, dsdrRef[i][j],
                  100.0*(dsdrEst-dsdrRef[i][j])/dsdrRef[i][j]);
        }

        rStep[j] = r[j];
    }
    NSLog(@"<><><> End of Test <><><>");
}

-(void)hmix:(int)mask
          t:(double)t
          p:(double)p
          r:(double [NR])r
       hmix:(double *)hmix // Enthalpy of mixing BINARY MASK: 1
{
    double s[NS], g, entropy;
    double rTotal = 0.0;
    int i;
    for (i=0; i<NR; i++) rTotal += r[i];

    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillxSpeciesWith:r andWith:s];
    [self fillnSpeciesWith:s];

    g = [self fillGwithT:t andwithP:p];
    entropy = [self fillSwithT:t andWithP:p];
    *hmix = g + t*entropy;
}

-(void)smix:(int)mask
          t:(double)t
          p:(double)p
          r:(double [NR])r
       smix:(double *)smix        // Entropy of mixing                  BINARY MASK: 001
         dx:(double [NR])dr       // (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010
        dx2:(double [NR][NR])dr2  // (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100
{
    double s[NS];
    double rTotal = 0.0;
    int i;
    for (i=0; i<NR; i++) rTotal += r[i];

    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillxSpeciesWith:r andWith:s];
    [self DxSpeciesDsWith:r andWith:s];
    [self DxSpeciesDrWith:r andWith:s];
    [self fillnSpeciesWith:s];

    if (mask & FIRST) {
        double entropy = [self fillSwithT:t andWithP:p];
        *smix = entropy;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS], dsdr[NS][NR], dsdt[NS];
        int k, l;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDRDTin:d2gdrdt withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];

        [self speciation:SECOND | THIRD t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            dr[i] = d2gdrdt[i];
            for (k=0; k<NS; k++) if(incOfxSpecies[k]) {
                dr[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
                for (l=0; l<NS; l++) if(incOfxSpecies[l]) dr[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
            }
            dr[i] *= -1.0;
        }
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdt[NS], d3gdr2ds[NR][NR][NS],
        d3gdr2dt[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
        d3gds2dt[NS][NS], dsdr[NS][NR], dsdt[NS],
        d2sdr2[NS][NR][NR], d2sdrdt[NS][NR];
        int i, j, k, l, m;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD3GDR2DSin:d3gdr2ds withT:t andWithP:p];
        [self fillD3GDR2DTin:d3gdr2dt withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDTin:d3gdrdsdt withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];

        [self speciation:SECOND | THIRD | FIFTH | SIXTH t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:NULL dr2:d2sdr2 drt:d2sdrdt drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dr2[i][j] = d3gdr2dt[i][j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k]) {
                    dr2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
                    + d3gdrdsdt[i][k]*dsdr[k][j]
                    + d3gdrdsdt[j][k]*dsdr[k][i]
                    + d2gdsdt[k]*d2sdr2[k][i][j]
                    + d2gdrds[i][k]*d2sdrdt[k][j]
                    + d2gdrds[j][k]*d2sdrdt[k][i];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l]) {
                        dr2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
                        + d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
                        + d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
                        + d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
                        + d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
                        + d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
                        for (m=0; m<NS; m++) if(incOfxSpecies[m])
                            dr2[i][j] += [self fillD3GDS3atI:k andL:l andM:m withT:t andWithP:p]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
                    }
                }
                dr2[i][j] *= -1.0;
            }
        }
    }
}

-(void)cpmix:(int)mask
           t:(double)t
           p:(double)p
           r:(double [NR])r
       cpmix:(double *)cpmix // Heat capacity of mixing         BINARY MASK: 001
          dt:(double *)dt    // d(cp)/d(t)                      BINARY MASK: 010
          dx:(double [NR])dr // d(cp)/d(x[])d(t)                BINARY MASK: 100
{
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    double rTotal = 0.0;
    int i, j, k, l;
    for (i=0; i<NR; i++) rTotal += r[i];

    [self speciation:FIRST | THIRD t:t p:p r:r s:s dr:NULL dt:dsdt dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillxSpeciesWith:r andWith:s];
    [self DxSpeciesDsWith:r andWith:s];
    [self DxSpeciesDrWith:r andWith:s];
    [self fillnSpeciesWith:s];
    [self DnSpeciesDsWith:s];
    [self D2nSpeciesDs2With:s];
    [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
    [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
    d2gdt2  = [self fillD2GDT2withT:t andWithP:p];

    if (mask & FIRST) {
        *cpmix = d2gdt2;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;
    }

    if(mask & SECOND) {
        double d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS], temp;
        double d3gdt3 = [self fillD3GDT3withT:t andWithP:p];

        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDSDT2in:d3gdsdt2 withT:t andWithP:p];

        [self speciation:EIGHTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:d2sdt2 dtp:NULL dp2:NULL];

        /* compute d2gdt2 */
        temp = d2gdt2;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            temp += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) temp += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }

        *dt = d3gdt3;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *dt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) {
                *dt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
                + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k])
                    *dt += [self fillD3GDS3atI:i andL:j andM:k withT:t andWithP:p]*dsdt[i]*dsdt[j]*dsdt[k];
            }
        }
        *dt = -t*(*dt) - temp;
    }

    if(mask & THIRD) {
        double d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
        d3gds2dt[NS][NS], d2gdrds[NR][NS], d3gdrdt2[NR], d3gdsdt2[NS],
        dsdr[NS][NR], d2sdrdt[NS][NR], d2sdt2[NS];

        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDTin:d3gdrdsdt withT:t andWithP:p];
        [self fillD3GDRDT2in:d3gdrdt2 withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDSDT2in:d3gdsdt2 withT:t andWithP:p];

        [self speciation:SECOND | SIXTH | EIGHTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:d2sdrdt drp:NULL dt2:d2sdt2 dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0,dr[i]=d3gdrdt2[i]; j<NS; j++) if(incOfxSpecies[j]) {
                dr[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
                2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k]) {
                    dr[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
                    2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
                    2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
                    d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l])
                        dr[i] += [self fillD3GDS3atI:j andL:k andM:l withT:t andWithP:p]*dsdr[j][i]*dsdt[k]*dsdt[l];
                }
            }
            dr[i] *= -t;
        }
    }
}

-(void)vmix:(int)mask
          t:(double)t
          p:(double)p
          r:(double [NR])r
       vmix:(double *)vmix       // Volume of mixing                BINARY MASK: 0000000001
         dx:(double [NR])dr      // pointer to dx[]) d(v)/d(x[])    BINARY MASK: 0000000010
        dx2:(double [NR][NR])dr2 // pointer to dx2[][]) d(v)/d(x[])2BINARY MASK: 0000000100
         dt:(double *)dt         // d(v)/d(t)                       BINARY MASK: 0000001000
         dp:(double *)dp         // d(v)/d(p)                       BINARY MASK: 0000010000
        dt2:(double *)dt2        // d2(v)/d(t)2                     BINARY MASK: 0000100000
       dtdp:(double *)dtdp       // d2(v)/d(t)d(p)                  BINARY MASK: 0001000000
        dp2:(double *)dp2        // d2(v)/d(p)2                     BINARY MASK: 0010000000
       dxdt:(double [NR])drdt    // d2(v)/d(x[])d(t)                BINARY MASK: 0100000000
       dxdp:(double [NR])drdp    // d2(v)/d(x[])d(p)                BINARY MASK: 1000000000
{
    double s[NS];
    double rTotal = 0.0;
    int i, j;
    for (i=0; i<NR; i++) rTotal += r[i];

    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    [self fillxSpeciesWith:r andWith:s];
    [self DxSpeciesDsWith:r andWith:s];
    [self DxSpeciesDrWith:r andWith:s];
    [self fillnSpeciesWith:s];

    if (mask & FIRST) {
        double v = [self fillVwithT:t andWithP:p];
        *vmix = v;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS], dsdr[NS][NR], dsdp[NS];
        int k;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDRDPin:d2gdrdp withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];

        [self speciation:SECOND | FOURTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            dr[i] = d2gdrdp[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) {
                dr[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NS; k++) if(incOfxSpecies[k]) dr[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdp[NS], d3gdr2ds[NR][NR][NS],
        d3gdr2dp[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdp[NR][NS],
        d3gds2dp[NS][NS], dsdr[NS][NR], dsdp[NS],
        d2sdr2[NS][NR][NR], d2sdrdp[NS][NR];
        int k, l, m;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDR2DSin:d3gdr2ds withT:t andWithP:p];
        [self fillD3GDR2DPin:d3gdr2dp withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDPin:d3gdrdsdp withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];

        [self speciation:SECOND | FOURTH | FIFTH | SEVENTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:d2sdr2 drt:NULL drp:d2sdrdp dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dr2[i][j] = d3gdr2dp[i][j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k]) {
                    dr2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
                    + d3gdrdsdp[i][k]*dsdr[k][j]
                    + d3gdrdsdp[j][k]*dsdr[k][i]
                    + d2gdsdp[k]*d2sdr2[k][i][j]
                    + d2gdrds[i][k]*d2sdrdp[k][j]
                    + d2gdrds[j][k]*d2sdrdp[k][i];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l]) {
                        dr2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
                        + d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
                        + d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
                        + d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
                        + d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
                        + d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
                        for (m=0; m<NS; m++) if(incOfxSpecies[m])
                            dr2[i][j] += [self fillD3GDS3atI:k andL:l andM:m withT:t andWithP:p]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
                    }
                }
            }
        }
    }

    if(mask & FOURTH) {
        double d2gdtdp = [self fillD2GDTDPwithT:t andWithP:p];
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS];

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];

        [self speciation:THIRD | FOURTH t:t p:p r:r s:NULL dr:NULL dt:dsdt dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        *dt = d2gdtdp;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }
    }

    if(mask & FIFTH) {
        double d2gdp2 = [self fillD2GDP2withT:t andWithP:p];
        double d2gds2[NS][NS], d2gdsdp[NS], dsdp[NS];

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];

        [self speciation:FOURTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        *dp = d2gdp2;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *dp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }
    }

    if(mask & SIXTH) {
        double d3gdt2dp = [self fillD2GDTDPwithT:t andWithP:p];
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS],
        d3gds2dp[NS][NS], d3gds2dt[NS][NS], d3gdsdtdp[NS], d3gdsdt2[NS],
        dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];
        int k;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];
        [self fillD3GDSDT2in:d3gdsdt2 withT:t andWithP:p];
        [self fillD3GDSDTDPin:d3gdsdtdp withT:t andWithP:p];

        [self speciation:THIRD | FOURTH | EIGHTH | NINTH t:t p:p r:r s:NULL dr:NULL dt:dsdt dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:d2sdt2 dtp:d2sdtdp dp2:NULL];

        *dt2 = d3gdt2dp;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *dt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
            + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) {
                *dt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k])
                    *dt2 += [self fillD3GDS3atI:i andL:j andM:k withT:t andWithP:p]*dsdt[i]*dsdt[j]*dsdp[k];
            }
        }
    }

    if(mask & SEVENTH) {
        double d3gdtdp2 = [self fillD3GDTDP2withT:t andWithP:p];
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS],
        d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS],
        dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];
        int k;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];
        [self fillD3GDSDTDPin:d3gdsdtdp withT:t andWithP:p];
        [self fillD3GDSDP2in:d3gdsdp2 withT:t andWithP:p];

        [self speciation:THIRD | FOURTH | NINTH | TENTH t:t p:p r:r s:NULL dr:NULL dt:dsdt dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:d2sdtdp dp2:d2sdp2];

        *dtdp = d3gdtdp2;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *dtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
            + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) {
                *dtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
                + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
                + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
                + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k])
                    *dtdp += [self fillD3GDS3atI:i andL:j andM:k withT:t andWithP:p]*dsdt[i]*dsdp[j]*dsdp[k];
            }
        }
    }

    if(mask & EIGHTH) {
        double d3gdp3 = [self fillD3GDP3withT:t andWithP:p];
        double d2gds2[NS][NS], d2gdsdp[NS], d3gds2dp[NS][NS],
        d3gdsdp2[NS], dsdp[NS], d2sdp2[NS];
        int k;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];
        [self fillD3GDSDP2in:d3gdsdp2 withT:t andWithP:p];

        [self speciation:FOURTH | TENTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:d2sdp2];

        *dp2 = d3gdp3;
        for (i=0; i<NS; i++) if(incOfxSpecies[i]) {
            *dp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
            for (j=0; j<NS; j++) if(incOfxSpecies[j]) {
                *dp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) if(incOfxSpecies[k])
                    *dp2 += [self fillD3GDS3atI:i andL:j andM:k withT:t andWithP:p]*dsdp[i]*dsdp[j]*dsdp[k];
            }
        }
    }

    if(mask & NINTH) {
        double d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
        d3gds2dp[NS][NS], d2gdrds[NR][NS], d3gdrdtdp[NR], d3gdsdtdp[NS],
        dsdt[NS], dsdp[NS], dsdr[NS][NR], d2sdrdt[NS][NR], d2sdrdp[NS][NR],
        d2gds2[NS][NS], d2gdsdt[NS], d3gdrdsdp[NR][NS], d2gdsdp[NS],
        d2sdtdp[NS], d3gds2dt[NS][NS];
        int k, l;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDTin:d2gdsdt withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDTin:d3gdrdsdt withT:t andWithP:p];
        [self fillD3GDRDSDPin:d3gdrdsdp withT:t andWithP:p];
        [self fillD3GDRDTDPin:d3gdrdtdp withT:t andWithP:p];
        [self fillD3GDS2DTin:d3gds2dt withT:t andWithP:p];
        [self fillD3GDSDTDPin:d3gdsdtdp withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];

        [self speciation:SECOND | THIRD | FOURTH | SIXTH | SEVENTH | NINTH t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:dsdp dr2:NULL drt:d2sdrdt drp:d2sdrdp dt2:NULL dtp:d2sdtdp dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0,drdt[i]=d3gdrdtdp[i]; j<NS; j++) if(incOfxSpecies[j]) {
                drdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
                d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
                d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
                for (k=0; k<NS; k++) if(incOfxSpecies[k]) {
                    drdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
                    d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
                    d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
                    d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
                    d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
                    d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l])
                        drdt[i] += [self fillD3GDS3atI:j andL:k andM:l withT:t andWithP:p]*dsdr[j][i]*dsdt[k]*dsdp[l];
                }
            }
        }
    }

    if(mask & TENTH) {
        double d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
        d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
        d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS];
        int k, l;

        [self DnSpeciesDsWith:s];
        [self D2nSpeciesDs2With:s];
        [self fillD2GDRDSin:d2gdrds withT:t andWithP:p];
        [self fillD2GDS2in:d2gds2 withT:t andWithP:p];
        [self fillD2GDSDPin:d2gdsdp withT:t andWithP:p];
        [self fillD3GDRDS2in:d3gdrds2 withT:t andWithP:p];
        [self fillD3GDRDSDPin:d3gdrdsdp withT:t andWithP:p];
        [self fillD3GDRDP2in:d3gdrdp2 withT:t andWithP:p];
        [self fillD3GDS2DPin:d3gds2dp withT:t andWithP:p];
        [self fillD3GDSDP2in:d3gdsdp2 withT:t andWithP:p];

        [self speciation:SECOND | FOURTH | SEVENTH | TENTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:NULL drt:NULL drp:d2sdrdp dt2:NULL dtp:NULL dp2:d2sdp2];

        for (i=0; i<NR; i++) {
            for (j=0,drdp[i]=d3gdrdp2[i]; j<NS; j++) if(incOfxSpecies[j]) {
                drdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
                2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
                d2gdsdp[j]*d2sdrdp[j][i];
                for (k=0; k<NS; k++) if(incOfxSpecies[k]) {
                    drdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
                    2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
                    2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
                    d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
                    for (l=0; l<NS; l++) if(incOfxSpecies[l])
                        drdp[i] += [self fillD3GDS3atI:j andL:k andM:l withT:t andWithP:p]*dsdr[j][i]*dsdp[k]*dsdp[l];
                }
            }
        }
    }
}

-(NSUInteger)numberOfSolutionSpecies {
    return NE;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
    id object = [endmembers objectAtIndex:index];
    if ([[object class] instancesRespondToSelector:@selector(phaseName)]) return [[endmembers objectAtIndex:index] phaseName];
    return [@"Species-" stringByAppendingFormat:@"%ld", index];
}

-(NSDictionary *)getSpeciesMoleFractionsForBulkComposition:(double *)m aT:(double)t andP:(double)p {
    double r[NR], s[NS];
    double rTotal = 0.0;
    int i;

    [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
    for (i=0; i<NR; i++) rTotal += r[i];
    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

    [self fillxSpeciesWith:r andWith:s];

    NSMutableDictionary *result = [NSMutableDictionary dictionaryWithCapacity:NE];
    for (i=0; i<NE; i++) [result setObject:[NSNumber numberWithDouble:xSpecies[i]] forKey:[self nameOfSolutionSpeciesAtIndex:i]];

    return [NSDictionary dictionaryWithDictionary:result];
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
    DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
    double *mComponents = [mComponentWrapper pointerToDouble];
    for (NSUInteger i=0; i<NA; i++) mComponents[i] = mSpecies[i];
    NSAssert(NO, @"DEWFluid:Function not yet implemented");
    return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *muSpeciesWrapper = [[DoubleVector alloc] initWithSize:NE];
    double *muSpecies = [muSpeciesWrapper pointerToDouble];

    double r[NR], s[NS];
    [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        // standard state (solvent and solute)
        muSpecies[i] = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
        // mixing properties
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        if (i == mIndH2O) {
            // solvent - osmotic term
            muSpecies[i] += dhMuExcess;
        } else {
            // solute - ideal term
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            muSpecies[i] += R*t*log(molality);
            // solute - actvity coefficient
            muSpecies[i] += dhMuExcess;
        }
    } else muSpecies[i] = 0.0;

    return muSpeciesWrapper;
}

-(DoubleVector *)activitiesOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *activitiesWrapper = [[DoubleVector alloc] initWithSize:NE];
    double *activities = [activitiesWrapper pointerToDouble];

    double r[NR], s[NS];
    [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
    [self speciation:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

    for (NSUInteger i=0; i<NE; i++) if (xSpecies[i] > 0.0) {
        activities[i] = 0.0;
        // mixing properties
        double dhMuExcess = [self fillDebyeHuckelExcessFreeEnergyForSpeciesIndex:i withT:t andWithP:p];
        if (i == mIndH2O) {
            // solvent - osmotic term
            activities[i] += dhMuExcess;
        } else {
            // solute - ideal term
            double molality = xSpecies[i]*CapGam/xSpecies[mIndH2O];
            activities[i] += R*t*log(molality);
            // solute - actvity coefficient
            activities[i] += dhMuExcess;
        }
        activities[i] = exp(activities[i]/R/t);
    } else activities[i] = 0.0;

    return activitiesWrapper;
}


-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
    DoubleVector *speciesElementArrayWrapper = [[endmembers objectAtIndex:index] formulaAsElementArray];
    return speciesElementArrayWrapper;
}

-(void)correctActivityCoefficients:(double [NA])gamma forComposition:(double [NA])x {
    for (NSUInteger i=0, j=0; i<NA; i++) if (x[i] != 0.0) {
        if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
        else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
        j++;
    }
}

#define NATOMS 1.0

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
    NSMutableArray *results = [NSMutableArray arrayWithCapacity:NA];
    double mu0[NA], deltaMu[NA], xNz[NA], x[NA], gamma[NA], xLast[NA], affinity = 0.0;
    NSUInteger i, j, nz = 0, index[NA];

    if (self.debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");

    // Compute solid -> liquid delta mus and deflate composition space
    for (i=0; i<NA; i++) {
        if (chemicalPotentials[i] != 0.0) {
            mu0[i] = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
            deltaMu[nz] = chemicalPotentials[i] - mu0[i];
            index[nz] = i;
            gamma[nz] = 1.0;
            nz++;
        }
        x[i] = 0.0;
        xLast[i] = 0.0;
    }

    // There are no non-zero chemical potentials, so the phase can never form.  Return an affinity of zero and zero the mole fractions
    if (nz == 0) {
        [results addObject:[NSNumber numberWithDouble:affinity]];
        for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];
        [results addObject:[NSNumber numberWithBool:YES]];          // convergence flag
        [results addObject:[NSNumber numberWithUnsignedInteger:0]]; // iteration count
        [results addObject:[NSNumber numberWithDouble:NATOMS]];     // number of atoms used to scale affinity
        [results addObject:[NSNumber numberWithDouble:0.0]];        // likely error in affinity

        if (self.debugV) NSLog(@"... Terminated. Trivial case.");
        return [NSArray arrayWithArray:results];
    }

    NSUInteger count = 0;
    BOOL converged = NO;
    double affinityLast, xReduced[NA];
    do {
        affinityLast = affinity;
        // Solve for mole ractions in the deflated composition space
        if (nz == 1) {
            xNz[0] = 1.0;
            affinity = -(deltaMu[0]-R*t*log(gamma[0]));
        } else {
            double sum = 1.0;
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                xNz[i] = exp(((deltaMu[i]-R*t*log(gamma[i]))-(deltaMu[nz-2]-R*t*log(gamma[nz-2])))/(R*t));
                sum += xNz[i];
            }
            xNz[nz-2] = exp(((deltaMu[nz-2]-R*t*log(gamma[nz-2]))-(deltaMu[nz-1]-R*t*log(gamma[nz-1])))/(R*t));

            xNz[nz-2] /= 1.0 + xNz[nz-2]*sum;
            xNz[nz-1] = 1.0 - xNz[nz-2];
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                xNz[i] *= xNz[nz-2];
                xNz[nz-1] -= xNz[i];
            }

            /* compute the chemical affinity (choice of mu[] is arbitrary) */
            affinity = -(deltaMu[0]-R*t*log(gamma[0])) + R*t*log(xNz[0]);
        }

        // Reinflate the solution
        for (i=0; i<nz; i++) x[index[i]] = xNz[i];

        // Determine activity coefficients
        double a[NA], mu[NA], r[NR];
        for (i=0; i<NA; i++) {
            NSString *qualifier = (x[i] < 0.0) ? @"Adjusted" : @"";
            xReduced[i] = (x[i] < 0.0) ? DBL_EPSILON : x[i];
            if (self.debugV && (xReduced[i] != 0.0)) NSLog(@"... X%2.2lu %@ = %g %@", i,
                [[self nameOfSolutionSpeciesAtIndex:i] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                xReduced[i], qualifier);
        }
        if (![self testPermissibleValuesOfComponents:xReduced]) {
            NSLog(@"Composition estimate is infeasible.");
#if TARGET_OS_IPHONE || TARGET_IPHONE_SIMULATOR
            for (i=0; i<NA; i++)
                NSLog(@"species X%2.2lu %@ = %g", (unsigned long)i,
                    [[self nameOfSolutionSpeciesAtIndex:i] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[i]);
#else
            for (i=0; i<NA; i++)
                NSLog(@"species X%2.2lu %@ = %g", i,
                    [[self nameOfSolutionSpeciesAtIndex:i] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[i]);
#endif
        }
        [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:xReduced r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
        [self activity:FIRST | SECOND t:t p:p r:r a:a mu:mu dx:NULL];

        for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];

        [self correctActivityCoefficients:gamma forComposition:x];
        converged = (fabs(affinity-affinityLast) < 0.1);
        count++;

    } while (count < 50 && !converged);

    for (i=0; i<NA; i++) if (xReduced[i] <= 1000.0*DBL_EPSILON) xReduced[i] = 0.0;

    if (self.debugS) {
#if TARGET_OS_IPHONE || TARGET_IPHONE_SIMULATOR
        NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
              converged ? @"YES" : @"NO", [self phaseName], (unsigned long)count, fabs(affinity-affinityLast), NATOMS);
#else
        NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
              converged ? @"YES" : @"NO", [self phaseName], count, fabs(affinity-affinityLast), NATOMS);
#endif
        for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
                                                         [[endmembers objectAtIndex:i] phaseName], gamma[j++], x[i]);
    }
    if (self.debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");


    [results addObject:[NSNumber numberWithDouble:affinity]];                         // affinity in J
    for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:xReduced[i]]]; // composition in mole fraction of endmembers
    [results addObject:[NSNumber numberWithBool:converged]];                          // convergence flag
    [results addObject:[NSNumber numberWithUnsignedInteger:count]];                   // iteration count
    [results addObject:[NSNumber numberWithDouble:NATOMS]];                           // number of atoms used to scale affinity
    [results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]];      // likely error in affinity

    return [NSArray arrayWithArray:results];
}

#pragma mark -
#pragma mark Solution Phase Protocol functions

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
        double *componentToElements = [[(HKFspeciesProperties *)[endmembers objectAtIndex:i] formulaAsElementArray] pointerToDouble];
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

    if (!computeMixingQuantities) return muWrapper;

    for (NSUInteger i=0; i<NA; i++) {
        id component = [endmembers objectAtIndex:i];
        mu[i] += [component getGibbsFreeEnergyFromT:t andP:p];
    }

    return muWrapper;
}

// --> SolutionPhaseProtocol public function
-(DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleMatrix *dadmWrapper = [[DoubleMatrix alloc] initWithRowSize:NA andWithColumnSize:NR];
    double **dadm = [dadmWrapper pointerToPointerToDouble];
    double activity[NA], dadr[NA][NR], r[NR], drdm[NR][NA];
    [self convert:SECOND outMask:THIRD | FIFTH t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:drdm d2m:NULL dr:NULL d3m:NULL];
    [self activity:FIRST | THIRD t:t p:p r:r a:activity mu:NULL dx:dadr];
    for (NSUInteger i=0; i<NA; i++) {
        [self convertIntensiveToExtensiveGradients:activity[i] dpMix:dadr[i] mToTal:[self totalMolesFromMolesOfComponents:m] drdm:drdm dp:dadm[i]];
    }

    return dadmWrapper;
}

// --> SolutionPhaseProtocol public function
-(double)getGibbsFreeEnergyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    double g, r[NR];
    [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
    [self gmix:FIRST t:t p:p r:r gmix:&g dx:NULL dx2:NULL dx3:NULL];
    g *= [self totalMolesFromMolesOfComponents:m];

    if (!computeMixingQuantities) return g;

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

    if (!computeMixingQuantities) return dgdmWrapper;

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

    if (!computeMixingQuantities) return h;

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

    if (!computeMixingQuantities) return s;

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

    if (!computeMixingQuantities) return dsdmWrapper;

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

    if (!computeMixingQuantities) return cp;

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

    if (!computeMixingQuantities) return dcpdt;

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

    if (!computeMixingQuantities) return dcpdmWrapper;

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

    if (!computeMixingQuantities) return v;

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

    if (!computeMixingQuantities) return dvdmWrapper;

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

    if (!computeMixingQuantities) return dvdt;

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

    if (!computeMixingQuantities) return dvdp;

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

    if (!computeMixingQuantities) return d2vdt2;

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

    if (!computeMixingQuantities) return d2vdtdp;

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

    if (!computeMixingQuantities) return d2vdp2;

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

    if (!computeMixingQuantities) return dvdmdtWrapper;

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

    if (!computeMixingQuantities) return dvdmdpWrapper;

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

@end
