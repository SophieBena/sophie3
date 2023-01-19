//
//  DEWlibraryInterface.m
//  DEWTestProgram
//
//  Created by Mark Ghiorso on 11/7/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "DEWFluid.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"

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

#define NR      16      // Number of independent mole fraction variables
#define NA      17      // Number of liquid components

@interface DEWFluidWrapper : NSObject
@property (strong) DEWFluid *dewFluid;
@property (strong) DoubleVector *molesWrapper;
@end

@implementation DEWFluidWrapper

-(instancetype)init {
    if ((self = [super init])) {
        _dewFluid = [[DEWFluid alloc] init];
        _molesWrapper = [[DoubleVector alloc] initWithSize:NA andInitialValue:0.0];
    }
    return self;
}

@end

static DEWFluidWrapper *dewFluidWrapper;
static BOOL dewFluidInitialized = NO;

/* returns TRUE if values are correct/within bounds, else returns FALSE */
int testDEW (int mask, double t, double p,
             int    na,       /* # of components in solution     BINARY MASK: 000001 */
             int    nr,       /* # of indep compos variables     BINARY MASK: 000010 */
             char **names,    /* names compon, expected order    BINARY MASK: 000100 */
             char **formulas, /* form of compon, expected order  BINARY MASK: 001000 */
             double *r,       /* indep compositional variables   BINARY MASK: 010000 */
             double *m        /* moles of endmember components   BINARY MASK: 100000 */
) {                           /* r[] and m[] are tested for bound constraints        */
    if (!dewFluidInitialized) { dewFluidWrapper = [[DEWFluidWrapper alloc] init]; dewFluidInitialized = YES; }
    return 0;
}

void conDEW (int inpMask, int outMask, double t, double p,
             double *e,     /* moles of elements               BINARY MASK: 00000001 */
             double *m,     /* moles of endmember components   BINARY MASK: 00000010 */
             double *r,     /* indep compositional variables   BINARY MASK: 00000100 */
             double *x,     /* mole fractions of endmember cmp BINARY MASK: 00001000 */
             double **dm,   /* matrix[i][j]: dr[i]/dm[j]       BINARY MASK: 00010000 */
             double ***d2m, /* cube[i][j][k]: d2r[i]/dm[j]dm[k]BINARY MASK: 00100000 */
             double **dr,   /* matrix[i][j]: dx[i]/dr[j]       BINARY MASK: 01000000 */
             double ****d3m /* 4d[i][j][k][l]: d3r[i]/dm[j]dm[k]dm[l] MASK: 10000000 */
) {
    if (!dewFluidInitialized) { dewFluidWrapper = [[DEWFluidWrapper alloc] init]; dewFluidInitialized = YES; }
}

void actDEW (int mask, double t, double p, double *x,
             double *a,   /* (pointer to a[]) activities           BINARY MASK: 0001 */
             double *mu,  /* (pointer to mu[]) chemical potentials BINARY MASK: 0010 */
             double **dx  /* (pointer to dx[][]) d(a[])/d(x[])     BINARY MASK: 0100 */
) {                       /* exclusion applied to activities if:   BINARY MASK: 1000 */
    double *moles = [dewFluidWrapper.molesWrapper pointerToDouble];
    conDEW(THIRD, SECOND, t, p, NULL, moles, x, NULL, NULL, NULL, NULL, NULL);
    if (mask & FIRST) {
        DoubleVector *activityWrapper = [dewFluidWrapper.dewFluid getActivityFromMolesOfComponents:moles andT:t andP:p];
        double *activity = [activityWrapper pointerToDouble];
        for (NSUInteger i=0; i<NA; i++) a[i] = activity[i];
    }
    if (mask & SECOND) {
        DoubleVector *chemicalPotentialWrapper = [dewFluidWrapper.dewFluid getChemicalPotentialFromMolesOfComponents:moles andT:t andP:p];
        double *chemicalPotential = [chemicalPotentialWrapper pointerToDouble];
        for (NSUInteger i=0; i<NA; i++) mu[i] = chemicalPotential[i];
    }
}

void gmixDEW (int mask, double t, double p, double *x,
              double *gmix,  /* Gibbs energy of mixing              BINARY MASK: 0001 */
              double *dx,    /* (pointer to dx[]) d(g)/d(x[])       BINARY MASK: 0010 */
              double **dx2,  /* (pointer to dx2[][]) d2(g)/d(x[])2  BINARY MASK: 0100 */
              double ***dx3  /* (pointer to dx3[][][]) d3(g)/d(x[])3BINARY MASK: 1000 */
) {
    if (!dewFluidInitialized) { dewFluidWrapper = [[DEWFluidWrapper alloc] init]; dewFluidInitialized = YES; }
}

void hmixDEW (int mask, double t, double p, double *x,
              double *hmix  /* Enthalpy of mixing                      BINARY MASK: 1 */
) {

}

void smixDEW (int mask, double t, double p, double *x,
              double *smix,  /* Entropy of mixing                    BINARY MASK: 001 */
              double *dx,    /* (pointer to dx[]) d(s)/d(x[])        BINARY MASK: 010 */
              double **dx2   /* (pointer to dx2[][]) d2(s)/d(x[])2   BINARY MASK: 100 */
) {
    if (!dewFluidInitialized) { dewFluidWrapper = [[DEWFluidWrapper alloc] init]; dewFluidInitialized = YES; }
}

void cpmixDEW (int mask, double t, double p, double *x,
               double *cpmix, /* Heat capacity of mixing              BINARY MASK: 001 */
               double *dt,    /* d(cp)/d(t)                           BINARY MASK: 010 */
               double *dx     /* d(cp)/d(x[])                         BINARY MASK: 100 */
) {

}

void vmixDEW (int mask, double t, double p, double *x,
              double *vmix, /* Volume of mixing               BINARY MASK: 0000000001 */
              double *dx,   /* (pointer to dx[]) d(v)/d(x[])  BINARY MASK: 0000000010 */
              double **dx2, /* (point dx2[][]) d(v)/d(x[])2   BINARY MASK: 0000000100 */
              double *dt,   /* d(v)/d(t)                      BINARY MASK: 0000001000 */
              double *dp,   /* d(v)/d(p)                      BINARY MASK: 0000010000 */
              double *dt2,  /* d2(v)/d(t)2                    BINARY MASK: 0000100000 */
              double *dtdp, /* d2(v)/d(t)d(p)                 BINARY MASK: 0001000000 */
              double *dp2,  /* d2(v)/d(p)2                    BINARY MASK: 0010000000 */
              double *dxdt, /* d2(v)/d(x[])d(t)               BINARY MASK: 0100000000 */
              double *dxdp  /* d2(v)/d(x[])d(p)               BINARY MASK: 1000000000 */
) {
    if (!dewFluidInitialized) { dewFluidWrapper = [[DEWFluidWrapper alloc] init]; dewFluidInitialized = YES; }
}

void dispDEW (int mask, double t, double p, double *x,
              char **formula /* Mineral formula for interface display  BINARY MASK: 1 */
) {
    if (!dewFluidInitialized) { dewFluidWrapper = [[DEWFluidWrapper alloc] init]; dewFluidInitialized = YES; }
}
