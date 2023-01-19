//
//  main.c
//  PureCtestProgram
//
//  Created by Mark Ghiorso on 11/7/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#include <stdio.h>

extern int testDEW (int mask, double t, double p, int na, int nr, char **names, char **formulas, double *r, double *m);
extern void conDEW (int inpMask, int outMask, double t, double p, double *e, double *m, double *r, double *x,
                    double **dm, double ***d2m, double **dr, double ****d3m);
extern void actDEW (int mask, double t, double p, double *x, double *a, double *mu, double **dx);
extern void gmixDEW (int mask, double t, double p, double *x, double *gmix, double *dx, double **dx2, double ***dx3);
extern void hmixDEW (int mask, double t, double p, double *x, double *hmix);
extern void smixDEW (int mask, double t, double p, double *x, double *smix, double *dx, double **dx2);
extern void cpmixDEW (int mask, double t, double p, double *x, double *cpmix, double *dt, double *dx);
extern void vmixDEW (int mask, double t, double p, double *x, double *vmix, double *dx, double **dx2,
                     double *dt, double *dp, double *dt2, double *dtdp, double *dp2, double *dxdt, double *dxdp);
extern void dispDEW (int mask, double t, double p, double *x, char **formula);

int main(int argc, const char * argv[]) {
    // insert code here...
    printf("Hello, World!\n");
    return 0;
}
