//
//  ZhangDuan2009.m
//  ThermoFit
//
//  Created by Mark Ghiorso on 4/18/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "ZhangDuan2009.h"

#pragma mark -
#pragma mark C preprocessor defines

#define R 8.3143

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#pragma mark -
#pragma mark static constants

/* pure EOS terms for 0.1 to 10 GPa */
static const double H2Oa1  =  2.95177298930e-002;
static const double H2Oa2  = -6.33756452413e+003;
static const double H2Oa3  = -2.75265428882e+005;
static const double H2Oa4  =  1.29128089283e-003;
static const double H2Oa5  = -1.45797416153e+002;
static const double H2Oa6  =  7.65938947237e+004;
static const double H2Oa7  =  2.58661493537e-006;
static const double H2Oa8  =  0.52126532146e+000;
static const double H2Oa9  = -1.39839523753e+002;
static const double H2Oa10 = -2.36335007175e-008;
static const double H2Oa11 =  5.35026383543e-003;
static const double H2Oa12 = -0.27110649951e+000;
static const double H2Oa13 =  2.50387836486e+004;
static const double H2Oa14 =  0.73226726041e+000;
static const double H2Oa15 =  1.54833359970e-002;
static const double H2Ogam =  H2Oa15;

/* H2O, critical constants, K, bars, J/bar */

static const double H2Oepsilon = 510.00;
static const double H2Osigma   =   2.88;

static const double idealCoeff[13] = {
    3.10409601236035e+01,
    -3.91422080460869e+01,
    3.79695277233575e+01,
    -2.18374910952284e+01,
    7.42251494566339e+00,
    -1.38178929609470e+00,
    1.08807067571454e-01,
    -1.20771176848589e+01,
    3.39105078851732e+00,
    -5.84520979955060e-01,
    5.89930846488082e-02,
    -3.12970001415882e-03,
    6.57460740981757e-05
};

#pragma mark -
#pragma mark Virial functions

static void BVcAndDerivative(double t, double *bv, double *dbvdt, double *d2bvdt2, double *d3bvdt3) {
    double H2OTr    = 154.0*t/H2Oepsilon;
    double dH2OTrdt = 154.0/H2Oepsilon;

    double bEnd = H2Oa1 + H2Oa2/H2OTr/H2OTr + H2Oa3/H2OTr/H2OTr/H2OTr;
    double dbEnddt = - 2.0*H2Oa2*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa3*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2bEnddt2 = 6.0*H2Oa2*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa3*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3bEnddt3 = - 24.0*H2Oa2*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr
    - 60.0*H2Oa3*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *bv = bEnd;
    *dbvdt   = dbEnddt;
    *d2bvdt2 = d2bEnddt2;
    *d3bvdt3 = d3bEnddt3;

    return;
}

static void CVcAndDerivative(double t, double *cv, double *dcvdt, double *d2cvdt2, double *d3cvdt3) {
    double H2OTr    = 154.0*t/H2Oepsilon;
    double dH2OTrdt = 154.0/H2Oepsilon;

    double cEnd = H2Oa4 + H2Oa5/H2OTr/H2OTr + H2Oa6/H2OTr/H2OTr/H2OTr;
    double dcEnddt = - 2.0*H2Oa5*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa6*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2cEnddt2 = 6.0*H2Oa5*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa6*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3cEnddt3 = -24.0*H2Oa5*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2Oa6*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *cv = cEnd;
    *dcvdt   = dcEnddt;
    *d2cvdt2 = d2cEnddt2;
    *d3cvdt3 = d3cEnddt3;

    return;
}

static void DVcAndDerivative(double t, double *dv, double *ddvdt, double *d2dvdt2, double *d3dvdt3) {
    double H2OTr    = 154.0*t/H2Oepsilon;
    double dH2OTrdt = 154.0/H2Oepsilon;

    double dEnd = H2Oa7 + H2Oa8/H2OTr/H2OTr + H2Oa9/H2OTr/H2OTr/H2OTr;
    double ddEnddt = - 2.0*H2Oa8*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa9*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2dEnddt2 = 6.0*H2Oa8*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa9*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3dEnddt3 = - 24.0*H2Oa8*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2Oa9*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *dv      = dEnd;
    *ddvdt   = ddEnddt;
    *d2dvdt2 = d2dEnddt2;
    *d3dvdt3 = d3dEnddt3;

    return;
}

static void EVcAndDerivative(double t, double *ev, double *devdt, double *d2evdt2, double *d3evdt3) {
    double H2OTr    = 154.0*t/H2Oepsilon;
    double dH2OTrdt = 154.0/H2Oepsilon;

    double eEnd = H2Oa10 + H2Oa11/H2OTr/H2OTr + H2Oa12/H2OTr/H2OTr/H2OTr;
    double deEnddt = - 2.0*H2Oa11*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa12*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2eEnddt2 = 6.0*H2Oa11*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa12*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3eEnddt3 = - 24.0*H2Oa11*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr
    - 60.0*H2Oa12*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *ev      = eEnd;
    *devdt   = deEnddt;
    *d2evdt2 = d2eEnddt2;
    *d3evdt3 = d3eEnddt3;

    return;
}

static void FVcAndDerivative(double t, double *fv, double *dfvdt, double *d2fvdt2, double *d3fvdt3) {
    double H2OTr    = 154.0*t/H2Oepsilon;
    double dH2OTrdt = 154.0/H2Oepsilon;

    double fEnd      = H2Oa13*H2Oa14/H2OTr/H2OTr/H2OTr;
    double dfEnddt   = -  3.0*H2Oa13*H2Oa14*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2fEnddt2 =   12.0*H2Oa13*H2Oa14*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3fEnddt3 = - 60.0*H2Oa13*H2Oa14*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *fv      = fEnd;
    *dfvdt   = dfEnddt;
    *d2fvdt2 = d2fEnddt2;
    *d3fvdt3 = d3fEnddt3;

    return;
}

static void GVcAndDerivative(double t, double *gv, double *dgvdt, double *d2gvdt2, double *d3gvdt3) {
    double H2OTr    = 154.0*t/H2Oepsilon;
    double dH2OTrdt = 154.0/H2Oepsilon;

    double gEnd      = H2Oa13*H2Oa15/H2OTr/H2OTr/H2OTr;
    double dgEnddt   = -  3.0*H2Oa13*H2Oa15*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2gEnddt2 =   12.0*H2Oa13*H2Oa15*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3gEnddt3 = - 60.0*H2Oa13*H2Oa15*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *gv      = gEnd;
    *dgvdt   = dgEnddt;
    *d2gvdt2 = d2gEnddt2;
    *d3gvdt3 = d3gEnddt3;

    return;
}

static void GammaVcAndDerivative(double t, double *gammav) {
    *gammav = H2Ogam;

    return;
}

#pragma mark -
#pragma mark Ideal gas functions

static void idealGasH2O(double t, double *cp, double *s0, double *h0, double *dcpdt) {
    int i;

    for (i=0, *cp=0.0; i<7; i++) *cp += idealCoeff[i]*pow(t/1000.0, (double) i);
    for (i=7; i<13; i++)         *cp += idealCoeff[i]/pow(t/1000.0, (double) (i-6));

    for (i=1, *dcpdt=0.0; i<7; i++) *dcpdt +=  ((double) i)  *idealCoeff[i]*pow(t/1000.0, (double) i-1);
    for (i=7; i<13; i++)            *dcpdt += -((double) i-6)*idealCoeff[i]/pow(t/1000.0, (double) (i+1-6));

    for (i=0, *h0=0.0; i<7; i++) *h0 += idealCoeff[i]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    *h0 += idealCoeff[7]*log(t/1000.0);
    for (i=8; i<13; i++)         *h0 += idealCoeff[i]/pow(t/1000.0, (double) (i-7))/((double) (7-i));

    *s0  = idealCoeff[0]*log(t/1000.0);
    for (i=1; i<7; i++)	       *s0 += idealCoeff[i]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++)         *s0 += idealCoeff[i]/pow(t/1000.0, (double) (i-6))/((double) (6-i));

    *cp    *= 8.31451;
    *h0    *= 8.31451*1000.0;
    *s0    *= 8.31451;
    *dcpdt *= 8.31451/1000.0;

    *h0 += - 355665.4136;
    *s0 +=   359.6505;
}

#pragma mark -
#pragma mark Duan driver routine (pure water)

static int ZD2009H2ODriver(double t, double p,
                         double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                         double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {
    double bv, cv, dv, ev, fv, gv, gammav, v, z = 1.0, dzdv, dzdt, d2zdv2, d2zdvdt, d2zdt2;
    double dbvdt, dcvdt, ddvdt, devdt, dfvdt, dgvdt;
    double d2bvdt2, d2cvdt2, d2dvdt2, d2evdt2, d2fvdt2, d2gvdt2;
    double d3bvdt3, d3cvdt3, d3dvdt3, d3evdt3, d3fvdt3, d3gvdt3;
    double lnPhiH2O, dlnPhiH2Odv, dlnPhiH2Odt, d2lnPhiH2Odv2, d2lnPhiH2Odt2, d2lnPhiH2Odvdt;

    BVcAndDerivative    (t, &bv,     &dbvdt,    &d2bvdt2,  &d3bvdt3);
    CVcAndDerivative    (t, &cv,     &dcvdt,    &d2cvdt2,  &d3cvdt3);
    DVcAndDerivative    (t, &dv,     &ddvdt,    &d2dvdt2,  &d3dvdt3);
    EVcAndDerivative    (t, &ev,     &devdt,    &d2evdt2,  &d3evdt3);
    FVcAndDerivative    (t, &fv,     &dfvdt,    &d2fvdt2,  &d3fvdt3);
    GVcAndDerivative    (t, &gv,     &dgvdt,    &d2gvdt2,  &d3gvdt3);
    GammaVcAndDerivative(t, &gammav);

    double pm, tm, vm;
    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0, vLowest = 0.017;

        v = 0.08314467*t/p;
        pm = 3.063*pow(H2Osigma, 3.0)*p/H2Oepsilon;
        tm = 154.0*t/H2Oepsilon;
        vm = v*pow(3.691/H2Osigma, 3.0);
        z = 1.0 + bv/vm + cv/vm/vm + dv/vm/vm/vm/vm + ev/vm/vm/vm/vm/vm + (fv/vm/vm + gv/vm/vm/vm/vm) * exp(-gammav/vm/vm);
        while ((z < 0.0) && (iter < 200)) {
            double dzdvm = - bv/vm/vm - 2.0*cv/vm/vm/vm + - 4.0*dv/vm/vm/vm/vm/vm - 5.0*ev/vm/vm/vm/vm/vm/vm
            - 2.0*(fv/vm/vm/vm) * exp(-gammav/vm/vm) + 2.0*(fv/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
            - 4.0*(gv/vm/vm/vm/vm/vm) * exp(-gammav/vm/vm) + 2.0*(gv/vm/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm);
            vm += (0.1-z)/dzdvm;
            z = 1.0 + bv/vm + cv/vm/vm + dv/vm/vm/vm/vm + ev/vm/vm/vm/vm/vm + (fv/vm/vm + gv/vm/vm/vm/vm) * exp(-gammav/vm/vm);
            iter++;
        }
        if (z < 0.0) return FALSE;

        iter = 0;
        while (iter < 200) {
            delv = z*0.08314467*tm/pm - vm;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < vm*100.0*DBL_EPSILON) ) break;
            vPrevious = vm;
            delvPrevious = delv;
            vm = (z*0.08314467*tm/pm + vm)/2.0;
            double increment = vm - vPrevious;
            while ((z = 1.0 + bv/vm + cv/vm/vm + dv/vm/vm/vm/vm + ev/vm/vm/vm/vm/vm + (fv/vm/vm + gv/vm/vm/vm/vm) * exp(-gammav/vm/vm)) < 0.0) {
                increment /= 2.0;
                if (fabs(increment) < 100.0*DBL_EPSILON) return FALSE;
                vm = vPrevious + increment;
                printf("Z&D2009: t = %.1f, p = %.1f, z = %g, vm = %g, vInc = %g\n", t, p, z, vm, increment);
            }
            if (vm < vLowest) vm = vLowest;
            iter++;
            if (iter == 200) {
                printf("Z&D2009: t = %.1f, p = %.1f, z = %g, vm = %g, vmPrev = %g, delvm = %g, delvmPrev = %g, vLow = %g\n",
                       t, p, z, vm, vPrevious, delv, delvPrevious, vLowest);
                if (fabs(vm-vLowest) <= 100.0*DBL_EPSILON) {
                    vLowest *= 0.90;
                    iter = 0;
                }
            }
        }
        if (fabs(delv) > vm*100.0*DBL_EPSILON) {
            double dx;
            double rtb = (delv < 0.0) ? (dx = vPrevious-vm,vm) : (dx = vm-vPrevious,vPrevious);
            iter = 0;
            while (iter < 200) {
                vm = rtb + (dx *= 0.5);
                z  = 1.0 + bv/vm + cv/vm/vm + dv/vm/vm/vm/vm + ev/vm/vm/vm/vm/vm + (fv/vm/vm + gv/vm/vm/vm/vm) * exp(-gammav/vm/vm);
                delv = z*0.08314467*tm/pm - vm;
                if (delv <= 0.0) rtb = vm;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
            if ( (iter == 200) || (fabs(dx) > 100.0*DBL_EPSILON) ) {
                printf("dewH2ODriver: z = %g, vm = %g, delvm = %g, iter = %d\n", z, vm, delv, iter);
                return FALSE;
            }
        } else if (iter == 200) return FALSE;

        v = 1000.0*vm*pow(H2Osigma/3.691, 3.0)/10.0;        // J/bar
    }
    double DvmDv = (1000.0)*pow(3.691/H2Osigma, 3.0)/10.0;  // J/bar

    lnPhiH2O  = 0.0;
    lnPhiH2O += -log(z) + z - 1.0;
    lnPhiH2O += bv/vm;
    lnPhiH2O += cv/2.0/vm/vm;
    lnPhiH2O += dv/4.0/vm/vm/vm/vm;
    lnPhiH2O += ev/5.0/vm/vm/vm/vm/vm;
    lnPhiH2O += (fv/2.0/gammav)*(1.0-exp(-gammav/vm/vm));
    lnPhiH2O += (gv/2.0/gammav/gammav)*(1.0 - (gammav/vm/vm+1.0)*exp(-gammav/vm/vm));

    // z = 1.0 + bv/vm + cv/vm/vm + dv/vm/vm/vm/vm + ev/vm/vm/vm/vm/vm + (fv/vm/vm + gv/vm/vm/vm/vm) * exp(-gammav/vm/vm);
    dzdv = - bv/vm/vm - 2.0*cv/vm/vm/vm + - 4.0*dv/vm/vm/vm/vm/vm - 5.0*ev/vm/vm/vm/vm/vm/vm
    - 2.0*(fv/vm/vm/vm) * exp(-gammav/vm/vm) + 2.0*(fv/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
    - 4.0*(gv/vm/vm/vm/vm/vm) * exp(-gammav/vm/vm) + 2.0*(gv/vm/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm);
    dzdt = dbvdt/vm + dcvdt/vm/vm + ddvdt/vm/vm/vm/vm + devdt/vm/vm/vm/vm/vm + (dfvdt/vm/vm + dgvdt/vm/vm/vm/vm) * exp(-gammav/vm/vm);

    dlnPhiH2Odv  = 0.0;
    dlnPhiH2Odv += -dzdv/z + dzdv;
    dlnPhiH2Odv += -bv/vm/vm;
    dlnPhiH2Odv += -cv/vm/vm/vm;
    dlnPhiH2Odv += -dv/vm/vm/vm/vm/vm;
    dlnPhiH2Odv += -ev/vm/vm/vm/vm/vm/vm;
    dlnPhiH2Odv += -(fv/vm/vm/vm)*exp(-gammav/vm/vm);
    dlnPhiH2Odv += -(gv/vm/vm/vm/vm/vm)*exp(-gammav/vm/vm);

    dlnPhiH2Odt  = 0.0;
    dlnPhiH2Odt += -dzdt/z + dzdt;
    dlnPhiH2Odt += dbvdt/vm;
    dlnPhiH2Odt += dcvdt/2.0/vm/vm;
    dlnPhiH2Odt += ddvdt/4.0/vm/vm/vm/vm;
    dlnPhiH2Odt += devdt/5.0/vm/vm/vm/vm/vm;
    dlnPhiH2Odt += (dfvdt/2.0/gammav)*(1.0-exp(-gammav/vm/vm));
    dlnPhiH2Odt += (dgvdt/2.0/gammav/gammav)*(1.0 - (gammav/vm/vm+1.0)*exp(-gammav/vm/vm));

    d2zdv2 =  2.0*bv/vm/vm/vm + 6.0*cv/vm/vm/vm/vm + 20.0*dv/vm/vm/vm/vm/vm/vm + 30.0*ev/vm/vm/vm/vm/vm/vm/vm
    + 6.0*(fv/vm/vm/vm/vm) * exp(-gammav/vm/vm) - 4.0*(fv/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
    - 4.0*(fv/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
    - 6.0*(fv/vm/vm) * (gammav/vm/vm/vm/vm) * exp(-gammav/vm/vm)
    + 4.0*(fv/vm/vm) * (gammav*gammav/vm/vm/vm/vm/vm/vm) * exp(-gammav/vm/vm)
    + 20.0*(gv/vm/vm/vm/vm/vm/vm) * exp(-gammav/vm/vm) - 8.0*(gv/vm/vm/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
    - 8.0*(gv/vm/vm/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
    - 6.0*(gv/vm/vm/vm/vm) * (gammav/vm/vm/vm/vm) * exp(-gammav/vm/vm)
    + 4.0*(gv/vm/vm/vm/vm) * (gammav/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm);

    d2zdvdt = - dbvdt/vm/vm - 2.0*dcvdt/vm/vm/vm + - 4.0*ddvdt/vm/vm/vm/vm/vm - 5.0*devdt/vm/vm/vm/vm/vm/vm
    - 2.0*(dfvdt/vm/vm/vm) * exp(-gammav/vm/vm) + 2.0*(dfvdt/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm)
    - 4.0*(dgvdt/vm/vm/vm/vm/vm) * exp(-gammav/vm/vm) + 2.0*(dgvdt/vm/vm/vm/vm) * (gammav/vm/vm/vm) * exp(-gammav/vm/vm);
    d2zdt2 = d2bvdt2/vm + d2cvdt2/vm/vm + d2dvdt2/vm/vm/vm/vm + d2evdt2/vm/vm/vm/vm/vm + (d2fvdt2/vm/vm + d2gvdt2/vm/vm/vm/vm) * exp(-gammav/vm/vm);

    d2lnPhiH2Odv2  = 0.0;
    d2lnPhiH2Odv2 += dzdv*dzdv/z/z - d2zdv2/z + d2zdv2;
    d2lnPhiH2Odv2 +=  2.0*bv/vm/vm/vm;
    d2lnPhiH2Odv2 +=  3.0*cv/vm/vm/vm/vm;
    d2lnPhiH2Odv2 +=  5.0*dv/vm/vm/vm/vm/vm/vm;
    d2lnPhiH2Odv2 +=  6.0*ev/vm/vm/vm/vm/vm/vm/vm;
    d2lnPhiH2Odv2 +=  3.0*(fv/vm/vm/vm/vm)*exp(-gammav/vm/vm) - 2.0*(fv/vm/vm/vm)*(gammav/vm/vm/vm)*exp(-gammav/vm/vm);
    d2lnPhiH2Odv2 +=  5.0*(gv/vm/vm/vm/vm/vm/vm)*exp(-gammav/vm/vm) - 2.0*(gv/vm/vm/vm/vm/vm)*(gammav/vm/vm/vm)*exp(-gammav/vm/vm);

    d2lnPhiH2Odvdt  = 0.0;
    d2lnPhiH2Odvdt += dzdv*dzdt/z/z -d2zdvdt/z + d2zdvdt;
    d2lnPhiH2Odvdt += -dbvdt/vm/vm;
    d2lnPhiH2Odvdt += -dcvdt/vm/vm/vm;
    d2lnPhiH2Odvdt += -ddvdt/vm/vm/vm/vm/vm;
    d2lnPhiH2Odvdt += -devdt/vm/vm/vm/vm/vm/vm;
    d2lnPhiH2Odvdt += -(dfvdt/vm/vm/vm)*exp(-gammav/vm/vm);
    d2lnPhiH2Odvdt += -(dgvdt/vm/vm/vm/vm/vm)*exp(-gammav/vm/vm);

    d2lnPhiH2Odt2  = 0.0;
    d2lnPhiH2Odt2 += dzdt*dzdt/z/z - d2zdt2/z + d2zdt2;
    d2lnPhiH2Odt2 += d2bvdt2/vm;
    d2lnPhiH2Odt2 += d2cvdt2/2.0/vm/vm;
    d2lnPhiH2Odt2 += d2dvdt2/4.0/vm/vm/vm/vm;
    d2lnPhiH2Odt2 += d2evdt2/5.0/vm/vm/vm/vm/vm;
    d2lnPhiH2Odt2 += (d2fvdt2/2.0/gammav)*(1.0-exp(-gammav/vm/vm));
    d2lnPhiH2Odt2 += (d2gvdt2/2.0/gammav/gammav)*(1.0 - (gammav/vm/vm+1.0)*exp(-gammav/vm/vm));

    dzdv           *= DvmDv;
    dlnPhiH2Odv    *= DvmDv;
    d2zdv2         *= DvmDv*DvmDv;
    d2zdvdt        *= DvmDv;
    d2lnPhiH2Odv2  *= DvmDv*DvmDv;
    d2lnPhiH2Odvdt *= DvmDv;

    *vPt      = v;
    *zPt      = z;
    *phi      = exp(lnPhiH2O);
    *dvdp     = 1.0/( p*(dzdv/z - 1.0/v) );
    *dvdt     = (1.0/t + dzdt/z)/(1.0/v - dzdv/z);
    *d2vdp2   = p*(1.0/v-dzdv/z)/v - dzdv/(*dvdp)/z + 1.0/(*dvdp)/(*dvdp)/p + p*d2zdv2/z;
    *d2vdp2  *= -(*dvdp)*(*dvdp)*(*dvdp);
    *d2vdtdp  = -(*dvdp)*(1.0/t + dzdt/z) - p*(*d2vdp2)*(1.0/t + dzdt/z) - p*(*dvdp)*(*dvdp)*(-dzdv*dzdt/z/z + d2zdvdt/z);
    *d2vdt2   = -p*(*d2vdtdp)*(1.0/t + dzdt/z) + p*(*dvdp)/t/t + p*(*dvdp)*dzdt*(dzdv*(*dvdt) + dzdt)/z/z
    - p*(*dvdp)*(*dvdt)*d2zdvdt/z - p*(*dvdp)*d2zdt2/z;

    *dlnphidt    = dlnPhiH2Odv*(*dvdt) + dlnPhiH2Odt;
    *dlnphidp    = dlnPhiH2Odv*(*dvdp);
    *d2lnphidt2  = d2lnPhiH2Odv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiH2Odvdt*(*dvdt) + dlnPhiH2Odv*(*d2vdt2) + d2lnPhiH2Odt2;
    *d2lnphidtdp = d2lnPhiH2Odv2*(*dvdt)*(*dvdp) + dlnPhiH2Odv*(*d2vdtdp) + d2lnPhiH2Odvdt*(*dvdp);
    *d2lnphidp2  = d2lnPhiH2Odv2*(*dvdp)*(*dvdp) + dlnPhiH2Odv*(*d2vdp2);

    return TRUE;
}

#pragma mark -
#pragma mark Duan pure H2O

static void propertiesOfPureH2O(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt,
                                double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2) {
    double z, phi, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2;

    idealGasH2O(t, cp, s, h, dcpdt);
    if (ZD2009H2ODriver(t, p, v, &z, &phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, &dlnphidt, &dlnphidp, &d2lnphidt2, &d2lnphidtdp, &d2lnphidp2)) {
        *g      = *h - t*(*s) + R*t*log(phi*p);
        *s     += - (R*log(phi*p) + R*t*dlnphidt);
        *h     += R*t*log(phi*p) - t*(R*log(phi*p) + R*t*dlnphidt);
        *cp    += -t*(2.0*R*dlnphidt + R*t*d2lnphidt2);
        {
            double zTemp, phiTemp, dlnphidtTemp, dlnphidpTemp, d2lnphidt2Temp, d2lnphidtdpTemp, d2lnphidp2Temp,
            vTemp, dvdtTemp, dvdpTemp, d2vdt2Temp, d2vdtdpTemp, d2vdp2Temp, d3lnphidt3;

            ZD2009H2ODriver(t*(1.0+sqrt(DBL_EPSILON)), p, &vTemp, &zTemp, &phiTemp, &dvdtTemp, &dvdpTemp, &d2vdt2Temp, &d2vdtdpTemp, &d2vdp2Temp,
                          &dlnphidtTemp, &dlnphidpTemp, &d2lnphidt2Temp, &d2lnphidtdpTemp, &d2lnphidp2Temp);

            d3lnphidt3 = (d2lnphidt2Temp - d2lnphidt2)/t/sqrt(DBL_EPSILON);
            *dcpdt += -(2.0*R*dlnphidt + R*t*d2lnphidt2) -t*(3.0*R*d2lnphidt2 + R*t*d3lnphidt3);
        }
    } else {
        *g       = 0.0;
        *h       = 0.0;
        *s       = 0.0;
        *cp      = 0.0;
        *dcpdt   = 0.0;
        *v       = 1.0;
        *dvdt    = 0.0;
        *dvdp    = 0.0;
        *d2vdt2  = 0.0;
        *d2vdtdp = 0.0;
        *d2vdp2  = 0.0;
    }
}

#pragma mark -
#pragma mark spline functions from Numerical Recipies

//static void spline(double x[], double y[], int n, double yp1, double ypn, double y2[], double u[]) {
//    int i,k;
//    double p, qn, sig, un;
//
//    if (yp1 > 0.99e30)
//        y2[0] = u[0] = 0.0;
//    else {
//        y2[0] = -0.5;
//        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
//    }
//
//    for (i=1; i<(n-1); i++) {
//        sig   = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
//        p     = sig*y2[i-1] + 2.0;
//        y2[i] = (sig-1.0)/p;
//        u[i]  = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
//        u[i]  = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
//    }
//
//    if (ypn > 0.99e30)
//        qn = un = 0.0;
//    else {
//        qn = 0.5;
//        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
//    }
//
//    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
//
//    for (k=n-2;k>=0;k--) y2[k] = y2[k]*y2[k+1]+u[k];
//}
//
//static double splint(double xa[], double ya[], double y2a[], int n, double x) {
//    int klo, khi, k;
//    double h, b, a;
//
//    klo = 0;
//    khi = n-1;
//    while (khi-klo > 1) {
//        k = (khi+klo) >> 1;
//        if (xa[k] > x) khi = k;
//        else klo = k;
//    }
//
//    h = xa[khi] - xa[klo];
//    if (h == 0.0) NSLog(@"Internal error in spline function.");
//
//    a = (xa[khi]-x)/h;
//    b = (x-xa[klo])/h;
//
//    return a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi])*(h*h)/6.0;
//}

static const unsigned int returnValueOfG       =  1;
static const unsigned int returnValueOfH       =  2;
static const unsigned int returnValueOfS       =  3;
static const unsigned int returnValueOfCP      =  4;
static const unsigned int returnValueOfDCPDT   =  5;
static const unsigned int returnValueOfV       =  6;
static const unsigned int returnValueOfdVdT    =  7;
static const unsigned int returnValueOfdVdP    =  8;
static const unsigned int returnValueOfd2VdT2  =  9;
static const unsigned int returnValueOfd2VdTdP = 10;
static const unsigned int returnValueOfd2VdP2  = 11;

#pragma mark -
#pragma mark Zhang & Duan H2O class

@interface ZhangDuan2009() {
//    double xSpline[nZandDcorrections];
//    double uSpline[nZandDcorrections];
//    double ySplineForG[nZandDcorrections];
//    double ySplineForH[nZandDcorrections];
//    double ySplineForS[nZandDcorrections];
//    double ySplineForCp[nZandDcorrections];
//    double y2SplineForG[nZandDcorrections];
//    double y2SplineForH[nZandDcorrections];
//    double y2SplineForS[nZandDcorrections];
//    double y2SplineForCp[nZandDcorrections];
}

@end

@implementation ZhangDuan2009

-(id)init {
    if ((self = [super init])) {
        [self setPhaseName:@"Water"];
        [self setPhaseFormula:@"H2O"];

//        for (NSInteger i=0; i<nZandDcorrections; i++) {
//            xSpline[i]      = zAndDcorrections[i].t;
//            ySplineForG[i]  = zAndDcorrections[i].g;
//            ySplineForH[i]  = zAndDcorrections[i].h;
//            ySplineForS[i]  = zAndDcorrections[i].s;
//            ySplineForCp[i] = zAndDcorrections[i].cp;
//        }
//        spline(xSpline, ySplineForG,  nZandDcorrections, 1.0e30, 1.0e30, y2SplineForG,  uSpline);
//        spline(xSpline, ySplineForH,  nZandDcorrections, 1.0e30, 1.0e30, y2SplineForH,  uSpline);
//        spline(xSpline, ySplineForS,  nZandDcorrections, 1.0e30, 1.0e30, y2SplineForS,  uSpline);
//        spline(xSpline, ySplineForCp, nZandDcorrections, 1.0e30, 1.0e30, y2SplineForCp, uSpline);
    }
    return self;
}

-(double)calculateWithT:(double)t andWithP:(double)p forCase:(NSUInteger)returnMode {
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
    double result = 0.0;

    propertiesOfPureH2O(t, p, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

    switch (returnMode) {
        case 1:
            result = gH2O; // - splint(xSpline, ySplineForG, y2SplineForG, nZandDcorrections, t);
            break;
        case 2:
            result = hH2O; // - splint(xSpline, ySplineForH, y2SplineForH, nZandDcorrections, t);
            break;
        case 3:
            result = sH2O; //  - splint(xSpline, ySplineForS, y2SplineForS, nZandDcorrections, t);
            break;
        case 4:
            result = cpH2O; //  - splint(xSpline, ySplineForCp, y2SplineForCp, nZandDcorrections, t);
            break;
        case 5:
            result = dcpdtH2O;
            break;
        case 6:
            result = vH2O;
            break;
        case 7:
            result = dvdtH2O;
            break;
        case 8:
            result = dvdpH2O;
            break;
        case 9:
            result = d2vdt2H2O;
            break;
        case 10:
            result = d2vdtdpH2O;
            break;
        case 11:
            result = d2vdp2H2O;
            break;
        default:
            break;
    }
    return result;
}


-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfG];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfH];
}

-(double)getEntropyFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfS];
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfCP];
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfDCPDT];
}

-(double)getVolumeFromT:(double)t andP:(double)p {
    return  [self calculateWithT:t andWithP:p forCase:returnValueOfV];
}

-(double)getDvDtFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfdVdT];
}

-(double)getDvDpFromT:(double)t andP:(double)p {
    return  [self calculateWithT:t andWithP:p forCase:returnValueOfdVdP];
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdT2];
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdTdP];
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
    return [self calculateWithT:t andWithP:p forCase:returnValueOfd2VdP2];
}

@end
