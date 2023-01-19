//
//  FluidDuan.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/25/12.
//  Copyright 2012 OFM Research Inc.. All rights reserved.
//

/*
 * H2O-CO2 model of Duan, Z, Zhang, Z (2006)
 * Equation of state of the H2O, CO2, and H2O-CO2 systems up to 10 GPa and 2573.15 K:
 * Molecular dynamics simulations with ab initio potential surface
 * Geochimica et Cosmochimica Acta, 70, 2311-2324
 */

#import "FluidDuan.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

#define DO_NOT_THROW_EXCEPTIONS

#pragma mark -
#pragma mark C preprocessor defines

#define R 8.3143

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

#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

#pragma mark -
#pragma mark static constants

/* pure EOSterms for 0 to 0.2 GPa */
static const double H2OLa1  =  4.38269941E-02;
static const double H2OLa2  = -1.68244362E-01;
static const double H2OLa3  = -2.36923373E-01;
static const double H2OLa4  =  1.13027462E-02;
static const double H2OLa5  = -7.67764181E-02;
static const double H2OLa6  =  9.71820593E-02;
static const double H2OLa7  =  6.62674916E-05;
static const double H2OLa8  =  1.06637349E-03;
static const double H2OLa9  = -1.23265258E-03;
static const double H2OLa10 = -8.93953948E-06;
static const double H2OLa11 = -3.88124606E-05;
static const double H2OLa12 =  5.61510206E-05;
static const double H2OLa   =  7.51274488E-03;
static const double H2OLb   =  2.51598931E+00;
static const double H2OLc   =  3.94000000E-02;

static const double CO2La1  =  1.14400435E-01;
static const double CO2La2  = -9.38526684E-01;
static const double CO2La3  =  7.21857006E-01;
static const double CO2La4  =  8.81072902E-03;
static const double CO2La5  =  6.36473911E-02;
static const double CO2La6  = -7.70822213E-02;
static const double CO2La7  =  9.01506064E-04;
static const double CO2La8  = -6.81834166E-03;
static const double CO2La9  =  7.32364258E-03;
static const double CO2La10 = -1.10288237E-04;
static const double CO2La11 =  1.26524193E-03;
static const double CO2La12 = -1.49730823E-03;
static const double CO2La   =  7.81940730E-03;
static const double CO2Lb   = -4.22918013E+00;
static const double CO2Lc   =  1.58500000E-01;

/* pure EOS terms for 0.2 to 10 GPa */
static const double H2OHa1  =  4.68071541E-02;
static const double H2OHa2  = -2.81275941E-01;
static const double H2OHa3  = -2.43926365E-01;
static const double H2OHa4  =  1.10016958E-02;
static const double H2OHa5  = -3.86603525E-02;
static const double H2OHa6  =  9.30095461E-02;
static const double H2OHa7  = -1.15747171E-05;
static const double H2OHa8  =  4.19873848E-04;
static const double H2OHa9  = -5.82739501E-04;
static const double H2OHa10 =  1.00936000E-06;
static const double H2OHa11 = -1.01713593E-05;
static const double H2OHa12 =  1.63934213E-05;
static const double H2OHa   = -4.49505919E-02;
static const double H2OHb   = -3.15028174E-01;
static const double H2OHc   =  1.25000000E-02;

static const double CO2Ha1  =  5.72573440E-03;
static const double CO2Ha2  =  7.94836769E+00;
static const double CO2Ha3  = -3.84236281E+01;
static const double CO2Ha4  =  3.71600369E-02;
static const double CO2Ha5  = -1.92888994E+00;
static const double CO2Ha6  =  6.64254770E+00;
static const double CO2Ha7  = -7.02203950E-06;
static const double CO2Ha8  =  1.77093234E-02;
static const double CO2Ha9  = -4.81892026E-02;
static const double CO2Ha10 =  3.88344869E-06;
static const double CO2Ha11 = -5.54833167E-04;
static const double CO2Ha12 =  1.70489748E-03;
static const double CO2Ha   = -4.13039220E-01;
static const double CO2Hb   = -8.47988634E+00;
static const double CO2Hc   =  2.80000000E-02;

/* H2O, critical constants, K, bars, J/bar */

static const double H2OTc = 647.25;
static const double H2OPc = 221.19;
#define H2OVc (8.314467*H2OTc/H2OPc)

/* CO2 */

static const double CO2Tc = 304.1282;
static const double CO2Pc = 73.773;
#define CO2Vc (8.314467*CO2Tc/CO2Pc)

#define H2O 0
#define CO2 1

static const double idealCoeff[13][2] = {
    {  3.10409601236035e+01, -0.18188731e+01 },
    { -3.91422080460869e+01,  0.12903022e+02 },
    {  3.79695277233575e+01, -0.96634864e+01 },
    { -2.18374910952284e+01,  0.42251879e+01 },
    {  7.42251494566339e+00, -0.10421640e+01 },
    { -1.38178929609470e+00,  0.12683515e+00 },
    {  1.08807067571454e-01, -0.49939675e-02 },
    { -1.20771176848589e+01,  0.24950242e+01 },
    {  3.39105078851732e+00, -0.82723750e+00 },
    { -5.84520979955060e-01,  0.15372481e+00 },
    {  5.89930846488082e-02, -0.15861243e-01 },
    { -3.12970001415882e-03,  0.86017150e-03 },
    {  6.57460740981757e-05, -0.19222165e-04 }
};

#pragma mark -
#pragma mark Math functions

static double powSum(double a, double fa, double b, double fb) {
    double sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    return (sum >= 0.0) ? pow(sum/(fa+fb), 3.0) : -pow(-sum/(fa+fb), 3.0);
}

static double DpowSum(double a, double fa, double b, double fb, double da, double db) {
    double sum = 0.0, dsum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    sum = (sum >= 0.0) ? pow(sum/(fa+fb), 2.0) : pow(-sum/(fa+fb), 2.0);
    dsum += (a >= 0.0) ? da*fa/pow(a, 2.0/3.0)/3.0 : da*fa/pow(-a, 2.0/3.0)/3.0;
    dsum += (b >= 0.0) ? db*fb/pow(b, 2.0/3.0)/3.0 : db*fb/pow(-b, 2.0/3.0)/3.0;
    return 3.0*sum*dsum/(fa+fb);
}

static double D2powSum(double a, double fa, double b, double fb, double da, double db, double d2a, double d2b) {
    double sum = 0.0, dsum = 0.0, d2sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    sum /= fa + fb;

    dsum += (a >= 0.0) ? da*fa/pow(a, 2.0/3.0)/3.0 : da*fa/pow(-a, 2.0/3.0)/3.0;
    dsum += (b >= 0.0) ? db*fb/pow(b, 2.0/3.0)/3.0 : db*fb/pow(-b, 2.0/3.0)/3.0;
    dsum /= fa + fb;

    d2sum += (a >= 0.0) ? d2a*fa/pow(a, 2.0/3.0)/3.0 : d2a*fa/pow(-a, 2.0/3.0)/3.0;
    d2sum += (b >= 0.0) ? d2b*fb/pow(b, 2.0/3.0)/3.0 : d2b*fb/pow(-b, 2.0/3.0)/3.0;
    d2sum += (a >= 0.0) ? -2.0*da*da*fa/pow(a, 5.0/3.0)/9.0 : 2.0*da*da*fa/pow(-a, 5.0/3.0)/9.0;
    d2sum += (b >= 0.0) ? -2.0*db*db*fb/pow(b, 5.0/3.0)/9.0 : 2.0*db*db*fb/pow(-b, 5.0/3.0)/9.0;
    d2sum /= fa + fb;

    return 6.0*sum*dsum*dsum + 3.0*sum*sum*d2sum;
}

static double D3powSum(double a, double fa, double b, double fb, double da, double db, double d2a, double d2b, double d3a, double d3b) {
    double sum = 0.0, dsum = 0.0, d2sum = 0.0, d3sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    sum /= fa + fb;

    dsum += (a >= 0.0) ? da*fa/pow(a, 2.0/3.0)/3.0 : da*fa/pow(-a, 2.0/3.0)/3.0;
    dsum += (b >= 0.0) ? db*fb/pow(b, 2.0/3.0)/3.0 : db*fb/pow(-b, 2.0/3.0)/3.0;
    dsum /= fa + fb;

    d2sum += (a >= 0.0) ? d2a*fa/pow(a, 2.0/3.0)/3.0 : d2a*fa/pow(-a, 2.0/3.0)/3.0;
    d2sum += (b >= 0.0) ? d2b*fb/pow(b, 2.0/3.0)/3.0 : d2b*fb/pow(-b, 2.0/3.0)/3.0;
    d2sum += (a >= 0.0) ? -2.0*da*da*fa/pow(a, 5.0/3.0)/9.0 : 2.0*da*da*fa/pow(-a, 5.0/3.0)/9.0;
    d2sum += (b >= 0.0) ? -2.0*db*db*fb/pow(b, 5.0/3.0)/9.0 : 2.0*db*db*fb/pow(-b, 5.0/3.0)/9.0;
    d2sum /= fa + fb;

    d3sum += (a >= 0.0) ? d3a*fa/pow(a, 2.0/3.0)/3.0 : d3a*fa/pow(-a, 2.0/3.0)/3.0;
    d3sum += (b >= 0.0) ? d3b*fb/pow(b, 2.0/3.0)/3.0 : d3b*fb/pow(-b, 2.0/3.0)/3.0;
    d3sum += (a >= 0.0) ? -6.0*da*d2a*fa/pow(a, 5.0/3.0)/9.0 : 6.0*da*d2a*fa/pow(-a, 5.0/3.0)/9.0;
    d3sum += (b >= 0.0) ? -6.0*db*d2b*fb/pow(b, 5.0/3.0)/9.0 : 6.0*db*d2b*fb/pow(-b, 5.0/3.0)/9.0;
    d3sum += (a >= 0.0) ? 10.0*da*da*da*fa/pow(a, 8.0/3.0)/27.0 : 10.0*da*da*da*fa/pow(-a, 8.0/3.0)/27.0;
    d3sum += (b >= 0.0) ? 10.0*db*db*db*fb/pow(b, 8.0/3.0)/27.0 : 10.0*db*db*db*fb/pow(-b, 8.0/3.0)/27.0;
    d3sum /= fa + fb;

    return 6.0*dsum*dsum*dsum + 18.0*sum*dsum*d2sum + 3.0*sum*sum*d3sum;
}

#pragma mark -
#pragma mark Virial functions

static void BVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *bv, double bvPrime[2], double *dbvdt, double *d2bvdt2, double *d3bvdt3, double dbvPrimedt[2], double d2bvPrimedt2[2], double d3bvPrimedt3[2],
                             double b2vPrime2[2][2]) {
    double bEnd[2], dbEnddt[2], d2bEnddt2[2], d3bEnddt3[2];
    double H2OTr    = t/H2OTc;
    double CO2Tr    = t/CO2Tc;
    double k1       = 0.0;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;
    double dk1dt    = 0.0;
    double d2k1dt2  = 0.0;
    double d3k1dt3  = 0.0;

    if (useLowPcoeff) {
        bEnd[H2O] = H2OLa1 + H2OLa2/H2OTr/H2OTr + H2OLa3/H2OTr/H2OTr/H2OTr;
        bEnd[CO2] = CO2La1 + CO2La2/CO2Tr/CO2Tr + CO2La3/CO2Tr/CO2Tr/CO2Tr;
        k1 = 3.131 - 5.0624e-3*t + 1.8641e-6*t*t - 31.409/t;
        dbEnddt[H2O] = - 2.0*H2OLa2*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa3*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dbEnddt[CO2] = - 2.0*CO2La2*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La3*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk1dt = - 5.0624e-3 + 2.0*1.8641e-6*t + 31.409/t/t;
        d2bEnddt2[H2O] = 6.0*H2OLa2*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa3*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2bEnddt2[CO2] = 6.0*CO2La2*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La3*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k1dt2 = 2.0*1.8641e-6 - 2.0*31.409/t/t/t;
        d3bEnddt3[H2O] = - 24.0*H2OLa2*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa3*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3bEnddt3[CO2] = - 24.0*CO2La2*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La3*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k1dt3 =  6.0*31.409/t/t/t/t;
    } else {
        bEnd[H2O] = H2OHa1 + H2OHa2/H2OTr/H2OTr + H2OHa3/H2OTr/H2OTr/H2OTr;
        bEnd[CO2] = CO2Ha1 + CO2Ha2/CO2Tr/CO2Tr + CO2Ha3/CO2Tr/CO2Tr/CO2Tr;
        k1 = 9.034 - 7.9212e-3*t + 2.3285e-6*t*t - 2.4221e3/t;
        dbEnddt[H2O] = - 2.0*H2OHa2*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa3*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dbEnddt[CO2] = - 2.0*CO2Ha2*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha3*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk1dt = - 7.9212e-3 + 2.0*2.3285e-6*t + 2.4221e3/t/t;
        d2bEnddt2[H2O] = 6.0*H2OHa2*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa3*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2bEnddt2[CO2] = 6.0*CO2Ha2*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha3*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k1dt2 = 2.0*2.3285e-6 - 2.0*2.4221e3/t/t/t;
        d3bEnddt3[H2O] = - 24.0*H2OHa2*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa3*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3bEnddt3[CO2] = - 24.0*CO2Ha2*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha3*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k1dt3 = 6.0*2.4221e3/t/t/t/t;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *bv = 0.0;
        bvPrime[H2O] = 0.0;
        bvPrime[CO2] = 0.0;

        b2vPrime2[H2O][H2O] = 0.0;
        b2vPrime2[H2O][CO2] = 0.0;
        b2vPrime2[CO2][CO2] = 0.0;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt   = 0.0;
        *d2bvdt2 = 0.0;
        *d3bvdt3 = 0.0;

        dbvPrimedt[H2O] = 0.0;
        dbvPrimedt[CO2] = 0.0;

        d2bvPrimedt2[H2O] = 0.0;
        d2bvPrimedt2[CO2] = 0.0;

        d3bvPrimedt3[H2O] = 0.0;
        d3bvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *bv = bEnd[H2O]*H2OVc;
        bvPrime[H2O] = 2.0*bEnd[H2O]*H2OVc;
        bvPrime[CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);

        b2vPrime2[H2O][H2O] = 2.0*bEnd[H2O]*H2OVc;
        b2vPrime2[H2O][CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        b2vPrime2[CO2][CO2] = 2.0*bEnd[CO2]*CO2Vc;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt   = dbEnddt[H2O]*H2OVc;
        *d2bvdt2 = d2bEnddt2[H2O]*H2OVc;
        *d3bvdt3 = d3bEnddt3[H2O]*H2OVc;

        dbvPrimedt[H2O] = 2.0*dbEnddt[H2O]*H2OVc;
        dbvPrimedt[CO2] = 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                               + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                               )*powSum(H2OVc, 1.0, CO2Vc, 1.0);

        d2bvPrimedt2[H2O] = 2.0*d2bEnddt2[H2O]*H2OVc;
        d2bvPrimedt2[CO2] = 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                 + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);

        d3bvPrimedt3[H2O] = 2.0*d3bEnddt3[H2O]*H2OVc;
        d3bvPrimedt3[CO2] = 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                 + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                 + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *bv = bEnd[CO2]*CO2Vc;
        bvPrime[H2O] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        bvPrime[CO2] = 2.0*bEnd[CO2]*CO2Vc;

        b2vPrime2[H2O][H2O] = 2.0*bEnd[H2O]*H2OVc;
        b2vPrime2[H2O][CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        b2vPrime2[CO2][CO2] = 2.0*bEnd[CO2]*CO2Vc;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt   = dbEnddt[CO2]*CO2Vc;
        *d2bvdt2 = d2bEnddt2[CO2]*CO2Vc;
        *d3bvdt3 = d3bEnddt3[CO2]*CO2Vc;

        dbvPrimedt[H2O] = 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                               + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                               )*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        dbvPrimedt[CO2] = 2.0*dbEnddt[CO2]*CO2Vc;

        d2bvPrimedt2[H2O] = 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                 + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        d2bvPrimedt2[CO2] = 2.0*d2bEnddt2[CO2]*CO2Vc;

        d3bvPrimedt3[H2O] = 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                 + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                 + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        d3bvPrimedt3[CO2] = 2.0*d3bEnddt3[CO2]*CO2Vc;

    } else {
        *bv  = 0.0;
        *bv += bEnd[H2O]*H2OVc*x[H2O]*x[H2O];
        *bv += 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *bv += bEnd[CO2]*CO2Vc*x[CO2]*x[CO2];

        bvPrime[H2O]  = 0.0;
        bvPrime[H2O] += 2.0*bEnd[H2O]*H2OVc*x[H2O];
        bvPrime[H2O] += 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        bvPrime[CO2]  = 0.0;
        bvPrime[CO2] += 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        bvPrime[CO2] += 2.0*bEnd[CO2]*CO2Vc*x[CO2];

        b2vPrime2[H2O][H2O] = 2.0*bEnd[H2O]*H2OVc;
        b2vPrime2[H2O][CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        b2vPrime2[CO2][CO2] = 2.0*bEnd[CO2]*CO2Vc;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt  = 0.0;
        *dbvdt += dbEnddt[H2O]*H2OVc*x[H2O]*x[H2O];
        *dbvdt += 2.0*(DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1 + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                       )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *dbvdt += dbEnddt[CO2]*CO2Vc*x[CO2]*x[CO2];

        *d2bvdt2  = 0.0;
        *d2bvdt2 += d2bEnddt2[H2O]*H2OVc*x[H2O]*x[H2O];
        *d2bvdt2 += 2.0*(  D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                         + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                         + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                         )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *d2bvdt2 += d2bEnddt2[CO2]*CO2Vc*x[CO2]*x[CO2];

        *d3bvdt3  = 0.0;
        *d3bvdt3 += d3bEnddt3[H2O]*H2OVc*x[H2O]*x[H2O];
        *d3bvdt3 += 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                         + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                         + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                         +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                         )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *d3bvdt3 += d3bEnddt3[CO2]*CO2Vc*x[CO2]*x[CO2];

        dbvPrimedt[H2O]  = 0.0;
        dbvPrimedt[H2O] += 2.0*dbEnddt[H2O]*H2OVc*x[H2O];
        dbvPrimedt[H2O] += 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                                + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                                )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        dbvPrimedt[CO2]  = 0.0;
        dbvPrimedt[CO2] += 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                                + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                                )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        dbvPrimedt[CO2] += 2.0*dbEnddt[CO2]*CO2Vc*x[CO2];

        d2bvPrimedt2[H2O]  = 0.0;
        d2bvPrimedt2[H2O] += 2.0*d2bEnddt2[H2O]*H2OVc*x[H2O];
        d2bvPrimedt2[H2O] += 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                  + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        d2bvPrimedt2[CO2]  = 0.0;
        d2bvPrimedt2[CO2] += 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                  + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        d2bvPrimedt2[CO2] += 2.0*d2bEnddt2[CO2]*CO2Vc*x[CO2];

        d3bvPrimedt3[H2O]  = 0.0;
        d3bvPrimedt3[H2O] += 2.0*d3bEnddt3[H2O]*H2OVc*x[H2O];
        d3bvPrimedt3[H2O] += 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                  + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                  + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        d3bvPrimedt3[CO2]  = 0.0;
        d3bvPrimedt3[CO2] += 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                  + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                  + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        d3bvPrimedt3[CO2] += 2.0*d3bEnddt3[CO2]*CO2Vc*x[CO2];

    }

    return;
}

static void CVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *cv, double cvPrime[2], double *dcvdt, double *d2cvdt2, double *d3cvdt3, double dcvPrimedt[2], double d2cvPrimedt2[2], double d3cvPrimedt3[2],
                             double c2vPrime2[2][2]) {
    double cEnd[2], dcEnddt[2], d2cEnddt2[2], d3cEnddt3[2];
    double H2OTr    = t/H2OTc;
    double CO2Tr    = t/CO2Tc;
    double k2       = 0.0;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;
    double dk2dt    = 0.0;
    double d2k2dt2  = 0.0;
    double d3k2dt3  = 0.0;

    if (useLowPcoeff) {
        cEnd[H2O] = H2OLa4 + H2OLa5/H2OTr/H2OTr + H2OLa6/H2OTr/H2OTr/H2OTr;
        cEnd[CO2] = CO2La4 + CO2La5/CO2Tr/CO2Tr + CO2La6/CO2Tr/CO2Tr/CO2Tr;
        k2 = -46.646 + 4.2877e-2*t - 1.0892e-5*t*t + 1.5782e4/t;
        dcEnddt[H2O] = - 2.0*H2OLa5*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa6*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dcEnddt[CO2] = - 2.0*CO2La5*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La6*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk2dt = 4.2877e-2 - 2.0*1.0892e-5*t - 1.5782e4/t/t;
        d2cEnddt2[H2O] = 6.0*H2OLa5*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa6*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2cEnddt2[CO2] = 6.0*CO2La5*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La6*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k2dt2 = - 2.0*1.0892e-5 + 2.0*1.5782e4/t/t/t;
        d3cEnddt3[H2O] = -24.0*H2OLa5*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa6*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3cEnddt3[CO2] = -24.0*CO2La5*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La6*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k2dt3 = - 3.0*2.0*1.5782e4/t/t/t/t;
    } else {
        cEnd[H2O] = H2OHa4 + H2OHa5/H2OTr/H2OTr + H2OHa6/H2OTr/H2OTr/H2OTr;
        cEnd[CO2] = CO2Ha4 + CO2Ha5/CO2Tr/CO2Tr + CO2Ha6/CO2Tr/CO2Tr/CO2Tr;
        k2 = -1.068 + 1.8756e-3*t - 4.9371e-7*t*t + 6.6180e2/t;
        dcEnddt[H2O] = - 2.0*H2OHa5*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa6*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dcEnddt[CO2] = - 2.0*CO2Ha5*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha6*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk2dt = 1.8756e-3 - 2.0*4.9371e-7*t - 6.6180e2/t/t;
        d2cEnddt2[H2O] = 6.0*H2OHa5*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa6*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2cEnddt2[CO2] = 6.0*CO2Ha5*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha6*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k2dt2 = - 2.0*4.9371e-7 + 2.0*6.6180e2/t/t/t;
        d3cEnddt3[H2O] = -24.0*H2OHa5*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa6*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3cEnddt3[CO2] = -24.0*CO2Ha5*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha6*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k2dt3 = - 3.0*2.0*6.6180e2/t/t/t/t;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *cv = 0.0;
        cvPrime[H2O] = 0.0;
        cvPrime[CO2] = 0.0;

        c2vPrime2[H2O][H2O] = 0.0;
        c2vPrime2[H2O][CO2] = 0.0;
        c2vPrime2[CO2][CO2] = 0.0;
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt   = 0.0;
        *d2cvdt2 = 0.0;
        *d3cvdt3 = 0.0;

        dcvPrimedt[H2O] = 0.0;
        dcvPrimedt[CO2] = 0.0;

        d2cvPrimedt2[H2O] = 0.0;
        d2cvPrimedt2[CO2] = 0.0;

        d3cvPrimedt3[H2O] = 0.0;
        d3cvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *cv = cEnd[H2O]*H2OVc*H2OVc;
        cvPrime[H2O] = 3.0*cEnd[H2O]*H2OVc*H2OVc;
        cvPrime[CO2] = 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        c2vPrime2[H2O][H2O] = 6.0*cEnd[H2O]*H2OVc*H2OVc;
        c2vPrime2[H2O][CO2] = 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        c2vPrime2[CO2][CO2] = 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt   = dcEnddt[H2O]*H2OVc*H2OVc;
        *d2cvdt2 = d2cEnddt2[H2O]*H2OVc*H2OVc;
        *d3cvdt3 = d3cEnddt3[H2O]*H2OVc*H2OVc;

        dcvPrimedt[H2O] = 3.0*dcEnddt[H2O]*H2OVc*H2OVc;
        dcvPrimedt[CO2] = 3.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        d2cvPrimedt2[H2O] = 3.0*d2cEnddt2[H2O]*H2OVc*H2OVc;
        d2cvPrimedt2[CO2] = 3.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 6.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        d3cvPrimedt3[H2O] = 3.0*d3cEnddt3[H2O]*H2OVc*H2OVc;
        d3cvPrimedt3[CO2] = 3.0*D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 9.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 9.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *cv = cEnd[CO2]*CO2Vc*CO2Vc;
        cvPrime[H2O] = 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        cvPrime[CO2] = 3.0*cEnd[CO2]*CO2Vc*CO2Vc;

        c2vPrime2[H2O][H2O] = 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        c2vPrime2[H2O][CO2] = 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        c2vPrime2[CO2][CO2] = 6.0*cEnd[CO2]*CO2Vc*CO2Vc;
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt   = dcEnddt[CO2]*CO2Vc*CO2Vc;
        *d2cvdt2 = d2cEnddt2[CO2]*CO2Vc*CO2Vc;
        *d3cvdt3 = d3cEnddt3[CO2]*CO2Vc*CO2Vc;

        dcvPrimedt[H2O] = 3.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        dcvPrimedt[CO2] = 3.0*dcEnddt[CO2]*CO2Vc*CO2Vc;

        d2cvPrimedt2[H2O] = 3.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 6.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        d2cvPrimedt2[CO2] = 3.0*d2cEnddt2[CO2]*CO2Vc*CO2Vc;

        d3cvPrimedt3[H2O] = 3.0*D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 9.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 9.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        d3cvPrimedt3[CO2] = 3.0*d3cEnddt3[CO2]*CO2Vc*CO2Vc;

    } else {
        *cv  = 0.0;
        *cv += cEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *cv += 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *cv += 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *cv += cEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        cvPrime[H2O]  = 0.0;
        cvPrime[H2O] += 3.0*cEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        cvPrime[H2O] += 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        cvPrime[H2O] += 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        cvPrime[CO2]  = 0.0;
        cvPrime[CO2] += 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        cvPrime[CO2] += 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        cvPrime[CO2] += 3.0*cEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        c2vPrime2[H2O][H2O] = 6.0*cEnd[H2O]*H2OVc*H2OVc*x[H2O]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[CO2];
        c2vPrime2[H2O][CO2] = 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2];
        c2vPrime2[CO2][CO2] = 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]
        + 6.0*cEnd[CO2]*CO2Vc*CO2Vc*x[CO2];
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt  = 0.0;
        *dcvdt += dcEnddt[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *dcvdt += 3.0*(DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2 + powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt
                       )*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *dcvdt += 3.0*(DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2 + powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt
                       )*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *dcvdt += dcEnddt[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        *d2cvdt2  = 0.0;
        *d2cvdt2 += d2cEnddt2[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *d2cvdt2 += 3.0*(  D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2
                         + 2.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt
                         + powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2
                         )*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *d2cvdt2 += 3.0*(  D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2
                         + 2.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt
                         + powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2
                         )*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *d2cvdt2 += d2cEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        *d3cvdt3  = 0.0;
        *d3cvdt3 += d3cEnddt3[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *d3cvdt3 += 3.0*(  D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2
                         + 3.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt
                         + 3.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2
                         + powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3
                         )*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *d3cvdt3 += 3.0*(  D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2
                         + 3.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt
                         + 3.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2
                         + powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3
                         )*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *d3cvdt3 += d3cEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        dcvPrimedt[H2O]  = 0.0;
        dcvPrimedt[H2O] += 3.0*dcEnddt[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        dcvPrimedt[H2O] += 6.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        dcvPrimedt[H2O] += 3.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        dcvPrimedt[CO2]  = 0.0;
        dcvPrimedt[CO2] += 3.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        dcvPrimedt[CO2] += 6.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        dcvPrimedt[CO2] += 3.0*dcEnddt[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        d2cvPrimedt2[H2O]  = 0.0;
        d2cvPrimedt2[H2O] += 3.0*d2cEnddt2[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        d2cvPrimedt2[H2O] += 6.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 12.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        d2cvPrimedt2[H2O] += 3.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 6.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        d2cvPrimedt2[CO2]  = 0.0;
        d2cvPrimedt2[CO2] += 3.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 6.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        d2cvPrimedt2[CO2] += 6.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 12.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        d2cvPrimedt2[CO2] += 3.0*d2cEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        d3cvPrimedt3[H2O]  = 0.0;
        d3cvPrimedt3[H2O] += 3.0*d3cEnddt3[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        d3cvPrimedt3[H2O] += 6.0*D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        d3cvPrimedt3[H2O] += 3.0*D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 9.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 9.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        d3cvPrimedt3[CO2]  = 0.0;
        d3cvPrimedt3[CO2] += 3.0*D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 9.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 9.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        d3cvPrimedt3[CO2] += 6.0*D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        d3cvPrimedt3[CO2] += 3.0*d3cEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

    }

    return;
}

static void DVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *dv, double dvPrime[2], double *ddvdt, double *d2dvdt2, double *d3dvdt3, double ddvPrimedt[2], double d2dvPrimedt2[2], double d3dvPrimedt3[2],
                             double d2vPrime2[2][2]) {
    double dEnd[2], ddEnddt[2], d2dEnddt2[2], d3dEnddt3[2];
    double H2OTr = t/H2OTc;
    double CO2Tr = t/CO2Tc;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;

    if (useLowPcoeff) {
        dEnd[H2O] = H2OLa7 + H2OLa8/H2OTr/H2OTr + H2OLa9/H2OTr/H2OTr/H2OTr;
        dEnd[CO2] = CO2La7 + CO2La8/CO2Tr/CO2Tr + CO2La9/CO2Tr/CO2Tr/CO2Tr;
        ddEnddt[H2O] = - 2.0*H2OLa8*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa9*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        ddEnddt[CO2] = - 2.0*CO2La8*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La9*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2dEnddt2[H2O] = 6.0*H2OLa8*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa9*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2dEnddt2[CO2] = 6.0*CO2La8*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La9*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3dEnddt3[H2O] = - 24.0*H2OLa8*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa9*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3dEnddt3[CO2] = - 24.0*CO2La8*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La9*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    } else {
        dEnd[H2O] = H2OHa7 + H2OHa8/H2OTr/H2OTr + H2OHa9/H2OTr/H2OTr/H2OTr;
        dEnd[CO2] = CO2Ha7 + CO2Ha8/CO2Tr/CO2Tr + CO2Ha9/CO2Tr/CO2Tr/CO2Tr;
        ddEnddt[H2O] = - 2.0*H2OHa8*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa9*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        ddEnddt[CO2] = - 2.0*CO2Ha8*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha9*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2dEnddt2[H2O] = 6.0*H2OHa8*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa9*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2dEnddt2[CO2] = 6.0*CO2Ha8*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha9*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3dEnddt3[H2O] = - 24.0*H2OHa8*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa9*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3dEnddt3[CO2] = - 24.0*CO2Ha8*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha9*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *dv = 0.0;
        dvPrime[H2O] = 0.0;
        dvPrime[CO2] = 0.0;

        d2vPrime2[H2O][H2O] = 0.0;
        d2vPrime2[H2O][CO2] = 0.0;
        d2vPrime2[CO2][CO2] = 0.0;
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt   = 0.0;
        *d2dvdt2 = 0.0;
        *d3dvdt3 = 0.0;

        ddvPrimedt[H2O] = 0.0;
        ddvPrimedt[CO2] = 0.0;

        d2dvPrimedt2[H2O] = 0.0;
        d2dvPrimedt2[CO2] = 0.0;

        d3dvPrimedt3[H2O] = 0.0;
        d3dvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *dv = dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        dvPrime[H2O] = 5.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        dvPrime[CO2] = 5.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

        d2vPrime2[H2O][H2O] = 20.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        d2vPrime2[H2O][CO2] = 20.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        d2vPrime2[CO2][CO2] = 20.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0);
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt   = ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        *d2dvdt2 = d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        *d3dvdt3 = d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;

        ddvPrimedt[H2O] = 5.0*ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        ddvPrimedt[CO2] = 5.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

        d2dvPrimedt2[H2O] = 5.0*d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        d2dvPrimedt2[CO2] = 5.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

        d3dvPrimedt3[H2O] = 5.0*d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        d3dvPrimedt3[CO2] = 5.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *dv = dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        dvPrime[H2O] = 5.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        dvPrime[CO2] = 5.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d2vPrime2[H2O][H2O] = 20.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0);
        d2vPrime2[H2O][CO2] = 20.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);
        d2vPrime2[CO2][CO2] = 20.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt   = ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d2dvdt2 = d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d3dvdt3 = d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        ddvPrimedt[H2O] = 5.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        ddvPrimedt[CO2] = 5.0*ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d2dvPrimedt2[H2O] = 5.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        d2dvPrimedt2[CO2] = 5.0*d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d3dvPrimedt3[H2O] = 5.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        d3dvPrimedt3[CO2] = 5.0*d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

    } else {
        *dv  = 0.0;
        *dv += dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *dv +=  5.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *dv += 10.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *dv += 10.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *dv +=  5.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *dv += dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        dvPrime[H2O]  =  0.0;
        dvPrime[H2O] +=  5.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        dvPrime[H2O] += 20.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        dvPrime[H2O] += 30.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        dvPrime[H2O] += 20.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        dvPrime[H2O] +=  5.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        dvPrime[CO2]  =  0.0;
        dvPrime[CO2] +=  5.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        dvPrime[CO2] += 20.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        dvPrime[CO2] += 30.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        dvPrime[CO2] += 20.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        dvPrime[CO2] += 5.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2vPrime2[H2O][H2O] = 20.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]
        + 60.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[CO2]
        + 60.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[CO2]*x[CO2]
        + 20.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[CO2]*x[CO2]*x[CO2];
        d2vPrime2[H2O][CO2] = 20.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]
        + 60.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]
        + 60.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]
        + 20.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2];
        d2vPrime2[CO2][CO2] = 20.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]
        + 60.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]
        + 60.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]
        + 20.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt  = 0.0;
        *ddvdt += ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *ddvdt +=  5.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *ddvdt += 10.0*DpowSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *ddvdt += 10.0*DpowSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *ddvdt +=  5.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *ddvdt += ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d2dvdt2  = 0.0;
        *d2dvdt2 += d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d2dvdt2 +=  5.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d2dvdt2 += 10.0*D2powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d2dvdt2 += 10.0*D2powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d2dvdt2 +=  5.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d2dvdt2 += d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d3dvdt3  = 0.0;
        *d3dvdt3 += d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d3dvdt3 +=  5.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d3dvdt3 += 10.0*D3powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d3dvdt3 += 10.0*D3powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d3dvdt3 +=  5.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d3dvdt3 += d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        ddvPrimedt[H2O]  =  0.0;
        ddvPrimedt[H2O] +=  5.0*ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        ddvPrimedt[H2O] += 20.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        ddvPrimedt[H2O] += 30.0*DpowSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        ddvPrimedt[H2O] += 20.0*DpowSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        ddvPrimedt[H2O] +=  5.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        ddvPrimedt[CO2]  =  0.0;
        ddvPrimedt[CO2] +=  5.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        ddvPrimedt[CO2] += 20.0*DpowSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        ddvPrimedt[CO2] += 30.0*DpowSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        ddvPrimedt[CO2] += 20.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        ddvPrimedt[CO2] += 5.0*ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2dvPrimedt2[H2O]  =  0.0;
        d2dvPrimedt2[H2O] +=  5.0*d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2dvPrimedt2[H2O] += 20.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2dvPrimedt2[H2O] += 30.0*D2powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2dvPrimedt2[H2O] += 20.0*D2powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2dvPrimedt2[H2O] +=  5.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2dvPrimedt2[CO2]  =  0.0;
        d2dvPrimedt2[CO2] +=  5.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2dvPrimedt2[CO2] += 20.0*D2powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2dvPrimedt2[CO2] += 30.0*D2powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2dvPrimedt2[CO2] += 20.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2dvPrimedt2[CO2] += 5.0*d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d3dvPrimedt3[H2O]  =  0.0;
        d3dvPrimedt3[H2O] +=  5.0*d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3dvPrimedt3[H2O] += 20.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3dvPrimedt3[H2O] += 30.0*D3powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3dvPrimedt3[H2O] += 20.0*D3powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3dvPrimedt3[H2O] +=  5.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d3dvPrimedt3[CO2]  =  0.0;
        d3dvPrimedt3[CO2] +=  5.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3dvPrimedt3[CO2] += 20.0*D3powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3dvPrimedt3[CO2] += 30.0*D3powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3dvPrimedt3[CO2] += 20.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3dvPrimedt3[CO2] += 5.0*d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

    }

    return;
}

static void EVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *ev, double evPrime[2], double *devdt, double *d2evdt2, double *d3evdt3, double devPrimedt[2], double d2evPrimedt2[2], double d3evPrimedt3[2],
                             double e2vPrime2[2][2]) {
    double eEnd[2], deEnddt[2], d2eEnddt2[2], d3eEnddt3[2];
    double H2OTr = t/H2OTc;
    double CO2Tr = t/CO2Tc;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;

    if (useLowPcoeff) {
        eEnd[H2O] = H2OLa10 + H2OLa11/H2OTr/H2OTr + H2OLa12/H2OTr/H2OTr/H2OTr;
        eEnd[CO2] = CO2La10 + CO2La11/CO2Tr/CO2Tr + CO2La12/CO2Tr/CO2Tr/CO2Tr;
        deEnddt[H2O] = - 2.0*H2OLa11*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa12*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        deEnddt[CO2] = - 2.0*CO2La11*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La12*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2eEnddt2[H2O] = 6.0*H2OLa11*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa12*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2eEnddt2[CO2] = 6.0*CO2La11*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La12*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3eEnddt3[H2O] = - 24.0*H2OLa11*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa12*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3eEnddt3[CO2] = - 24.0*CO2La11*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La12*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    } else {
        eEnd[H2O] = H2OHa10 + H2OHa11/H2OTr/H2OTr + H2OHa12/H2OTr/H2OTr/H2OTr;
        eEnd[CO2] = CO2Ha10 + CO2Ha11/CO2Tr/CO2Tr + CO2Ha12/CO2Tr/CO2Tr/CO2Tr;
        deEnddt[H2O] = - 2.0*H2OHa11*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa12*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        deEnddt[CO2] = - 2.0*CO2Ha11*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha12*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2eEnddt2[H2O] = 6.0*H2OHa11*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa12*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2eEnddt2[CO2] = 6.0*CO2Ha11*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha12*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3eEnddt3[H2O] = - 24.0*H2OHa11*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa12*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3eEnddt3[CO2] = - 24.0*CO2Ha11*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha12*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *ev = 0.0;
        evPrime[H2O] = 0.0;
        evPrime[CO2] = 0.0;

        e2vPrime2[H2O][H2O] = 0.0;
        e2vPrime2[H2O][CO2] = 0.0;
        e2vPrime2[CO2][CO2] = 0.0;
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt   = 0.0;
        *d2evdt2 = 0.0;
        *d3evdt3 = 0.0;

        devPrimedt[H2O] = 0.0;
        devPrimedt[CO2] = 0.0;

        d2evPrimedt2[H2O] = 0.0;
        d2evPrimedt2[CO2] = 0.0;

        d3evPrimedt3[H2O] = 0.0;
        d3evPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *ev = eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        evPrime[H2O] = 6.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        evPrime[CO2] = 6.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

        e2vPrime2[H2O][H2O] =  30.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        e2vPrime2[H2O][CO2] =  30.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);
        e2vPrime2[CO2][CO2] =  30.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0);
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt = deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        *d2evdt2 = d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        *d3evdt3 = d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;

        devPrimedt[H2O] = 6.0*deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        devPrimedt[CO2] = 6.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

        d2evPrimedt2[H2O] = 6.0*d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        d2evPrimedt2[CO2] = 6.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

        d3evPrimedt3[H2O] = 6.0*d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        d3evPrimedt3[CO2] = 6.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *ev = eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        evPrime[H2O] = 6.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        evPrime[CO2] = 6.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        e2vPrime2[H2O][H2O] = 30.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0);
        e2vPrime2[H2O][CO2] = 30.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        e2vPrime2[CO2][CO2] = 30.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt = deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d2evdt2 = d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d3evdt3 = d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        devPrimedt[H2O] = 6.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        devPrimedt[CO2] = 6.0*deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d2evPrimedt2[H2O] = 6.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        d2evPrimedt2[CO2] = 6.0*d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d3evPrimedt3[H2O] = 6.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        d3evPrimedt3[CO2] = 6.0*d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

    } else {
        *ev  = 0.0;
        *ev += eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *ev +=  6.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *ev += 15.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *ev += 20.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *ev += 15.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *ev +=  6.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *ev += eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        evPrime[H2O]  =  0.0;
        evPrime[H2O] +=  6.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        evPrime[H2O] += 30.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        evPrime[H2O] += 60.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        evPrime[H2O] += 60.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        evPrime[H2O] += 30.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        evPrime[H2O] +=  6.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        evPrime[CO2]  =  0.0;
        evPrime[CO2] +=  6.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        evPrime[CO2] += 30.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        evPrime[CO2] += 60.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        evPrime[CO2] += 60.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        evPrime[CO2] += 30.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        evPrime[CO2] +=  6.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        e2vPrime2[H2O][H2O] =  30.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]
        + 120.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]
        + 180.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]
        + 120.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]
        +  30.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        e2vPrime2[H2O][CO2] =  30.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]
        + 120.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]
        + 180.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]
        + 120.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]
        +  30.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        e2vPrime2[CO2][CO2] =  30.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]
        + 120.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]
        + 180.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]
        + 120.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]
        +  30.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt  = 0.0;
        *devdt += deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *devdt +=  6.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *devdt += 15.0*DpowSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *devdt += 20.0*DpowSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *devdt += 15.0*DpowSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *devdt +=  6.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *devdt += deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d2evdt2  = 0.0;
        *d2evdt2 += d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d2evdt2 +=  6.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d2evdt2 += 15.0*D2powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d2evdt2 += 20.0*D2powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d2evdt2 += 15.0*D2powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d2evdt2 +=  6.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d2evdt2 += d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d3evdt3  = 0.0;
        *d3evdt3 += d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d3evdt3 +=  6.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d3evdt3 += 15.0*D3powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d3evdt3 += 20.0*D3powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d3evdt3 += 15.0*D3powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d3evdt3 +=  6.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d3evdt3 += d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        devPrimedt[H2O]  =  0.0;
        devPrimedt[H2O] +=  6.0*deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        devPrimedt[H2O] += 30.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        devPrimedt[H2O] += 60.0*DpowSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        devPrimedt[H2O] += 60.0*DpowSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[H2O] += 30.0*DpowSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[H2O] +=  6.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];;

        devPrimedt[CO2]  =  0.0;
        devPrimedt[CO2] +=  6.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        devPrimedt[CO2] += 30.0*DpowSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        devPrimedt[CO2] += 60.0*DpowSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        devPrimedt[CO2] += 60.0*DpowSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[CO2] += 30.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[CO2] +=  6.0*deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2evPrimedt2[H2O]  =  0.0;
        d2evPrimedt2[H2O] +=  6.0*d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2evPrimedt2[H2O] += 30.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2evPrimedt2[H2O] += 60.0*D2powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2evPrimedt2[H2O] += 60.0*D2powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[H2O] += 30.0*D2powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[H2O] +=  6.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];;

        d2evPrimedt2[CO2]  =  0.0;
        d2evPrimedt2[CO2] +=  6.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2evPrimedt2[CO2] += 30.0*D2powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2evPrimedt2[CO2] += 60.0*D2powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2evPrimedt2[CO2] += 60.0*D2powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[CO2] += 30.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[CO2] +=  6.0*d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d3evPrimedt3[H2O]  =  0.0;
        d3evPrimedt3[H2O] +=  6.0*d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3evPrimedt3[H2O] += 30.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3evPrimedt3[H2O] += 60.0*D3powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3evPrimedt3[H2O] += 60.0*D3powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[H2O] += 30.0*D3powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[H2O] +=  6.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];;

        d3evPrimedt3[CO2]  =  0.0;
        d3evPrimedt3[CO2] +=  6.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3evPrimedt3[CO2] += 30.0*D3powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3evPrimedt3[CO2] += 60.0*D3powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3evPrimedt3[CO2] += 60.0*D3powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[CO2] += 30.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[CO2] +=  6.0*d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

    }

    return;
}

static void FVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *fv, double fvPrime[2], double *dfvdt, double *d2fvdt2, double *d3fvdt3, double dfvPrimedt[2], double d2fvPrimedt2[2], double d3fvPrimedt3[2],
                             double f2vPrime2[2][2]) {
    double fEnd[2], dfEnddt[2], d2fEnddt2[2], d3fEnddt3[2];
    double H2OTr = t/H2OTc;
    double CO2Tr = t/CO2Tc;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;

    if (useLowPcoeff) {
        fEnd[H2O] = H2OLa/H2OTr/H2OTr/H2OTr;
        fEnd[CO2] = CO2La/CO2Tr/CO2Tr/CO2Tr;
        dfEnddt[H2O] = - 3.0*H2OLa*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dfEnddt[CO2] = - 3.0*CO2La*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2fEnddt2[H2O] = 12.0*H2OLa*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2fEnddt2[CO2] = 12.0*CO2La*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3fEnddt3[H2O] = - 60.0*H2OLa*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3fEnddt3[CO2] = - 60.0*CO2La*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    } else {
        fEnd[H2O] = H2OHa/H2OTr/H2OTr/H2OTr;
        fEnd[CO2] = CO2Ha/CO2Tr/CO2Tr/CO2Tr;
        dfEnddt[H2O] = - 3.0*H2OHa*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dfEnddt[CO2] = - 3.0*CO2Ha*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2fEnddt2[H2O] = 12.0*H2OHa*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2fEnddt2[CO2] = 12.0*CO2Ha*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3fEnddt3[H2O] = - 60.0*H2OHa*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3fEnddt3[CO2] = - 60.0*CO2Ha*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *fv = 0.0;
        fvPrime[H2O] = 0.0;
        fvPrime[CO2] = 0.0;

        f2vPrime2[H2O][H2O] = 0.0;
        f2vPrime2[H2O][CO2] = 0.0;
        f2vPrime2[CO2][CO2] = 0.0;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt   = 0.0;
        *d2fvdt2 = 0.0;
        *d3fvdt3 = 0.0;

        dfvPrimedt[H2O] = 0.0;
        dfvPrimedt[CO2] = 0.0;

        d2fvPrimedt2[H2O] = 0.0;
        d2fvPrimedt2[CO2] = 0.0;

        d3fvPrimedt3[H2O] = 0.0;
        d3fvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *fv = fEnd[H2O]*H2OVc*H2OVc;
        fvPrime[H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        fvPrime[CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

        f2vPrime2[H2O][H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        f2vPrime2[H2O][CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        f2vPrime2[CO2][CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt   = dfEnddt[H2O]*H2OVc*H2OVc;
        *d2fvdt2 = d2fEnddt2[H2O]*H2OVc*H2OVc;
        *d3fvdt3 = d3fEnddt3[H2O]*H2OVc*H2OVc;

        dfvPrimedt[H2O] = 2.0*dfEnddt[H2O]*H2OVc*H2OVc;
        dfvPrimedt[CO2] = 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

        d2fvPrimedt2[H2O] = 2.0*d2fEnddt2[H2O]*H2OVc*H2OVc;
        d2fvPrimedt2[CO2] = 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

        d3fvPrimedt3[H2O] = 2.0*d3fEnddt3[H2O]*H2OVc*H2OVc;
        d3fvPrimedt3[CO2] = 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *fv = fEnd[CO2]*CO2Vc*CO2Vc;
        fvPrime[H2O] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        fvPrime[CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;

        f2vPrime2[H2O][H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        f2vPrime2[H2O][CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        f2vPrime2[CO2][CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt   = dfEnddt[CO2]*CO2Vc*CO2Vc;
        *d2fvdt2 = d2fEnddt2[CO2]*CO2Vc*CO2Vc;
        *d3fvdt3 = d3fEnddt3[CO2]*CO2Vc*CO2Vc;

        dfvPrimedt[H2O] = 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        dfvPrimedt[CO2] = 2.0*dfEnddt[CO2]*CO2Vc*CO2Vc;

        d2fvPrimedt2[H2O] = 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        d2fvPrimedt2[CO2] = 2.0*d2fEnddt2[CO2]*CO2Vc*CO2Vc;

        d3fvPrimedt3[H2O] = 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        d3fvPrimedt3[CO2] = 2.0*d3fEnddt3[CO2]*CO2Vc*CO2Vc;

    } else {
        *fv  = 0.0;
        *fv += fEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *fv += 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *fv += fEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        fvPrime[H2O]  = 0.0;
        fvPrime[H2O] += 2.0*fEnd[H2O]*H2OVc*H2OVc*x[H2O];
        fvPrime[H2O] += 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        fvPrime[CO2]  = 0.0;
        fvPrime[CO2] += 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        fvPrime[CO2] += 2.0*fEnd[CO2]*CO2Vc*CO2Vc*x[CO2];

        f2vPrime2[H2O][H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        f2vPrime2[H2O][CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        f2vPrime2[CO2][CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt  = 0.0;
        *dfvdt += dfEnddt[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *dfvdt += 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *dfvdt += dfEnddt[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        *d2fvdt2  = 0.0;
        *d2fvdt2 += d2fEnddt2[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *d2fvdt2 += 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *d2fvdt2 += d2fEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        *d3fvdt3  = 0.0;
        *d3fvdt3 += d3fEnddt3[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *d3fvdt3 += 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *d3fvdt3 += d3fEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        dfvPrimedt[H2O]  = 0.0;
        dfvPrimedt[H2O] += 2.0*dfEnddt[H2O]*H2OVc*H2OVc*x[H2O];
        dfvPrimedt[H2O] += 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        dfvPrimedt[CO2]  = 0.0;
        dfvPrimedt[CO2] += 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        dfvPrimedt[CO2] += 2.0*dfEnddt[CO2]*CO2Vc*CO2Vc*x[CO2];

        d2fvPrimedt2[H2O]  = 0.0;
        d2fvPrimedt2[H2O] += 2.0*d2fEnddt2[H2O]*H2OVc*H2OVc*x[H2O];
        d2fvPrimedt2[H2O] += 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        d2fvPrimedt2[CO2]  = 0.0;
        d2fvPrimedt2[CO2] += 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        d2fvPrimedt2[CO2] += 2.0*d2fEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2];

        d3fvPrimedt3[H2O]  = 0.0;
        d3fvPrimedt3[H2O] += 2.0*d3fEnddt3[H2O]*H2OVc*H2OVc*x[H2O];
        d3fvPrimedt3[H2O] += 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        d3fvPrimedt3[CO2]  = 0.0;
        d3fvPrimedt3[CO2] += 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        d3fvPrimedt3[CO2] += 2.0*d3fEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2];

    }

    return;
}

static void BetaAndDerivative(int useLowPcoeff, double t, double x[2], double *beta, double betaPrime[2]) {
    double betaEnd[2];

    if (useLowPcoeff) {
        betaEnd[H2O] = H2OLb;
        betaEnd[CO2] = CO2Lb;
    } else {
        betaEnd[H2O] = H2OHb;
        betaEnd[CO2] = CO2Hb;
    }

    *beta  = betaEnd[H2O]*x[H2O] + betaEnd[CO2]*x[CO2];

    betaPrime[H2O]  = betaEnd[H2O];
    betaPrime[CO2]  = betaEnd[CO2];

    return;
}

static void GammaVcAndDerivative(int useLowPcoeff, double t, double x[2], double *gammav, double gammavPrime[2], double gamma2vPrime2[2][2]) {
    double gammaEnd[2];
    double k3 = 0.0;

    if (useLowPcoeff) {
        gammaEnd[H2O] = H2OLc;
        gammaEnd[CO2] = CO2Lc;
        k3 = 0.9;
    } else {
        gammaEnd[H2O] = H2OHc;
        gammaEnd[CO2] = CO2Hc;
        k3 = 1.0;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *gammav = 0.0;
        gammavPrime[H2O] = 0.0;
        gammavPrime[CO2] = 0.0;

        gamma2vPrime2[H2O][H2O] = 0.0;
        gamma2vPrime2[H2O][CO2] = 0.0;
        gamma2vPrime2[CO2][CO2] = 0.0;
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *gammav = gammaEnd[H2O]*H2OVc*H2OVc;
        gammavPrime[H2O] = 3.0*gammaEnd[H2O]*H2OVc*H2OVc;
        gammavPrime[CO2] = 3.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        gamma2vPrime2[H2O][H2O] = 6.0*gammaEnd[H2O]*H2OVc*H2OVc;
        gamma2vPrime2[H2O][CO2] = 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        gamma2vPrime2[CO2][CO2] = 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *gammav = gammaEnd[CO2]*CO2Vc*CO2Vc;
        gammavPrime[H2O] = 3.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        gammavPrime[CO2] = 3.0*gammaEnd[CO2]*CO2Vc*CO2Vc;

        gamma2vPrime2[H2O][H2O] = 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        gamma2vPrime2[H2O][CO2] = 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        gamma2vPrime2[CO2][CO2] = 6.0*gammaEnd[CO2]*CO2Vc*CO2Vc;
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    } else {
        *gammav  = 0.0;
        *gammav += gammaEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *gammav += 3.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *gammav += 3.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *gammav += gammaEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        gammavPrime[H2O]  = 0.0;
        gammavPrime[H2O] += 3.0*gammaEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        gammavPrime[H2O] += 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        gammavPrime[H2O] += 3.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        gammavPrime[CO2]  = 0.0;
        gammavPrime[CO2] += 3.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        gammavPrime[CO2] += 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        gammavPrime[CO2] += 3.0*gammaEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        gamma2vPrime2[H2O][H2O] = 6.0*gammaEnd[H2O]*H2OVc*H2OVc*x[H2O]
        + 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[CO2];
        gamma2vPrime2[H2O][CO2] = 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]
        + 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2];
        gamma2vPrime2[CO2][CO2] = 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]
        + 6.0*gammaEnd[CO2]*CO2Vc*CO2Vc*x[CO2];
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    }

    return;
}

#pragma mark -
#pragma mark Ideal gas functions

static void idealGasH2O(double t, double *cp, double *s0, double *h0, double *dcpdt) {
    int i;

    for (i=0, *cp=0.0; i<7; i++) *cp += idealCoeff[i][H2O]*pow(t/1000.0, (double) i);
    for (i=7; i<13; i++)         *cp += idealCoeff[i][H2O]/pow(t/1000.0, (double) (i-6));

    for (i=1, *dcpdt=0.0; i<7; i++) *dcpdt +=  ((double) i)  *idealCoeff[i][H2O]*pow(t/1000.0, (double) i-1);
    for (i=7; i<13; i++)            *dcpdt += -((double) i-6)*idealCoeff[i][H2O]/pow(t/1000.0, (double) (i+1-6));

    for (i=0, *h0=0.0; i<7; i++) *h0 += idealCoeff[i][H2O]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    *h0 += idealCoeff[7][H2O]*log(t/1000.0);
    for (i=8; i<13; i++)         *h0 += idealCoeff[i][H2O]/pow(t/1000.0, (double) (i-7))/((double) (7-i));

    *s0  = idealCoeff[0][H2O]*log(t/1000.0);
    for (i=1; i<7; i++)	       *s0 += idealCoeff[i][H2O]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++)         *s0 += idealCoeff[i][H2O]/pow(t/1000.0, (double) (i-6))/((double) (6-i));

    *cp    *= 8.31451;
    *h0    *= 8.31451*1000.0;
    *s0    *= 8.31451;
    *dcpdt *= 8.31451/1000.0;

    *h0 += - 355665.4136;
    *s0 +=   359.6505;
}

static void idealGasCO2(double t, double *cp, double *s0, double *h0, double *dcpdt) {
    int i;

    for (i=0, *cp=0.0; i<7; i++) *cp += idealCoeff[i][CO2]*pow(t/1000.0, (double) i);
    for (i=7; i<13; i++)         *cp += idealCoeff[i][CO2]/pow(t/1000.0, (double) (i-6));

    for (i=1, *dcpdt=0.0; i<7; i++) *dcpdt +=  ((double) i)  *idealCoeff[i][CO2]*pow(t/1000.0, (double) i-1);
    for (i=7; i<13; i++)            *dcpdt += -((double) i-6)*idealCoeff[i][CO2]/pow(t/1000.0, (double) (i+1-6));

    for (i=0, *h0=0.0; i<7; i++) *h0 += idealCoeff[i][CO2]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    *h0 += idealCoeff[7][CO2]*log(t/1000.0);
    for (i=8; i<13; i++)         *h0 += idealCoeff[i][CO2]/pow(t/1000.0, (double) (i-7))/((double) (7-i));

    *s0  = idealCoeff[0][CO2]*log(t/1000.0);
    for (i=1; i<7; i++)	       *s0 += idealCoeff[i][CO2]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++)         *s0 += idealCoeff[i][CO2]/pow(t/1000.0, (double) (i-6))/((double) (6-i));

    *cp    *= 8.31451;
    *h0    *= 8.31451*1000.0;
    *s0    *= 8.31451;
    *dcpdt *= 8.31451/1000.0;

    *h0 += - 385358.2260;
    *s0 +=   210.0304;
}

#pragma mark -
#pragma mark Duan driver routine (mixing functions)

static void duanDriver(int useLowPcoeff, double t, double p, double x[2],
                       double *vPt, double *zPt, double phi[2], double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                       double dlnphidt[2], double dlnphidp[2], double d2lnphidt2[2], double d2lnphidtdp[2], double d2lnphidp2[2], double dlnphidr[2]) {

    double bv, cv, dv, ev, fv, beta, gammav, v, z = 1.0, dzdv, dzdt, d2zdv2, d2zdvdt, d2zdt2, dzdn[2];
    double dbvdt, dcvdt, ddvdt, devdt, dfvdt, d2bvdt2, d2cvdt2, d2dvdt2, d2evdt2, d2fvdt2, d3bvdt3, d3cvdt3, d3dvdt3, d3evdt3, d3fvdt3;
    double bvPrime[2], cvPrime[2], dvPrime[2], evPrime[2], fvPrime[2], betaPrime[2], gammavPrime[2];
    double b2vPrime2[2][2], c2vPrime2[2][2], d2vPrime2[2][2], e2vPrime2[2][2], f2vPrime2[2][2], gamma2vPrime2[2][2];
    double dbvPrimedt[2], d2bvPrimedt2[2], d3bvPrimedt3[2], dcvPrimedt[2], d2cvPrimedt2[2], d3cvPrimedt3[2],
    ddvPrimedt[2], d2dvPrimedt2[2], d3dvPrimedt3[2], devPrimedt[2], d2evPrimedt2[2], d3evPrimedt3[2],
    dfvPrimedt[2], d2fvPrimedt2[2], d3fvPrimedt3[2];
    double lnPhiH2O, lnPhiCO2, dlnPhiH2Odv, dlnPhiCO2dv, dlnPhiH2Odt, dlnPhiCO2dt;
    double d2lnPhiH2Odv2, d2lnPhiCO2dv2, d2lnPhiH2Odt2, d2lnPhiCO2dt2, d2lnPhiH2Odvdt, d2lnPhiCO2dvdt;
    double dlnPhiH2OdnH2O, dlnPhiH2OdnCO2, dlnPhiCO2dnH2O, dlnPhiCO2dnCO2;

    BVcAndDerivative    (useLowPcoeff, t, x, &bv,     bvPrime,    &dbvdt,    &d2bvdt2,  &d3bvdt3, dbvPrimedt, d2bvPrimedt2, d3bvPrimedt3, b2vPrime2);
    CVcAndDerivative    (useLowPcoeff, t, x, &cv,     cvPrime,    &dcvdt,    &d2cvdt2,  &d3cvdt3, dcvPrimedt, d2cvPrimedt2, d3cvPrimedt3, c2vPrime2);
    DVcAndDerivative    (useLowPcoeff, t, x, &dv,     dvPrime,	&ddvdt,    &d2dvdt2,  &d3dvdt3, ddvPrimedt, d2dvPrimedt2, d3dvPrimedt3, d2vPrime2);
    EVcAndDerivative    (useLowPcoeff, t, x, &ev,     evPrime,	&devdt,    &d2evdt2,  &d3evdt3, devPrimedt, d2evPrimedt2, d3evPrimedt3, e2vPrime2);
    FVcAndDerivative    (useLowPcoeff, t, x, &fv,     fvPrime,	&dfvdt,    &d2fvdt2,  &d3fvdt3, dfvPrimedt, d2fvPrimedt2, d3fvPrimedt3, f2vPrime2);
    BetaAndDerivative   (useLowPcoeff, t, x, &beta,   betaPrime);
    GammaVcAndDerivative(useLowPcoeff, t, x, &gammav, gammavPrime, gamma2vPrime2);

    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0;

        v = 8.314467*t/p;
        while (iter < 200) {
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            if (isnan(z)) {
                printf("duanDriver NAN: t = %g, p = %g, z = %g, v = %g, bv = %g, cv = %g, dv = %g, ev = %g, fv = %g, beta = %g, gammav = %g, iter = %d\n", t, p, z, v, bv, cv, dv, ev, fv, beta, gammav, iter);
                break;
            }
            delv = z*8.314467*t/p - v;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
            vPrevious = v;
            delvPrevious = delv;
            v = (z*8.314467*t/p + v)/2.0;
            iter++;
        }
        if (iter == 200) printf("duanDriver-1: t = %g, p = %g, z = %g, v = %g, vPrev = %g, delv = %g, delvPrev = %g, iter = %d\n", t, p, z, v, vPrevious, delv, delvPrevious, iter);
        else if (fabs(delv) > v*100.0*DBL_EPSILON) {
            double dx;
            double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
            iter = 0;
            while (iter < 200) {
                v = rtb + (dx *= 0.5);
                z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
                delv = z*8.314467*t/p - v;
                if (delv <= 0.0) rtb = v;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
            if ( (iter == 200) || (fabs(dx) > 100.0*DBL_EPSILON) ) printf("duanDriver-2: t = %g, p = %g, z = %g, v = %g, delv = %g, dx = %g, iter = %d\n", t, p, z, v, delv, dx, iter);
        }
    }

    lnPhiH2O  = 0.0;
    lnPhiH2O += -log(z);
    lnPhiH2O += bvPrime[H2O]/v;
    lnPhiH2O += cvPrime[H2O]/2.0/v/v;
    lnPhiH2O += dvPrime[H2O]/4.0/v/v/v/v;
    lnPhiH2O += evPrime[H2O]/5.0/v/v/v/v/v;
    lnPhiH2O += ((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiH2O += ((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)
    *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiH2O += ((gammavPrime[H2O]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    lnPhiCO2  = 0.0;
    lnPhiCO2 += -log(z);
    lnPhiCO2 += bvPrime[CO2]/v;
    lnPhiCO2 += cvPrime[CO2]/2.0/v/v;
    lnPhiCO2 += dvPrime[CO2]/4.0/v/v/v/v;
    lnPhiCO2 += evPrime[CO2]/5.0/v/v/v/v/v;
    lnPhiCO2 += ((fvPrime[CO2]*beta + betaPrime[CO2]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiCO2 += ((fvPrime[CO2]*gammav+gammavPrime[CO2]*fv-fv*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav)
    *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiCO2 += ((gammavPrime[CO2]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    dzdv = - bv/v/v - 2.0*cv/v/v/v + - 4.0*dv/v/v/v/v/v - 5.0*ev/v/v/v/v/v/v
    - 2.0*(fv/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(fv/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(fv/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    dzdt = dbvdt/v + dcvdt/v/v + ddvdt/v/v/v/v + devdt/v/v/v/v/v + (dfvdt/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    dzdn[H2O] = 1.0 - z + bvPrime[H2O]/v + cvPrime[H2O]/v/v + dvPrime[H2O]/v/v/v/v + evPrime[H2O]/v/v/v/v/v
    + (fvPrime[H2O]/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    + (fv/v/v) * (betaPrime[H2O] + gammavPrime[H2O]/v/v) * exp(-gammav/v/v)
    - (fv/v/v/v/v) * (gammavPrime[H2O]-gammav) * (beta + gammav/v/v) * exp(-gammav/v/v);
    dzdn[CO2] = 1.0 - z + bvPrime[CO2]/v + cvPrime[CO2]/v/v + dvPrime[CO2]/v/v/v/v + evPrime[CO2]/v/v/v/v/v
    + (fvPrime[CO2]/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    + (fv/v/v) * (betaPrime[CO2] + gammavPrime[CO2]/v/v) * exp(-gammav/v/v)
    - (fv/v/v/v/v) * (gammavPrime[CO2]-gammav) * (beta + gammav/v/v) * exp(-gammav/v/v);

    { /* Second derivative block */
        double H1, H2, H3, H4, H5, H6;
        double dbvPrimedn, dcvPrimedn, ddvPrimedn, devPrimedn, dfvPrimedn, dgammavPrimedn, dfvdn, dbetadn, dgammavdn;
        double dH1dn, dH2dn, dH3dn, dH4dn, dH5dn, dH6dn;

        /* H2O - H2O */
        H1 = (fvPrime[H2O]*beta + betaPrime[H2O]*fv)/2.0/gammav;
        H2 = exp(-gammav/v/v);
        H3 = (fvPrime[H2O]*gammav + gammavPrime[H2O]*fv - fv*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav;
        H4 = gammav/v/v + 1.0;
        H5 = (gammavPrime[H2O]-gammav)*fv/2.0/gammav/gammav;
        H6 = gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0;

        dbvPrimedn     =         -bvPrime[H2O] +     b2vPrime2[H2O][H2O];
        dcvPrimedn     =     -2.0*cvPrime[H2O] +     c2vPrime2[H2O][H2O];
        ddvPrimedn     =     -4.0*dvPrime[H2O] +     d2vPrime2[H2O][H2O];
        devPrimedn     =     -5.0*evPrime[H2O] +     e2vPrime2[H2O][H2O];
        dfvPrimedn     =         -fvPrime[H2O] +     f2vPrime2[H2O][H2O];
        dgammavPrimedn = -2.0*gammavPrime[H2O] + gamma2vPrime2[H2O][H2O];

        dfvdn     =     fvPrime[H2O] - 2.0*fv;
        dbetadn   =   betaPrime[H2O] - beta;
        dgammavdn = gammavPrime[H2O] - 3.0*gammav;

        dH1dn = (dfvPrimedn*beta + fvPrime[H2O]*dbetadn + betaPrime[H2O]*dfvdn)/2.0/gammav
        - (fvPrime[H2O]*beta + betaPrime[H2O]*fv)*dgammavdn/2.0/gammav/gammav;
        dH2dn = (- 2.0*gammav/v/v - dgammavdn/v/v)*exp(-gammav/v/v);
        dH3dn = (dfvPrimedn*gammav + fvPrime[H2O]*dgammavdn)/2.0/gammav/gammav
        + (dgammavPrimedn*fv + gammavPrime[H2O]*dfvdn)/2.0/gammav/gammav
        - (dfvdn*beta*(gammavPrime[H2O]-gammav) + fv*dbetadn*(gammavPrime[H2O]-gammav) + fv*beta*(dgammavPrimedn-dgammavdn))/2.0/gammav/gammav
        - 2.0*(fvPrime[H2O]*gammav + gammavPrime[H2O]*fv - fv*beta*(gammavPrime[H2O]-gammav))*dgammavdn/2.0/gammav/gammav/gammav;
        dH4dn = 2.0*gammav/v/v + dgammavdn/v/v;
        dH5dn = (dgammavPrimedn-dgammavdn)*fv/2.0/gammav/gammav
        + (gammavPrime[H2O]-gammav)*dfvdn/2.0/gammav/gammav
        - 2.0*(gammavPrime[H2O]-gammav)*fv*dgammavdn/2.0/gammav/gammav/gammav;
        dH6dn = 4.0*gammav*gammav/v/v/v/v + 2.0*gammav*dgammavdn/v/v/v/v + 4.0*gammav/v/v + 2.0*dgammavdn/v/v;

        dlnPhiH2OdnH2O  = 0.0;
        dlnPhiH2OdnH2O += -dzdn[H2O]/z;
        dlnPhiH2OdnH2O += bvPrime[H2O]/v         + dbvPrimedn/v;
        dlnPhiH2OdnH2O += cvPrime[H2O]/v/v       + dcvPrimedn/2.0/v/v;
        dlnPhiH2OdnH2O += dvPrime[H2O]/v/v/v/v   + ddvPrimedn/4.0/v/v/v/v;
        dlnPhiH2OdnH2O += evPrime[H2O]/v/v/v/v/v + devPrimedn/5.0/v/v/v/v/v;
        dlnPhiH2OdnH2O += dH1dn*(1.0-H2) - H1*dH2dn;
        dlnPhiH2OdnH2O += dH3dn*(1.0-H4*H2) - H3*(dH4dn*H2 + H4*dH2dn);
        dlnPhiH2OdnH2O += dH5dn*(H6*H2-2.0) + H5*(dH6dn*H2 + H6*dH2dn);

        /* H2O - CO2 */
        dbvPrimedn     =         -bvPrime[H2O] +     b2vPrime2[H2O][CO2];
        dcvPrimedn     =     -2.0*cvPrime[H2O] +     c2vPrime2[H2O][CO2];
        ddvPrimedn     =     -4.0*dvPrime[H2O] +     d2vPrime2[H2O][CO2];
        devPrimedn     =     -5.0*evPrime[H2O] +     e2vPrime2[H2O][CO2];
        dfvPrimedn     =         -fvPrime[H2O] +     f2vPrime2[H2O][CO2];
        dgammavPrimedn = -2.0*gammavPrime[H2O] + gamma2vPrime2[H2O][CO2];

        dfvdn     =     fvPrime[CO2] - 2.0*fv;
        dbetadn   =   betaPrime[CO2] - beta;
        dgammavdn = gammavPrime[CO2] - 3.0*gammav;

        dH1dn = (dfvPrimedn*beta + fvPrime[H2O]*dbetadn + betaPrime[H2O]*dfvdn)/2.0/gammav
        - (fvPrime[H2O]*beta + betaPrime[H2O]*fv)*dgammavdn/2.0/gammav/gammav;
        dH2dn = (- 2.0*gammav/v/v - dgammavdn/v/v)*exp(-gammav/v/v);
        dH3dn = (dfvPrimedn*gammav + fvPrime[H2O]*dgammavdn)/2.0/gammav/gammav
        + (dgammavPrimedn*fv + gammavPrime[H2O]*dfvdn)/2.0/gammav/gammav
        - (dfvdn*beta*(gammavPrime[H2O]-gammav) + fv*dbetadn*(gammavPrime[H2O]-gammav) + fv*beta*(dgammavPrimedn-dgammavdn))/2.0/gammav/gammav
        - 2.0*(fvPrime[H2O]*gammav + gammavPrime[H2O]*fv - fv*beta*(gammavPrime[H2O]-gammav))*dgammavdn/2.0/gammav/gammav/gammav;
        dH4dn = 2.0*gammav/v/v + dgammavdn/v/v;
        dH5dn = (dgammavPrimedn-dgammavdn)*fv/2.0/gammav/gammav
        + (gammavPrime[H2O]-gammav)*dfvdn/2.0/gammav/gammav
        - 2.0*(gammavPrime[H2O]-gammav)*fv*dgammavdn/2.0/gammav/gammav/gammav;
        dH6dn = 4.0*gammav*gammav/v/v/v/v + 2.0*gammav*dgammavdn/v/v/v/v + 4.0*gammav/v/v + 2.0*dgammavdn/v/v;

        dlnPhiH2OdnCO2  = 0.0;
        dlnPhiH2OdnCO2 += -dzdn[CO2]/z;
        dlnPhiH2OdnCO2 += bvPrime[H2O]/v         + dbvPrimedn/v;
        dlnPhiH2OdnCO2 += cvPrime[H2O]/v/v       + dcvPrimedn/2.0/v/v;
        dlnPhiH2OdnCO2 += dvPrime[H2O]/v/v/v/v   + ddvPrimedn/4.0/v/v/v/v;
        dlnPhiH2OdnCO2 += evPrime[H2O]/v/v/v/v/v + devPrimedn/5.0/v/v/v/v/v;
        dlnPhiH2OdnCO2 += dH1dn*(1.0-H2) - H1*dH2dn;
        dlnPhiH2OdnCO2 += dH3dn*(1.0-H4*H2) - H3*(dH4dn*H2 + H4*dH2dn);
        dlnPhiH2OdnCO2 += dH5dn*(H6*H2-2.0) + H5*(dH6dn*H2 + H6*dH2dn);

        /* CO2 - H2O */
        H1 = (fvPrime[CO2]*beta + betaPrime[CO2]*fv)/2.0/gammav;
        H2 = exp(-gammav/v/v);
        H3 = (fvPrime[CO2]*gammav + gammavPrime[CO2]*fv - fv*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav;
        H4 = gammav/v/v + 1.0;
        H5 = (gammavPrime[CO2]-gammav)*fv/2.0/gammav/gammav;
        H6 = gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0;

        dbvPrimedn     =         -bvPrime[CO2] +     b2vPrime2[CO2][H2O];
        dcvPrimedn     =     -2.0*cvPrime[CO2] +     c2vPrime2[CO2][H2O];
        ddvPrimedn     =     -4.0*dvPrime[CO2] +     d2vPrime2[CO2][H2O];
        devPrimedn     =     -5.0*evPrime[CO2] +     e2vPrime2[CO2][H2O];
        dfvPrimedn     =         -fvPrime[CO2] +     f2vPrime2[CO2][H2O];
        dgammavPrimedn = -2.0*gammavPrime[CO2] + gamma2vPrime2[CO2][H2O];

        dfvdn     =     fvPrime[H2O] - 2.0*fv;
        dbetadn   =   betaPrime[H2O] - beta;
        dgammavdn = gammavPrime[H2O] - 3.0*gammav;

        dH1dn = (dfvPrimedn*beta + fvPrime[CO2]*dbetadn + betaPrime[CO2]*dfvdn)/2.0/gammav
        - (fvPrime[CO2]*beta + betaPrime[CO2]*fv)*dgammavdn/2.0/gammav/gammav;
        dH2dn = (- 2.0*gammav/v/v - dgammavdn/v/v)*exp(-gammav/v/v);
        dH3dn = (dfvPrimedn*gammav + fvPrime[CO2]*dgammavdn)/2.0/gammav/gammav
        + (dgammavPrimedn*fv + gammavPrime[CO2]*dfvdn)/2.0/gammav/gammav
        - (dfvdn*beta*(gammavPrime[CO2]-gammav) + fv*dbetadn*(gammavPrime[CO2]-gammav) + fv*beta*(dgammavPrimedn-dgammavdn))/2.0/gammav/gammav
        - 2.0*(fvPrime[CO2]*gammav + gammavPrime[CO2]*fv - fv*beta*(gammavPrime[CO2]-gammav))*dgammavdn/2.0/gammav/gammav/gammav;
        dH4dn = 2.0*gammav/v/v + dgammavdn/v/v;
        dH5dn = (dgammavPrimedn-dgammavdn)*fv/2.0/gammav/gammav
        + (gammavPrime[CO2]-gammav)*dfvdn/2.0/gammav/gammav
        - 2.0*(gammavPrime[CO2]-gammav)*fv*dgammavdn/2.0/gammav/gammav/gammav;
        dH6dn = 4.0*gammav*gammav/v/v/v/v + 2.0*gammav*dgammavdn/v/v/v/v + 4.0*gammav/v/v + 2.0*dgammavdn/v/v;

        dlnPhiCO2dnH2O  = 0.0;
        dlnPhiCO2dnH2O += -dzdn[H2O]/z;
        dlnPhiCO2dnH2O += bvPrime[CO2]/v         + dbvPrimedn/v;
        dlnPhiCO2dnH2O += cvPrime[CO2]/v/v       + dcvPrimedn/2.0/v/v;
        dlnPhiCO2dnH2O += dvPrime[CO2]/v/v/v/v   + ddvPrimedn/4.0/v/v/v/v;
        dlnPhiCO2dnH2O += evPrime[CO2]/v/v/v/v/v + devPrimedn/5.0/v/v/v/v/v;
        dlnPhiCO2dnH2O += dH1dn*(1.0-H2) - H1*dH2dn;
        dlnPhiCO2dnH2O += dH3dn*(1.0-H4*H2) - H3*(dH4dn*H2 + H4*dH2dn);
        dlnPhiCO2dnH2O += dH5dn*(H6*H2-2.0) + H5*(dH6dn*H2 + H6*dH2dn);

        /* CO2 - CO2 */
        dbvPrimedn     =         -bvPrime[CO2] +     b2vPrime2[CO2][CO2];
        dcvPrimedn     =     -2.0*cvPrime[CO2] +     c2vPrime2[CO2][CO2];
        ddvPrimedn     =     -4.0*dvPrime[CO2] +     d2vPrime2[CO2][CO2];
        devPrimedn     =     -5.0*evPrime[CO2] +     e2vPrime2[CO2][CO2];
        dfvPrimedn     =         -fvPrime[CO2] +     f2vPrime2[CO2][CO2];
        dgammavPrimedn = -2.0*gammavPrime[CO2] + gamma2vPrime2[CO2][CO2];

        dfvdn     =     fvPrime[CO2] - 2.0*fv;
        dbetadn   =   betaPrime[CO2] - beta;
        dgammavdn = gammavPrime[CO2] - 3.0*gammav;

        dH1dn = (dfvPrimedn*beta + fvPrime[CO2]*dbetadn + betaPrime[CO2]*dfvdn)/2.0/gammav
        - (fvPrime[CO2]*beta + betaPrime[CO2]*fv)*dgammavdn/2.0/gammav/gammav;
        dH2dn = (- 2.0*gammav/v/v - dgammavdn/v/v)*exp(-gammav/v/v);
        dH3dn = (dfvPrimedn*gammav + fvPrime[CO2]*dgammavdn)/2.0/gammav/gammav
        + (dgammavPrimedn*fv + gammavPrime[CO2]*dfvdn)/2.0/gammav/gammav
        - (dfvdn*beta*(gammavPrime[CO2]-gammav) + fv*dbetadn*(gammavPrime[CO2]-gammav) + fv*beta*(dgammavPrimedn-dgammavdn))/2.0/gammav/gammav
        - 2.0*(fvPrime[CO2]*gammav + gammavPrime[CO2]*fv - fv*beta*(gammavPrime[CO2]-gammav))*dgammavdn/2.0/gammav/gammav/gammav;
        dH4dn = 2.0*gammav/v/v + dgammavdn/v/v;
        dH5dn = (dgammavPrimedn-dgammavdn)*fv/2.0/gammav/gammav
        + (gammavPrime[CO2]-gammav)*dfvdn/2.0/gammav/gammav
        - 2.0*(gammavPrime[CO2]-gammav)*fv*dgammavdn/2.0/gammav/gammav/gammav;
        dH6dn = 4.0*gammav*gammav/v/v/v/v + 2.0*gammav*dgammavdn/v/v/v/v + 4.0*gammav/v/v + 2.0*dgammavdn/v/v;

        dlnPhiCO2dnCO2  = 0.0;
        dlnPhiCO2dnCO2 += -dzdn[CO2]/z;
        dlnPhiCO2dnCO2 += bvPrime[CO2]/v         + dbvPrimedn/v;
        dlnPhiCO2dnCO2 += cvPrime[CO2]/v/v       + dcvPrimedn/2.0/v/v;
        dlnPhiCO2dnCO2 += dvPrime[CO2]/v/v/v/v   + ddvPrimedn/4.0/v/v/v/v;
        dlnPhiCO2dnCO2 += evPrime[CO2]/v/v/v/v/v + devPrimedn/5.0/v/v/v/v/v;
        dlnPhiCO2dnCO2 += dH1dn*(1.0-H2) - H1*dH2dn;
        dlnPhiCO2dnCO2 += dH3dn*(1.0-H4*H2) - H3*(dH4dn*H2 + H4*dH2dn);
        dlnPhiCO2dnCO2 += dH5dn*(H6*H2-2.0) + H5*(dH6dn*H2 + H6*dH2dn);
    }

    dlnPhiH2Odv  = 0.0;
    dlnPhiH2Odv += -dzdv/z;
    dlnPhiH2Odv += -bvPrime[H2O]/v/v;
    dlnPhiH2Odv += -cvPrime[H2O]/v/v/v;
    dlnPhiH2Odv += -dvPrime[H2O]/v/v/v/v/v;
    dlnPhiH2Odv += -evPrime[H2O]/v/v/v/v/v/v;
    dlnPhiH2Odv += -((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv += -((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv +=  ((gammavPrime[H2O]-gammav)*fv*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    dlnPhiCO2dv  = 0.0;
    dlnPhiCO2dv += -dzdv/z;
    dlnPhiCO2dv += -bvPrime[CO2]/v/v;
    dlnPhiCO2dv += -cvPrime[CO2]/v/v/v;
    dlnPhiCO2dv += -dvPrime[CO2]/v/v/v/v/v;
    dlnPhiCO2dv += -evPrime[CO2]/v/v/v/v/v/v;
    dlnPhiCO2dv += -((fvPrime[CO2]*beta + betaPrime[CO2]*fv)/v/v/v)*exp(-gammav/v/v);
    dlnPhiCO2dv += -((fvPrime[CO2]*gammav+gammavPrime[CO2]*fv-fv*beta*(gammavPrime[CO2]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    dlnPhiCO2dv +=  ((gammavPrime[CO2]-gammav)*fv*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    dlnPhiH2Odt  = 0.0;
    dlnPhiH2Odt += -dzdt/z;
    dlnPhiH2Odt += dbvPrimedt[H2O]/v;
    dlnPhiH2Odt += dcvPrimedt[H2O]/2.0/v/v;
    dlnPhiH2Odt += ddvPrimedt[H2O]/4.0/v/v/v/v;
    dlnPhiH2Odt += devPrimedt[H2O]/5.0/v/v/v/v/v;
    dlnPhiH2Odt += ((dfvPrimedt[H2O]*beta + betaPrime[H2O]*dfvdt)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    dlnPhiH2Odt += ((dfvPrimedt[H2O]*gammav+gammavPrime[H2O]*dfvdt-dfvdt*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    dlnPhiH2Odt += ((gammavPrime[H2O]-gammav)*dfvdt/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    dlnPhiCO2dt  = 0.0;
    dlnPhiCO2dt += -dzdt/z;
    dlnPhiCO2dt += dbvPrimedt[CO2]/v;
    dlnPhiCO2dt += dcvPrimedt[CO2]/2.0/v/v;
    dlnPhiCO2dt += ddvPrimedt[CO2]/4.0/v/v/v/v;
    dlnPhiCO2dt += devPrimedt[CO2]/5.0/v/v/v/v/v;
    dlnPhiCO2dt += ((dfvPrimedt[CO2]*beta + betaPrime[CO2]*dfvdt)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    dlnPhiCO2dt += ((dfvPrimedt[CO2]*gammav+gammavPrime[CO2]*dfvdt-dfvdt*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    dlnPhiCO2dt += ((gammavPrime[CO2]-gammav)*dfvdt/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    d2zdv2 =  2.0*bv/v/v/v + 6.0*cv/v/v/v/v + 20.0*dv/v/v/v/v/v/v + 30.0*ev/v/v/v/v/v/v/v
    +  6.0*(fv/v/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    + 14.0*(fv/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    - 14.0*(fv/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    -  8.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    +  4.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
    d2zdvdt = - dbvdt/v/v - 2.0*dcvdt/v/v/v + - 4.0*ddvdt/v/v/v/v/v - 5.0*devdt/v/v/v/v/v/v
    - 2.0*(dfvdt/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(dfvdt/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(dfvdt/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    d2zdt2 = d2bvdt2/v + d2cvdt2/v/v + d2dvdt2/v/v/v/v + d2evdt2/v/v/v/v/v + (d2fvdt2/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    d2lnPhiH2Odv2  = 0.0;
    d2lnPhiH2Odv2 += dzdv*dzdv/z/z - d2zdv2/z;
    d2lnPhiH2Odv2 +=  2.0*bvPrime[H2O]/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*cvPrime[H2O]/v/v/v/v;
    d2lnPhiH2Odv2 +=  5.0*dvPrime[H2O]/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  6.0*evPrime[H2O]/v/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -2.0*((fvPrime[H2O]*beta*gammav + betaPrime[H2O]*fv*gammav)/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  5.0*((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -2.0*((fvPrime[H2O]*gammav*gammav+gammavPrime[H2O]*fv*gammav-fv*beta*gammav*(gammavPrime[H2O]-gammav))/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -7.0*((gammavPrime[H2O]-gammav)*fv*gammav/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  2.0*((gammavPrime[H2O]-gammav)*fv*gammav*gammav/v/v/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiCO2dv2  = 0.0;
    d2lnPhiCO2dv2 += dzdv*dzdv/z/z - d2zdv2/z;
    d2lnPhiCO2dv2 +=  2.0*bvPrime[CO2]/v/v/v;
    d2lnPhiCO2dv2 +=  3.0*cvPrime[CO2]/v/v/v/v;
    d2lnPhiCO2dv2 +=  5.0*dvPrime[CO2]/v/v/v/v/v/v;
    d2lnPhiCO2dv2 +=  6.0*evPrime[CO2]/v/v/v/v/v/v/v;
    d2lnPhiCO2dv2 +=  3.0*((fvPrime[CO2]*beta + betaPrime[CO2]*fv)/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 += -2.0*((fvPrime[CO2]*beta*gammav + betaPrime[CO2]*fv*gammav)/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 +=  5.0*((fvPrime[CO2]*gammav+gammavPrime[CO2]*fv-fv*beta*(gammavPrime[CO2]-gammav))/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 += -2.0*((fvPrime[CO2]*gammav*gammav+gammavPrime[CO2]*fv*gammav-fv*beta*gammav*(gammavPrime[CO2]-gammav))/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 += -7.0*((gammavPrime[CO2]-gammav)*fv*gammav/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 +=  2.0*((gammavPrime[CO2]-gammav)*fv*gammav*gammav/v/v/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odvdt  = 0.0;
    d2lnPhiH2Odvdt += dzdv*dzdt/z/z -d2zdvdt/z;
    d2lnPhiH2Odvdt += -dbvPrimedt[H2O]/v/v;
    d2lnPhiH2Odvdt += -dcvPrimedt[H2O]/v/v/v;
    d2lnPhiH2Odvdt += -ddvPrimedt[H2O]/v/v/v/v/v;
    d2lnPhiH2Odvdt += -devPrimedt[H2O]/v/v/v/v/v/v;
    d2lnPhiH2Odvdt += -((dfvPrimedt[H2O]*beta + betaPrime[H2O]*dfvdt)/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt += -((dfvPrimedt[H2O]*gammav+gammavPrime[H2O]*dfvdt-dfvdt*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt +=  ((gammavPrime[H2O]-gammav)*dfvdt*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiCO2dvdt  = 0.0;
    d2lnPhiCO2dvdt += dzdv*dzdt/z/z -d2zdvdt/z;
    d2lnPhiCO2dvdt += -dbvPrimedt[CO2]/v/v;
    d2lnPhiCO2dvdt += -dcvPrimedt[CO2]/v/v/v;
    d2lnPhiCO2dvdt += -ddvPrimedt[CO2]/v/v/v/v/v;
    d2lnPhiCO2dvdt += -devPrimedt[CO2]/v/v/v/v/v/v;
    d2lnPhiCO2dvdt += -((dfvPrimedt[CO2]*beta + betaPrime[CO2]*dfvdt)/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dvdt += -((dfvPrimedt[CO2]*gammav+gammavPrime[CO2]*dfvdt-dfvdt*beta*(gammavPrime[CO2]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dvdt +=  ((gammavPrime[CO2]-gammav)*dfvdt*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odt2  = 0.0;
    d2lnPhiH2Odt2 += dzdt*dzdt/z/z - d2zdt2/z;
    d2lnPhiH2Odt2 += d2bvPrimedt2[H2O]/v;
    d2lnPhiH2Odt2 += d2cvPrimedt2[H2O]/2.0/v/v;
    d2lnPhiH2Odt2 += d2dvPrimedt2[H2O]/4.0/v/v/v/v;
    d2lnPhiH2Odt2 += d2evPrimedt2[H2O]/5.0/v/v/v/v/v;
    d2lnPhiH2Odt2 += ((d2fvPrimedt2[H2O]*beta + betaPrime[H2O]*d2fvdt2)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    d2lnPhiH2Odt2 += ((d2fvPrimedt2[H2O]*gammav+gammavPrime[H2O]*d2fvdt2-d2fvdt2*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    d2lnPhiH2Odt2 += ((gammavPrime[H2O]-gammav)*d2fvdt2/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    d2lnPhiCO2dt2  = 0.0;
    d2lnPhiCO2dt2 += dzdt*dzdt/z/z - d2zdt2/z;
    d2lnPhiCO2dt2 += d2bvPrimedt2[CO2]/v;
    d2lnPhiCO2dt2 += d2cvPrimedt2[CO2]/2.0/v/v;
    d2lnPhiCO2dt2 += d2dvPrimedt2[CO2]/4.0/v/v/v/v;
    d2lnPhiCO2dt2 += d2evPrimedt2[CO2]/5.0/v/v/v/v/v;
    d2lnPhiCO2dt2 += ((d2fvPrimedt2[CO2]*beta + betaPrime[CO2]*d2fvdt2)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    d2lnPhiCO2dt2 += ((d2fvPrimedt2[CO2]*gammav+gammavPrime[CO2]*d2fvdt2-d2fvdt2*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    d2lnPhiCO2dt2 += ((gammavPrime[CO2]-gammav)*d2fvdt2/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    *vPt = v;
    *zPt = z;
    phi[H2O] = exp(lnPhiH2O);
    phi[CO2] = exp(lnPhiCO2);
    *dvdp = 1.0/( p*(dzdv/z - 1.0/v) );
    *dvdt = (1.0/t + dzdt/z)/(1.0/v - dzdv/z);
    *d2vdp2  = p*(1.0/v-dzdv/z)/v - dzdv/(*dvdp)/z + 1.0/(*dvdp)/(*dvdp)/p + p*d2zdv2/z;
    *d2vdp2 *= -(*dvdp)*(*dvdp)*(*dvdp);
    *d2vdtdp = -(*dvdp)*(1.0/t + dzdt/z) - p*(*d2vdp2)*(1.0/t + dzdt/z) - p*(*dvdp)*(*dvdp)*(-dzdv*dzdt/z/z + d2zdvdt/z);
    *d2vdt2  = -p*(*d2vdtdp)*(1.0/t + dzdt/z) + p*(*dvdp)/t/t + p*(*dvdp)*dzdt*(dzdv*(*dvdt) + dzdt)/z/z - p*(*dvdp)*(*dvdt)*d2zdvdt/z - p*(*dvdp)*d2zdt2/z;

    dlnphidt[H2O] = dlnPhiH2Odv*(*dvdt) + dlnPhiH2Odt;
    dlnphidt[CO2] = dlnPhiCO2dv*(*dvdt) + dlnPhiCO2dt;

    dlnphidp[H2O] = dlnPhiH2Odv*(*dvdp);
    dlnphidp[CO2] = dlnPhiCO2dv*(*dvdp);

    d2lnphidt2[H2O] = d2lnPhiH2Odv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiH2Odvdt*(*dvdt) + dlnPhiH2Odv*(*d2vdt2) + d2lnPhiH2Odt2;
    d2lnphidt2[CO2] = d2lnPhiCO2dv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiCO2dvdt*(*dvdt) + dlnPhiCO2dv*(*d2vdt2) + d2lnPhiCO2dt2;

    d2lnphidtdp[H2O] = d2lnPhiH2Odv2*(*dvdt)*(*dvdp) + dlnPhiH2Odv*(*d2vdtdp) + d2lnPhiH2Odvdt*(*dvdp);
    d2lnphidtdp[CO2] = d2lnPhiCO2dv2*(*dvdt)*(*dvdp) + dlnPhiCO2dv*(*d2vdtdp) + d2lnPhiCO2dvdt*(*dvdp);

    d2lnphidp2[H2O] = d2lnPhiH2Odv2*(*dvdp)*(*dvdp) + dlnPhiH2Odv*(*d2vdp2);
    d2lnphidp2[CO2] = d2lnPhiCO2dv2*(*dvdp)*(*dvdp) + dlnPhiCO2dv*(*d2vdp2);

    /* dV/dr = R*t*(dlnphidp[CO2] - dlnphidp[H2O]) */
    dlnphidr[H2O] = dlnPhiH2Odv*R*t*(dlnphidp[CO2] - dlnphidp[H2O]) - dlnPhiH2OdnH2O + dlnPhiH2OdnCO2;
    dlnphidr[CO2] = dlnPhiCO2dv*R*t*(dlnphidp[CO2] - dlnphidp[H2O]) - dlnPhiCO2dnH2O + dlnPhiCO2dnCO2;
}

#pragma mark -
#pragma mark Duan entry function (mixing properties)

static void duan(int useLowPcoeff, double t, double p, double x[2],
                 double *vPt, double *zPt, double phi[2], double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                 double dlnphidt[2], double dlnphidp[2], double d2lnphidt2[2], double d2lnphidtdp[2], double d2lnphidp2[2], double dlnphidr[2]) {

    duanDriver(useLowPcoeff, t, p, x, vPt, zPt, phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);

    if (!useLowPcoeff) {
        double vRef, zRef, phiRef[2], dvdtRef, dvdpRef, d2vdt2Ref, d2vdtdpRef, d2vdp2Ref, dlnphidtRef[2], dlnphidpRef[2], d2lnphidt2Ref[2], d2lnphidtdpRef[2], d2lnphidp2Ref[2], dlnphidrRef[2];

        duanDriver(1, t, 2000.0, x, &vRef, &zRef, phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, dlnphidtRef, dlnphidpRef, d2lnphidt2Ref, d2lnphidtdpRef, d2lnphidp2Ref, dlnphidrRef);
        phi[H2O]         *= phiRef[H2O];
        phi[CO2]         *= phiRef[CO2];
        dlnphidt[H2O]    += dlnphidtRef[H2O];
        dlnphidt[CO2]    += dlnphidtRef[CO2];
        d2lnphidt2[H2O]  += d2lnphidt2Ref[H2O];
        d2lnphidt2[CO2]  += d2lnphidt2Ref[CO2];
        dlnphidr[H2O]    += dlnphidrRef[H2O];
        dlnphidr[CO2]    += dlnphidrRef[CO2];

        duanDriver(0, t, 2000.0, x, &vRef, &zRef, phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, dlnphidtRef, dlnphidpRef, d2lnphidt2Ref, d2lnphidtdpRef, d2lnphidp2Ref, dlnphidrRef);
        phi[H2O]         /= phiRef[H2O];
        phi[CO2]         /= phiRef[CO2];
        dlnphidt[H2O]    -= dlnphidtRef[H2O];
        dlnphidt[CO2]    -= dlnphidtRef[CO2];
        d2lnphidt2[H2O]  -= d2lnphidt2Ref[H2O];
        d2lnphidt2[CO2]  -= d2lnphidt2Ref[CO2];
        dlnphidr[H2O]    -= dlnphidrRef[H2O];
        dlnphidr[CO2]    -= dlnphidrRef[CO2];
    }

}

#pragma mark -
#pragma mark Duan driver routine (pure water)

static void duanH2ODriver(int useLowPcoeff, double t, double p, double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                          double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {
    double bv, cv, dv, ev, fv, beta, gammav, v, z = 1.0, dzdv, dzdt, d2zdv2, d2zdvdt, d2zdt2;
    double dbvdt, dcvdt, ddvdt, devdt, dfvdt, d2bvdt2, d2cvdt2, d2dvdt2, d2evdt2, d2fvdt2, d3bvdt3, d3cvdt3, d3dvdt3, d3evdt3, d3fvdt3;
    double bvPrime[2], cvPrime[2], dvPrime[2], evPrime[2], fvPrime[2], betaPrime[2], gammavPrime[2];
    double b2vPrime2[2][2], c2vPrime2[2][2], d2vPrime2[2][2], e2vPrime2[2][2], f2vPrime2[2][2], gamma2vPrime2[2][2];
    double dbvPrimedt[2], d2bvPrimedt2[2], d3bvPrimedt3[2], dcvPrimedt[2], d2cvPrimedt2[2], d3cvPrimedt3[2],
    ddvPrimedt[2], d2dvPrimedt2[2], d3dvPrimedt3[2], devPrimedt[2], d2evPrimedt2[2], d3evPrimedt3[2],
    dfvPrimedt[2], d2fvPrimedt2[2], d3fvPrimedt3[2];
    double lnPhiH2O, dlnPhiH2Odv, dlnPhiH2Odt, d2lnPhiH2Odv2, d2lnPhiH2Odt2, d2lnPhiH2Odvdt;
    double x[2] = { 1.0, 0.0 }; /* H2O */

    BVcAndDerivative    (useLowPcoeff, t, x, &bv,     bvPrime,    &dbvdt,    &d2bvdt2,  &d3bvdt3, dbvPrimedt, d2bvPrimedt2, d3bvPrimedt3, b2vPrime2);
    CVcAndDerivative    (useLowPcoeff, t, x, &cv,     cvPrime,    &dcvdt,    &d2cvdt2,  &d3cvdt3, dcvPrimedt, d2cvPrimedt2, d3cvPrimedt3, c2vPrime2);
    DVcAndDerivative    (useLowPcoeff, t, x, &dv,     dvPrime,	&ddvdt,    &d2dvdt2,  &d3dvdt3, ddvPrimedt, d2dvPrimedt2, d3dvPrimedt3, d2vPrime2);
    EVcAndDerivative    (useLowPcoeff, t, x, &ev,     evPrime,	&devdt,    &d2evdt2,  &d3evdt3, devPrimedt, d2evPrimedt2, d3evPrimedt3, e2vPrime2);
    FVcAndDerivative    (useLowPcoeff, t, x, &fv,     fvPrime,	&dfvdt,    &d2fvdt2,  &d3fvdt3, dfvPrimedt, d2fvPrimedt2, d3fvPrimedt3, f2vPrime2);
    BetaAndDerivative   (useLowPcoeff, t, x, &beta,   betaPrime);
    GammaVcAndDerivative(useLowPcoeff, t, x, &gammav, gammavPrime, gamma2vPrime2);

    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0, vLowest = 1.7, vHighest = 1.0e6; // 8/21/18 - MSG

        v = 8.314467*t/p;
        // if block added 8/21/18 - MSG
        if ((p > 1.0) && (p < 220.0) && (t > 273.15) && (t < 646.15)) {
            static const double psat_coeff[5] = { 1.44021565e+00, -2.75944904e-02,  3.50602876e-04, -2.44834016e-06,
                                                  1.57085668e-08 };
            static const double vsat_coeff[6] = { 1.84252342e+01, -3.06586710e-02,  5.65750627e-04, -2.69937313e-06,
                                                  4.67555414e-09, -8.89632469e+00 };
            double x = t - 273.15;
            double psat = psat_coeff[0] + psat_coeff[1]*x + psat_coeff[2]*x*x + psat_coeff[3]*x*x*x + psat_coeff[4]*x*x*x*x;
            if (p > psat) vHighest = (vsat_coeff[0] + vsat_coeff[1]*x + vsat_coeff[2]*x*x + vsat_coeff[3]*x*x*x
                + vsat_coeff[4]*x*x*x*x)/10.0 + vsat_coeff[5]/(x-374.0);
        }
        while (iter < 200) {
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            delv = z*8.314467*t/p - v;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
            vPrevious = v;
            delvPrevious = delv;
            v = (z*8.314467*t/p + v)/2.0;
            if (v < 1.0) v = vLowest;
            if (v > vHighest) v = vHighest; // 8/21/18 - MSG
            iter++;
            if (iter == 200) {
                printf("duanH2ODriver(a): t = %g, p = %g, z = %g, v = %g, vPrev = %g, delv = %g, delvPrev = %g, iter = %d\n", t, p, z, v, vPrevious, delv, delvPrevious, iter);
                if (fabs(v-vLowest) <= 100.0*DBL_EPSILON) {
                    vLowest *= 0.90;
                    v = 8.314467*t/p;
                    iter = 0;
                }
            }
        }
        if (fabs(delv) > v*100.0*DBL_EPSILON) {
            double dx;
            double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
            iter = 0;
            while (iter < 200) {
                v = rtb + (dx *= 0.5);
                z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
                delv = z*8.314467*t/p - v;
                if (delv <= 0.0) rtb = v;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
            if ( (iter == 200) || (fabs(dx) > 100.0*DBL_EPSILON) ) printf("duanH2ODriver(b): t = %g, p = %g, z = %g, v = %g, delv = %g, dx = %g, iter = %d\n", t, p, z, v, delv, dx, iter);
        }
    }

    lnPhiH2O  = 0.0;
    lnPhiH2O += -log(z);
    lnPhiH2O += bvPrime[H2O]/v;
    lnPhiH2O += cvPrime[H2O]/2.0/v/v;
    lnPhiH2O += dvPrime[H2O]/4.0/v/v/v/v;
    lnPhiH2O += evPrime[H2O]/5.0/v/v/v/v/v;
    lnPhiH2O += ((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiH2O += ((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)
    *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiH2O += ((gammavPrime[H2O]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    dzdv = - bv/v/v - 2.0*cv/v/v/v + - 4.0*dv/v/v/v/v/v - 5.0*ev/v/v/v/v/v/v
    - 2.0*(fv/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(fv/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(fv/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    dzdt = dbvdt/v + dcvdt/v/v + ddvdt/v/v/v/v + devdt/v/v/v/v/v + (dfvdt/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    dlnPhiH2Odv  = 0.0;
    dlnPhiH2Odv += -dzdv/z;
    dlnPhiH2Odv += -bvPrime[H2O]/v/v;
    dlnPhiH2Odv += -cvPrime[H2O]/v/v/v;
    dlnPhiH2Odv += -dvPrime[H2O]/v/v/v/v/v;
    dlnPhiH2Odv += -evPrime[H2O]/v/v/v/v/v/v;
    dlnPhiH2Odv += -((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv += -((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv += +((gammavPrime[H2O]-gammav)*fv*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    dlnPhiH2Odt  = 0.0;
    dlnPhiH2Odt += -dzdt/z;
    dlnPhiH2Odt += dbvPrimedt[H2O]/v;
    dlnPhiH2Odt += dcvPrimedt[H2O]/2.0/v/v;
    dlnPhiH2Odt += ddvPrimedt[H2O]/4.0/v/v/v/v;
    dlnPhiH2Odt += devPrimedt[H2O]/5.0/v/v/v/v/v;
    dlnPhiH2Odt += ((dfvPrimedt[H2O]*beta + betaPrime[H2O]*dfvdt)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    dlnPhiH2Odt += ((dfvPrimedt[H2O]*gammav+gammavPrime[H2O]*dfvdt-dfvdt*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    dlnPhiH2Odt += ((gammavPrime[H2O]-gammav)*dfvdt/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    d2zdv2 =  2.0*bv/v/v/v + 6.0*cv/v/v/v/v + 20.0*dv/v/v/v/v/v/v + 30.0*ev/v/v/v/v/v/v/v
    +  6.0*(fv/v/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    + 14.0*(fv/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    - 14.0*(fv/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    -  8.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    +  4.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
    d2zdvdt = - dbvdt/v/v - 2.0*dcvdt/v/v/v + - 4.0*ddvdt/v/v/v/v/v - 5.0*devdt/v/v/v/v/v/v
    - 2.0*(dfvdt/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(dfvdt/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(dfvdt/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    d2zdt2 = d2bvdt2/v + d2cvdt2/v/v + d2dvdt2/v/v/v/v + d2evdt2/v/v/v/v/v + (d2fvdt2/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    d2lnPhiH2Odv2  = 0.0;
    d2lnPhiH2Odv2 += dzdv*dzdv/z/z - d2zdv2/z;
    d2lnPhiH2Odv2 +=  2.0*bvPrime[H2O]/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*cvPrime[H2O]/v/v/v/v;
    d2lnPhiH2Odv2 +=  5.0*dvPrime[H2O]/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  6.0*evPrime[H2O]/v/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -2.0*((fvPrime[H2O]*beta*gammav + betaPrime[H2O]*fv*gammav)/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  5.0*((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -2.0*((fvPrime[H2O]*gammav*gammav+gammavPrime[H2O]*fv*gammav-fv*beta*gammav*(gammavPrime[H2O]-gammav))/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -7.0*((gammavPrime[H2O]-gammav)*fv*gammav/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  2.0*((gammavPrime[H2O]-gammav)*fv*gammav*gammav/v/v/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odvdt  = 0.0;
    d2lnPhiH2Odvdt += dzdv*dzdt/z/z -d2zdvdt/z;
    d2lnPhiH2Odvdt += -dbvPrimedt[H2O]/v/v;
    d2lnPhiH2Odvdt += -dcvPrimedt[H2O]/v/v/v;
    d2lnPhiH2Odvdt += -ddvPrimedt[H2O]/v/v/v/v/v;
    d2lnPhiH2Odvdt += -devPrimedt[H2O]/v/v/v/v/v/v;
    d2lnPhiH2Odvdt += -((dfvPrimedt[H2O]*beta + betaPrime[H2O]*dfvdt)/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt += -((dfvPrimedt[H2O]*gammav+gammavPrime[H2O]*dfvdt-dfvdt*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt +=  ((gammavPrime[H2O]-gammav)*dfvdt*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odt2  = 0.0;
    d2lnPhiH2Odt2 += dzdt*dzdt/z/z - d2zdt2/z;
    d2lnPhiH2Odt2 += d2bvPrimedt2[H2O]/v;
    d2lnPhiH2Odt2 += d2cvPrimedt2[H2O]/2.0/v/v;
    d2lnPhiH2Odt2 += d2dvPrimedt2[H2O]/4.0/v/v/v/v;
    d2lnPhiH2Odt2 += d2evPrimedt2[H2O]/5.0/v/v/v/v/v;
    d2lnPhiH2Odt2 += ((d2fvPrimedt2[H2O]*beta + betaPrime[H2O]*d2fvdt2)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    d2lnPhiH2Odt2 += ((d2fvPrimedt2[H2O]*gammav+gammavPrime[H2O]*d2fvdt2-d2fvdt2*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    d2lnPhiH2Odt2 += ((gammavPrime[H2O]-gammav)*d2fvdt2/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    *vPt = v;
    *zPt = z;
    *phi = exp(lnPhiH2O);
    *dvdp = 1.0/( p*(dzdv/z - 1.0/v) );
    *dvdt = (1.0/t + dzdt/z)/(1.0/v - dzdv/z);
    *d2vdp2  = p*(1.0/v-dzdv/z)/v - dzdv/(*dvdp)/z + 1.0/(*dvdp)/(*dvdp)/p + p*d2zdv2/z;
    *d2vdp2 *= -(*dvdp)*(*dvdp)*(*dvdp);
    *d2vdtdp = -(*dvdp)*(1.0/t + dzdt/z) - p*(*d2vdp2)*(1.0/t + dzdt/z) - p*(*dvdp)*(*dvdp)*(-dzdv*dzdt/z/z + d2zdvdt/z);
    *d2vdt2  = -p*(*d2vdtdp)*(1.0/t + dzdt/z) + p*(*dvdp)/t/t + p*(*dvdp)*dzdt*(dzdv*(*dvdt) + dzdt)/z/z - p*(*dvdp)*(*dvdt)*d2zdvdt/z - p*(*dvdp)*d2zdt2/z;

    *dlnphidt    = dlnPhiH2Odv*(*dvdt) + dlnPhiH2Odt;
    *dlnphidp    = dlnPhiH2Odv*(*dvdp);
    *d2lnphidt2  = d2lnPhiH2Odv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiH2Odvdt*(*dvdt) + dlnPhiH2Odv*(*d2vdt2) + d2lnPhiH2Odt2;
    *d2lnphidtdp = d2lnPhiH2Odv2*(*dvdt)*(*dvdp) + dlnPhiH2Odv*(*d2vdtdp) + d2lnPhiH2Odvdt*(*dvdp);
    *d2lnphidp2  = d2lnPhiH2Odv2*(*dvdp)*(*dvdp) + dlnPhiH2Odv*(*d2vdp2);
}

#pragma mark -
#pragma mark Duan entry function (pure H2O)

static void duanH2O(int useLowPcoeff, double t, double p, double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                    double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {

    duanH2ODriver(useLowPcoeff, t, p, vPt, zPt, phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2);

    if (!useLowPcoeff) {
        double vRef, zRef, phiRef, dvdtRef, dvdpRef, d2vdt2Ref, d2vdtdpRef, d2vdp2Ref, dlnphidtRef, dlnphidpRef, d2lnphidt2Ref, d2lnphidtdpRef, d2lnphidp2Ref;

        duanH2ODriver(1, t, 2000.0, &vRef, &zRef, &phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, &dlnphidtRef, &dlnphidpRef, &d2lnphidt2Ref, &d2lnphidtdpRef, &d2lnphidp2Ref);
        *phi         *= phiRef;
        *dlnphidt    += dlnphidtRef;
        *d2lnphidt2  += d2lnphidt2Ref;

        duanH2ODriver(0, t, 2000.0, &vRef, &zRef, &phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, &dlnphidtRef, &dlnphidpRef, &d2lnphidt2Ref, &d2lnphidtdpRef, &d2lnphidp2Ref);
        *phi         /= phiRef;
        *dlnphidt    -= dlnphidtRef;
        *d2lnphidt2  -= d2lnphidt2Ref;
    }

}

#pragma mark -
#pragma mark Duan driver routine (pure CO2)

static void duanCO2Driver(int useLowPcoeff, double t, double p, double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                          double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {
    double bv, cv, dv, ev, fv, beta, gammav, v, z = 1.0, dzdv, dzdt, d2zdv2, d2zdvdt, d2zdt2;
    double dbvdt, dcvdt, ddvdt, devdt, dfvdt, d2bvdt2, d2cvdt2, d2dvdt2, d2evdt2, d2fvdt2, d3bvdt3, d3cvdt3, d3dvdt3, d3evdt3, d3fvdt3;
    double bvPrime[2], cvPrime[2], dvPrime[2], evPrime[2], fvPrime[2], betaPrime[2], gammavPrime[2];
    double b2vPrime2[2][2], c2vPrime2[2][2], d2vPrime2[2][2], e2vPrime2[2][2], f2vPrime2[2][2], gamma2vPrime2[2][2];
    double dbvPrimedt[2], d2bvPrimedt2[2], d3bvPrimedt3[2], dcvPrimedt[2], d2cvPrimedt2[2], d3cvPrimedt3[2],
    ddvPrimedt[2], d2dvPrimedt2[2], d3dvPrimedt3[2], devPrimedt[2], d2evPrimedt2[2], d3evPrimedt3[2],
    dfvPrimedt[2], d2fvPrimedt2[2], d3fvPrimedt3[2];
    double lnPhiCO2, dlnPhiCO2dv, dlnPhiCO2dt, d2lnPhiCO2dv2, d2lnPhiCO2dt2, d2lnPhiCO2dvdt;
    double x[2] = { 0.0, 1.0 }; /* CO2 */

    BVcAndDerivative    (useLowPcoeff, t, x, &bv,     bvPrime,    &dbvdt,    &d2bvdt2,  &d3bvdt3, dbvPrimedt, d2bvPrimedt2, d3bvPrimedt3, b2vPrime2);
    CVcAndDerivative    (useLowPcoeff, t, x, &cv,     cvPrime,    &dcvdt,    &d2cvdt2,  &d3cvdt3, dcvPrimedt, d2cvPrimedt2, d3cvPrimedt3, c2vPrime2);
    DVcAndDerivative    (useLowPcoeff, t, x, &dv,     dvPrime,	&ddvdt,    &d2dvdt2,  &d3dvdt3, ddvPrimedt, d2dvPrimedt2, d3dvPrimedt3, d2vPrime2);
    EVcAndDerivative    (useLowPcoeff, t, x, &ev,     evPrime,	&devdt,    &d2evdt2,  &d3evdt3, devPrimedt, d2evPrimedt2, d3evPrimedt3, e2vPrime2);
    FVcAndDerivative    (useLowPcoeff, t, x, &fv,     fvPrime,	&dfvdt,    &d2fvdt2,  &d3fvdt3, dfvPrimedt, d2fvPrimedt2, d3fvPrimedt3, f2vPrime2);
    BetaAndDerivative   (useLowPcoeff, t, x, &beta,   betaPrime);
    GammaVcAndDerivative(useLowPcoeff, t, x, &gammav, gammavPrime, gamma2vPrime2);

    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0;

        v = 8.314467*t/p;
        while (iter < 200) {
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            delv = z*8.314467*t/p - v;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
            vPrevious = v;
            delvPrevious = delv;
            v = (z*8.314467*t/p + v)/2.0;
            iter++;
        }
        if (iter == 200) printf("duanCO2Driver(a): t = %g, p = %g, z = %g, v = %g, vPrev = %g, delv = %g, delvPrev = %g, iter = %d\n", t, p, z, v, vPrevious, delv, delvPrevious, iter);
        else if (fabs(delv) > v*100.0*DBL_EPSILON) {
            double dx;
            double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
            iter = 0;
            while (iter < 200) {
                v = rtb + (dx *= 0.5);
                z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
                delv = z*8.314467*t/p - v;
                if (delv <= 0.0) rtb = v;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
            if ( (iter == 200) || (fabs(dx) > 100.0*DBL_EPSILON) ) printf("duanCO2Driver(b): t = %g, p = %g, z = %g, v = %g, delv = %g, dx = %g, iter = %d\n", t, p, z, v, delv, dx, iter);
        }
    }

    lnPhiCO2  = 0.0;
    lnPhiCO2 += -log(z);
    lnPhiCO2 += bvPrime[CO2]/v;
    lnPhiCO2 += cvPrime[CO2]/2.0/v/v;
    lnPhiCO2 += dvPrime[CO2]/4.0/v/v/v/v;
    lnPhiCO2 += evPrime[CO2]/5.0/v/v/v/v/v;
    lnPhiCO2 += ((fvPrime[CO2]*beta + betaPrime[CO2]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiCO2 += ((fvPrime[CO2]*gammav+gammavPrime[CO2]*fv-fv*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav)
    *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiCO2 += ((gammavPrime[CO2]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    dzdv = - bv/v/v - 2.0*cv/v/v/v + - 4.0*dv/v/v/v/v/v - 5.0*ev/v/v/v/v/v/v
    - 2.0*(fv/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(fv/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(fv/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    dzdt = dbvdt/v + dcvdt/v/v + ddvdt/v/v/v/v + devdt/v/v/v/v/v + (dfvdt/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    dlnPhiCO2dv  = 0.0;
    dlnPhiCO2dv += -dzdv/z;
    dlnPhiCO2dv += -bvPrime[CO2]/v/v;
    dlnPhiCO2dv += -cvPrime[CO2]/v/v/v;
    dlnPhiCO2dv += -dvPrime[CO2]/v/v/v/v/v;
    dlnPhiCO2dv += -evPrime[CO2]/v/v/v/v/v/v;
    dlnPhiCO2dv += -((fvPrime[CO2]*beta + betaPrime[CO2]*fv)/v/v/v)*exp(-gammav/v/v);
    dlnPhiCO2dv += -((fvPrime[CO2]*gammav+gammavPrime[CO2]*fv-fv*beta*(gammavPrime[CO2]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    dlnPhiCO2dv += +((gammavPrime[CO2]-gammav)*fv*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    dlnPhiCO2dt  = 0.0;
    dlnPhiCO2dt += -dzdt/z;
    dlnPhiCO2dt += dbvPrimedt[CO2]/v;
    dlnPhiCO2dt += dcvPrimedt[CO2]/2.0/v/v;
    dlnPhiCO2dt += ddvPrimedt[CO2]/4.0/v/v/v/v;
    dlnPhiCO2dt += devPrimedt[CO2]/5.0/v/v/v/v/v;
    dlnPhiCO2dt += ((dfvPrimedt[CO2]*beta + betaPrime[CO2]*dfvdt)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    dlnPhiCO2dt += ((dfvPrimedt[CO2]*gammav+gammavPrime[CO2]*dfvdt-dfvdt*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    dlnPhiCO2dt += ((gammavPrime[CO2]-gammav)*dfvdt/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    d2zdv2 =  2.0*bv/v/v/v + 6.0*cv/v/v/v/v + 20.0*dv/v/v/v/v/v/v + 30.0*ev/v/v/v/v/v/v/v
    +  6.0*(fv/v/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    + 14.0*(fv/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    - 14.0*(fv/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    -  8.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    +  4.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
    d2zdvdt = - dbvdt/v/v - 2.0*dcvdt/v/v/v + - 4.0*ddvdt/v/v/v/v/v - 5.0*devdt/v/v/v/v/v/v
    - 2.0*(dfvdt/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(dfvdt/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(dfvdt/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    d2zdt2 = d2bvdt2/v + d2cvdt2/v/v + d2dvdt2/v/v/v/v + d2evdt2/v/v/v/v/v + (d2fvdt2/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    d2lnPhiCO2dv2  = 0.0;
    d2lnPhiCO2dv2 += dzdv*dzdv/z/z - d2zdv2/z;
    d2lnPhiCO2dv2 +=  2.0*bvPrime[CO2]/v/v/v;
    d2lnPhiCO2dv2 +=  3.0*cvPrime[CO2]/v/v/v/v;
    d2lnPhiCO2dv2 +=  5.0*dvPrime[CO2]/v/v/v/v/v/v;
    d2lnPhiCO2dv2 +=  6.0*evPrime[CO2]/v/v/v/v/v/v/v;
    d2lnPhiCO2dv2 +=  3.0*((fvPrime[CO2]*beta + betaPrime[CO2]*fv)/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 += -2.0*((fvPrime[CO2]*beta*gammav + betaPrime[CO2]*fv*gammav)/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 +=  5.0*((fvPrime[CO2]*gammav+gammavPrime[CO2]*fv-fv*beta*(gammavPrime[CO2]-gammav))/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 += -2.0*((fvPrime[CO2]*gammav*gammav+gammavPrime[CO2]*fv*gammav-fv*beta*gammav*(gammavPrime[CO2]-gammav))/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 += -7.0*((gammavPrime[CO2]-gammav)*fv*gammav/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dv2 +=  2.0*((gammavPrime[CO2]-gammav)*fv*gammav*gammav/v/v/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiCO2dvdt  = 0.0;
    d2lnPhiCO2dvdt += dzdv*dzdt/z/z -d2zdvdt/z;
    d2lnPhiCO2dvdt += -dbvPrimedt[CO2]/v/v;
    d2lnPhiCO2dvdt += -dcvPrimedt[CO2]/v/v/v;
    d2lnPhiCO2dvdt += -ddvPrimedt[CO2]/v/v/v/v/v;
    d2lnPhiCO2dvdt += -devPrimedt[CO2]/v/v/v/v/v/v;
    d2lnPhiCO2dvdt += -((dfvPrimedt[CO2]*beta + betaPrime[CO2]*dfvdt)/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dvdt += -((dfvPrimedt[CO2]*gammav+gammavPrime[CO2]*dfvdt-dfvdt*beta*(gammavPrime[CO2]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiCO2dvdt +=  ((gammavPrime[CO2]-gammav)*dfvdt*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiCO2dt2  = 0.0;
    d2lnPhiCO2dt2 += dzdt*dzdt/z/z - d2zdt2/z;
    d2lnPhiCO2dt2 += d2bvPrimedt2[CO2]/v;
    d2lnPhiCO2dt2 += d2cvPrimedt2[CO2]/2.0/v/v;
    d2lnPhiCO2dt2 += d2dvPrimedt2[CO2]/4.0/v/v/v/v;
    d2lnPhiCO2dt2 += d2evPrimedt2[CO2]/5.0/v/v/v/v/v;
    d2lnPhiCO2dt2 += ((d2fvPrimedt2[CO2]*beta + betaPrime[CO2]*d2fvdt2)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    d2lnPhiCO2dt2 += ((d2fvPrimedt2[CO2]*gammav+gammavPrime[CO2]*d2fvdt2-d2fvdt2*beta*(gammavPrime[CO2]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    d2lnPhiCO2dt2 += ((gammavPrime[CO2]-gammav)*d2fvdt2/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    *vPt = v;
    *zPt = z;
    *phi = exp(lnPhiCO2);
    *dvdp = 1.0/( p*(dzdv/z - 1.0/v) );
    *dvdt = (1.0/t + dzdt/z)/(1.0/v - dzdv/z);
    *d2vdp2  = p*(1.0/v-dzdv/z)/v - dzdv/(*dvdp)/z + 1.0/(*dvdp)/(*dvdp)/p + p*d2zdv2/z;
    *d2vdp2 *= -(*dvdp)*(*dvdp)*(*dvdp);
    *d2vdtdp = -(*dvdp)*(1.0/t + dzdt/z) - p*(*d2vdp2)*(1.0/t + dzdt/z) - p*(*dvdp)*(*dvdp)*(-dzdv*dzdt/z/z + d2zdvdt/z);
    *d2vdt2  = -p*(*d2vdtdp)*(1.0/t + dzdt/z) + p*(*dvdp)/t/t + p*(*dvdp)*dzdt*(dzdv*(*dvdt) + dzdt)/z/z - p*(*dvdp)*(*dvdt)*d2zdvdt/z - p*(*dvdp)*d2zdt2/z;

    *dlnphidt    = dlnPhiCO2dv*(*dvdt) + dlnPhiCO2dt;
    *dlnphidp    = dlnPhiCO2dv*(*dvdp);
    *d2lnphidt2  = d2lnPhiCO2dv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiCO2dvdt*(*dvdt) + dlnPhiCO2dv*(*d2vdt2) + d2lnPhiCO2dt2;
    *d2lnphidtdp = d2lnPhiCO2dv2*(*dvdt)*(*dvdp) + dlnPhiCO2dv*(*d2vdtdp) + d2lnPhiCO2dvdt*(*dvdp);
    *d2lnphidp2  = d2lnPhiCO2dv2*(*dvdp)*(*dvdp) + dlnPhiCO2dv*(*d2vdp2);
}

#pragma mark -
#pragma mark Duan entry function (pure CO2)

static void duanCO2(int useLowPcoeff, double t, double p, double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                    double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {

    duanCO2Driver(useLowPcoeff, t, p, vPt, zPt, phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2);

    if (!useLowPcoeff) {
        double vRef, zRef, phiRef, dvdtRef, dvdpRef, d2vdt2Ref, d2vdtdpRef, d2vdp2Ref, dlnphidtRef, dlnphidpRef, d2lnphidt2Ref, d2lnphidtdpRef, d2lnphidp2Ref;

        duanCO2Driver(1, t, 2000.0, &vRef, &zRef, &phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, &dlnphidtRef, &dlnphidpRef, &d2lnphidt2Ref, &d2lnphidtdpRef, &d2lnphidp2Ref);
        *phi         *= phiRef;
        *dlnphidt    += dlnphidtRef;
        *d2lnphidt2  += d2lnphidt2Ref;

        duanCO2Driver(0, t, 2000.0, &vRef, &zRef, &phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, &dlnphidtRef, &dlnphidpRef, &d2lnphidt2Ref, &d2lnphidtdpRef, &d2lnphidp2Ref);
        *phi         /= phiRef;
        *dlnphidt    -= dlnphidtRef;
        *d2lnphidt2  -= d2lnphidt2Ref;
    }

}

#pragma mark -
#pragma mark Duan pure H2O

static void propertiesOfPureH2O(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt,
                                double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2) {
    double z, phi, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2;

    idealGasH2O(t, cp, s, h, dcpdt);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, v, &z, &phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, &dlnphidt, &dlnphidp, &d2lnphidt2, &d2lnphidtdp, &d2lnphidp2);

    *g      = *h - t*(*s) + R*t*log(phi*p);
    *s     += - (R*log(phi*p) + R*t*dlnphidt);
    *h     += R*t*log(phi*p) - t*(R*log(phi*p) + R*t*dlnphidt);
    *cp    += -t*(2.0*R*dlnphidt + R*t*d2lnphidt2);
    {
        double zTemp, phiTemp, dlnphidtTemp, dlnphidpTemp, d2lnphidt2Temp, d2lnphidtdpTemp, d2lnphidp2Temp, vTemp, dvdtTemp, dvdpTemp, d2vdt2Temp, d2vdtdpTemp, d2vdp2Temp, d3lnphidt3;

        duanH2O((p <= 2000.0) ? 1 : 0, t*(1.0+sqrt(DBL_EPSILON)), p, &vTemp, &zTemp, &phiTemp, &dvdtTemp, &dvdpTemp, &d2vdt2Temp, &d2vdtdpTemp, &d2vdp2Temp,
                &dlnphidtTemp, &dlnphidpTemp, &d2lnphidt2Temp, &d2lnphidtdpTemp, &d2lnphidp2Temp);

        d3lnphidt3 = (d2lnphidt2Temp - d2lnphidt2)/t/sqrt(DBL_EPSILON);
        *dcpdt += -(2.0*R*dlnphidt + R*t*d2lnphidt2) -t*(3.0*R*d2lnphidt2 + R*t*d3lnphidt3);
    }
}

#pragma mark -
#pragma mark Duan pure CO2

static void propertiesOfPureCO2(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt,
                                double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2) {
    double z, phi, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2;

    idealGasCO2(t, cp, s, h, dcpdt);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, v, &z, &phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, &dlnphidt, &dlnphidp, &d2lnphidt2, &d2lnphidtdp, &d2lnphidp2);

    *g = *h - t*(*s) + R*t*log(phi*p);
    *s     += - (R*log(phi*p) + R*t*dlnphidt);
    *h     += R*t*log(phi*p) - t*(R*log(phi*p) + R*t*dlnphidt);
    *cp    += -t*(2.0*R*dlnphidt + R*t*d2lnphidt2);
    {
        double zTemp, phiTemp, dlnphidtTemp, dlnphidpTemp, d2lnphidt2Temp, d2lnphidtdpTemp, d2lnphidp2Temp, vTemp, dvdtTemp, dvdpTemp, d2vdt2Temp, d2vdtdpTemp, d2vdp2Temp, d3lnphidt3;

        duanCO2((p <= 2000.0) ? 1 : 0, t*(1.0+sqrt(DBL_EPSILON)), p, &vTemp, &zTemp, &phiTemp, &dvdtTemp, &dvdpTemp, &d2vdt2Temp, &d2vdtdpTemp, &d2vdp2Temp,
                &dlnphidtTemp, &dlnphidpTemp, &d2lnphidt2Temp, &d2lnphidtdpTemp, &d2lnphidp2Temp);

        d3lnphidt3 = (d2lnphidt2Temp - d2lnphidt2)/t/sqrt(DBL_EPSILON);
        *dcpdt += -(2.0*R*dlnphidt + R*t*d2lnphidt2) -t*(3.0*R*d2lnphidt2 + R*t*d3lnphidt3);
    }
}

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
#pragma mark Duan H2O class

@implementation DuanH2O

-(id)init {
	if ((self = [super init])) {
		[self setPhaseName:@"Water"];
		[self setPhaseFormula:@"H2O"];
	}
	return self;
}

-(double)calculateWithT:(double)t andWithP:(double)p forCase:(NSUInteger)returnMode {
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
	double result = 0.0;

    propertiesOfPureH2O(t, p, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

	switch (returnMode) {
		case 1:
			result = gH2O;
			break;
		case 2:
			result = hH2O;
			break;
		case 3:
			result = sH2O;
			break;
		case 4:
			result = cpH2O;
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

#pragma mark -
#pragma mark Duan CO2 class

@implementation DuanCO2

-(id)init {
	if ((self = [super init])) {
		[self setPhaseName:@"Carbon Dioxide"];
		[self setPhaseFormula:@"CO2"];
	}
	return self;
}

-(double)calculateWithT:(double)t andWithP:(double)p forCase:(NSUInteger)returnMode {
    double gCO2, hCO2, sCO2, cpCO2, dcpdtCO2, vCO2, dvdtCO2, dvdpCO2, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2;
	double result = 0.0;

    propertiesOfPureCO2(t, p, &gCO2, &hCO2, &sCO2, &cpCO2, &dcpdtCO2, &vCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2);

	switch (returnMode) {
		case 1:
			result = gCO2;
			break;
		case 2:
			result = hCO2;
			break;
		case 3:
			result = sCO2;
			break;
		case 4:
			result = cpCO2;
			break;
		case 5:
			result = dcpdtCO2;
			break;
		case 6:
			result = vCO2;
			break;
		case 7:
			result = dvdtCO2;
			break;
		case 8:
			result = dvdpCO2;
			break;
		case 9:
			result = d2vdt2CO2;
			break;
		case 10:
			result = d2vdtdpCO2;
			break;
		case 11:
			result = d2vdp2CO2;
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

#pragma mark -
#pragma mark Duan Fluid class

@implementation FluidDuan

static NSArray *endmembers;

#define NR      1    // Number of independent mole fraction variables
#define NA      2    // Number of liquid components
#define NATOMS  3.0  // Average number of atoms in the formula unit (guess)

+(void)initialize {
	if (self == [FluidDuan class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(FluidDuan) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

        PhaseBase *h2o = [[DuanH2O alloc] init];
		[mutableEndmembers addObject:h2o];
		if (debug) NSLog(@"... allocated h2o ...");

        PhaseBase *co2 = [[DuanCO2 alloc] init];
		[mutableEndmembers addObject:co2];
		if (debug) NSLog(@"... allocated co2 ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(FluidDuan) ... entry ...");
		computeMixingQuantities = NO;
		[self setPhaseName:@"Fluid"];
		if (debug) NSLog(@"... exiting.");
	}
	return self;
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
	const char *NAMES[NA]    = { "Water", "Carbon dioxide" };
	const char *FORMULAS[NA] = { "H2O",   "CO2" };
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

		m[0] = e[H]/2.0; /* moles of H2O */
        m[1] = e[C];     /* Moles of CO2 */

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

			r[0] = (sum != 0.0) ? m[1]/sum : 0.0;
		}

		if (outMask & FOURTH) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of mole fractions of endmember components                            */

			for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
		}

		if (outMask & FIFTH) {
			/* Calculates the matrix dr[i]/dm[j] using m[] as input                 */

			if (sum == 0.0) {
				for (j=0; j<NA; j++) dm[0][j] = 0.0;
			} else {
				for (j=0; j<NA; j++) dm[0][j] = (1 == j) ? (1.0-m[1]/sum)/sum : - m[1]/SQUARE(sum);
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
            x[0] = 1.0 - r[0];
            x[1] = r[0];
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
		NSLog(@"Illegal call to convert with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double h2o = 1.0 - r[0];
	double co2 =       r[0];
	double sum = h2o + co2;
	if (sum == 0.0) sum = 1.0;

	return [NSString stringWithFormat:@"X:H2O %5.3f CO2 %5.3f", h2o/sum, co2/sum];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dr // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
    double x[2], v, z, phi[2], dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt[2], dlnphidp[2], d2lnphidt2[2], d2lnphidtdp[2], d2lnphidp2[2], dlnphidr[2];
    double vH2O, vCO2, zH2O, zCO2, phiH2O, phiCO2, dvdtH2O, dvdpH2O, dvdtCO2, dvdpCO2, d2vdt2H2O,
    d2vdtdpH2O, d2vdp2H2O, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2, dlnphidtH2O, dlnphidtCO2, dlnphidpH2O, dlnphidpCO2,
    d2lnphidt2H2O, d2lnphidtdpH2O, d2lnphidp2H2O, d2lnphidt2CO2, d2lnphidtdpCO2, d2lnphidp2CO2;
    int i, j;

    x[H2O] = 1.0 - r[0];
    x[CO2] = r[0];

    if        (fabs(x[CO2]) < 100.0*DBL_EPSILON) {
        if (mask & FIRST)  { a[H2O]   = 1.0;          a[CO2]  = 0.0; }
        if (mask & SECOND) { mu[H2O]  = DBL_EPSILON; mu[CO2]  = 0.0; }
        if (mask & THIRD)  { dr[0][0] = 0.0;         dr[1][0] = 0.0; }
        if (mask & FOURTH) { ; }
        return;
    } else if (fabs(x[H2O]) < 100.0*DBL_EPSILON) {
        if (mask & FIRST)  { a[H2O]   = 0.0; a[CO2]   = 1.0;         }
        if (mask & SECOND) { mu[H2O]  = 0.0; mu[CO2]  = DBL_EPSILON; }
        if (mask & THIRD)  { dr[0][0] = 0.0; dr[1][0] = 0.0;         }
        if (mask & FOURTH) { ; }
        return;
    }

    duan((p <= 2000.0) ? 1 : 0, t, p, x, &v, &z, phi, &dvdt, &dvdp, &d2vdt2, &d2vdtdp, &d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, &vH2O, &zH2O, &phiH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, &dlnphidtH2O, &dlnphidpH2O, &d2lnphidt2H2O, &d2lnphidtdpH2O, &d2lnphidp2H2O);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, &vCO2, &zCO2, &phiCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2, &dlnphidtCO2, &dlnphidpCO2, &d2lnphidt2CO2, &d2lnphidtdpCO2, &d2lnphidp2CO2);

    if (mask & FIRST) {
        a[H2O] = (phiH2O != 0.0) ? x[H2O]*phi[H2O]/phiH2O : 0.0;
        a[CO2] = (phiCO2 != 0.0) ? x[CO2]*phi[CO2]/phiCO2 : 0.0;
    }

    if (mask & SECOND) {
        mu[H2O] = ((phi[H2O] != 0.0) || (phiH2O != 0.0)) ? R*t*log(x[H2O]*phi[H2O]/phiH2O) : 0.0;
        mu[CO2] = ((phi[CO2] != 0.0) || (phiCO2 != 0.0)) ? R*t*log(x[CO2]*phi[CO2]/phiCO2) : 0.0;
    }

    if (mask & THIRD) {
        double g, dgdr[NR], fr[NA][NR];
        double d2gdr2[NR][NR], dfrdr[NA][NR], sum;
        int k;

        fr[0][0] =     - r[0];
        fr[1][0] = 1.0 - r[0];

        g            = ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) ? R*t*(x[H2O]*log(x[H2O]*phi[H2O]/phiH2O) + x[CO2]*log(x[CO2]*phi[CO2]/phiCO2)) : 0.0;
        dgdr[0]      = ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) ? R*t*(log(x[CO2]*phi[CO2]/phiCO2) - log(x[H2O]*phi[H2O]/phiH2O))               : 0.0;
        d2gdr2[0][0] = ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) ? R*t*(1.0/x[CO2] + 1.0/x[H2O] - dlnphidr[H2O] + dlnphidr[CO2])                 : 0.0;

        for(i=0; i<NA; i++) dfrdr[i][0] = - 1.0;

        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                for (dr[i][k]=g, j=0; j<NR; j++) dr[i][k] += fr[i][j]*dgdr[j];
                dr[i][k] = exp(dr[i][k]/(R*t));
                sum = (1.0+dfrdr[i][k])*dgdr[k];
                for (j=0; j<NR; j++) sum += fr[i][j]*d2gdr2[j][k];
                dr[i][k] *= sum/(R*t);
            }
        }
    }
}

-(void)gmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   gmix:(double *)gmix              // Gibbs energy of mixing             BINARY MASK: 0001
		 dx:(double [NR])dr          // (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010
		dx2:(double [NR][NR])dr2     // (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100
		dx3:(double [NR][NR][NR])dr3 // (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000
{
    double x[2], v, z, phi[2], dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt[2], dlnphidp[2], d2lnphidt2[2], d2lnphidtdp[2], d2lnphidp2[2], dlnphidr[2];
    double vH2O, vCO2, zH2O, zCO2, phiH2O, phiCO2, dvdtH2O, dvdpH2O, dvdtCO2, dvdpCO2, d2vdt2H2O,
    d2vdtdpH2O, d2vdp2H2O, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2, dlnphidtH2O, dlnphidtCO2, dlnphidpH2O, dlnphidpCO2,
    d2lnphidt2H2O, d2lnphidtdpH2O, d2lnphidp2H2O, d2lnphidt2CO2, d2lnphidtdpCO2, d2lnphidp2CO2;

    x[H2O] = 1.0 - r[0];
    x[CO2] = r[0];

    if ((fabs(x[CO2]) < 100.0*DBL_EPSILON) || (fabs(x[H2O]) < 100.0*DBL_EPSILON)) {
        if (mask & FIRST)  { *gmix           = 0.0; }
        if (mask & SECOND) { dr[0]        = 0.0; }
        if (mask & THIRD)  { dr2[0][0]    = 0.0; }
        if (mask & FOURTH) { dr3[0][0][0] = 0.0; }
        return;
    }

    duan((p <= 2000.0) ? 1 : 0, t, p, x, &v, &z, phi, &dvdt, &dvdp, &d2vdt2, &d2vdtdp, &d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, &vH2O, &zH2O, &phiH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, &dlnphidtH2O, &dlnphidpH2O, &d2lnphidt2H2O, &d2lnphidtdpH2O, &d2lnphidp2H2O);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, &vCO2, &zCO2, &phiCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2, &dlnphidtCO2, &dlnphidpCO2, &d2lnphidt2CO2, &d2lnphidtdpCO2, &d2lnphidp2CO2);

    if (mask & FIRST) {
        *gmix = ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) ? R*t*(x[H2O]*log(x[H2O]*phi[H2O]/phiH2O) + x[CO2]*log(x[CO2]*phi[CO2]/phiCO2)) : 0.0;
    }

    if(mask & SECOND) {
        dr[0] = ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) ? R*t*(log(x[CO2]*phi[CO2]/phiCO2) - log(x[H2O]*phi[H2O]/phiH2O)) : 0.0;
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR];
        int i, j;

        d2gdr2[0][0] = R*t*(1.0/x[CO2] + 1.0/x[H2O] - dlnphidr[H2O] + dlnphidr[CO2]);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) dr2[i][j] = d2gdr2[i][j];
        }
    }

    if (mask & FOURTH) {
#ifdef DO_NOT_THROW_EXCEPTIONS
        int i,j,k;
        for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) dr3[i][j][k] = 0.0;
#else
        NSException *duanException = [NSException exceptionWithName:@"Duan fluid Internal Error"
                                                             reason:@"d2g(mix)dr2 not yet implemented."
                                                           userInfo:nil];
		@throw duanException;
#endif
    }
}

-(void)hmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   hmix:(double *)hmix // Enthalpy of mixing BINARY MASK: 1
{
    double x[2], v, z, phi[2], dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt[2], dlnphidp[2], d2lnphidt2[2], d2lnphidtdp[2], d2lnphidp2[2], dlnphidr[2];
    double vH2O, vCO2, zH2O, zCO2, phiH2O, phiCO2, dvdtH2O, dvdpH2O, dvdtCO2, dvdpCO2, d2vdt2H2O,
    d2vdtdpH2O, d2vdp2H2O, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2, dlnphidtH2O, dlnphidtCO2, dlnphidpH2O, dlnphidpCO2,
    d2lnphidt2H2O, d2lnphidtdpH2O, d2lnphidp2H2O, d2lnphidt2CO2, d2lnphidtdpCO2, d2lnphidp2CO2;

    x[H2O] = 1.0 - r[0];
    x[CO2] = r[0];

    if ((fabs(x[CO2]) < 100.0*DBL_EPSILON) || (fabs(x[H2O]) < 100.0*DBL_EPSILON)) {
        if (mask & FIRST)  { *hmix        = 0.0; }
        return;
    }

    duan((p <= 2000.0) ? 1 : 0, t, p, x, &v, &z, phi, &dvdt, &dvdp, &d2vdt2, &d2vdtdp, &d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, &vH2O, &zH2O, &phiH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, &dlnphidtH2O, &dlnphidpH2O, &d2lnphidt2H2O, &d2lnphidtdpH2O, &d2lnphidp2H2O);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, &vCO2, &zCO2, &phiCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2, &dlnphidtCO2, &dlnphidpCO2, &d2lnphidt2CO2, &d2lnphidtdpCO2, &d2lnphidp2CO2);

    if ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) {
        *hmix  = R*(x[H2O]*log(x[H2O]) + x[H2O]*log(phi[H2O]) - x[H2O]*log(phiH2O) + x[CO2]*log(x[CO2]) + x[CO2]*log(phi[CO2]) - x[CO2]*log(phiCO2))
        + R*t*(x[H2O]*dlnphidt[H2O] - x[H2O]*dlnphidtH2O + x[CO2]*dlnphidt[CO2] - x[CO2]*dlnphidtCO2);
        *hmix *= -t;
        *hmix +=  R*t*(x[H2O]*log(x[H2O]*phi[H2O]/phiH2O) + x[CO2]*log(x[CO2]*phi[CO2]/phiCO2));
    } else *hmix = 0.0;
}

-(void)smix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   smix:(double *)smix        // Entropy of mixing                  BINARY MASK: 001
		 dx:(double [NR])dr       // (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010
		dx2:(double [NR][NR])dr2  // (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100
{
    double x[2], v, z, phi[2], dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt[2], dlnphidp[2], d2lnphidt2[2], d2lnphidtdp[2], d2lnphidp2[2], dlnphidr[2];
    double vH2O, vCO2, zH2O, zCO2, phiH2O, phiCO2, dvdtH2O, dvdpH2O, dvdtCO2, dvdpCO2, d2vdt2H2O,
    d2vdtdpH2O, d2vdp2H2O, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2, dlnphidtH2O, dlnphidtCO2, dlnphidpH2O, dlnphidpCO2,
    d2lnphidt2H2O, d2lnphidtdpH2O, d2lnphidp2H2O, d2lnphidt2CO2, d2lnphidtdpCO2, d2lnphidp2CO2;

    x[H2O] = 1.0 - r[0];
    x[CO2] = r[0];

    if ((fabs(x[CO2]) < 100.0*DBL_EPSILON) || (fabs(x[H2O]) < 100.0*DBL_EPSILON)) {
        if (mask & FIRST)  { *smix        = 0.0; }
        if (mask & SECOND) { dr[0]        = 0.0; }
        if (mask & THIRD)  { dr2[0][0]    = 0.0; }
        return;
    }

    duan((p <= 2000.0) ? 1 : 0, t, p, x, &v, &z, phi, &dvdt, &dvdp, &d2vdt2, &d2vdtdp, &d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, &vH2O, &zH2O, &phiH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, &dlnphidtH2O, &dlnphidpH2O, &d2lnphidt2H2O, &d2lnphidtdpH2O, &d2lnphidp2H2O);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, &vCO2, &zCO2, &phiCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2, &dlnphidtCO2, &dlnphidpCO2, &d2lnphidt2CO2, &d2lnphidtdpCO2, &d2lnphidp2CO2);

    if (mask & FIRST) {
        *smix = 0.0;
        if ((phi[H2O] != 0.0) && (phi[CO2] != 0.0)) {
            *smix = R*(x[H2O]*log(x[H2O]) + x[H2O]*log(phi[H2O]) - x[H2O]*log(phiH2O) + x[CO2]*log(x[CO2]) + x[CO2]*log(phi[CO2]) - x[CO2]*log(phiCO2))
            + R*t*(x[H2O]*dlnphidt[H2O] - x[H2O]*dlnphidtH2O + x[CO2]*dlnphidt[CO2] - x[CO2]*dlnphidtCO2);
            *smix *= -1.0;
        }
    }

    if(mask & SECOND) {
        double d2gdrdt[NR];
        int i;

        d2gdrdt[0] = R*(- log(x[H2O]) - log(phi[H2O]) + log(phiH2O) + log(x[CO2]) + log(phi[CO2]) - log(phiCO2))
        + R*t*(- dlnphidt[H2O] + dlnphidtH2O + dlnphidt[CO2] - dlnphidtCO2);

        for (i=0; i<NR; i++) dr[i] = - d2gdrdt[i];
    }

    if(mask & THIRD) {
#ifdef DO_NOT_THROW_EXCEPTIONS
        int i, j;
        for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr2[i][j] = 0.0;
#else
        NSException *duanException = [NSException exceptionWithName:@"Duan fluid Internal Error"
                                                             reason:@"d2s(mix)dr2 not yet implemented."
                                                           userInfo:nil];
		@throw duanException;
#endif
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
    double x[2], v, z, phi[2], dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt[2], dlnphidp[2], d2lnphidt2[2], d2lnphidtdp[2], d2lnphidp2[2], dlnphidr[2];
    double vH2O, vCO2, zH2O, zCO2, phiH2O, phiCO2, dvdtH2O, dvdpH2O, dvdtCO2, dvdpCO2, d2vdt2H2O,
    d2vdtdpH2O, d2vdp2H2O, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2, dlnphidtH2O, dlnphidtCO2, dlnphidpH2O, dlnphidpCO2,
    d2lnphidt2H2O, d2lnphidtdpH2O, d2lnphidp2H2O, d2lnphidt2CO2, d2lnphidtdpCO2, d2lnphidp2CO2;
    double d2gdt2 = 0.0;

    x[H2O] = 1.0 - r[0];
    x[CO2] = r[0];

    if ((fabs(x[CO2]) < 100.0*DBL_EPSILON) || (fabs(x[H2O]) < 100.0*DBL_EPSILON)) {
        if (mask & FIRST)  { *cpmix = 0.0; }
        if (mask & SECOND) { *dt    = 0.0; }
        if (mask & THIRD)  { dr[0]  = 0.0; }
        return;
    }

    duan((p <= 2000.0) ? 1 : 0, t, p, x, &v, &z, phi, &dvdt, &dvdp, &d2vdt2, &d2vdtdp, &d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, &vH2O, &zH2O, &phiH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, &dlnphidtH2O, &dlnphidpH2O, &d2lnphidt2H2O, &d2lnphidtdpH2O, &d2lnphidp2H2O);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, &vCO2, &zCO2, &phiCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2, &dlnphidtCO2, &dlnphidpCO2, &d2lnphidt2CO2, &d2lnphidtdpCO2, &d2lnphidp2CO2);

    d2gdt2 = 2.0*R*(x[H2O]*dlnphidt[H2O] - x[H2O]*dlnphidtH2O + x[CO2]*dlnphidt[CO2] - x[CO2]*dlnphidtCO2)
    + R*t*(x[H2O]*d2lnphidt2[H2O] - x[H2O]*d2lnphidt2H2O + x[CO2]*d2lnphidt2[CO2] - x[CO2]*d2lnphidt2CO2);

    if (mask & FIRST) {
        *cpmix = - t*d2gdt2;
    }

    if(mask & SECOND) {
#ifdef DO_NOT_THROW_EXCEPTIONS
        *dt = 0.0;
#else
        NSException *duanException = [NSException exceptionWithName:@"Duan fluid Internal Error"
                                                             reason:@"dcp(mix)dt not yet implemented."
                                                           userInfo:nil];
		@throw duanException;
#endif
    }

    if(mask & THIRD) {
        double d3gdrdt2[NR];
        int i;

        d3gdrdt2[0] = 2.0*R*(- dlnphidt[H2O] + dlnphidtH2O  + dlnphidt[CO2] - dlnphidtCO2)
        + R*t*(- d2lnphidt2[H2O] + d2lnphidt2H2O + d2lnphidt2[CO2] - d2lnphidt2CO2);

        for (i=0; i<NR; i++) dr[i] = -t*d3gdrdt2[i];
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
    double x[2], v, z, phi[2], dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt[2], dlnphidp[2], d2lnphidt2[2], d2lnphidtdp[2], d2lnphidp2[2], dlnphidr[2];
    double vH2O, vCO2, zH2O, zCO2, phiH2O, phiCO2, dvdtH2O, dvdpH2O, dvdtCO2, dvdpCO2, d2vdt2H2O,
    d2vdtdpH2O, d2vdp2H2O, d2vdt2CO2, d2vdtdpCO2, d2vdp2CO2, dlnphidtH2O, dlnphidtCO2, dlnphidpH2O, dlnphidpCO2,
    d2lnphidt2H2O, d2lnphidtdpH2O, d2lnphidp2H2O, d2lnphidt2CO2, d2lnphidtdpCO2, d2lnphidp2CO2;

    x[H2O] = 1.0 - r[0];
    x[CO2] = r[0];

    if ((fabs(x[CO2]) < 100.0*DBL_EPSILON) || (fabs(x[H2O]) < 100.0*DBL_EPSILON)) {
        if (mask & FIRST)   { *vmix     = 0.0; }
        if (mask & SECOND)  { dr[0]     = 0.0; }
        if (mask & THIRD)   { dr2[0][0] = 0.0; }
        if (mask & FOURTH)  { *dt       = 0.0; }
        if (mask & FIFTH)   { *dp       = 0.0; }
        if (mask & SIXTH)   { *dt2      = 0.0; }
        if (mask & SEVENTH) { *dtdp     = 0.0; }
        if (mask & EIGHTH)  { *dp2      = 0.0; }
        if (mask & NINTH)   { *drdt     = 0.0; }
        if (mask & TENTH)   { *drdp     = 0.0; }
        return;
    }

    duan((p <= 2000.0) ? 1 : 0, t, p, x, &v, &z, phi, &dvdt, &dvdp, &d2vdt2, &d2vdtdp, &d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2, dlnphidr);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, &vH2O, &zH2O, &phiH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, &dlnphidtH2O, &dlnphidpH2O, &d2lnphidt2H2O, &d2lnphidtdpH2O, &d2lnphidp2H2O);
    duanCO2((p <= 2000.0) ? 1 : 0, t, p, &vCO2, &zCO2, &phiCO2, &dvdtCO2, &dvdpCO2, &d2vdt2CO2, &d2vdtdpCO2, &d2vdp2CO2, &dlnphidtCO2, &dlnphidpCO2, &d2lnphidt2CO2, &d2lnphidtdpCO2, &d2lnphidp2CO2);

    if (mask & FIRST) {
        *vmix = v - x[0]*vH2O - x[1]*vCO2;
    }

    if(mask & SECOND) {
        double d2gdrdp[NR];
        int i;

        d2gdrdp[0] = R*t*(dlnphidp[CO2] - dlnphidp[H2O]) - R*t*(dlnphidpCO2 - dlnphidpH2O);

        for (i=0; i<NR; i++) dr[i] = d2gdrdp[i];
    }

    if(mask & THIRD) {
#ifdef DO_NOT_THROW_EXCEPTIONS
        int i,j;
        for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr2[i][j] = 0.0;
#else
        NSException *duanException = [NSException exceptionWithName:@"Duan fluid Internal Error"
                                                             reason:@"d2v(mix)dr2 not yet implemented."
                                                           userInfo:nil];
		@throw duanException;
#endif
    }

    if(mask & FOURTH) {
        *dt = dvdt - x[0]*dvdtH2O - x[1]*dvdtCO2;
    }

    if(mask & FIFTH) {
        *dp = dvdp - x[0]*dvdpH2O - x[1]*dvdpCO2;
    }

    if(mask & SIXTH) {
        *dt2 = d2vdt2 - x[0]*d2vdt2H2O - x[1]*d2vdt2CO2;
    }

    if(mask & SEVENTH) {
        *dtdp = d2vdtdp - x[0]*d2vdtdpH2O - x[1]*d2vdtdpCO2;
    }

    if(mask & EIGHTH) {
        *dp2 = d2vdp2 - x[0]*d2vdp2H2O - x[1]*d2vdp2CO2;
    }

    if(mask & NINTH) {
        double d3gdrdtdp[NR];
        int i;

        d3gdrdtdp[0] = R*(- dlnphidp[H2O] + dlnphidpH2O + dlnphidp[CO2] - dlnphidpCO2)
        + R*t*(- d2lnphidtdp[H2O] + d2lnphidtdpH2O + d2lnphidtdp[CO2] - d2lnphidtdpCO2);

        for (i=0; i<NR; i++) drdt[i] = d3gdrdtdp[i];
    }

    if(mask & TENTH) {
        double d3gdrdp2[NR];
        int i;

        d3gdrdp2[0] = R*t*(- d2lnphidp2[H2O] + d2lnphidp2H2O + d2lnphidp2[CO2] - d2lnphidp2CO2);

        for (i=0; i<NR; i++) drdp[i] = d3gdrdp2[i];
    }
}

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"H2O";
			break;
		case 1:
			return @"CO2";
			break;
		default:
			return @"";
			break;
	}
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
	DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
	double *mComponents = [mComponentWrapper pointerToDouble];
	for (NSUInteger i=0; i<NA; i++) mComponents[i] = mSpecies[i];
	return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	return [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
	return [[endmembers objectAtIndex:index] formulaAsElementArray];
}

-(void)correctActivityCoefficients:(double [NA])gamma forComposition:(double [NA])x {
	for (NSUInteger i=0, j=0; i<NA; i++) if (x[i] != 0.0) {
		if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
		else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
		j++;
	}
}

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
	NSMutableArray *results = [NSMutableArray arrayWithCapacity:NA+1];
	double deltaMu[NA], xNz[NA], x[NA], gamma[NA], xLast[NA], affinity = 0.0;
	NSUInteger i, j, nz = 0, index[NA];

	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");

	// Compute solid -> liquid delta mus and deflate composition space
    for (i=0; i<NA; i++) { deltaMu[i] = 0.0; index[i] = 0.0; gamma[i] = 0.0; }
	for (i=0; i<NA; i++) {
		if (chemicalPotentials[i] != 0.0) {
			deltaMu[nz] = chemicalPotentials[i] - [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
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

		if (debugV) NSLog(@"... Terminated. Trivial case.");
		return [NSArray arrayWithArray:results];
	}

	NSUInteger count = 0;
	BOOL converged = NO;
	double affinityLast;
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
		double a[NA], r[NR];
		[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:x r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
		[self activity:FIRST t:t p:p r:r a:a mu:NULL dx:NULL];

		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];
		[self correctActivityCoefficients:gamma forComposition:x];
		converged = (fabs(affinity-affinityLast) < 0.1);
		count++;

	} while (count < 50 && !converged);

	if (debugS) {
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
	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");


	[results addObject:[NSNumber numberWithDouble:affinity]];                    // affinity in J
	for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];   // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                     // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];              // iteration count
	[results addObject:[NSNumber numberWithDouble:NATOMS]];                      // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]]; // likely error in affinity

	return [NSArray arrayWithArray:results];
}

#pragma mark -
#pragma mark Solution Phase protocol functions


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
		double *componentToElements = [[[endmembers objectAtIndex:i] formulaAsElementArray] pointerToDouble];
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

@end
