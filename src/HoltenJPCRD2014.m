//
//  HoltenJPCRD2014.m
//  ThermoFit
//
//	Program implementing the equation of state in the article
//	"Equation of state for supercooled water at pressures up to 400 MPa",
//	Journal of Physical and Chemical Reference Data,
//	Vincent Holten, Jan V. Sengers, and Mikhail A. Anisimov,
//	Institute for Physical Science and Technology and
//	Department of Chemical and Biomolecular Engineering,
//	University of Maryland, College Park, Maryland 20742, U.S.A.
//
//	August 2014
//
//  Created by Mark Ghiorso on 2/29/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#import "HoltenJPCRD2014.h"

#pragma mark -
#pragma mark Model parameters

static const double Tc     =  228.2;
static const double Pc     =    0;
static const double rho0   = 1081.6482;
static const double R      =  461.523087;
static const double omega0 =    0.52122690;
static const double L0     =    0.76317954;
static const double k0     =    0.072158686;
static const double k1     =  - 0.31569232;
static const double k2     =    5.2992608;

static const double MolesPerKg = 55.508435;

#pragma mark -
#pragma mark Model background coefficients

static const int n = 20;

static const double c[n] = {-8.1570681381655, 1.2875032e+000, 7.0901673598012,
    -3.2779161e-002, 7.3703949e-001, -2.1628622e-001, -5.1782479e+000,
    4.2293517e-004, 2.3592109e-002, 4.3773754e+000, -2.9967770e-003,
    -9.6558018e-001, 3.7595286e+000, 1.2632441e+000, 2.8542697e-001,
    -8.5994947e-001, -3.2916153e-001, 9.0019616e-002, 8.1149726e-002,
    -3.2788213e+000};
static const double a[n] = {0, 0, 1, -0.2555, 1.5762, 1.64, 3.6385, -0.3828,
    1.6219, 4.3287, 3.4763, 5.1556, -0.3593, 5.0361, 2.9786, 6.2373,
    4.046, 5.3558, 9.0157, 1.2194};
static const double b[n] = {0, 1, 0, 2.1051, 1.1422, 0.951, 0, 3.6402,
    2.076, -0.0016, 2.2769, 0.0008, 0.3706, -0.3975, 2.973, -0.318,
    2.9805, 2.9265, 0.4456, 0.1298};
static const double d[n] = {0, 0, 0, -0.0016, 0.6894, 0.013, 0.0002, 0.0435,
    0.05, 0.0004, 0.0528, 0.0147, 0.8584, 0.9924, 1.0041, 1.0961,
    1.0228, 1.0303, 1.618, 0.5213};

@interface HoltenJPCRD2014() {
    double tOld, pOld;
    double G, S, rho, Kap, Alp, CP, CV, U;
    double EPS, EPS2;
}

@end

@implementation HoltenJPCRD2014

- (instancetype)init {
    self = [super init];
    if (self) {
        self.phaseName    = @"Supercooled Water";
        self.phaseFormula = @"H2O";
        tOld = -1.0;
        pOld = -1.0;
        EPS = sqrt(DBL_EPSILON);
        EPS2 = sqrt(sqrt(DBL_EPSILON));
    }
    return self;
}

#pragma mark -
#pragma mark Model Background derivatives

double B(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * pow(tau, a[i]) * pow(pi, b[i]) * exp(-d[i]*pi);
    return sum;
}

double Bp(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * pow(tau, a[i]) * pow(pi, b[i]-1) * (b[i]-d[i]*pi) * exp(-d[i]*pi);
    return sum;
}

double Bt(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * a[i] * pow(tau,a[i]-1) * pow(pi, b[i]) * exp(-d[i]*pi);
    return sum;
}

double Bpp(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * pow(tau,a[i]) * pow(pi, b[i]-2) * (pow(d[i]*pi - b[i], 2) - b[i]) * exp(-d[i]*pi);
    return sum;
}

double Btp(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * a[i] * pow(tau,a[i]-1) * pow(pi, b[i]-1) * (b[i]-d[i]*pi) *  exp(-d[i]*pi);
    return sum;
}

double Btt(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * a[i] * (a[i]-1) * pow(tau,a[i]-2) * pow(pi, b[i]) * exp(-d[i]*pi);
    return sum;
}

#pragma mark -
#pragma mark Functions for calculating the equilibrium fraction xe

double xefun(double x, double L, double W) {
    return L + log(x/(1-x)) + W*(1-2*x);
}

- (double)findxe:(double)L and:(double)W {
    double x0, x1;
    bool flip = false;

    if (L < 0) {
        flip = true;
        L = -L;
    }

    // Find starting values for x
    if (W < 1.1111111*(2.944439 - L)) { // xe = 0.05 isoline, W = (10/9) * (ln(19) - L)
        x0 = 0.049;
        x1 = 0.5;
    } else if (W < 1.0204081*(4.595119 - L)) { // xe = 0.01 isoline, W = (50/49) * (ln(99) - L)
        x0 = 0.0099;
        x1 = 0.051;
    } else {
        x0 = 0.99 * exp(-1.0204081 * L - W);
        x1 = 1.01 * 1.087 * exp(-L - W);
        if (x1 > 0.0101) x1 = 0.0101;
    }

    double y0 = xefun(x0, L, W),
    y1 = xefun(x1, L, W);
    NSAssert(y0*y1 < 0, @"Error in findxe(): starting values for x are incorrect.");

    // Bisection algorithm
    // This could be replaced with a root-finding function from a numerical library
    int N = 0;
    while (fabs(x1 - x0) > 10.0*DBL_EPSILON) {
        double x = (x0 + x1) / 2.0;
        double y = xefun(x, L, W);
        if (y0 * y >= 0) { //same sign
            x0 = x;
            y0 = y;
        } else {
            x1 = x;
//            y1 = y;
        }
        NSAssert(N++ < 50, @"Error in findxe(): bisection does not converge.");
    }

    double x = (x0 + x1) / 2.0;
    return (flip) ? 1.0 - x : x;
}

#pragma mark -
#pragma mark Main property evaluation function

- (void)evalAtTinK:(double)T andPinPa:(double)P {
    if ((fabs(T - tOld) < 10.0*DBL_EPSILON) && (fabs(P - pOld) < 10.0*DBL_EPSILON)) return;

    const double P0 = -300e6;

    // Dimensionless temperature and pressure
    double t   =  (T - Tc)/Tc;
    double p   =  (P - Pc)/(rho0*R*Tc);
    double tau = T/Tc;
    double pi  = (P - P0)/(rho0*R*Tc);

    // Field L and its derivatives
    double K1  = sqrt(pow(1.0+k0*k2+k1*(p-k2*t), 2.0) - 4.0*k0*k1*k2*(p-k2*t));
    double K3  = pow(K1, 3.0);
    double K2  = sqrt(1.0 + k2*k2);
    double L   = L0 * K2 * (1.0 - K1 + k0*k2 + k1*(p + k2*t)) / (2.0*k1*k2);
    double Lt  = L0 * 0.5 * K2 * (1 + (1 - k0*k2 + k1*(p - k2*t))/K1);
    double Lp  = L0 * K2 * (K1 + k0*k2 - k1*p + k1*k2*t - 1) / (2.0*k2*K1);
    double Ltt = -2.0*L0*K2*k0*k1*k2*k2 / K3;
    double Ltp =  2.0*L0*K2*k0*k1*k2 / K3;
    double Lpp = -2.0*L0*K2*k0*k1 / K3;

    // Interaction parameter omega
    double omega = 2.0 + omega0 * p;

    // Calculate equilibrium fraction xe
    double x = [self findxe:L and:omega];

    // Order parameter f and susceptibility chi
    double f   = 2.0 * x - 1.0;
    double f2  = f*f;
    double chi = 1.0 / (2.0 / (1.0 - f2) - omega);

    // Dimensionless properties
    double g0  = x*L + x*log(x) + (1.0-x)*log(1.0-x) + omega*x*(1-x); // g0 = (g - gA)/tau
    double g   = B(tau,pi) + tau*g0;
    double s   = -0.5*(f+1)*Lt*tau - g0 - Bt(tau,pi);
    double v   = 0.5 * tau * (omega0/2.0*(1.0-f2) + Lp*(f+1)) + Bp(tau,pi);
    double kap = (1.0/v) * (tau/2.0 * (chi * pow(Lp - omega0*f, 2.0) - (f+1.0)*Lpp) - Bpp(tau,pi));
    double alp = (1.0/v) * (Ltp/2.0 * tau*(f+1.0) + (omega0/2.0*(1.0-f2) + Lp*(f+1.0))/2.0 - tau*Lt/2.0 * chi*(Lp - omega0*f) + Btp(tau,pi));
    double cp  = tau * ( -Lt * (f+1) + tau*(Lt*Lt*chi - Ltt*(f+1)) / 2 - Btt(tau,pi));

    // Properties in SI units
    G   = R * Tc * g;					  // Specific Gibbs energy       J/kg
    S   = R * s;						  // Specific entropy            J/kg-K
    rho = rho0 / v;						  // Density                     kg/m^3
    Kap = kap / (rho0 * R * Tc);		  // Isothermal compressibility  Pa
    Alp = alp / Tc;						  // Expansion coefficient       1/K
    CP  = R * cp;						  // Isobaric heat capacity      J/kg-K
    CV  = CP - T*Alp*Alp / (rho * Kap);   // Isochoric heat capacity     J/kg-K
    U   = 1/sqrt(rho*Kap - T*Alp*Alp/CP); // Speed of sound              m/s

    tOld = T;
    pOld = P;

//    printf("Temperature   \t%f\tK\n"
//           "Pressure      \t%f\tMPa\n"
//           "Fraction xe   \t%.8f\t\n"
//           "Gibbs energy  \t%f\tJ/kg\n"
//           "Entropy       \t%f\tJ/(kg K)\n"
//           "Density       \t%.5f\tkg/m3\n"
//           "Isoth. comp.  \t%e\t1/MPa\n"
//           "Expansivity   \t%e\t1/K\n"
//           "CP            \t%.4f\tJ/(kg K)\n"
//           "CV            \t%.4f\tJ/(kg K)\n"
//           "Speed of sound\t%.4f\tm/s\n\n",
//           T, 1e-6*P, x, G, S, rho, 1e6 * Kap, Alp, CP, CV, U);
}

#pragma mark -
#pragma mark Additional public methods

- (double)homogeneousIceNucleationTemperatureForPressureInBars:(double)p {
    double p1 = p/10.0; // MPa
    double result = 172.82 + 0.03718*p1 + 3.403e-5*p1*p1 - 1.573e-8*p1*p1*p1;
    return (result > 0.0) ? result : 1.0;
}

- (double)homogeneousIceNucleationPressureForTemperatureInK:(double)t {
    double p0 = 0.1;    // MPa
    double t0 = 235.15; // K
    double theta = t/t0;
    double result = 1.0 + 2282.7*(1.0-pow(theta, 6.243)) + 157.24*(1.0-pow(theta, 79.81));
    return result*p0*10.0; //bars
}

#pragma mark -
#pragma mark Stoichiometric Phase Protocol

- (double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return G/MolesPerKg; // J/mol
}

- (double)getEnthalpyFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return (G+t*S)/MolesPerKg; // J/mol
}

- (double)getEntropyFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return S/MolesPerKg; // J/mol-K
}

- (double)getHeatCapacityFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return CP/MolesPerKg; // J/mol-K
}

- (double)getDcpDtFromT:(double)t andP:(double)p {
    [self evalAtTinK:t*(1.0+EPS) andPinPa:p*1.0e5];
    double result = CP;
    [self evalAtTinK:t*(1.0-EPS) andPinPa:p*1.0e5];
    result -= CP;
    return result/2.0/t/EPS/MolesPerKg;
}

- (double)getVolumeFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return 1.0e5/rho/MolesPerKg; // J/mol-bar
}

- (double)getDvDtFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return Alp*1.0e5/rho/MolesPerKg; // J/mol-bar-K
}

- (double)getDvDpFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5];
    return  -1.0e10*Kap/rho/MolesPerKg; // J/mol-bar^2
}

- (double)getD2vDt2FromT:(double)t andP:(double)p {
    [self evalAtTinK:t*(1.0+EPS2) andPinPa:p*1.0e5];
    double result = Alp/rho;
    [self evalAtTinK:t*(1.0-EPS2) andPinPa:p*1.0e5];
    result -= Alp/rho;
    return result*1.0e5/2.0/t/EPS2/MolesPerKg;
}

- (double)getD2vDtDpFromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5*(1.0-2.0*EPS2)];
    double result = Alp/rho/12.0;
    [self evalAtTinK:t andPinPa:p*1.0e5*(1.0-EPS2)];
    result -= 2.0*Alp/rho/3.0;
    [self evalAtTinK:t andPinPa:p*1.0e5*(1.0+EPS2)];
    result += 2.0*Alp/rho/3.0;
    [self evalAtTinK:t andPinPa:p*1.0e5*(1.0+2.0*EPS2)];
    result -= Alp/rho/12.0;
    return result*1.0e5/p/EPS2/MolesPerKg; // *1.0e5
}

- (double)getD2vDp2FromT:(double)t andP:(double)p {
    [self evalAtTinK:t andPinPa:p*1.0e5*(1.0+EPS2)];
    double result = Kap/rho;
    [self evalAtTinK:t andPinPa:p*1.0e5*(1.0-EPS2)];
    result -= Kap/rho;
    return -result*1.0e10/2.0/p/EPS2/MolesPerKg;
}

@end
