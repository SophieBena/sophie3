/* 
	Program implementing the equation of state in the article
	"Equation of state for supercooled water at pressures up to 400 MPa",
	Journal of Physical and Chemical Reference Data,
	Vincent Holten, Jan V. Sengers, and Mikhail A. Anisimov,
	Institute for Physical Science and Technology and
	Department of Chemical and Biomolecular Engineering,
	University of Maryland, College Park, Maryland 20742, U.S.A.

	August 2014
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Parameters

const double Tc = 228.2, Pc = 0,
	rho0 = 1081.6482,
	R = 461.523087,
	omega0 = 0.52122690,
	L0 = 0.76317954,
	k0 = 0.072158686,
	k1 = -0.31569232,
	k2 = 5.2992608;

// Background coefficients

const int n = 20;

const double c[n] = {-8.1570681381655, 1.2875032e+000, 7.0901673598012,
	-3.2779161e-002, 7.3703949e-001, -2.1628622e-001, -5.1782479e+000,
	4.2293517e-004, 2.3592109e-002, 4.3773754e+000, -2.9967770e-003,
	-9.6558018e-001, 3.7595286e+000, 1.2632441e+000, 2.8542697e-001,
	-8.5994947e-001, -3.2916153e-001, 9.0019616e-002, 8.1149726e-002,
	-3.2788213e+000};
const double a[n] = {0, 0, 1, -0.2555, 1.5762, 1.64, 3.6385, -0.3828,
	1.6219, 4.3287, 3.4763, 5.1556, -0.3593, 5.0361, 2.9786, 6.2373,
	4.046, 5.3558, 9.0157, 1.2194};
const double b[n] = {0, 1, 0, 2.1051, 1.1422, 0.951, 0, 3.6402,
	2.076, -0.0016, 2.2769, 0.0008, 0.3706, -0.3975, 2.973, -0.318,
	2.9805, 2.9265, 0.4456, 0.1298};
const double d[n] = {0, 0, 0, -0.0016, 0.6894, 0.013, 0.0002, 0.0435,
	0.05, 0.0004, 0.0528, 0.0147, 0.8584, 0.9924, 1.0041, 1.0961,
	1.0228, 1.0303, 1.618, 0.5213};

// Background derivatives

double B(double tau, double pi)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += c[i] * pow(tau, a[i]) * pow(pi, b[i]) * exp(-d[i]*pi);
	return sum;
}

double Bp(double tau, double pi)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += c[i] * pow(tau, a[i]) * pow(pi, b[i]-1) * (b[i]-d[i]*pi) * exp(-d[i]*pi);
	return sum;
}

double Bt(double tau, double pi)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += c[i] * a[i] * pow(tau,a[i]-1) * pow(pi, b[i]) * exp(-d[i]*pi);
	return sum;
}

double Bpp(double tau, double pi)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += c[i] * pow(tau,a[i]) * pow(pi, b[i]-2) * (pow(d[i]*pi - b[i], 2) - b[i]) * exp(-d[i]*pi);
	return sum;
}

double Btp(double tau, double pi)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += c[i] * a[i] * pow(tau,a[i]-1) * pow(pi, b[i]-1) * (b[i]-d[i]*pi) *  exp(-d[i]*pi);
	return sum;
}

double Btt(double tau, double pi)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += c[i] * a[i] * (a[i]-1) * pow(tau,a[i]-2) * pow(pi, b[i]) * exp(-d[i]*pi);
	return sum;
}

// Functions for calculating the equilibrium fraction xe

double xefun(double x, double L, double W)
{
	return L + log(x/(1-x)) + W*(1-2*x);
}

double findxe(double L, double W)
{
	double x0, x1;
	bool flip = false;

	if (L < 0)
	{
		flip = true;
		L = -L;
	}

	// Find starting values for x
	if (W < 1.1111111*(2.944439 - L)) // xe = 0.05 isoline, W = (10/9) * (ln(19) - L)
	{
		x0 = 0.049;
		x1 = 0.5;
	}
	else if (W < 1.0204081*(4.595119 - L)) // xe = 0.01 isoline, W = (50/49) * (ln(99) - L)
	{
		x0 = 0.0099;
		x1 = 0.051;
	}
	else
	{
		x0 = 0.99 * exp(-1.0204081 * L - W);
		x1 = 1.01 * 1.087 * exp(-L - W);
		if (x1 > 0.0101) x1 = 0.0101;
	}

	double y0 = xefun(x0, L, W),
		   y1 = xefun(x1, L, W);

	if (y0 * y1 >= 0)  // same sign, no zero crossing
	{
		printf("Error in findxe(): starting values for x are incorrect.\n");
		exit(1);
	}

	// Bisection algorithm
	// This could be replaced with a root-finding function from a numerical library
	int N = 0;
	while (fabs(x1 - x0) > 1e-10)
	{
		double x = (x0 + x1) / 2;
		double y = xefun(x, L, W);
		if (y0 * y >= 0) //same sign
		{
			x0 = x;
			y0 = y;
		}
		else
		{
			x1 = x;
			y1 = y;
		}
		if (N++ > 50)
		{
			printf("Error in findxe(): bisection does not converge.\n");
			exit(1);
		}
	}

	double x = (x0 + x1) / 2;
	return flip? 1 - x : x;
}

// Main property evaluation function

void eval(double T, double P) // T in K, P in Pa
{
	const double P0 = -300e6;

	// Dimensionless temperature and pressure
	double t =  (T - Tc) / Tc,
		   p =  (P - Pc) / (rho0 * R * Tc),
		   tau = T / Tc,
		   pi = (P - P0) / (rho0 * R * Tc);

	// Field L and its derivatives
	double K1 = sqrt(pow(1+k0*k2+k1*(p-k2*t), 2) - 4*k0*k1*k2*(p-k2*t));
	double
		K3 = pow(K1, 3),
		K2 = sqrt(1 + k2*k2),
		L = L0 * K2 * (1 - K1 + k0*k2 + k1*(p + k2*t)) / (2*k1*k2),
		Lt = L0 * 0.5 * K2 * (1 + (1 - k0*k2 + k1*(p - k2*t))/K1),
		Lp = L0 * K2 * (K1 + k0*k2 - k1*p + k1*k2*t - 1) / (2*k2*K1),
		Ltt = -2*L0*K2*k0*k1*k2*k2 / K3,
		Ltp =  2*L0*K2*k0*k1*k2 / K3,
		Lpp = -2*L0*K2*k0*k1 / K3;

	// Interaction parameter omega
	double omega = 2 + omega0 * p;

	// Calculate equilibrium fraction xe
	double x = findxe(L, omega);

	// Order parameter f and susceptibility chi
	double f = 2 * x - 1,
		   f2 = f*f,
		   chi = 1 / (2 / (1 - f2) - omega);

	// Dimensionless properties
	double g0 = x*L + x*log(x) + (1-x)*log(1-x) + omega*x*(1-x), // g0 = (g - gA)/tau
		   g = B(tau,pi) + tau*g0,
		   s = -0.5*(f+1)*Lt*tau - g0 - Bt(tau,pi),
		   v = 0.5 * tau * (omega0/2*(1-f2) + Lp*(f+1)) + Bp(tau,pi),
		   kap = (1/v) * (tau/2 * (chi * pow(Lp - omega0 * f, 2) - (f+1)*Lpp) - Bpp(tau,pi)),
		   alp = (1/v) * (Ltp/2 * tau*(f+1) + (omega0/2*(1-f2) + Lp*(f+1))/2
				- tau*Lt/2 * chi*(Lp - omega0*f) + Btp(tau,pi)),
		   cp = tau * ( -Lt * (f+1) + tau*(Lt*Lt*chi - Ltt*(f+1)) / 2 - Btt(tau,pi));

	// Properties in SI units
	double G = R * Tc * g,						// Specific Gibbs energy
		   S = R * s,							// Specific entropy
		   rho = rho0 / v,						// Density
	       Kap = kap / (rho0 * R * Tc),			// Isothermal compressibility
		   Alp = alp / Tc,						// Expansion coefficient
		   CP = R * cp,							// Isobaric heat capacity
		   CV = CP - T*Alp*Alp / (rho * Kap),	// Isochoric heat capacity
		   U = 1/sqrt(rho*Kap - T*Alp*Alp/CP);	// Speed of sound

	printf("Temperature   \t%f\tK\n"
		   "Pressure      \t%f\tMPa\n"
		   "Fraction xe   \t%.8f\t\n"
		   "Gibbs energy  \t%f\tJ/kg\n"
		   "Entropy       \t%f\tJ/(kg K)\n"
		   "Density       \t%.5f\tkg/m3\n"
		   "Isoth. comp.  \t%e\t1/MPa\n"
		   "Expansivity   \t%e\t1/K\n"
		   "CP            \t%.4f\tJ/(kg K)\n"
		   "CV            \t%.4f\tJ/(kg K)\n"
		   "Speed of sound\t%.4f\tm/s\n\n",
		   T, 1e-6*P, x, G, S, rho, 1e6 * Kap, Alp, CP, CV, U);
}

int main()
{
	// Calculate test values to compare with table in article
	const int ntest = 6;
	double Tlist[ntest] = {  273.15,   235.15, 250, 200, 250, 273.16};
	double Plist[ntest] = {0.101325, 0.101325, 200, 400, 400, 611.65477100789442e-6};

	for (int j = 0; j < ntest; j++)
	{
		eval(Tlist[j], 1e6*Plist[j]);
	}
}
