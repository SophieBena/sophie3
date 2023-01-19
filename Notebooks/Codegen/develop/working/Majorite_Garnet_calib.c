
static char *identifier = "Thu Jan 31 09:09:40 2019";

static double T_r = 298.15;
static double P_r = 1.0;
static double H_TrPr = -6286547.62;
static double S_TrPr = 266.359;
static double k0 = 640.71997;
static double k1 = -4542.07;
static double k2 = -4701900;
static double k3 = 0;
static double V_TrPr = 11.316;
static double v1 = -5.762e-07;
static double v2 = 0;
static double v3 = 2.2518700000000002e-05;
static double v4 = 3.7e-09;
static double l1 = 0.0;
static double l2 = 0.0;
static double k_lambda = 0.0;
static double T_lambda_Pr = 0.0;
static double T_lambda_ref = 0.0;
static double H_t = 0.0;
static double d0 = 0.0;
static double d1 = 0.0;
static double d2 = 0.0;
static double d3 = 0.0;
static double d4 = 0.0;
static double d5 = 0.0;
static double T_D = 0.0;
static double T_D_ref = 0.0;


#include <math.h>

static double Garnet_g(double T, double P) {
    double result = 0.0;
    result += H_TrPr + 2*sqrt(T)*k1 + T*k0 - T*(S_TrPr + k0*log(T) - k0*log(T_r) + (1.0/2.0)*k2/((T_r)*(T_r)) + (1.0/3.0)*k3/((T_r)*(T_r)*(T_r)) + 2*k1/sqrt(T_r) - 1.0/2.0*k2/((T)*(T)) - 1.0/3.0*k3/((T)*(T)*(T)) - 2*k1/sqrt(T)) - 2*sqrt(T_r)*k1 - T_r*k0 + k2/T_r + (1.0/2.0)*k3/((T_r)*(T_r)) - k2/T - 1.0/2.0*k3/((T)*(T));
    result += (1.0/3.0)*((P)*(P)*(P))*V_TrPr*v2 + ((P)*(P))*(-P_r*V_TrPr*v2 + (1.0/2.0)*V_TrPr*v1) + P*(((P_r)*(P_r))*V_TrPr*v2 - P_r*V_TrPr*v1 + ((T)*(T))*V_TrPr*v4 - 2*T*T_r*V_TrPr*v4 + T*V_TrPr*v3 + ((T_r)*(T_r))*V_TrPr*v4 - T_r*V_TrPr*v3 + V_TrPr) - 1.0/3.0*((P_r)*(P_r)*(P_r))*V_TrPr*v2 - ((P_r)*(P_r))*(-P_r*V_TrPr*v2 + (1.0/2.0)*V_TrPr*v1) - P_r*(((P_r)*(P_r))*V_TrPr*v2 - P_r*V_TrPr*v1 + ((T)*(T))*V_TrPr*v4 - 2*T*T_r*V_TrPr*v4 + T*V_TrPr*v3 + ((T_r)*(T_r))*V_TrPr*v4 - T_r*V_TrPr*v3 + V_TrPr);
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += (1.0/4.0)*((T)*(T)*(T)*(T))*((l2)*(l2)) + ((T)*(T)*(T))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) + ((T)*(T))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - T*((1.0/3.0)*((T)*(T)*(T))*((l2)*(l2)) + ((T)*(T))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + T*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)))) + T*(-((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1))) - 1.0/4.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - (T_lambda_ref + k_lambda*(P - P_r))*(-((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1)));
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += (-T + T_lambda_Pr + k_lambda*(P - P_r))*(((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += H_t*(-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += H_t*(-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r));
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += 2*sqrt(T)*d1 + (1.0/3.0)*((T)*(T)*(T))*d4 + (1.0/2.0)*((T)*(T))*d3 + T*d0 - T*((1.0/2.0)*((T)*(T))*d4 + T*d3 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d4 - T_D_ref*d3 + 2*d0*log(sqrt(T)) - 2*d0*log(sqrt(T_D_ref)) + (1.0/2.0)*d2/((T_D_ref)*(T_D_ref)) + 2*d1/sqrt(T_D_ref) - 1.0/2.0*d2/((T)*(T)) - 2*d1/sqrt(T)) - 2*sqrt(T_D_ref)*d1 - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref))*d4 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d3 - T_D_ref*d0 + (P - P_r)*(2*sqrt(T)*d1 + (1.0/3.0)*((T)*(T)*(T))*d4 + (1.0/2.0)*((T)*(T))*d3 + T*d0 - 2*sqrt(T_D_ref)*d1 - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref))*d4 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d3 - T_D_ref*d0 + d2/T_D_ref - d2/T)/d5 + d2/T_D_ref - d2/T;
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D))) {
             result += (-T + T_D)*(-1.0/2.0*((T_D)*(T_D))*d4 + T_D*(T_D*d4 + d3 + d0/T_D + d2/((T_D)*(T_D)*(T_D)) + d1/pow(T_D, 3.0/2.0)) - 1.0/2.0*((T_D_ref)*(T_D_ref))*d4 - T_D_ref*d3 + 2*d0*log(sqrt(T_D)) - 2*d0*log(sqrt(T_D_ref)) - d0 - (P - P_r)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5 + (1.0/2.0)*d2/((T_D_ref)*(T_D_ref)) + 2*d1/sqrt(T_D_ref) - 3.0/2.0*d2/((T_D)*(T_D)) - 3*d1/sqrt(T_D));
        }
    }
    return result;
}

static double Garnet_dgdt(double T, double P) {
    double result = 0.0;
    result += -S_TrPr - T*(k0/T + k2/((T)*(T)*(T)) + k3/((T)*(T)*(T)*(T)) + k1/pow(T, 3.0/2.0)) - k0*log(T) + k0*log(T_r) + k0 - 1.0/2.0*k2/((T_r)*(T_r)) - 1.0/3.0*k3/((T_r)*(T_r)*(T_r)) - 2*k1/sqrt(T_r) + (3.0/2.0)*k2/((T)*(T)) + (4.0/3.0)*k3/((T)*(T)*(T)) + 3*k1/sqrt(T);
    result += P*(2*T*V_TrPr*v4 - 2*T_r*V_TrPr*v4 + V_TrPr*v3) - P_r*(2*T*V_TrPr*v4 - 2*T_r*V_TrPr*v4 + V_TrPr*v3);
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1)) + (2.0/3.0)*((T)*(T)*(T))*((l2)*(l2)) - ((T)*(T))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + 3*((T)*(T))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) + 2*T*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - T*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - T*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((T)*(T))*((l2)*(l2)) + 2*T*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + ((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + (1.0/3.0)*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)));
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + (2.0/3.0)*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (1.0/3.0)*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) + ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -H_t/(T_lambda_Pr + k_lambda*(P - P_r));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -H_t/(T_lambda_Pr + k_lambda*(P - P_r));
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += (1.0/2.0)*((T)*(T))*d4 - T*(T*d4 + d3 + d0/T + d2/((T)*(T)*(T)) + d1/pow(T, 3.0/2.0)) + (1.0/2.0)*((T_D_ref)*(T_D_ref))*d4 + T_D_ref*d3 - 2*d0*log(sqrt(T)) + 2*d0*log(sqrt(T_D_ref)) + d0 + (P - P_r)*(((T)*(T))*d4 + T*d3 + d0 + d2/((T)*(T)) + d1/sqrt(T))/d5 - 1.0/2.0*d2/((T_D_ref)*(T_D_ref)) - 2*d1/sqrt(T_D_ref) + (3.0/2.0)*d2/((T)*(T)) + 3*d1/sqrt(T);
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D))) {
             result += (1.0/2.0)*((T_D)*(T_D))*d4 - T_D*(T_D*d4 + d3 + d0/T_D + d2/((T_D)*(T_D)*(T_D)) + d1/pow(T_D, 3.0/2.0)) + (1.0/2.0)*((T_D_ref)*(T_D_ref))*d4 + T_D_ref*d3 - 2*d0*log(sqrt(T_D)) + 2*d0*log(sqrt(T_D_ref)) + d0 + (P - P_r)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5 - 1.0/2.0*d2/((T_D_ref)*(T_D_ref)) - 2*d1/sqrt(T_D_ref) + (3.0/2.0)*d2/((T_D)*(T_D)) + 3*d1/sqrt(T_D);
        }
    }
    return result;
}

static double Garnet_dgdp(double T, double P) {
    double result = 0.0;
    result += ((P)*(P))*V_TrPr*v2 + 2*P*(-P_r*V_TrPr*v2 + (1.0/2.0)*V_TrPr*v1) + ((P_r)*(P_r))*V_TrPr*v2 - P_r*V_TrPr*v1 + ((T)*(T))*V_TrPr*v4 - 2*T*T_r*V_TrPr*v4 + T*V_TrPr*v3 + ((T_r)*(T_r))*V_TrPr*v4 - T_r*V_TrPr*v3 + V_TrPr;
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -((T)*(T)*(T))*k_lambda*((l2)*(l2)) + ((T)*(T))*(3*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*k_lambda*l1*l2) + T*(-3*((P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*((k_lambda)*(k_lambda))*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*((k_lambda)*(k_lambda))*l1*l2 - k_lambda*((l1)*(l1))) - T*(-3.0/2.0*((T)*(T))*k_lambda*((l2)*(l2)) + T*(6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*k_lambda*l1*l2) + 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + (1.0/2.0)*k_lambda*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) + k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (-T_lambda_ref - k_lambda*(P - P_r))*(6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*k_lambda*l1*l2)) - 3*k_lambda*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - k_lambda*(-((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1))) + (-T_lambda_ref - k_lambda*(P - P_r))*(-3*((P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*((k_lambda)*(k_lambda))*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*((k_lambda)*(k_lambda))*l1*l2 - k_lambda*((l1)*(l1))) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(3*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*k_lambda*l1*l2);
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += k_lambda*(((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*((k_lambda)*(k_lambda))*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*((k_lambda)*(k_lambda))*l1*l2 + 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((l1)*(l1)) - 1.0/2.0*k_lambda*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (1.0/2.0)*k_lambda*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*k_lambda*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 6*k_lambda*(T_lambda_Pr + k_lambda*(P - P_r))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + k_lambda*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) + (-2*T_lambda_Pr - 2*k_lambda*(P - P_r))*(3*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*k_lambda*l1*l2) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*k_lambda*l1*l2) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*k_lambda*l1*l2 + 2*k_lambda*((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) - 3.0/2.0*k_lambda*((l2)*(l2))*(2*T_lambda_Pr + 2*k_lambda*(P - P_r)) + 2*k_lambda*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*k_lambda*l1*l2));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += H_t*k_lambda/(T_lambda_Pr + k_lambda*(P - P_r)) - H_t*k_lambda*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += H_t*k_lambda/(T_lambda_Pr + k_lambda*(P - P_r)) - H_t*k_lambda*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += (2*sqrt(T)*d1 + (1.0/3.0)*((T)*(T)*(T))*d4 + (1.0/2.0)*((T)*(T))*d3 + T*d0 - 2*sqrt(T_D_ref)*d1 - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref))*d4 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d3 - T_D_ref*d0 + d2/T_D_ref - d2/T)/d5;
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D))) {
             result += -(-T + T_D)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5;
        }
    }
    return result;
}

static double Garnet_d2gdt2(double T, double P) {
    double result = 0.0;
    result += (1.0/2.0)*T*(2*k0/((T)*(T)) + 6*k2/((T)*(T)*(T)*(T)) + 8*k3/pow(T, 5) + 3*k1/pow(T, 5.0/2.0)) - 2*k0/T - 4*k2/((T)*(T)*(T)) - 5*k3/((T)*(T)*(T)*(T)) - 5.0/2.0*k1/pow(T, 3.0/2.0);
    result += 2*V_TrPr*v4*(P - P_r);
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 + ((T)*(T))*((l2)*(l2)) - T*(-3*P*k_lambda*((l2)*(l2)) + 3*P_r*k_lambda*((l2)*(l2)) + 2*T*((l2)*(l2)) + 2*l1*l2 + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T))) - ((l1)*(l1)) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T;
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += (1.0/2.0)*T*(-2*d4 + 2*d0/((T)*(T)) + 6*d2/((T)*(T)*(T)*(T)) + 3*d1/pow(T, 5.0/2.0)) - d3 + (1.0/2.0)*(P - P_r)*(4*T*d4 + 2*d3 - 4*d2/((T)*(T)*(T)) - d1/pow(T, 3.0/2.0))/d5 - 2*d0/T - 4*d2/((T)*(T)*(T)) - 5.0/2.0*d1/pow(T, 3.0/2.0);
        }
    }
    return result;
}

static double Garnet_d2gdtdp(double T, double P) {
    double result = 0.0;
    result += V_TrPr*(2*T*v4 - 2*T_r*v4 + v3);
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += k_lambda*(-3.0/2.0*((T)*(T))*((l2)*(l2)) + T*(-6*P*k_lambda*((l2)*(l2)) + 6*P_r*k_lambda*((l2)*(l2)) + 3*T*((l2)*(l2)) + 4*l1*l2 - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 1.0/2.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += k_lambda*(-3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((l1)*(l1)) - 1.0/2.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/2.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += H_t*k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += H_t*k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += (((T)*(T))*d4 + T*d3 + d0 + d2/((T)*(T)) + d1/sqrt(T))/d5;
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D))) {
             result += (((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5;
        }
    }
    return result;
}

static double Garnet_d2gdp2(double T, double P) {
    double result = 0.0;
    result += V_TrPr*(2*P*v2 - 2*P_r*v2 + v1);
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += ((k_lambda)*(k_lambda))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + 3*((T)*(T))*((l2)*(l2)) + 2*T*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + T*(-6*T*((l2)*(l2)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 3*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r))) + ((l1)*(l1)));
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += ((k_lambda)*(k_lambda))*(6*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 12*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*k_lambda*l1*l2 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P_r*k_lambda*l1*l2 + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-18*P*k_lambda*((l2)*(l2)) + 18*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*l1*l2 + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 4*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += 2*H_t*((k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += 2*H_t*((k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    return result;
}

static double Garnet_d3gdt3(double T, double P) {
    double result = 0.0;
    result += -1.0/4.0*T*(8*k0/((T)*(T)*(T)) + 48*k2/pow(T, 5) + 80*k3/pow(T, 6) + 15*k1/pow(T, 7.0/2.0)) + 3*k0/((T)*(T)) + 15*k2/((T)*(T)*(T)*(T)) + 24*k3/pow(T, 5) + (21.0/4.0)*k1/pow(T, 5.0/2.0);
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += 3*P*k_lambda*((l2)*(l2)) - 3*P_r*k_lambda*((l2)*(l2)) - 2*T*(((l2)*(l2)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T)*(T))) - 2*l1*l2 - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T));
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += -1.0/4.0*T*(8*d0/((T)*(T)*(T)) + 48*d2/pow(T, 5) + 15*d1/pow(T, 7.0/2.0)) - d4 + (1.0/4.0)*(P - P_r)*(8*d4 + 24*d2/((T)*(T)*(T)*(T)) + 3*d1/pow(T, 5.0/2.0))/d5 + 3*d0/((T)*(T)) + 15*d2/((T)*(T)*(T)*(T)) + (21.0/4.0)*d1/pow(T, 5.0/2.0);
        }
    }
    return result;
}

static double Garnet_d3gdt2dp(double T, double P) {
    double result = 0.0;
    result += 2*V_TrPr*v4;
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += k_lambda*(-6*P*k_lambda*((l2)*(l2)) + 6*P_r*k_lambda*((l2)*(l2)) + T*(3*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T))) + 4*l1*l2 - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T);
        }
    }
    if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
        if((T > (T_D_ref)) && (T <= (T_D))) {
             result += (2*T*d4 + d3 - 2*d2/((T)*(T)*(T)) - 1.0/2.0*d1/pow(T, 3.0/2.0))/d5;
        }
    }
    return result;
}

static double Garnet_d3gdtdp2(double T, double P) {
    double result = 0.0;
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += ((k_lambda)*(k_lambda))*(-6*P*k_lambda*((l2)*(l2)) + 6*P_r*k_lambda*((l2)*(l2)) - 2*T*l2*(3*l2 - k_lambda*l2*(P - P_r)/T + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 4*l1*l2 + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 3*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += ((k_lambda)*(k_lambda))*(-18*P*k_lambda*((l2)*(l2)) + 18*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*l1*l2 + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 4*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -2*H_t*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += -2*H_t*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    return result;
}

static double Garnet_d3gdp3(double T, double P) {
    double result = 0.0;
    result += 2*V_TrPr*v2;
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += ((k_lambda)*(k_lambda)*(k_lambda))*(-6*T*((l2)*(l2)) - T*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 11*((l2)*(l2)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))) - 2*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1));
        }
    }
    if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += ((k_lambda)*(k_lambda)*(k_lambda))*(54*P*k_lambda*((l2)*(l2)) - 54*P_r*k_lambda*((l2)*(l2)) - 6*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 24*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 9*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 36*l1*l2 - 6*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + 12*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 6*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 6*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 12*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += 6*H_t*((k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    if(((H_t) != 0.0)) {
        if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
             result += 6*H_t*((k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
        }
    }
    return result;
}


static double Garnet_s(double T, double P) {
    double result = -Garnet_dgdt(T, P);
    return result;
}

static double Garnet_v(double T, double P) {
    double result = Garnet_dgdp(T, P);
    return result;
}

static double Garnet_cv(double T, double P) {
    double result = -T*Garnet_d2gdt2(T, P);
    double dvdt = Garnet_d2gdtdp(T, P);
    double dvdp = Garnet_d2gdp2(T, P);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double Garnet_cp(double T, double P) {
    double result = -T*Garnet_d2gdt2(T, P);
    return result;
}

static double Garnet_dcpdt(double T, double P) {
    double result = -T*Garnet_d3gdt3(T, P) - Garnet_d2gdt2(T, P);
    return result;
}

static double Garnet_alpha(double T, double P) {
    double result = Garnet_d2gdtdp(T, P)/Garnet_dgdp(T, P);
    return result;
}

static double Garnet_beta(double T, double P) {
    double result = -Garnet_d2gdp2(T, P)/Garnet_dgdp(T, P);
    return result;
}

static double Garnet_K(double T, double P) {
    double result = -Garnet_dgdp(T, P)/Garnet_d2gdp2(T, P);
    return result;
}

static double Garnet_Kp(double T, double P) {
    double result = Garnet_dgdp(T, P);
    result *= Garnet_d3gdp3(T, P);
    result /= pow(Garnet_d2gdp2(T, P), 2.0);
    return result - 1.0;
}


#include <math.h>

static double Garnet_dparam_g(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += -T*(-k0/T_r - k2/((T_r)*(T_r)*(T_r)) - k3/((T_r)*(T_r)*(T_r)*(T_r)) - k1/pow(T_r, 3.0/2.0)) - k0 - k2/((T_r)*(T_r)) - k3/((T_r)*(T_r)*(T_r)) - k1/sqrt(T_r);
        result += P*(-2*T*V_TrPr*v4 + 2*T_r*V_TrPr*v4 - V_TrPr*v3) - P_r*(-2*T*V_TrPr*v4 + 2*T_r*V_TrPr*v4 - V_TrPr*v3);
        break;
    case 1: /* P_r */ 
        result += -((P)*(P))*V_TrPr*v2 + P*(2*P_r*V_TrPr*v2 - V_TrPr*v1) - ((P_r)*(P_r))*V_TrPr*v2 + P_r*V_TrPr*v1 - 2*P_r*(-P_r*V_TrPr*v2 + (1.0/2.0)*V_TrPr*v1) - P_r*(2*P_r*V_TrPr*v2 - V_TrPr*v1) - ((T)*(T))*V_TrPr*v4 + 2*T*T_r*V_TrPr*v4 - T*V_TrPr*v3 - ((T_r)*(T_r))*V_TrPr*v4 + T_r*V_TrPr*v3 - V_TrPr;
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((T)*(T)*(T))*k_lambda*((l2)*(l2)) + ((T)*(T))*(-3*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 3*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*k_lambda*l1*l2) + T*(3*((P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*((k_lambda)*(k_lambda))*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*((k_lambda)*(k_lambda))*l1*l2 + k_lambda*((l1)*(l1))) - T*((3.0/2.0)*((T)*(T))*k_lambda*((l2)*(l2)) + T*(-6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*k_lambda*l1*l2) - 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 1.0/2.0*k_lambda*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) - k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (-T_lambda_ref - k_lambda*(P - P_r))*(-6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*k_lambda*l1*l2)) + 3*k_lambda*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) + 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + k_lambda*(-((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1))) + (-T_lambda_ref - k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*((k_lambda)*(k_lambda))*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*((k_lambda)*(k_lambda))*l1*l2 + k_lambda*((l1)*(l1))) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 3*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*k_lambda*l1*l2);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -k_lambda*(((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-3*((P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*((k_lambda)*(k_lambda))*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*((k_lambda)*(k_lambda))*l1*l2 - 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*((l1)*(l1)) + (1.0/2.0)*k_lambda*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/2.0*k_lambda*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + 6*k_lambda*(T_lambda_Pr + k_lambda*(P - P_r))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) + 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - k_lambda*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) + (-2*T_lambda_Pr - 2*k_lambda*(P - P_r))*(-3*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 3*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*k_lambda*l1*l2) + (T_lambda_Pr + k_lambda*(P - P_r))*(-6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*k_lambda*l1*l2) + (T_lambda_Pr + k_lambda*(P - P_r))*(-6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((k_lambda)*(k_lambda))*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*l1*l2 - 2*k_lambda*((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + (3.0/2.0)*k_lambda*((l2)*(l2))*(2*T_lambda_Pr + 2*k_lambda*(P - P_r)) - 2*k_lambda*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + k_lambda*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(-6*P*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*k_lambda*l1*l2));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -H_t*k_lambda/(T_lambda_Pr + k_lambda*(P - P_r)) + H_t*k_lambda*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -H_t*k_lambda/(T_lambda_Pr + k_lambda*(P - P_r)) + H_t*k_lambda*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(2*sqrt(T)*d1 + (1.0/3.0)*((T)*(T)*(T))*d4 + (1.0/2.0)*((T)*(T))*d3 + T*d0 - 2*sqrt(T_D_ref)*d1 - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref))*d4 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d3 - T_D_ref*d0 + d2/T_D_ref - d2/T)/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5;
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 1;
        break;
    case 3: /* S_TrPr */ 
        result += -T;
        break;
    case 4: /* k0 */ 
        result += -T*(log(T) - log(T_r)) + T - T_r;
        break;
    case 5: /* k1 */ 
        result += 2*sqrt(T) - T*(2/sqrt(T_r) - 2/sqrt(T)) - 2*sqrt(T_r);
        break;
    case 6: /* k2 */ 
        result += -T*((1.0/2.0)/((T_r)*(T_r)) - 1.0/2.0/((T)*(T))) + 1.0/T_r - 1/T;
        break;
    case 7: /* k3 */ 
        result += -T*((1.0/3.0)/((T_r)*(T_r)*(T_r)) - 1.0/3.0/((T)*(T)*(T))) + (1.0/2.0)/((T_r)*(T_r)) - 1.0/2.0/((T)*(T));
        break;
    case 8: /* V_TrPr */ 
        result += (1.0/3.0)*((P)*(P)*(P))*v2 + ((P)*(P))*(-P_r*v2 + (1.0/2.0)*v1) + P*(((P_r)*(P_r))*v2 - P_r*v1 + ((T)*(T))*v4 - 2*T*T_r*v4 + T*v3 + ((T_r)*(T_r))*v4 - T_r*v3 + 1) - 1.0/3.0*((P_r)*(P_r)*(P_r))*v2 - ((P_r)*(P_r))*(-P_r*v2 + (1.0/2.0)*v1) - P_r*(((P_r)*(P_r))*v2 - P_r*v1 + ((T)*(T))*v4 - 2*T*T_r*v4 + T*v3 + ((T_r)*(T_r))*v4 - T_r*v3 + 1);
        break;
    case 9: /* v1 */ 
        result += (1.0/2.0)*((P)*(P))*V_TrPr - P*P_r*V_TrPr + (1.0/2.0)*((P_r)*(P_r))*V_TrPr;
        break;
    case 10: /* v2 */ 
        result += (1.0/3.0)*((P)*(P)*(P))*V_TrPr - ((P)*(P))*P_r*V_TrPr + P*((P_r)*(P_r))*V_TrPr - 1.0/3.0*((P_r)*(P_r)*(P_r))*V_TrPr;
        break;
    case 11: /* v3 */ 
        result += P*(T*V_TrPr - T_r*V_TrPr) - P_r*(T*V_TrPr - T_r*V_TrPr);
        break;
    case 12: /* v4 */ 
        result += P*(((T)*(T))*V_TrPr - 2*T*T_r*V_TrPr + ((T_r)*(T_r))*V_TrPr) - P_r*(((T)*(T))*V_TrPr - 2*T*T_r*V_TrPr + ((T_r)*(T_r))*V_TrPr);
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (2.0/3.0)*((T)*(T)*(T))*l2 + ((T)*(T))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1) + T*(2*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1) - T*(((T)*(T))*l2 + T*(-4*P*k_lambda*l2 + 4*P_r*k_lambda*l2 + 2*l1) - k_lambda*(P - P_r)*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*l1)*log(T) + k_lambda*(P - P_r)*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(-4*P*k_lambda*l2 + 4*P_r*k_lambda*l2 + 2*l1)) - 2.0/3.0*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(2*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-2*((P)*(P))*((k_lambda)*(k_lambda))*l2 + 4*P*P_r*((k_lambda)*(k_lambda))*l2 + 2*P*k_lambda*l1 - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 - 2*P_r*k_lambda*l1 - k_lambda*(P - P_r)*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-2*T_lambda_Pr - 2*k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(-4*P*k_lambda*l2 + 4*P_r*k_lambda*l2 + 2*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(-4*P*k_lambda*l2 + 4*P_r*k_lambda*l2 - k_lambda*(P - P_r)*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l1 + l2*(2*T_lambda_Pr + 2*k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(-4*P*k_lambda*l2 + 4*P_r*k_lambda*l2 + 2*l1));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (1.0/2.0)*((T)*(T)*(T)*(T))*l2 + ((T)*(T)*(T))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + (2.0/3.0)*l1) + ((T)*(T))*(3*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 6*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1) + T*(-2*((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 6*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1 - 6*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l1 + 2*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1) - T*((2.0/3.0)*((T)*(T)*(T))*l2 + ((T)*(T))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + T*(6*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 12*P*P_r*((k_lambda)*(k_lambda))*l2 - 4*P*k_lambda*l1 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 4*P_r*k_lambda*l1) - k_lambda*(P - P_r)*(-2*P*k_lambda + 2*P_r*k_lambda)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + k_lambda*(P - P_r)*(-2*P*k_lambda + 2*P_r*k_lambda)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(6*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 12*P*P_r*((k_lambda)*(k_lambda))*l2 - 4*P*k_lambda*l1 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 4*P_r*k_lambda*l1) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1)) - 1.0/2.0*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-T_lambda_ref - k_lambda*(P - P_r))*(-2*((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 6*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1 - 6*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l1 + 2*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + (2.0/3.0)*l1) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(3*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 6*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (-T + T_lambda_Pr + k_lambda*(P - P_r))*(2*((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 6*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1 + 6*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 4*P*P_r*((k_lambda)*(k_lambda))*l1 - 2*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1 - k_lambda*(P - P_r)*(-2*P*k_lambda + 2*P_r*k_lambda)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-2*P*k_lambda + 2*P_r*k_lambda)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4.0/3.0*l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2.0/3.0*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-2*T_lambda_Pr - 2*k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 6*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + (2.0/3.0)*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 12*P*P_r*((k_lambda)*(k_lambda))*l2 - 4*P*k_lambda*l1 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 4*P_r*k_lambda*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 12*P*P_r*((k_lambda)*(k_lambda))*l2 - 4*P*k_lambda*l1 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 4*P_r*k_lambda*l1 - k_lambda*(P - P_r)*(-2*P*k_lambda + 2*P_r*k_lambda)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (2*T_lambda_Pr + 2*k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1)) + (-T_lambda_ref - k_lambda*(P - P_r))*(6*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 12*P*P_r*((k_lambda)*(k_lambda))*l2 - 4*P*k_lambda*l1 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 4*P_r*k_lambda*l1) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((T)*(T)*(T))*(-P*((l2)*(l2)) + P_r*((l2)*(l2))) + ((T)*(T))*(3*((P)*(P))*k_lambda*((l2)*(l2)) - 6*P*P_r*k_lambda*((l2)*(l2)) - 2*P*l1*l2 + 3*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 2*P_r*l1*l2) + T*(-3*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 9*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P)*(P))*k_lambda*l1*l2 - 9*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*P_r*k_lambda*l1*l2 - P*((l1)*(l1)) + 3*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P_r)*(P_r))*k_lambda*l1*l2 + P_r*((l1)*(l1))) - T*(((T)*(T))*(-3.0/2.0*P*((l2)*(l2)) + (3.0/2.0)*P_r*((l2)*(l2))) + T*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*(-2*P*l2 + 2*P_r*l2)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + k_lambda*(P - P_r)*(-2*P*l2 + 2*P_r*l2)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 1.0/3.0*((l2)*(l2))*(3*P - 3*P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-P + P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) + (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - (2*P - 2*P_r)*(T_lambda_ref + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (-T_lambda_ref - k_lambda*(P - P_r))*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*((l2)*(l2)) + (3.0/2.0)*P_r*((l2)*(l2)))) - 1.0/4.0*((l2)*(l2))*(4*P - 4*P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-P + P_r)*(-((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 - 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 - P*k_lambda*((l1)*(l1)) + ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 + P_r*k_lambda*((l1)*(l1))) - (2*P - 2*P_r)*(T_lambda_ref + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) - (3*P - 3*P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) + (-T_lambda_ref - k_lambda*(P - P_r))*(-3*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 9*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P)*(P))*k_lambda*l1*l2 - 9*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*P_r*k_lambda*l1*l2 - P*((l1)*(l1)) + 3*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P_r)*(P_r))*k_lambda*l1*l2 + P_r*((l1)*(l1))) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-P*((l2)*(l2)) + P_r*((l2)*(l2))) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(3*((P)*(P))*k_lambda*((l2)*(l2)) - 6*P*P_r*k_lambda*((l2)*(l2)) - 2*P*l1*l2 + 3*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 2*P_r*l1*l2);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (P - P_r)*(((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 9*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*((P)*(P))*k_lambda*l1*l2 + 9*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P*P_r*k_lambda*l1*l2 + P*((l1)*(l1)) - 3*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*((P_r)*(P_r))*k_lambda*l1*l2 - P_r*((l1)*(l1)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*(-2*P*l2 + 2*P_r*l2)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-2*P*l2 + 2*P_r*l2)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*((l2)*(l2))*(3*P - 3*P_r)*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*(3*P - 3*P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (-2*P + 2*P_r)*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + (-P + P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + (P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) + (2*P - 2*P_r)*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 3*(2*P - 2*P_r)*(T_lambda_Pr + k_lambda*(P - P_r))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - (2*P - 2*P_r)*(T_lambda_ref + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (-2*T_lambda_Pr - 2*k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*((l2)*(l2)) - 6*P*P_r*k_lambda*((l2)*(l2)) - 2*P*l1*l2 + 3*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 2*P_r*l1*l2) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*((l2)*(l2)) + (3.0/2.0)*P_r*((l2)*(l2))) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*((l2)*(l2)) + P_r*((l2)*(l2))) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 - k_lambda*(-P + P_r)*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - k_lambda*(P - P_r)*(-2*P*l2 + 2*P_r*l2)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l2)*(l2))*(2*P - 2*P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + (2*P - 2*P_r)*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (2*T_lambda_Pr + 2*k_lambda*(P - P_r))*(-3.0/2.0*P*((l2)*(l2)) + (3.0/2.0)*P_r*((l2)*(l2)))) + (-T_lambda_ref - k_lambda*(P - P_r))*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*((l2)*(l2)) + (3.0/2.0)*P_r*((l2)*(l2))));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(-P + P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + H_t*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(-P + P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + H_t*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2.0/3.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - 3*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-P*k_lambda*((l2)*(l2)) + P_r*k_lambda*((l2)*(l2)) + (2.0/3.0)*l1*l2) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*((3.0/2.0)*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P*k_lambda*l1*l2 + (3.0/2.0)*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P_r*k_lambda*l1*l2 + (1.0/2.0)*((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2)) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) - ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*((l2)*(l2)) + 3*P_r*k_lambda*((l2)*(l2)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + ((l2)*(l2))*(2*T_lambda_Pr + 2*k_lambda*(P - P_r))) + (2*T_lambda_Pr + 2*k_lambda*(P - P_r))*(-3.0/2.0*P*k_lambda*((l2)*(l2)) + (3.0/2.0)*P_r*k_lambda*((l2)*(l2)) + l1*l2) + (2*T_lambda_Pr + 2*k_lambda*(P - P_r))*(3*P*k_lambda*((l2)*(l2)) - 3*P_r*k_lambda*((l2)*(l2)) - 2*l1*l2));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t/(T_lambda_Pr + k_lambda*(P - P_r)) - H_t*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t/(T_lambda_Pr + k_lambda*(P - P_r)) - H_t*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - T*(-3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - ((l1)*(l1)) - ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (2*T_lambda_ref + 2*k_lambda*(P - P_r))*((3.0/2.0)*P*k_lambda*((l2)*(l2)) - 3.0/2.0*P_r*k_lambda*((l2)*(l2)) - l1*l2)) - ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(P*k_lambda*((l2)*(l2)) - P_r*k_lambda*((l2)*(l2)) - 2.0/3.0*l1*l2) + (2*T_lambda_ref + 2*k_lambda*(P - P_r))*(-3.0/2.0*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 3*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 2*P*k_lambda*l1*l2 - 3.0/2.0*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 2*P_r*k_lambda*l1*l2 - 1.0/2.0*((l1)*(l1)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - ((l1)*(l1)) - ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + (2*T_lambda_ref + 2*k_lambda*(P - P_r))*((3.0/2.0)*P*k_lambda*((l2)*(l2)) - 3.0/2.0*P_r*k_lambda*((l2)*(l2)) - l1*l2));
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        break;
    case 19: /* d0 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -T*(2*log(sqrt(T)) - 2*log(sqrt(T_D_ref))) + T - T_D_ref + (P - P_r)*(T - T_D_ref)/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*(2*log(sqrt(T_D)) - 2*log(sqrt(T_D_ref)) - (P - P_r)/d5);
            }
        }
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 2*sqrt(T) - T*(2/sqrt(T_D_ref) - 2/sqrt(T)) - 2*sqrt(T_D_ref) + (P - P_r)*(2*sqrt(T) - 2*sqrt(T_D_ref))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*(2/sqrt(T_D_ref) - 2/sqrt(T_D) - (P - P_r)/(sqrt(T_D)*d5));
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -T*((1.0/2.0)/((T_D_ref)*(T_D_ref)) - 1.0/2.0/((T)*(T))) + (P - P_r)*(1.0/T_D_ref - 1/T)/d5 + 1.0/T_D_ref - 1/T;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*((1.0/2.0)/((T_D_ref)*(T_D_ref)) - 1.0/2.0/((T_D)*(T_D)) - (P - P_r)/(((T_D)*(T_D))*d5));
            }
        }
        break;
    case 22: /* d3 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (1.0/2.0)*((T)*(T)) - T*(T - T_D_ref) - 1.0/2.0*((T_D_ref)*(T_D_ref)) + (P - P_r)*((1.0/2.0)*((T)*(T)) - 1.0/2.0*((T_D_ref)*(T_D_ref)))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*(T_D - T_D*(P - P_r)/d5 - T_D_ref);
            }
        }
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (1.0/3.0)*((T)*(T)*(T)) - T*((1.0/2.0)*((T)*(T)) - 1.0/2.0*((T_D_ref)*(T_D_ref))) - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref)) + (P - P_r)*((1.0/3.0)*((T)*(T)*(T)) - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref)))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*((1.0/2.0)*((T_D)*(T_D)) - ((T_D)*(T_D))*(P - P_r)/d5 - 1.0/2.0*((T_D_ref)*(T_D_ref)));
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(P - P_r)*(2*sqrt(T)*d1 + (1.0/3.0)*((T)*(T)*(T))*d4 + (1.0/2.0)*((T)*(T))*d3 + T*d0 - 2*sqrt(T_D_ref)*d1 - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref))*d4 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d3 - T_D_ref*d0 + d2/T_D_ref - d2/T)/((d5)*(d5));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (P - P_r)*(-T + T_D)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -1.0/2.0*((T_D)*(T_D))*d4 + T_D*(T_D*d4 + d3 + d0/T_D + d2/((T_D)*(T_D)*(T_D)) + d1/pow(T_D, 3.0/2.0)) - 1.0/2.0*((T_D_ref)*(T_D_ref))*d4 - T_D_ref*d3 + 2*d0*log(sqrt(T_D)) - 2*d0*log(sqrt(T_D_ref)) - d0 + (-T + T_D)*(T_D*(d4 - d0/((T_D)*(T_D)) - 3*d2/((T_D)*(T_D)*(T_D)*(T_D)) - 3.0/2.0*d1/pow(T_D, 5.0/2.0)) + d3 - (P - P_r)*(2*T_D*d4 + d3 - 2*d2/((T_D)*(T_D)*(T_D)) - 1.0/2.0*d1/pow(T_D, 3.0/2.0))/d5 + 2*d0/T_D + 4*d2/((T_D)*(T_D)*(T_D)) + (5.0/2.0)*d1/pow(T_D, 3.0/2.0)) - (P - P_r)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5 + (1.0/2.0)*d2/((T_D_ref)*(T_D_ref)) + 2*d1/sqrt(T_D_ref) - 3.0/2.0*d2/((T_D)*(T_D)) - 3*d1/sqrt(T_D);
            }
        }
        break;
    case 26: /* T_D_ref */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -T*(-T_D_ref*d4 - d3 - d0/T_D_ref - d2/((T_D_ref)*(T_D_ref)*(T_D_ref)) - d1/pow(T_D_ref, 3.0/2.0)) - ((T_D_ref)*(T_D_ref))*d4 - T_D_ref*d3 - d0 + (P - P_r)*(-((T_D_ref)*(T_D_ref))*d4 - T_D_ref*d3 - d0 - d2/((T_D_ref)*(T_D_ref)) - d1/sqrt(T_D_ref))/d5 - d2/((T_D_ref)*(T_D_ref)) - d1/sqrt(T_D_ref);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-T + T_D)*(-T_D_ref*d4 - d3 - d0/T_D_ref - d2/((T_D_ref)*(T_D_ref)*(T_D_ref)) - d1/pow(T_D_ref, 3.0/2.0));
            }
        }
        break;
    }
    return result;
}

static double Garnet_dparam_dgdt(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += k0/T_r + k2/((T_r)*(T_r)*(T_r)) + k3/((T_r)*(T_r)*(T_r)*(T_r)) + k1/pow(T_r, 3.0/2.0);
        result += 2*V_TrPr*v4*(-P + P_r);
        break;
    case 1: /* P_r */ 
        result += V_TrPr*(-2*T*v4 + 2*T_r*v4 - v3);
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*((3.0/2.0)*((T)*(T))*((l2)*(l2)) - T*(-6*P*k_lambda*((l2)*(l2)) + 6*P_r*k_lambda*((l2)*(l2)) + 3*T*((l2)*(l2)) + 4*l1*l2 - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + (1.0/2.0)*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + (1.0/2.0)*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (1.0/2.0)*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -H_t*k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -H_t*k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(((T)*(T))*d4 + T*d3 + d0 + d2/((T)*(T)) + d1/sqrt(T))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/d5;
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += -1;
        break;
    case 4: /* k0 */ 
        result += -log(T) + log(T_r);
        break;
    case 5: /* k1 */ 
        result += 2*(-1/sqrt(T_r) + pow(T, -1.0/2.0));
        break;
    case 6: /* k2 */ 
        result += (1.0/2.0)*(-1/((T_r)*(T_r)) + pow(T, -2));
        break;
    case 7: /* k3 */ 
        result += (1.0/3.0)*(-1/((T_r)*(T_r)*(T_r)) + pow(T, -3));
        break;
    case 8: /* V_TrPr */ 
        result += (P - P_r)*(2*T*v4 - 2*T_r*v4 + v3);
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 0.0;
        break;
    case 11: /* v3 */ 
        result += V_TrPr*(P - P_r);
        break;
    case 12: /* v4 */ 
        result += 2*V_TrPr*(P - P_r)*(T - T_r);
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1 + ((T)*(T))*l2 - 2*T*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + T*l2 + l1 - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*(T_lambda_ref + k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1 + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + l1 + l2*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_ref + k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1);
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 6*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1 - 6*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l1 + 2*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1 + (4.0/3.0)*((T)*(T)*(T))*l2 + ((T)*(T))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - 2*T*(3*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 6*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1 + ((T)*(T))*l2 + T*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + (2.0/3.0)*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 6*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P)*(P))*((k_lambda)*(k_lambda))*l1 - 6*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 4*P*P_r*((k_lambda)*(k_lambda))*l1 + 2*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1 - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + (4.0/3.0)*l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (2.0/3.0)*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 6*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1 + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1)) + ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1);
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -3*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 9*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P)*(P))*k_lambda*l1*l2 - 9*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*P_r*k_lambda*l1*l2 - P*((l1)*(l1)) + 3*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P_r)*(P_r))*k_lambda*l1*l2 + P_r*((l1)*(l1)) - 3.0/2.0*((T)*(T))*((l2)*(l2))*(P - P_r) - T*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 - 3*T*((l2)*(l2))*(P - P_r) + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T) - 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 1.0/2.0*((l2)*(l2))*(P - P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + (P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -3*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 9*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P)*(P))*k_lambda*l1*l2 - 9*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*P_r*k_lambda*l1*l2 - P*((l1)*(l1)) + 3*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P_r)*(P_r))*k_lambda*l1*l2 + P_r*((l1)*(l1)) - 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + (1.0/2.0)*((l2)*(l2))*(P - P_r)*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/2.0*((l2)*(l2))*(P - P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + l2*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + (P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - (P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + l2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1)) - (T_lambda_Pr + k_lambda*(P - P_r))*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(P - P_r)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(P - P_r)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - (T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*((l2)*(l2)) + 3*P_r*k_lambda*((l2)*(l2)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + 2*((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1);
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -1/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -1/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        break;
    case 19: /* d0 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -2*log(sqrt(T)) + 2*log(sqrt(T_D_ref)) + (P - P_r)/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -2*log(sqrt(T_D)) + 2*log(sqrt(T_D_ref)) + (P - P_r)/d5;
            }
        }
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -2/sqrt(T_D_ref) + 2/sqrt(T) + (P - P_r)/(sqrt(T)*d5);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -2/sqrt(T_D_ref) + 2/sqrt(T_D) + (P - P_r)/(sqrt(T_D)*d5);
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -1.0/2.0/((T_D_ref)*(T_D_ref)) + (1.0/2.0)/((T)*(T)) + (P - P_r)/(((T)*(T))*d5);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -1.0/2.0/((T_D_ref)*(T_D_ref)) + (1.0/2.0)/((T_D)*(T_D)) + (P - P_r)/(((T_D)*(T_D))*d5);
            }
        }
        break;
    case 22: /* d3 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -T + T*(P - P_r)/d5 + T_D_ref;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -T_D + T_D*(P - P_r)/d5 + T_D_ref;
            }
        }
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -1.0/2.0*((T)*(T)) + ((T)*(T))*(P - P_r)/d5 + (1.0/2.0)*((T_D_ref)*(T_D_ref));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -1.0/2.0*((T_D)*(T_D)) + ((T_D)*(T_D))*(P - P_r)/d5 + (1.0/2.0)*((T_D_ref)*(T_D_ref));
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(P - P_r)*(((T)*(T))*d4 + T*d3 + d0 + d2/((T)*(T)) + d1/sqrt(T))/((d5)*(d5));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -(P - P_r)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (1.0/2.0)*T_D*(-2*d4 + 2*d0/((T_D)*(T_D)) + 6*d2/((T_D)*(T_D)*(T_D)*(T_D)) + 3*d1/pow(T_D, 5.0/2.0)) - d3 + (1.0/2.0)*(P - P_r)*(4*T_D*d4 + 2*d3 - 4*d2/((T_D)*(T_D)*(T_D)) - d1/pow(T_D, 3.0/2.0))/d5 - 2*d0/T_D - 4*d2/((T_D)*(T_D)*(T_D)) - 5.0/2.0*d1/pow(T_D, 3.0/2.0);
            }
        }
        break;
    case 26: /* T_D_ref */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += T_D_ref*d4 + d3 + d0/T_D_ref + d2/((T_D_ref)*(T_D_ref)*(T_D_ref)) + d1/pow(T_D_ref, 3.0/2.0);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += T_D_ref*d4 + d3 + d0/T_D_ref + d2/((T_D_ref)*(T_D_ref)*(T_D_ref)) + d1/pow(T_D_ref, 3.0/2.0);
            }
        }
        break;
    }
    return result;
}

static double Garnet_dparam_dgdp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += V_TrPr*(-2*T*v4 + 2*T_r*v4 - v3);
        break;
    case 1: /* P_r */ 
        result += V_TrPr*(-2*P*v2 + 2*P_r*v2 - v1);
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(-3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 - 3*((T)*(T))*((l2)*(l2)) - 2*T*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - T*(-6*T*((l2)*(l2)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 3*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r))) - ((l1)*(l1)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(-6*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 12*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P*k_lambda*l1*l2 - 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P_r*k_lambda*l1*l2 - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((l1)*(l1)) - ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-18*P*k_lambda*((l2)*(l2)) + 18*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*l1*l2 + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 4*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*((k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*((k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += 0.0;
        break;
    case 5: /* k1 */ 
        result += 0.0;
        break;
    case 6: /* k2 */ 
        result += 0.0;
        break;
    case 7: /* k3 */ 
        result += 0.0;
        break;
    case 8: /* V_TrPr */ 
        result += ((P)*(P))*v2 - P*(2*P_r*v2 - v1) + ((P_r)*(P_r))*v2 - P_r*v1 + ((T)*(T))*v4 - 2*T*T_r*v4 + T*v3 + ((T_r)*(T_r))*v4 - T_r*v3 + 1;
        break;
    case 9: /* v1 */ 
        result += V_TrPr*(P - P_r);
        break;
    case 10: /* v2 */ 
        result += V_TrPr*(((P)*(P)) - 2*P*P_r + ((P_r)*(P_r)));
        break;
    case 11: /* v3 */ 
        result += V_TrPr*(T - T_r);
        break;
    case 12: /* v4 */ 
        result += V_TrPr*(((T)*(T)) - 2*T*T_r + ((T_r)*(T_r)));
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*k_lambda*(-((T)*(T))*l2 - T*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1) + T*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*T*l2 - k_lambda*l2*(P - P_r)*log(T) + k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + l1 - l2*(T_lambda_ref + k_lambda*(P - P_r)) + (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r))) - k_lambda*(((P)*(P))*k_lambda*l2 - 2*P*P_r*k_lambda*l2 - P*l1 + ((P_r)*(P_r))*k_lambda*l2 + P_r*l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-2*((P)*(P))*((k_lambda)*(k_lambda))*l2 + 4*P*P_r*((k_lambda)*(k_lambda))*l2 + 2*P*k_lambda*l1 - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 - 2*P_r*k_lambda*l1 - 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + l1 + l2*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*(T_lambda_ref + k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + l1) + 2*(-T + T_lambda_Pr + k_lambda*(P - P_r))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + k_lambda*l2*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + l1 + l2*(T_lambda_ref + k_lambda*(P - P_r)) + (T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2 - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-2*((T)*(T)*(T))*l2 - 2*((T)*(T))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - 2*T*k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + T*(3*((T)*(T))*l2 + 4*T*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T) - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1)) - 2*((k_lambda)*(k_lambda))*(-((P)*(P)*(P))*k_lambda*l2 + 3*((P)*(P))*P_r*k_lambda*l2 + ((P)*(P))*l1 - 3*P*((P_r)*(P_r))*k_lambda*l2 - 2*P*P_r*l1 + ((P_r)*(P_r)*(P_r))*k_lambda*l2 + ((P_r)*(P_r))*l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(2*((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 6*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1 + 6*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 + 4*P*P_r*((k_lambda)*(k_lambda))*l1 - 2*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*l2 - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1 + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - 4.0/3.0*l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2.0/3.0*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 6*P*P_r*((k_lambda)*(k_lambda))*l2 - 2*P*k_lambda*l1 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 2*P_r*k_lambda*l1 + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1)) - ((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(12*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 24*P*P_r*((k_lambda)*(k_lambda))*l2 - 8*P*k_lambda*l1 + 12*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 8*P_r*k_lambda*l1 - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*P*k_lambda*l2 - 3*P_r*k_lambda*l2 - ((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - l1 - l2*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1)));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - ((T)*(T)*(T))*((l2)*(l2)) - 2*((T)*(T))*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - T*(9*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 18*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*k_lambda*l1*l2 + 9*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (1.0/2.0)*T*(6*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 12*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*k_lambda*l1*l2 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P_r*k_lambda*l1*l2 + 3*((T)*(T))*((l2)*(l2)) + 8*T*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T) - 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) + 8*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 8*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - 6*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 2*((l1)*(l1)) - ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 8*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r))) + 3*k_lambda*((l2)*(l2))*(P - P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*l2*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - k_lambda*(-3*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 9*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P)*(P))*k_lambda*l1*l2 - 9*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*P_r*k_lambda*l1*l2 - P*((l1)*(l1)) + 3*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*((P_r)*(P_r))*k_lambda*l1*l2 + P_r*((l1)*(l1))) + 2*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (T_lambda_ref + k_lambda*(P - P_r))*(9*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 18*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*k_lambda*l1*l2 + 9*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P_r*k_lambda*l1*l2 + ((l1)*(l1)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((P)*(P)*(P))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 3*((P)*(P))*P_r*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P)*(P))*((k_lambda)*(k_lambda))*l1*l2 + 3*P*((P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*P_r*((k_lambda)*(k_lambda))*l1*l2 + P*k_lambda*((l1)*(l1)) - ((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda)*(k_lambda))*((l2)*(l2)) - 2*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l1*l2 - P_r*k_lambda*((l1)*(l1)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + (1.0/2.0)*k_lambda*(P - P_r)*(6*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 12*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*k_lambda*l1*l2 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P_r*k_lambda*l1*l2 + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r))) - 1.0/2.0*k_lambda*(-6*((P)*(P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 18*((P)*(P))*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*((P)*(P))*k_lambda*l1*l2 - 18*P*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 16*P*P_r*k_lambda*l1*l2 - 2*P*((l1)*(l1)) + 6*((P_r)*(P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*((P_r)*(P_r))*k_lambda*l1*l2 + 2*P_r*((l1)*(l1)) - 4*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l2)*(l2))*(P - P_r)*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*(P - P_r)*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*l2*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 2*l2*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 4*l2*(T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + 2*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*(P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) - 2*(P - P_r)*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + l2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1)) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)))) - 2.0/3.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/3.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 1.0/2.0*l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 1.0/2.0*l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + l2*(T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1)) - (T_lambda_ref + k_lambda*(P - P_r))*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 + ((l1)*(l1))) + (1.0/2.0)*(-T + T_lambda_Pr + k_lambda*(P - P_r))*(18*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 36*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 16*P*k_lambda*l1*l2 + 18*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 16*P_r*k_lambda*l1*l2 - 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) - 8*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 8*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + 8*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r)) + 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*l2*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + 2*k_lambda*(P - P_r)*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + 6*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((l1)*(l1)) + ((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 8*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-9*P*k_lambda*((l2)*(l2)) + 9*P_r*k_lambda*((l2)*(l2)) + 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*((l2)*(l2))*(P - P_r) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(-2*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(-2*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(6*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 12*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P*k_lambda*l1*l2 + 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P_r*k_lambda*l1*l2 + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((l1)*(l1)) - 1.0/2.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + (1.0/2.0)*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(-3*P*k_lambda*((l2)*(l2)) + 3*P_r*k_lambda*((l2)*(l2)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + 2*((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r))) + (T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(3*P*k_lambda*((l2)*(l2)) - 3*P_r*k_lambda*((l2)*(l2)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l1*l2 - (T_lambda_Pr + k_lambda*(P - P_r))*(2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + ((l2)*(l2)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*k_lambda*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*k_lambda*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += T*k_lambda*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -k_lambda*(3*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 6*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P*k_lambda*l1*l2 + 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P_r*k_lambda*l1*l2 - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + ((l1)*(l1)) + ((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r));
            }
        }
        break;
    case 19: /* d0 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (T - T_D_ref)/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (T - T_D)/d5;
            }
        }
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 2*(sqrt(T) - sqrt(T_D_ref))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (T - T_D)/(sqrt(T_D)*d5);
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (1.0/T_D_ref - 1/T)/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (T - T_D)/(((T_D)*(T_D))*d5);
            }
        }
        break;
    case 22: /* d3 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (1.0/2.0)*(((T)*(T)) - ((T_D_ref)*(T_D_ref)))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += T_D*(T - T_D)/d5;
            }
        }
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (1.0/3.0)*(((T)*(T)*(T)) - ((T_D_ref)*(T_D_ref)*(T_D_ref)))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += ((T_D)*(T_D))*(T - T_D)/d5;
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(2*sqrt(T)*d1 + (1.0/3.0)*((T)*(T)*(T))*d4 + (1.0/2.0)*((T)*(T))*d3 + T*d0 - 2*sqrt(T_D_ref)*d1 - 1.0/3.0*((T_D_ref)*(T_D_ref)*(T_D_ref))*d4 - 1.0/2.0*((T_D_ref)*(T_D_ref))*d3 - T_D_ref*d0 + d2/T_D_ref - d2/T)/((d5)*(d5));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -(T - T_D)*(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (-((T_D)*(T_D))*d4 - T_D*d3 - d0 + (1.0/2.0)*(T - T_D)*(4*T_D*d4 + 2*d3 - 4*d2/((T_D)*(T_D)*(T_D)) - d1/pow(T_D, 3.0/2.0)) - d2/((T_D)*(T_D)) - d1/sqrt(T_D))/d5;
            }
        }
        break;
    case 26: /* T_D_ref */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(((T_D_ref)*(T_D_ref))*d4 + T_D_ref*d3 + d0 + d2/((T_D_ref)*(T_D_ref)) + d1/sqrt(T_D_ref))/d5;
            }
        }
        break;
    }
    return result;
}

static double Garnet_dparam_d2gdt2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += -2*V_TrPr*v4;
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) - T*(3*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T))) - 4*l1*l2 + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(2*T*d4 + d3 - 2*d2/((T)*(T)*(T)) - 1.0/2.0*d1/pow(T, 3.0/2.0))/d5;
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += -1/T;
        break;
    case 5: /* k1 */ 
        result += -1/pow(T, 3.0/2.0);
        break;
    case 6: /* k2 */ 
        result += -1/((T)*(T)*(T));
        break;
    case 7: /* k3 */ 
        result += -1/((T)*(T)*(T)*(T));
        break;
    case 8: /* V_TrPr */ 
        result += 2*v4*(P - P_r);
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 0.0;
        break;
    case 11: /* v3 */ 
        result += 0.0;
        break;
    case 12: /* v4 */ 
        result += 2*V_TrPr*(P - P_r);
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*(2*P*k_lambda*l2 - 2*P_r*k_lambda*l2 - T*(l2 + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T))) - l1 + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T);
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*(-3*((P)*(P))*((k_lambda)*(k_lambda))*l2 + 6*P*P_r*((k_lambda)*(k_lambda))*l2 + 2*P*k_lambda*l1 - 3*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 - 2*P_r*k_lambda*l1 + ((T)*(T))*l2 - T*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*T*l2 + l1 - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T))) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T);
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -6*((P)*(P))*k_lambda*((l2)*(l2)) + 12*P*P_r*k_lambda*((l2)*(l2)) + 4*P*l1*l2 - 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) - 4*P_r*l1*l2 + T*(3*P*((l2)*(l2)) - 3*P_r*((l2)*(l2)) + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T))) - 4*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + 2*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T;
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        result += 0.0;
        break;
    case 17: /* T_lambda_ref */ 
        result += 0.0;
        break;
    case 18: /* H_t */ 
        result += 0.0;
        break;
    case 19: /* d0 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -1/T;
            }
        }
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(1 + (1.0/2.0)*(P - P_r)/d5)/pow(T, 3.0/2.0);
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(1 + 2*(P - P_r)/d5)/((T)*(T)*(T));
            }
        }
        break;
    case 22: /* d3 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -1 + (P - P_r)/d5;
            }
        }
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += T*(-1 + 2*(P - P_r)/d5);
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(P - P_r)*(2*T*d4 + d3 - 2*d2/((T)*(T)*(T)) - 1.0/2.0*d1/pow(T, 3.0/2.0))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        result += 0.0;
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static double Garnet_dparam_d2gdtdp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += -2*V_TrPr*v4;
        break;
    case 1: /* P_r */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*T*l2*(3*l2 - k_lambda*l2*(P - P_r)/T + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l1*l2 - 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + 3*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(18*P*k_lambda*((l2)*(l2)) - 18*P_r*k_lambda*((l2)*(l2)) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 12*l1*l2 - 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + 4*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += 0.0;
        break;
    case 5: /* k1 */ 
        result += 0.0;
        break;
    case 6: /* k2 */ 
        result += 0.0;
        break;
    case 7: /* k3 */ 
        result += 0.0;
        break;
    case 8: /* V_TrPr */ 
        result += 2*T*v4 - 2*T_r*v4 + v3;
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 0.0;
        break;
    case 11: /* v3 */ 
        result += V_TrPr;
        break;
    case 12: /* v4 */ 
        result += 2*V_TrPr*(T - T_r);
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*k_lambda*(T*(2*l2 - k_lambda*l2*(P - P_r)/T + (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) - k_lambda*l2*(P - P_r)*log(T) + k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - l2*(T_lambda_ref + k_lambda*(P - P_r)) + (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*k_lambda*(2*P*k_lambda*l2 - 2*P_r*k_lambda*l2 - k_lambda*l2*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - l1 - l2*(T_lambda_ref + k_lambda*(P - P_r)) - (T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2 - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) + (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-6*((P)*(P))*((k_lambda)*(k_lambda))*l2 + 12*P*P_r*((k_lambda)*(k_lambda))*l2 + 4*P*k_lambda*l1 - 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 - 4*P_r*k_lambda*l1 - 3*((T)*(T))*l2 + 2*T*(-6*P*k_lambda*l2 + 6*P_r*k_lambda*l2 + 3*T*l2 + 2*l1 + ((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/T - 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T) - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-12*((P)*(P))*((k_lambda)*(k_lambda))*l2 + 24*P*P_r*((k_lambda)*(k_lambda))*l2 + 8*P*k_lambda*l1 - 12*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 - 8*P_r*k_lambda*l1 + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 4*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*P*k_lambda*l2 - 3*P_r*k_lambda*l2 - ((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - l1 - l2*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -6*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 12*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 4*P*k_lambda*l1*l2 - 6*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 4*P_r*k_lambda*l1*l2 - 3.0/2.0*((T)*(T))*((l2)*(l2)) + T*(-12*P*k_lambda*((l2)*(l2)) + 12*P_r*k_lambda*((l2)*(l2)) + 3*T*((l2)*(l2)) + 4*l1*l2 + 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/T - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T) + 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T) - 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*l2*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 1.0/2.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -9*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 18*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) + 8*P*k_lambda*l1*l2 - 9*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 8*P_r*k_lambda*l1*l2 + 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r)) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*l2*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) - k_lambda*(P - P_r)*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - ((l1)*(l1)) - 1.0/2.0*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 1.0/2.0*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + (T_lambda_Pr + k_lambda*(P - P_r))*(-9*P*k_lambda*((l2)*(l2)) + 9*P_r*k_lambda*((l2)*(l2)) + 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*((l2)*(l2))*(P - P_r) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(-2*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += H_t*(-2*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-3*P*k_lambda*((l2)*(l2)) + 3*P_r*k_lambda*((l2)*(l2)) - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + (T_lambda_Pr + k_lambda*(P - P_r))*(2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + ((l2)*(l2)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*H_t*k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*H_t*k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 19: /* d0 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 1.0/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += 1.0/d5;
            }
        }
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 1/(sqrt(T)*d5);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += 1/(sqrt(T_D)*d5);
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 1/(((T)*(T))*d5);
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += 1/(((T_D)*(T_D))*d5);
            }
        }
        break;
    case 22: /* d3 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += T/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += T_D/d5;
            }
        }
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += ((T)*(T))/d5;
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += ((T_D)*(T_D))/d5;
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(((T)*(T))*d4 + T*d3 + d0 + d2/((T)*(T)) + d1/sqrt(T))/((d5)*(d5));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += -(((T_D)*(T_D))*d4 + T_D*d3 + d0 + d2/((T_D)*(T_D)) + d1/sqrt(T_D))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D))) {
                 result += (2*T_D*d4 + d3 - 2*d2/((T_D)*(T_D)*(T_D)) - 1.0/2.0*d1/pow(T_D, 3.0/2.0))/d5;
            }
        }
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static double Garnet_dparam_d2gdp2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += -2*V_TrPr*v2;
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda)*(k_lambda))*(6*T*((l2)*(l2)) + T*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 11*((l2)*(l2)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))) + 2*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda)*(k_lambda))*(-54*P*k_lambda*((l2)*(l2)) + 54*P_r*k_lambda*((l2)*(l2)) + 6*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 6*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 24*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 9*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 36*l1*l2 + 6*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 12*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 6*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - (-T + T_lambda_Pr + k_lambda*(P - P_r))*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 6*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 12*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += 0.0;
        break;
    case 5: /* k1 */ 
        result += 0.0;
        break;
    case 6: /* k2 */ 
        result += 0.0;
        break;
    case 7: /* k3 */ 
        result += 0.0;
        break;
    case 8: /* V_TrPr */ 
        result += 2*P*v2 - 2*P_r*v2 + v1;
        break;
    case 9: /* v1 */ 
        result += V_TrPr;
        break;
    case 10: /* v2 */ 
        result += 2*V_TrPr*(P - P_r);
        break;
    case 11: /* v3 */ 
        result += 0.0;
        break;
    case 12: /* v4 */ 
        result += 0.0;
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-2*P*k_lambda*l2 + 2*P_r*k_lambda*l2 + 2*T*l2 + T*(2*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*l2*log(T) + 2*l2*log(T_lambda_ref + k_lambda*(P - P_r)) - 3*l2 - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r))) + l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-4*P*k_lambda*l2 + 4*P_r*k_lambda*l2 + 2*k_lambda*l2*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l1 + 2*l2*(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2 - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) - (-T + T_lambda_Pr + k_lambda*(P - P_r))*(2*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l2*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l2*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(3*((T)*(T))*l2 + 2*T*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - T*(-9*P*k_lambda*l2 + 9*P_r*k_lambda*l2 + 6*T*l2 + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*k_lambda*l2*(P - P_r)*log(T) + 4*k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 3*l1 - 2*l2*(T_lambda_ref + k_lambda*(P - P_r)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r))) + k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(12*((P)*(P))*((k_lambda)*(k_lambda))*l2 - 24*P*P_r*((k_lambda)*(k_lambda))*l2 - 8*P*k_lambda*l1 + 12*((P_r)*(P_r))*((k_lambda)*(k_lambda))*l2 + 8*P_r*k_lambda*l1 - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + l2*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + l2*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*P*k_lambda*l2 - 3*P_r*k_lambda*l2 - ((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - l1 - l2*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-6*P*k_lambda*l2 + 6*P_r*k_lambda*l2 - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*l2*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l1 + 2*l2*(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2 - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(18*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 36*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 16*P*k_lambda*l1*l2 + 18*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 16*P_r*k_lambda*l1*l2 + 6*((T)*(T))*((l2)*(l2)) + 2*T*l2*(-9*P*k_lambda*l2 + 9*P_r*k_lambda*l2 + 4*l1) - T*(-12*P*k_lambda*((l2)*(l2)) + 12*P_r*k_lambda*((l2)*(l2)) + 12*T*((l2)*(l2)) + 6*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) + 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T) + 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*((l2)*(l2))*(P - P_r) - 20*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 8*l1*l2 - 4*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + 8*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) - 2*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) - 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r))) + 6*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*l2*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + 2*((l1)*(l1)) - 2*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-9*P*k_lambda*l2 + 9*P_r*k_lambda*l2 + 4*l1) + 8*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(24*((P)*(P))*((k_lambda)*(k_lambda))*((l2)*(l2)) - 48*P*P_r*((k_lambda)*(k_lambda))*((l2)*(l2)) - 24*P*k_lambda*l1*l2 + 24*((P_r)*(P_r))*((k_lambda)*(k_lambda))*((l2)*(l2)) + 24*P_r*k_lambda*l1*l2 - 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))*log(T_lambda_ref + k_lambda*(P - P_r)) - 8*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 8*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + 8*k_lambda*((l2)*(l2))*(P - P_r)*(T_lambda_ref + k_lambda*(P - P_r)) + 16*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 16*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*l2*(3*((P)*(P))*k_lambda*l2 - 6*P*P_r*k_lambda*l2 - 2*P*l1 + 3*((P_r)*(P_r))*k_lambda*l2 + 2*P_r*l1) + 2*k_lambda*(P - P_r)*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - k_lambda*(P - P_r)*(-18*P*k_lambda*((l2)*(l2)) + 18*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*l1*l2 + 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 4*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + 8*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 12*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(6*((P)*(P))*k_lambda*((l2)*(l2)) - 12*P*P_r*k_lambda*((l2)*(l2)) - 4*P*l1*l2 + 6*((P_r)*(P_r))*k_lambda*((l2)*(l2)) + 4*P_r*l1*l2 + 2*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2))*(P - P_r)*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(P - P_r)*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - (P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) + 4*((l1)*(l1)) + 2*((l2)*(l2))*((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((l2)*(l2))*((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 8*l2*(T_lambda_ref + k_lambda*(P - P_r))*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(6*P*k_lambda*((l2)*(l2)) - 6*P_r*k_lambda*((l2)*(l2)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*l1*l2 - ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-9*P*k_lambda*((l2)*(l2)) + 9*P_r*k_lambda*((l2)*(l2)) + 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*((l2)*(l2))*(P - P_r) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*l1*l2 + ((l2)*(l2))*(T_lambda_Pr + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-36*P*k_lambda*((l2)*(l2)) + 36*P_r*k_lambda*((l2)*(l2)) - 6*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) + 10*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 14*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*((l2)*(l2))*(P - P_r) + 20*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 36*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 13*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*l1*l2 + 4*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 8*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 5*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 10*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((l2)*(l2)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 8*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r))) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*k_lambda*(3*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2 + 2*(-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*k_lambda*(3*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2 + 2*(-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(24*P*k_lambda*((l2)*(l2)) - 24*P_r*k_lambda*((l2)*(l2)) - 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 16*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 9*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 16*l1*l2 - 2*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + 4*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + ((l2)*(l2)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 2*(-T + T_lambda_Pr + k_lambda*(P - P_r))*(2*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - (k_lambda*((l2)*(l2))*(P - P_r) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 8*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(3*P*k_lambda*((l2)*(l2)) - 3*P_r*k_lambda*((l2)*(l2)) - T*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))) - 2*l1*l2 + l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-6*P*k_lambda*((l2)*(l2)) + 6*P_r*k_lambda*((l2)*(l2)) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 4*l1*l2 + ((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 19: /* d0 */ 
        result += 0.0;
        break;
    case 20: /* d1 */ 
        result += 0.0;
        break;
    case 21: /* d2 */ 
        result += 0.0;
        break;
    case 22: /* d3 */ 
        result += 0.0;
        break;
    case 23: /* d4 */ 
        result += 0.0;
        break;
    case 24: /* d5 */ 
        result += 0.0;
        break;
    case 25: /* T_D */ 
        result += 0.0;
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static double Garnet_dparam_d3gdt3(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-3*((l2)*(l2)) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)) + 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(P*k_lambda*l2 - P_r*k_lambda*l2 + 2*k_lambda*l2*(P - P_r) - l1)/((T)*(T)));
            }
        }
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(2*d4 + 6*d2/((T)*(T)*(T)*(T)) + (3.0/4.0)*d1/pow(T, 5.0/2.0))/d5;
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += pow(T, -2);
        break;
    case 5: /* k1 */ 
        result += (3.0/2.0)/pow(T, 5.0/2.0);
        break;
    case 6: /* k2 */ 
        result += 3/((T)*(T)*(T)*(T));
        break;
    case 7: /* k3 */ 
        result += 4/pow(T, 5);
        break;
    case 8: /* V_TrPr */ 
        result += 0.0;
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 0.0;
        break;
    case 11: /* v3 */ 
        result += 0.0;
        break;
    case 12: /* v4 */ 
        result += 0.0;
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*(l2 + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*(3*P*k_lambda*l2 - 3*P_r*k_lambda*l2 - 2*T*(l2 + ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)*(T))) - l1 + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 3*P*((l2)*(l2)) - 3*P_r*((l2)*(l2)) + 6*k_lambda*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)) - 3*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T)) - 2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(P*k_lambda*l2 - P_r*k_lambda*l2 + 2*k_lambda*l2*(P - P_r) - l1)/((T)*(T));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        result += 0.0;
        break;
    case 17: /* T_lambda_ref */ 
        result += 0.0;
        break;
    case 18: /* H_t */ 
        result += 0.0;
        break;
    case 19: /* d0 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += pow(T, -2);
            }
        }
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += (3.0/4.0)*(2 + (P - P_r)/d5)/pow(T, 5.0/2.0);
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 3*(1 + 2*(P - P_r)/d5)/((T)*(T)*(T)*(T));
            }
        }
        break;
    case 22: /* d3 */ 
        result += 0.0;
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -1 + 2*(P - P_r)/d5;
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(P - P_r)*(2*d4 + 6*d2/((T)*(T)*(T)*(T)) + (3.0/4.0)*d1/pow(T, 5.0/2.0))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        result += 0.0;
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static double Garnet_dparam_d3gdt2dp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*l2*(3*l2 - 2*k_lambda*l2*(P - P_r)/T + 4*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + (2*P*k_lambda*l2 - 2*P_r*k_lambda*l2 + k_lambda*l2*(P - P_r) - 2*l1)/T);
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += 0.0;
        break;
    case 5: /* k1 */ 
        result += 0.0;
        break;
    case 6: /* k2 */ 
        result += 0.0;
        break;
    case 7: /* k3 */ 
        result += 0.0;
        break;
    case 8: /* V_TrPr */ 
        result += 2*v4;
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 0.0;
        break;
    case 11: /* v3 */ 
        result += 0.0;
        break;
    case 12: /* v4 */ 
        result += 2*V_TrPr;
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*k_lambda*(2*l2 - 2*k_lambda*l2*(P - P_r)/T + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + (P*k_lambda*l2 - P_r*k_lambda*l2 + k_lambda*l2*(P - P_r) - l1)/T);
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*k_lambda*(-6*P*k_lambda*l2 + 6*P_r*k_lambda*l2 + T*(3*l2 - ((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T)*(T)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T))) + 2*l1 + 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/T - 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T);
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -12*P*k_lambda*((l2)*(l2)) + 12*P_r*k_lambda*((l2)*(l2)) + T*(3*((l2)*(l2)) - 2*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T)*(T)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T)*(T)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T)*(T))) + 4*l1*l2 + 4*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/T - 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/T;
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        result += 0.0;
        break;
    case 17: /* T_lambda_ref */ 
        result += 0.0;
        break;
    case 18: /* H_t */ 
        result += 0.0;
        break;
    case 19: /* d0 */ 
        result += 0.0;
        break;
    case 20: /* d1 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -1.0/2.0/(pow(T, 3.0/2.0)*d5);
            }
        }
        break;
    case 21: /* d2 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -2/(((T)*(T)*(T))*d5);
            }
        }
        break;
    case 22: /* d3 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 1.0/d5;
            }
        }
        break;
    case 23: /* d4 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += 2*T/d5;
            }
        }
        break;
    case 24: /* d5 */ 
        if(((d0) != 0.0) || ((d1) != 0.0) || ((d2) != 0.0) || ((d3) != 0.0) || ((d4) != 0.0) || ((d5) != 0.0) || ((T_D) != 0.0) || ((T_D_ref) != 0.0)) {
            if((T > (T_D_ref)) && (T <= (T_D))) {
                 result += -(2*T*d4 + d3 - 2*d2/((T)*(T)*(T)) - 1.0/2.0*d1/pow(T, 3.0/2.0))/((d5)*(d5));
            }
        }
        break;
    case 25: /* T_D */ 
        result += 0.0;
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static double Garnet_dparam_d3gdtdp2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda)*(k_lambda))*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 11*((l2)*(l2)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda)*(k_lambda))*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -6*H_t*((k_lambda)*(k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -6*H_t*((k_lambda)*(k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += 0.0;
        break;
    case 5: /* k1 */ 
        result += 0.0;
        break;
    case 6: /* k2 */ 
        result += 0.0;
        break;
    case 7: /* k3 */ 
        result += 0.0;
        break;
    case 8: /* V_TrPr */ 
        result += 0.0;
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 0.0;
        break;
    case 11: /* v3 */ 
        result += 0.0;
        break;
    case 12: /* v4 */ 
        result += 0.0;
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(2*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*l2*log(T) + 2*l2*log(T_lambda_ref + k_lambda*(P - P_r)) - 3*l2 - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(2*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l2*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l2*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(3*P*k_lambda*l2 - 3*P_r*k_lambda*l2 - 2*T*(3*l2 - 2*k_lambda*l2*(P - P_r)/T + (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 4*k_lambda*l2*(P - P_r)*log(T) - 4*k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - l1 + 2*l2*(T_lambda_ref + k_lambda*(P - P_r)) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-6*P*k_lambda*l2 + 6*P_r*k_lambda*l2 - 2*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) + 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*l2*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 4*k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 8*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*l1 + 2*l2*(T_lambda_ref + k_lambda*(P - P_r)) + 2*(T_lambda_Pr + k_lambda*(P - P_r))*(-((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2 - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) - 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-6*P*k_lambda*((l2)*(l2)) + 6*P_r*k_lambda*((l2)*(l2)) - 2*T*l2*(6*l2 - 5*k_lambda*l2*(P - P_r)/T + 4*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/T) - 6*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) - 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T) - 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - k_lambda*((l2)*(l2))*(P - P_r) + 20*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 4*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 8*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T) + 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += k_lambda*(-36*P*k_lambda*((l2)*(l2)) + 36*P_r*k_lambda*((l2)*(l2)) - 6*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) + 10*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) - 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 14*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 10*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*((l2)*(l2))*(P - P_r) + 20*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 36*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 13*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*l1*l2 + 4*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) - 8*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 2*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) - 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 8*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 2*(T_lambda_Pr + k_lambda*(P - P_r))*(3*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 5*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 10*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((l2)*(l2)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) + 8*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*k_lambda*(3*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 2)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*H_t*k_lambda*(3*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 2)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(2*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - (k_lambda*((l2)*(l2))*(P - P_r) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((l2)*(l2)) + 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda))*(-k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + ((l2)*(l2)) + 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*((k_lambda)*(k_lambda))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 19: /* d0 */ 
        result += 0.0;
        break;
    case 20: /* d1 */ 
        result += 0.0;
        break;
    case 21: /* d2 */ 
        result += 0.0;
        break;
    case 22: /* d3 */ 
        result += 0.0;
        break;
    case 23: /* d4 */ 
        result += 0.0;
        break;
    case 24: /* d5 */ 
        result += 0.0;
        break;
    case 25: /* T_D */ 
        result += 0.0;
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static double Garnet_dparam_d3gdp3(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*((k_lambda)*(k_lambda)*(k_lambda)*(k_lambda))*(T*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 12*((l2)*(l2)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))))/(T_lambda_ref + k_lambda*(P - P_r)) + 3*((l2)*(l2)));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda)*(k_lambda)*(k_lambda))*(-12*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 12*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 12*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 4*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) - 12*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) + 24*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 24*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(6*k_lambda*((l2)*(l2))*(P - P_r)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*((l2)*(l2))*(P - P_r)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 8*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 12*((l2)*(l2))/(T_lambda_ref + k_lambda*(P - P_r)) + 12*((l2)*(l2))/(T_lambda_Pr + k_lambda*(P - P_r)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 4*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 6*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 24*H_t*((k_lambda)*(k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 24*H_t*((k_lambda)*(k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 2: /* H_TrPr */ 
        result += 0.0;
        break;
    case 3: /* S_TrPr */ 
        result += 0.0;
        break;
    case 4: /* k0 */ 
        result += 0.0;
        break;
    case 5: /* k1 */ 
        result += 0.0;
        break;
    case 6: /* k2 */ 
        result += 0.0;
        break;
    case 7: /* k3 */ 
        result += 0.0;
        break;
    case 8: /* V_TrPr */ 
        result += 2*v2;
        break;
    case 9: /* v1 */ 
        result += 0.0;
        break;
    case 10: /* v2 */ 
        result += 2*V_TrPr;
        break;
    case 11: /* v3 */ 
        result += 0.0;
        break;
    case 12: /* v4 */ 
        result += 0.0;
        break;
    case 13: /* l1 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += -2*((k_lambda)*(k_lambda)*(k_lambda))*(T*(3*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*l2 - 3*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)))/(T_lambda_ref + k_lambda*(P - P_r)) + 2*l2);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda)*(k_lambda))*(-6*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 6*l2*log(T_lambda_Pr + k_lambda*(P - P_r)) - 6*l2*log(T_lambda_ref + k_lambda*(P - P_r)) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(3*k_lambda*l2*(P - P_r)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 3*k_lambda*l2*(P - P_r)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*l2/(T_lambda_ref + k_lambda*(P - P_r)) + 6*l2/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 14: /* l2 */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda)*(k_lambda))*(6*P*k_lambda*l2 - 6*P_r*k_lambda*l2 - 6*T*l2 - T*(-3*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 12*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 6*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*l2*log(T) + 6*l2*log(T_lambda_ref + k_lambda*(P - P_r)) - 11*l2 - 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r))) - 2*l1);
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 2*((k_lambda)*(k_lambda)*(k_lambda))*(18*P*k_lambda*l2 - 18*P_r*k_lambda*l2 + 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) - 12*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 9*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 12*k_lambda*l2*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 12*k_lambda*l2*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 12*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 24*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 6*l1 - 6*l2*(T_lambda_ref + k_lambda*(P - P_r)) - 6*(T_lambda_Pr + k_lambda*(P - P_r))*(-((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - l2 - (-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) + (-T + T_lambda_Pr + k_lambda*(P - P_r))*(-3*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 2*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 12*k_lambda*l2*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 12*k_lambda*l2*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*l2*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*l2*log(T_lambda_ref + k_lambda*(P - P_r)) - 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r))) + 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 6*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)));
            }
        }
        break;
    case 15: /* k_lambda */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(18*P*k_lambda*((l2)*(l2)) - 18*P_r*k_lambda*((l2)*(l2)) - 18*T*((l2)*(l2)) + T*(12*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 16*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 42*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 42*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 14*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 18*((l2)*(l2))*log(T) - 18*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) + 33*((l2)*(l2)) + 36*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 9*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r)))) + 6*k_lambda*((l2)*(l2))*(P - P_r) - 12*l1*l2 - 6*l2*(-9*P*k_lambda*l2 + 9*P_r*k_lambda*l2 + 4*l1) + 12*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 6*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda))*(162*P*k_lambda*((l2)*(l2)) - 162*P_r*k_lambda*((l2)*(l2)) + 18*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)) - 30*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)) + 18*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 42*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 18*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 36*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_Pr + k_lambda*(P - P_r)) + 36*k_lambda*((l2)*(l2))*(P - P_r)*log(T_lambda_ref + k_lambda*(P - P_r)) - 6*k_lambda*((l2)*(l2))*(P - P_r) - 72*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 132*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 6*k_lambda*(P - P_r)*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + k_lambda*(P - P_r)*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - 18*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 48*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 72*l1*l2 - 18*((l2)*(l2))*(T_lambda_ref + k_lambda*(P - P_r)) + 24*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + l1) + 6*l2*(-3*P*k_lambda*l2 + 3*P_r*k_lambda*l2 + 2*l1) + 36*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_Pr + k_lambda*(P - P_r)) - 36*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*log(T_lambda_ref + k_lambda*(P - P_r)) - 6*(T_lambda_Pr + k_lambda*(P - P_r))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 6*(T_lambda_Pr + k_lambda*(P - P_r))*(3*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 3*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 5*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 10*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 5*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*((l2)*(l2)) + 4*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) - (-T + T_lambda_Pr + k_lambda*(P - P_r))*(12*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*((k_lambda)*(k_lambda))*((l2)*(l2))*((P - P_r)*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 16*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 4*((k_lambda)*(k_lambda))*l2*((P - P_r)*(P - P_r))*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 6*((k_lambda)*(k_lambda))*((P - P_r)*(P - P_r))*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 42*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 36*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) - 42*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 30*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 14*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 6*k_lambda*(P - P_r)*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))))/(T_lambda_Pr + k_lambda*(P - P_r)) + 8*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 18*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) - 18*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) + 36*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - 36*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 9*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 9*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)))) + 18*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_ref + k_lambda*(P - P_r)) - 36*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda))*(-4*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 3 - 3*(-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*H_t*((k_lambda)*(k_lambda))*(-4*k_lambda*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 4*k_lambda*(P - P_r)*(-T + T_lambda_Pr + k_lambda*(P - P_r))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 3 - 3*(-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 16: /* T_lambda_Pr */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += ((k_lambda)*(k_lambda)*(k_lambda))*(6*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) - 18*k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 42*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) + 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - 26*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 6*((l2)*(l2))*log(T_lambda_Pr + k_lambda*(P - P_r)) + 6*((l2)*(l2))*log(T_lambda_ref + k_lambda*(P - P_r)) + 6*((l2)*(l2)) - 12*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) + 36*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) - 3*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + 6*(-T + T_lambda_Pr + k_lambda*(P - P_r))*(-2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - ((l2)*(l2)) + (k_lambda*((l2)*(l2))*(P - P_r) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r)) + ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))))/(T_lambda_Pr + k_lambda*(P - P_r)) + 6*(k_lambda*((l2)*(l2))*(P - P_r) + 4*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_Pr + k_lambda*(P - P_r)) + 3*k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1) - 2*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/(T_lambda_Pr + k_lambda*(P - P_r)))/(T_lambda_Pr + k_lambda*(P - P_r)) + 21*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 24*H_t*((k_lambda)*(k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 24*H_t*((k_lambda)*(k_lambda)*(k_lambda))*(-1 + (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 17: /* T_lambda_ref */ 
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_ref + k_lambda*(P - P_r))) && (P > (P_r)) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*T*((k_lambda)*(k_lambda)*(k_lambda))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))))/(T_lambda_ref + k_lambda*(P - P_r));
            }
        }
        if(((l1) != 0.0) || ((l2) != 0.0) || ((k_lambda) != 0.0) || ((T_lambda_Pr) != 0.0) || ((T_lambda_ref) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*((k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_ref + k_lambda*(P - P_r)))*(k_lambda*((l2)*(l2))*(P - P_r)/(T_lambda_ref + k_lambda*(P - P_r)) + 2*k_lambda*l2*(P - P_r)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) + k_lambda*(P - P_r)*((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))) - ((l2)*(l2)) - 2*l2*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)/(T_lambda_ref + k_lambda*(P - P_r)) - ((-P*k_lambda*l2 + P_r*k_lambda*l2 + l1)*(-P*k_lambda*l2 + P_r*k_lambda*l2 + l1))/((T_lambda_ref + k_lambda*(P - P_r))*(T_lambda_ref + k_lambda*(P - P_r))));
            }
        }
        break;
    case 18: /* H_t */ 
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r))) && (T <= (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*((k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        if(((H_t) != 0.0)) {
            if((T > (T_lambda_Pr + k_lambda*(P - P_r)))) {
                 result += 6*((k_lambda)*(k_lambda)*(k_lambda))*(1 - (-T + T_lambda_Pr + k_lambda*(P - P_r))/(T_lambda_Pr + k_lambda*(P - P_r)))/((T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r))*(T_lambda_Pr + k_lambda*(P - P_r)));
            }
        }
        break;
    case 19: /* d0 */ 
        result += 0.0;
        break;
    case 20: /* d1 */ 
        result += 0.0;
        break;
    case 21: /* d2 */ 
        result += 0.0;
        break;
    case 22: /* d3 */ 
        result += 0.0;
        break;
    case 23: /* d4 */ 
        result += 0.0;
        break;
    case 24: /* d5 */ 
        result += 0.0;
        break;
    case 25: /* T_D */ 
        result += 0.0;
        break;
    case 26: /* T_D_ref */ 
        result += 0.0;
        break;
    }
    return result;
}

static int Garnet_get_param_number(void) {
    return 27;
}

static const char *paramNames[27] = { "T_r", "P_r", "H_TrPr", "S_TrPr", "k0", "k1", "k2", "k3", "V_TrPr", "v1", "v2", "v3", "v4", "l1", "l2", "k_lambda", "T_lambda_Pr", "T_lambda_ref", "H_t", "d0", "d1", "d2", "d3", "d4", "d5", "T_D", "T_D_ref"  };

static const char *paramUnits[27] = { "K", "bar", "J", "J/K", "J/K-m", "J-K^(1/2)-m", "J-K/m", "J-K^2", "J/bar-m", "1/bar", "1/bar^2", "1/K", "1/K^2", "(J/m)^(1/2)-K", "(J/m)^(1/2)/K^2", "K/bar", "K", "K", "J/m", "J/K-m", "J/K^(1/2)-m", "J-K/m", "J/K^2-m", "J/K^3-m", "bar", "K", "K"  };

static const char **Garnet_get_param_names(void) {
    return paramNames;
}

static const char **Garnet_get_param_units(void) {
    return paramUnits;
}

static void Garnet_get_param_values(double **values) {
    (*values)[0] = T_r;
    (*values)[1] = P_r;
    (*values)[2] = H_TrPr;
    (*values)[3] = S_TrPr;
    (*values)[4] = k0;
    (*values)[5] = k1;
    (*values)[6] = k2;
    (*values)[7] = k3;
    (*values)[8] = V_TrPr;
    (*values)[9] = v1;
    (*values)[10] = v2;
    (*values)[11] = v3;
    (*values)[12] = v4;
    (*values)[13] = l1;
    (*values)[14] = l2;
    (*values)[15] = k_lambda;
    (*values)[16] = T_lambda_Pr;
    (*values)[17] = T_lambda_ref;
    (*values)[18] = H_t;
    (*values)[19] = d0;
    (*values)[20] = d1;
    (*values)[21] = d2;
    (*values)[22] = d3;
    (*values)[23] = d4;
    (*values)[24] = d5;
    (*values)[25] = T_D;
    (*values)[26] = T_D_ref;
}

static int Garnet_set_param_values(double *values) {
    T_r= values[0];
    P_r= values[1];
    H_TrPr= values[2];
    S_TrPr= values[3];
    k0= values[4];
    k1= values[5];
    k2= values[6];
    k3= values[7];
    V_TrPr= values[8];
    v1= values[9];
    v2= values[10];
    v3= values[11];
    v4= values[12];
    l1= values[13];
    l2= values[14];
    k_lambda= values[15];
    T_lambda_Pr= values[16];
    T_lambda_ref= values[17];
    H_t= values[18];
    d0= values[19];
    d1= values[20];
    d2= values[21];
    d3= values[22];
    d4= values[23];
    d5= values[24];
    T_D= values[25];
    T_D_ref= values[26];
    return 1;
}

static double Garnet_get_param_value(int index) {
    double result = 0.0;
    switch (index) {
    case 0:
        result = T_r;
        break;
    case 1:
        result = P_r;
        break;
    case 2:
        result = H_TrPr;
        break;
    case 3:
        result = S_TrPr;
        break;
    case 4:
        result = k0;
        break;
    case 5:
        result = k1;
        break;
    case 6:
        result = k2;
        break;
    case 7:
        result = k3;
        break;
    case 8:
        result = V_TrPr;
        break;
    case 9:
        result = v1;
        break;
    case 10:
        result = v2;
        break;
    case 11:
        result = v3;
        break;
    case 12:
        result = v4;
        break;
    case 13:
        result = l1;
        break;
    case 14:
        result = l2;
        break;
    case 15:
        result = k_lambda;
        break;
    case 16:
        result = T_lambda_Pr;
        break;
    case 17:
        result = T_lambda_ref;
        break;
    case 18:
        result = H_t;
        break;
    case 19:
        result = d0;
        break;
    case 20:
        result = d1;
        break;
    case 21:
        result = d2;
        break;
    case 22:
        result = d3;
        break;
    case 23:
        result = d4;
        break;
    case 24:
        result = d5;
        break;
    case 25:
        result = T_D;
        break;
    case 26:
        result = T_D_ref;
        break;
     default:
         break;
    }
    return result;
}

static int Garnet_set_param_value(int index, double value) {
    int result = 1;
    switch (index) {
    case 0:
        T_r = value;
        break;
    case 1:
        P_r = value;
        break;
    case 2:
        H_TrPr = value;
        break;
    case 3:
        S_TrPr = value;
        break;
    case 4:
        k0 = value;
        break;
    case 5:
        k1 = value;
        break;
    case 6:
        k2 = value;
        break;
    case 7:
        k3 = value;
        break;
    case 8:
        V_TrPr = value;
        break;
    case 9:
        v1 = value;
        break;
    case 10:
        v2 = value;
        break;
    case 11:
        v3 = value;
        break;
    case 12:
        v4 = value;
        break;
    case 13:
        l1 = value;
        break;
    case 14:
        l2 = value;
        break;
    case 15:
        k_lambda = value;
        break;
    case 16:
        T_lambda_Pr = value;
        break;
    case 17:
        T_lambda_ref = value;
        break;
    case 18:
        H_t = value;
        break;
    case 19:
        d0 = value;
        break;
    case 20:
        d1 = value;
        break;
    case 21:
        d2 = value;
        break;
    case 22:
        d3 = value;
        break;
    case 23:
        d4 = value;
        break;
    case 24:
        d5 = value;
        break;
    case 25:
        T_D = value;
        break;
    case 26:
        T_D_ref = value;
        break;
     default:
         break;
    }
    return result;
}



const char *Majorite_Garnet_calib_identifier(void) {
    return identifier;
}

const char *Majorite_Garnet_calib_name(void) {
    return "Majorite";
}

const char *Majorite_Garnet_calib_formula(void) {
    return "Mg3MgSiSi3O12";
}

const double Majorite_Garnet_calib_mw(void) {
    return 401.5548;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,12.0,0.0,0.0,0.0,
        1.0,0.0,3.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0
    };

const double *Majorite_Garnet_calib_elements(void) {
    return elmformula;
}

double Majorite_Garnet_calib_g(double T, double P) {
    return Garnet_g(T, P);
}

double Majorite_Garnet_calib_dgdt(double T, double P) {
    return Garnet_dgdt(T, P);
}

double Majorite_Garnet_calib_dgdp(double T, double P) {
    return Garnet_dgdp(T, P);
}

double Majorite_Garnet_calib_d2gdt2(double T, double P) {
    return Garnet_d2gdt2(T, P);
}

double Majorite_Garnet_calib_d2gdtdp(double T, double P) {
    return Garnet_d2gdtdp(T, P);
}

double Majorite_Garnet_calib_d2gdp2(double T, double P) {
    return Garnet_d2gdp2(T, P);
}

double Majorite_Garnet_calib_d3gdt3(double T, double P) {
    return Garnet_d3gdt3(T, P);
}

double Majorite_Garnet_calib_d3gdt2dp(double T, double P) {
    return Garnet_d3gdt2dp(T, P);
}

double Majorite_Garnet_calib_d3gdtdp2(double T, double P) {
    return Garnet_d3gdtdp2(T, P);
}

double Majorite_Garnet_calib_d3gdp3(double T, double P) {
    return Garnet_d3gdp3(T, P);
}

double Majorite_Garnet_calib_s(double T, double P) {
    return Garnet_s(T, P);
}

double Majorite_Garnet_calib_v(double T, double P) {
    return Garnet_v(T, P);
}

double Majorite_Garnet_calib_cv(double T, double P) {
    return Garnet_cv(T, P);
}

double Majorite_Garnet_calib_cp(double T, double P) {
    return Garnet_cp(T, P);
}

double Majorite_Garnet_calib_dcpdt(double T, double P) {
    return Garnet_dcpdt(T, P);
}

double Majorite_Garnet_calib_alpha(double T, double P) {
    return Garnet_alpha(T, P);
}

double Majorite_Garnet_calib_beta(double T, double P) {
    return Garnet_beta(T, P);
}

double Majorite_Garnet_calib_K(double T, double P) {
    return Garnet_K(T, P);
}

double Majorite_Garnet_calib_Kp(double T, double P) {
    return Garnet_Kp(T, P);
}

int Majorite_Garnet_get_param_number(void) {
    return Garnet_get_param_number();
}

const char **Majorite_Garnet_get_param_names(void) {
    return Garnet_get_param_names();
}

const char **Majorite_Garnet_get_param_units(void) {
    return Garnet_get_param_units();
}

void Majorite_Garnet_get_param_values(double **values) {
    Garnet_get_param_values(values);
}

int Majorite_Garnet_set_param_values(double *values) {
    return Garnet_set_param_values(values);
}

double Majorite_Garnet_get_param_value(int index) {
    return Garnet_get_param_value(index);
}

int Majorite_Garnet_set_param_value(int index, double value) {
    return Garnet_set_param_value(index, value);
}

double Majorite_Garnet_dparam_g(double T, double P, int index) {
    return Garnet_dparam_g(T, P, index);
}

double Majorite_Garnet_dparam_dgdt(double T, double P, int index) {
    return Garnet_dparam_dgdt(T, P, index);
}

double Majorite_Garnet_dparam_dgdp(double T, double P, int index) {
    return Garnet_dparam_dgdp(T, P, index);
}

double Majorite_Garnet_dparam_d2gdt2(double T, double P, int index) {
    return Garnet_dparam_d2gdt2(T, P, index);
}

double Majorite_Garnet_dparam_d2gdtdp(double T, double P, int index) {
    return Garnet_dparam_d2gdtdp(T, P, index);
}

double Majorite_Garnet_dparam_d2gdp2(double T, double P, int index) {
    return Garnet_dparam_d2gdp2(T, P, index);
}

double Majorite_Garnet_dparam_d3gdt3(double T, double P, int index) {
    return Garnet_dparam_d3gdt3(T, P, index);
}

double Majorite_Garnet_dparam_d3gdt2dp(double T, double P, int index) {
    return Garnet_dparam_d3gdt2dp(T, P, index);
}

double Majorite_Garnet_dparam_d3gdtdp2(double T, double P, int index) {
    return Garnet_dparam_d3gdtdp2(T, P, index);
}

double Majorite_Garnet_dparam_d3gdp3(double T, double P, int index) {
    return Garnet_dparam_d3gdp3(T, P, index);
}

