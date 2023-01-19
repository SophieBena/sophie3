
#include <stdlib.h>
#include <stdio.h>

#include "Almandine_Garnet_calib.h"
#include "Grossular_Garnet_calib.h"
#include "Pyrope_Garnet_calib.h"


typedef struct _endmembers {
  const char *(*name) (void);
  const char *(*formula) (void);
  const double (*mw) (void);
  const double *(*elements) (void);
  double (*mu0) (double t, double p);
  double (*dmu0dT) (double t, double p);
  double (*dmu0dP) (double t, double p);
  double (*d2mu0dT2) (double t, double p);
  double (*d2mu0dTdP) (double t, double p);
  double (*d2mu0dP2) (double t, double p);
  double (*d3mu0dT3) (double t, double p);
  double (*d3mu0dT2dP) (double t, double p);
  double (*d3mu0dTdP2) (double t, double p);
  double (*d3mu0dP3) (double t, double p);
} Endmembers;

static Endmembers endmember[] = {
  {
    Almandine_Garnet_calib_name,
    Almandine_Garnet_calib_formula,
    Almandine_Garnet_calib_mw,
    Almandine_Garnet_calib_elements,
    Almandine_Garnet_calib_g,
    Almandine_Garnet_calib_dgdt,
    Almandine_Garnet_calib_dgdp,
    Almandine_Garnet_calib_d2gdt2,
    Almandine_Garnet_calib_d2gdtdp,
    Almandine_Garnet_calib_d2gdp2,
    Almandine_Garnet_calib_d3gdt3,
    Almandine_Garnet_calib_d3gdt2dp,
    Almandine_Garnet_calib_d3gdtdp2,
    Almandine_Garnet_calib_d3gdp3
  },
  {
    Grossular_Garnet_calib_name,
    Grossular_Garnet_calib_formula,
    Grossular_Garnet_calib_mw,
    Grossular_Garnet_calib_elements,
    Grossular_Garnet_calib_g,
    Grossular_Garnet_calib_dgdt,
    Grossular_Garnet_calib_dgdp,
    Grossular_Garnet_calib_d2gdt2,
    Grossular_Garnet_calib_d2gdtdp,
    Grossular_Garnet_calib_d2gdp2,
    Grossular_Garnet_calib_d3gdt3,
    Grossular_Garnet_calib_d3gdt2dp,
    Grossular_Garnet_calib_d3gdtdp2,
    Grossular_Garnet_calib_d3gdp3
  },
  {
    Pyrope_Garnet_calib_name,
    Pyrope_Garnet_calib_formula,
    Pyrope_Garnet_calib_mw,
    Pyrope_Garnet_calib_elements,
    Pyrope_Garnet_calib_g,
    Pyrope_Garnet_calib_dgdt,
    Pyrope_Garnet_calib_dgdp,
    Pyrope_Garnet_calib_d2gdt2,
    Pyrope_Garnet_calib_d2gdtdp,
    Pyrope_Garnet_calib_d2gdp2,
    Pyrope_Garnet_calib_d3gdt3,
    Pyrope_Garnet_calib_d3gdt2dp,
    Pyrope_Garnet_calib_d3gdtdp2,
    Pyrope_Garnet_calib_d3gdp3
  },

};
static int nc = (sizeof endmember / sizeof(struct _endmembers));

static const double R=8.3143;


static char *identifier = "Thu Jul 18 15:30:16 2019";
static double T_r = 298.15;
static double P_r = 1.0;
static double Wh12 = 11470.0;
static double Ws12 = 5.08;
static double Wv12 = 0.13;
static double dWh12 = -8850.0;
static double dWs12 = 0.0;
static double dWv12 = -0.04000000000000001;
static double Wh13 = 1975.0;
static double Ws13 = 0.0;
static double Wv13 = 0.034999999999999996;
static double dWh13 = 1745.0;
static double dWs13 = 0.0;
static double dWv13 = 0.024999999999999998;
static double Wh23 = 45380.0;
static double Ws23 = 18.79;
static double Wv23 = 0.1;
static double dWh23 = -23820.0;
static double dWs23 = 0.0;
static double dWv23 = 0.0;
static double Wh123 = 0;
static double Ws123 = 0;
static double Wv123 = 0;


#include "Berman_91_Garnet_Solution_calc.h"
#include "Berman_91_Garnet_Solution_calib.h"

const char *Garnet_Berman_91_Garnet_Solution_calib_identifier(void) {
    return identifier;
}

const char *Garnet_Berman_91_Garnet_Solution_calib_name(void) {
    return "Garnet";
}

char *Garnet_Berman_91_Garnet_Solution_calib_formula(double T, double P, double n[3]) {
    double sum, elm[3];
    const double *end0 = (*endmember[0].elements)();
    const double *end1 = (*endmember[1].elements)();
    const double *end2 = (*endmember[2].elements)();
    int i;
    const char *fmt = "(Fe%5.3fCa%5.3fMg%5.3f)Al2Si3O12";
    char *result = (char *) malloc(33*sizeof(char));
    for (i=0, sum=0.0; i<nc; i++) sum += n[i];
    if (sum == 0.0) return result;
    elm[0] = end0[26]*n[0]/sum + end1[26]*n[1]/sum + end2[26]*n[2]/sum;
    elm[1] = end0[20]*n[0]/sum + end1[20]*n[1]/sum + end2[20]*n[2]/sum;
    elm[2] = end0[12]*n[0]/sum + end1[12]*n[1]/sum + end2[12]*n[2]/sum;
    sprintf(result, fmt, elm[0], elm[1], elm[2]);
    return result;

}

double *Garnet_Berman_91_Garnet_Solution_calib_conv_elm_to_moles(double *e) {
    double *n = (double *) malloc(3*sizeof(double));
    n[0]=e[26];
    n[1]=e[20];
    n[2]=e[12];
    return n;

}

int Garnet_Berman_91_Garnet_Solution_calib_test_moles(double *n) {
    int result = 1;
    result &= (n[0] > 0.0);
    result &= (n[1] > 0.0);
    result &= (n[2] > 0.0);
    return result;

}

const char *Garnet_Berman_91_Garnet_Solution_calib_endmember_name(int index) {
    return (*endmember[index].name)();
}

const char *Garnet_Berman_91_Garnet_Solution_calib_endmember_formula(int index) {
    return (*endmember[index].formula)();
}

const double Garnet_Berman_91_Garnet_Solution_calib_endmember_mw(int index) {
    return (*endmember[index].mw)();
}

const double *Garnet_Berman_91_Garnet_Solution_calib_endmember_elements(int index) {
    return (*endmember[index].elements)();
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_mu0(int index, double t, double p) {
    return (*endmember[index].mu0)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_dmu0dT(int index, double t, double p) {
    return (*endmember[index].dmu0dT)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_dmu0dP(int index, double t, double p) {
    return (*endmember[index].dmu0dP)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d2mu0dT2(int index, double t, double p) {
    return (*endmember[index].d2mu0dT2)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d2mu0dTdP(int index, double t, double p) {
    return (*endmember[index].d2mu0dTdP)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d2mu0dP2(int index, double t, double p) {
    return (*endmember[index].d2mu0dP2)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d3mu0dT3(int index, double t, double p) {
    return (*endmember[index].d3mu0dT3)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d3mu0dT2dP(int index, double t, double p) {
    return (*endmember[index].d3mu0dT2dP)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d3mu0dTdP2(int index, double t, double p) {
    return (*endmember[index].d3mu0dTdP2)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_endmember_d3mu0dP3(int index, double t, double p) {
    return (*endmember[index].d3mu0dP3)(t, p);
}

double Garnet_Berman_91_Garnet_Solution_calib_g(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_g(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_dgdn(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_dgdn(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d2gdn2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d2gdn2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d3gdn3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d3gdn3(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_dgdt(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_dgdt(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d2gdndt(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d2gdndt(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d3gdn2dt(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d3gdn2dt(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdn3dt(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdn3dt(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_dgdp(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_dgdp(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d2gdndp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d2gdndp(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d3gdn2dp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d3gdn2dp(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdn3dp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdn3dp(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d2gdt2(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d2gdt2(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d3gdndt2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d3gdndt2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdn2dt2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdn2dt2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn3dt2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn3dt2(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d2gdtdp(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d2gdtdp(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d3gdndtdp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d3gdndtdp(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdn2dtdp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdn2dtdp(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn3dtdp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn3dtdp(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d2gdp2(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d2gdp2(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d3gdndp2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d3gdndp2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdn2dp2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdn2dp2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn3dp2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn3dp2(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d3gdt3(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d3gdt3(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdndt3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdndt3(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn2dt3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn2dt3(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d6gdn3dt3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d6gdn3dt3(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d3gdt2dp(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d3gdt2dp(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdndt2dp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdndt2dp(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn2dt2dp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn2dt2dp(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d6gdn3dt2dp(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d6gdn3dt2dp(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d3gdtdp2(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d3gdtdp2(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdndtdp2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdndtdp2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn2dtdp2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn2dtdp2(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d6gdn3dtdp2(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d6gdn3dtdp2(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_d3gdp3(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_d3gdp3(T, P, n);
}

void Garnet_Berman_91_Garnet_Solution_calib_d4gdndp3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d4gdndp3(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d5gdn2dp3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d5gdn2dp3(T, P, n, result);
}

void Garnet_Berman_91_Garnet_Solution_calib_d6gdn3dp3(double T, double P, double n[3], double result[3]) {
    Berman_91_Garnet_Solution_d6gdn3dp3(T, P, n, result);
}

double Garnet_Berman_91_Garnet_Solution_calib_s(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_s(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_v(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_v(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_cv(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_cv(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_cp(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_cp(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_dcpdt(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_dcpdt(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_alpha(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_alpha(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_beta(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_beta(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_K(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_K(T, P, n);
}

double Garnet_Berman_91_Garnet_Solution_calib_Kp(double T, double P, double n[3]) {
    return Berman_91_Garnet_Solution_Kp(T, P, n);
}

int Garnet_Berman_91_Garnet_Solution_get_param_number(void) {
    return Berman_91_Garnet_Solution_get_param_number();
}

const char **Garnet_Berman_91_Garnet_Solution_get_param_names(void) {
    return Berman_91_Garnet_Solution_get_param_names();
}

const char **Garnet_Berman_91_Garnet_Solution_get_param_units(void) {
    return Berman_91_Garnet_Solution_get_param_units();
}

void Garnet_Berman_91_Garnet_Solution_get_param_values(double **values) {
    Berman_91_Garnet_Solution_get_param_values(values);
}

int Garnet_Berman_91_Garnet_Solution_set_param_values(double *values) {
    return Berman_91_Garnet_Solution_set_param_values(values);
}

double Garnet_Berman_91_Garnet_Solution_get_param_value(int index) {
    return Berman_91_Garnet_Solution_get_param_value(index);
}

int Garnet_Berman_91_Garnet_Solution_set_param_value(int index, double value) {
    return Berman_91_Garnet_Solution_set_param_value(index, value);
}

double Garnet_Berman_91_Garnet_Solution_dparam_g(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_g(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_dgdt(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_dgdt(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_dgdp(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_dgdp(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d2gdt2(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d2gdt2(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d2gdtdp(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d2gdtdp(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d2gdp2(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d2gdp2(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d3gdt3(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d3gdt3(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d3gdt2dp(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d3gdt2dp(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d3gdtdp2(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d3gdtdp2(T, P, n, index);
}

double Garnet_Berman_91_Garnet_Solution_dparam_d3gdp3(double T, double P, double n[3], int index) {
    return Berman_91_Garnet_Solution_dparam_d3gdp3(T, P, n, index);
}

void Garnet_Berman_91_Garnet_Solution_dparam_dgdn(double T, double P, double n[3], int index, double result[3]) {
    Berman_91_Garnet_Solution_dparam_dgdn(T, P, n, index, result);
}

