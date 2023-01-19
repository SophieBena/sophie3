
const char *Garnet_Garnet_Solution_ideal_calib_identifier(void);
const char *Garnet_Garnet_Solution_ideal_calib_name(void);
char *Garnet_Garnet_Solution_ideal_calib_formula(double T, double P, double n[3]);
double *Garnet_Garnet_Solution_ideal_calib_conv_elm_to_moles(double *e);
int Garnet_Garnet_Solution_ideal_calib_test_moles(double *n);

const char *Garnet_Garnet_Solution_ideal_calib_endmember_name(int index);
const char *Garnet_Garnet_Solution_ideal_calib_endmember_formula(int index);
const double Garnet_Garnet_Solution_ideal_calib_endmember_mw(int index);
const double *Garnet_Garnet_Solution_ideal_calib_endmember_elements(int index);
double Garnet_Garnet_Solution_ideal_calib_endmember_mu0(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_dmu0dT(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_dmu0dP(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d2mu0dT2(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d2mu0dTdP(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d2mu0dP2(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d3mu0dT3(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d3mu0dT2dP(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d3mu0dTdP2(int index, double t, double p);
double Garnet_Garnet_Solution_ideal_calib_endmember_d3mu0dP3(int index, double t, double p);

double Garnet_Garnet_Solution_ideal_calib_g(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_dgdt(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_dgdp(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d2gdt2(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d2gdtdp(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d2gdp2(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d3gdt3(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d3gdt2dp(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d3gdtdp2(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_d3gdp3(double T, double P, double n[3]);

void Garnet_Garnet_Solution_ideal_calib_dgdn(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d2gdndt(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d2gdndp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d3gdndt2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d3gdndtdp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d3gdndp2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdndt3(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdndt2dp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdndtdp2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdndp3(double T, double P, double n[3], double result[3]);

void Garnet_Garnet_Solution_ideal_calib_d2gdn2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d3gdn2dt(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d3gdn2dp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdn2dt2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdn2dtdp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdn2dp2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn2dt3(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn2dt2dp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn2dtdp2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn2dp3(double T, double P, double n[3], double result[3]);

void Garnet_Garnet_Solution_ideal_calib_d3gdn3(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdn3dt(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d4gdn3dp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn3dt2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn3dtdp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d5gdn3dp2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d6gdn3dt3(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d6gdn3dt2dp(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d6gdn3dtdp2(double T, double P, double n[3], double result[3]);
void Garnet_Garnet_Solution_ideal_calib_d6gdn3dp3(double T, double P, double n[3], double result[3]);

double Garnet_Garnet_Solution_ideal_calib_s(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_v(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_cv(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_cp(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_dcpdt(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_alpha(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_beta(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_K(double T, double P, double n[3]);
double Garnet_Garnet_Solution_ideal_calib_Kp(double T, double P, double n[3]);


int Garnet_Garnet_Solution_ideal_get_param_number(void);
const char **Garnet_Garnet_Solution_ideal_get_param_names(void);
const char **Garnet_Garnet_Solution_ideal_get_param_units(void);
void Garnet_Garnet_Solution_ideal_get_param_values(double **values);
int Garnet_Garnet_Solution_ideal_set_param_values(double *values);
double Garnet_Garnet_Solution_ideal_get_param_value(int index);
int Garnet_Garnet_Solution_ideal_set_param_value(int index, double value);

double Garnet_Garnet_Solution_ideal_dparam_g(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_dgdt(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_dgdp(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d2gdt2(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d2gdtdp(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d2gdp2(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d3gdt3(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d3gdt2dp(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d3gdtdp2(double T, double P, double n[3], int index);
double Garnet_Garnet_Solution_ideal_dparam_d3gdp3(double T, double P, double n[3], int index);

void Garnet_Garnet_Solution_ideal_dparam_dgdn(double T, double P, double n[3], int index, double result[3]);

