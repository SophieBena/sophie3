
const char *Majorite_Garnet_calib_identifier(void);
const char *Majorite_Garnet_calib_name(void);
const char *Majorite_Garnet_calib_formula(void);
const double Majorite_Garnet_calib_mw(void);
const double *Majorite_Garnet_calib_elements(void);

double Majorite_Garnet_calib_g(double T, double P);
double Majorite_Garnet_calib_dgdt(double T, double P);
double Majorite_Garnet_calib_dgdp(double T, double P);
double Majorite_Garnet_calib_d2gdt2(double T, double P);
double Majorite_Garnet_calib_d2gdtdp(double T, double P);
double Majorite_Garnet_calib_d2gdp2(double T, double P);
double Majorite_Garnet_calib_d3gdt3(double T, double P);
double Majorite_Garnet_calib_d3gdt2dp(double T, double P);
double Majorite_Garnet_calib_d3gdtdp2(double T, double P);
double Majorite_Garnet_calib_d3gdp3(double T, double P);

double Majorite_Garnet_calib_s(double T, double P);
double Majorite_Garnet_calib_v(double T, double P);
double Majorite_Garnet_calib_cv(double T, double P);
double Majorite_Garnet_calib_cp(double T, double P);
double Majorite_Garnet_calib_dcpdt(double T, double P);
double Majorite_Garnet_calib_alpha(double T, double P);
double Majorite_Garnet_calib_beta(double T, double P);
double Majorite_Garnet_calib_K(double T, double P);
double Majorite_Garnet_calib_Kp(double T, double P);

int Majorite_Garnet_get_param_number(void);
const char **Majorite_Garnet_get_param_names(void);
const char **Majorite_Garnet_get_param_units(void);
void Majorite_Garnet_get_param_values(double **values);
int Majorite_Garnet_set_param_values(double *values);
double Majorite_Garnet_get_param_value(int index);
int Majorite_Garnet_set_param_value(int index, double value);

double Majorite_Garnet_dparam_g(double T, double P, int index);
double Majorite_Garnet_dparam_dgdt(double T, double P, int index);
double Majorite_Garnet_dparam_dgdp(double T, double P, int index);
double Majorite_Garnet_dparam_d2gdt2(double T, double P, int index);
double Majorite_Garnet_dparam_d2gdtdp(double T, double P, int index);
double Majorite_Garnet_dparam_d2gdp2(double T, double P, int index);
double Majorite_Garnet_dparam_d3gdt3(double T, double P, int index);
double Majorite_Garnet_dparam_d3gdt2dp(double T, double P, int index);
double Majorite_Garnet_dparam_d3gdtdp2(double T, double P, int index);
double Majorite_Garnet_dparam_d3gdp3(double T, double P, int index);

