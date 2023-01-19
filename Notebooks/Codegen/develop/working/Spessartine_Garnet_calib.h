
const char *Spessartine_Garnet_calib_identifier(void);
const char *Spessartine_Garnet_calib_name(void);
const char *Spessartine_Garnet_calib_formula(void);
const double Spessartine_Garnet_calib_mw(void);
const double *Spessartine_Garnet_calib_elements(void);

double Spessartine_Garnet_calib_g(double T, double P);
double Spessartine_Garnet_calib_dgdt(double T, double P);
double Spessartine_Garnet_calib_dgdp(double T, double P);
double Spessartine_Garnet_calib_d2gdt2(double T, double P);
double Spessartine_Garnet_calib_d2gdtdp(double T, double P);
double Spessartine_Garnet_calib_d2gdp2(double T, double P);
double Spessartine_Garnet_calib_d3gdt3(double T, double P);
double Spessartine_Garnet_calib_d3gdt2dp(double T, double P);
double Spessartine_Garnet_calib_d3gdtdp2(double T, double P);
double Spessartine_Garnet_calib_d3gdp3(double T, double P);

double Spessartine_Garnet_calib_s(double T, double P);
double Spessartine_Garnet_calib_v(double T, double P);
double Spessartine_Garnet_calib_cv(double T, double P);
double Spessartine_Garnet_calib_cp(double T, double P);
double Spessartine_Garnet_calib_dcpdt(double T, double P);
double Spessartine_Garnet_calib_alpha(double T, double P);
double Spessartine_Garnet_calib_beta(double T, double P);
double Spessartine_Garnet_calib_K(double T, double P);
double Spessartine_Garnet_calib_Kp(double T, double P);

int Spessartine_Garnet_get_param_number(void);
const char **Spessartine_Garnet_get_param_names(void);
const char **Spessartine_Garnet_get_param_units(void);
void Spessartine_Garnet_get_param_values(double **values);
int Spessartine_Garnet_set_param_values(double *values);
double Spessartine_Garnet_get_param_value(int index);
int Spessartine_Garnet_set_param_value(int index, double value);

double Spessartine_Garnet_dparam_g(double T, double P, int index);
double Spessartine_Garnet_dparam_dgdt(double T, double P, int index);
double Spessartine_Garnet_dparam_dgdp(double T, double P, int index);
double Spessartine_Garnet_dparam_d2gdt2(double T, double P, int index);
double Spessartine_Garnet_dparam_d2gdtdp(double T, double P, int index);
double Spessartine_Garnet_dparam_d2gdp2(double T, double P, int index);
double Spessartine_Garnet_dparam_d3gdt3(double T, double P, int index);
double Spessartine_Garnet_dparam_d3gdt2dp(double T, double P, int index);
double Spessartine_Garnet_dparam_d3gdtdp2(double T, double P, int index);
double Spessartine_Garnet_dparam_d3gdp3(double T, double P, int index);

