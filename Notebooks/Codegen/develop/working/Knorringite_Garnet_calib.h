
const char *Knorringite_Garnet_calib_identifier(void);
const char *Knorringite_Garnet_calib_name(void);
const char *Knorringite_Garnet_calib_formula(void);
const double Knorringite_Garnet_calib_mw(void);
const double *Knorringite_Garnet_calib_elements(void);

double Knorringite_Garnet_calib_g(double T, double P);
double Knorringite_Garnet_calib_dgdt(double T, double P);
double Knorringite_Garnet_calib_dgdp(double T, double P);
double Knorringite_Garnet_calib_d2gdt2(double T, double P);
double Knorringite_Garnet_calib_d2gdtdp(double T, double P);
double Knorringite_Garnet_calib_d2gdp2(double T, double P);
double Knorringite_Garnet_calib_d3gdt3(double T, double P);
double Knorringite_Garnet_calib_d3gdt2dp(double T, double P);
double Knorringite_Garnet_calib_d3gdtdp2(double T, double P);
double Knorringite_Garnet_calib_d3gdp3(double T, double P);

double Knorringite_Garnet_calib_s(double T, double P);
double Knorringite_Garnet_calib_v(double T, double P);
double Knorringite_Garnet_calib_cv(double T, double P);
double Knorringite_Garnet_calib_cp(double T, double P);
double Knorringite_Garnet_calib_dcpdt(double T, double P);
double Knorringite_Garnet_calib_alpha(double T, double P);
double Knorringite_Garnet_calib_beta(double T, double P);
double Knorringite_Garnet_calib_K(double T, double P);
double Knorringite_Garnet_calib_Kp(double T, double P);

int Knorringite_Garnet_get_param_number(void);
const char **Knorringite_Garnet_get_param_names(void);
const char **Knorringite_Garnet_get_param_units(void);
void Knorringite_Garnet_get_param_values(double **values);
int Knorringite_Garnet_set_param_values(double *values);
double Knorringite_Garnet_get_param_value(int index);
int Knorringite_Garnet_set_param_value(int index, double value);

double Knorringite_Garnet_dparam_g(double T, double P, int index);
double Knorringite_Garnet_dparam_dgdt(double T, double P, int index);
double Knorringite_Garnet_dparam_dgdp(double T, double P, int index);
double Knorringite_Garnet_dparam_d2gdt2(double T, double P, int index);
double Knorringite_Garnet_dparam_d2gdtdp(double T, double P, int index);
double Knorringite_Garnet_dparam_d2gdp2(double T, double P, int index);
double Knorringite_Garnet_dparam_d3gdt3(double T, double P, int index);
double Knorringite_Garnet_dparam_d3gdt2dp(double T, double P, int index);
double Knorringite_Garnet_dparam_d3gdtdp2(double T, double P, int index);
double Knorringite_Garnet_dparam_d3gdp3(double T, double P, int index);

