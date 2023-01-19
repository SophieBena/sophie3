# Cython numpy wrapper code for arrays is taken from:
# http://gael-varoquaux.info/programming/cython-example-of-exposing-c-computed-arrays-in-python-without-data-copies.html
# Author: Gael Varoquaux, BSD license

# Declare the prototype of the C functions
cdef extern from "Pyrope_Garnet_calib.h":
    const char *Pyrope_Garnet_calib_identifier();
    const char *Pyrope_Garnet_calib_name();
    const char *Pyrope_Garnet_calib_formula();
    const double Pyrope_Garnet_calib_mw();
    const double *Pyrope_Garnet_calib_elements();
    double Pyrope_Garnet_calib_g(double t, double p);
    double Pyrope_Garnet_calib_dgdt(double t, double p);
    double Pyrope_Garnet_calib_dgdp(double t, double p);
    double Pyrope_Garnet_calib_d2gdt2(double t, double p);
    double Pyrope_Garnet_calib_d2gdtdp(double t, double p);
    double Pyrope_Garnet_calib_d2gdp2(double t, double p);
    double Pyrope_Garnet_calib_d3gdt3(double t, double p);
    double Pyrope_Garnet_calib_d3gdt2dp(double t, double p);
    double Pyrope_Garnet_calib_d3gdtdp2(double t, double p);
    double Pyrope_Garnet_calib_d3gdp3(double t, double p);
    double Pyrope_Garnet_calib_s(double t, double p);
    double Pyrope_Garnet_calib_v(double t, double p);
    double Pyrope_Garnet_calib_cv(double t, double p);
    double Pyrope_Garnet_calib_cp(double t, double p);
    double Pyrope_Garnet_calib_dcpdt(double t, double p);
    double Pyrope_Garnet_calib_alpha(double t, double p);
    double Pyrope_Garnet_calib_beta(double t, double p);
    double Pyrope_Garnet_calib_K(double t, double p);
    double Pyrope_Garnet_calib_Kp(double t, double p);
    int Pyrope_Garnet_get_param_number();
    const char **Pyrope_Garnet_get_param_names();
    const char **Pyrope_Garnet_get_param_units();
    void Pyrope_Garnet_get_param_values(double **values);
    int Pyrope_Garnet_set_param_values(double *values);
    double Pyrope_Garnet_get_param_value(int index);
    int Pyrope_Garnet_set_param_value(int index, double value);
    double Pyrope_Garnet_dparam_g(double t, double p, int index);
    double Pyrope_Garnet_dparam_dgdt(double t, double p, int index);
    double Pyrope_Garnet_dparam_dgdp(double t, double p, int index);
    double Pyrope_Garnet_dparam_d2gdt2(double t, double p, int index);
    double Pyrope_Garnet_dparam_d2gdtdp(double t, double p, int index);
    double Pyrope_Garnet_dparam_d2gdp2(double t, double p, int index);
    double Pyrope_Garnet_dparam_d3gdt3(double t, double p, int index);
    double Pyrope_Garnet_dparam_d3gdt2dp(double t, double p, int index);
    double Pyrope_Garnet_dparam_d3gdtdp2(double t, double p, int index);
    double Pyrope_Garnet_dparam_d3gdp3(double t, double p, int index);

from libc.stdlib cimport malloc, free
from cpython cimport PyObject, Py_INCREF
import ctypes

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# here is the "wrapper" signature
def cy_Pyrope_Garnet_calib_identifier():
    result = <bytes> Pyrope_Garnet_calib_identifier()
    return result.decode('UTF-8')
def cy_Pyrope_Garnet_calib_name():
    result = <bytes> Pyrope_Garnet_calib_name()
    return result.decode('UTF-8')
def cy_Pyrope_Garnet_calib_formula():
    result = <bytes> Pyrope_Garnet_calib_formula()
    return result.decode('UTF-8')
def cy_Pyrope_Garnet_calib_mw():
    result = Pyrope_Garnet_calib_mw()
    return result
def cy_Pyrope_Garnet_calib_elements():
    cdef const double *e = Pyrope_Garnet_calib_elements()
    np_array = np.zeros(106)
    for i in range(0,106):
        np_array[i] = e[i]
    return np_array
def cy_Pyrope_Garnet_calib_g(double t, double p):
    result = Pyrope_Garnet_calib_g(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_dgdt(double t, double p):
    result = Pyrope_Garnet_calib_dgdt(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_dgdp(double t, double p):
    result = Pyrope_Garnet_calib_dgdp(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d2gdt2(double t, double p):
    result = Pyrope_Garnet_calib_d2gdt2(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d2gdtdp(double t, double p):
    result = Pyrope_Garnet_calib_d2gdtdp(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d2gdp2(double t, double p):
    result = Pyrope_Garnet_calib_d2gdp2(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d3gdt3(double t, double p):
    result = Pyrope_Garnet_calib_d3gdt3(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d3gdt2dp(double t, double p):
    result = Pyrope_Garnet_calib_d3gdt2dp(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d3gdtdp2(double t, double p):
    result = Pyrope_Garnet_calib_d3gdtdp2(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_d3gdp3(double t, double p):
    result = Pyrope_Garnet_calib_d3gdp3(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_s(double t, double p):
    result = Pyrope_Garnet_calib_s(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_v(double t, double p):
    result = Pyrope_Garnet_calib_v(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_cv(double t, double p):
    result = Pyrope_Garnet_calib_cv(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_cp(double t, double p):
    result = Pyrope_Garnet_calib_cp(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_dcpdt(double t, double p):
    result = Pyrope_Garnet_calib_dcpdt(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_alpha(double t, double p):
    result = Pyrope_Garnet_calib_alpha(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_beta(double t, double p):
    result = Pyrope_Garnet_calib_beta(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_K(double t, double p):
    result = Pyrope_Garnet_calib_K(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_calib_Kp(double t, double p):
    result = Pyrope_Garnet_calib_Kp(<double> t, <double> p)
    return result
def cy_Pyrope_Garnet_get_param_number():
    result = Pyrope_Garnet_get_param_number()
    return result
def cy_Pyrope_Garnet_get_param_names():
    cdef const char **names = Pyrope_Garnet_get_param_names()
    n = Pyrope_Garnet_get_param_number()
    result = []
    for i in range(0,n):
        entry = <bytes> names[i]
        result.append(entry.decode('UTF-8'))
    return result
def cy_Pyrope_Garnet_get_param_units():
    cdef const char **units = Pyrope_Garnet_get_param_units()
    n = Pyrope_Garnet_get_param_number()
    result = []
    for i in range(0,n):
        entry = <bytes> units[i]
        result.append(entry.decode('UTF-8'))
    return result
def cy_Pyrope_Garnet_get_param_values():
    n = Pyrope_Garnet_get_param_number()
    cdef double *m = <double *>malloc(n*sizeof(double))
    Pyrope_Garnet_get_param_values(&m)
    np_array = np.zeros(n)
    for i in range(n):
        np_array[i] = m[i]
    free(m)
    return np_array
def cy_Pyrope_Garnet_set_param_values(np_array):
    n = len(np_array)
    cdef double *m = <double *>malloc(n*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = Pyrope_Garnet_set_param_values(m);
    free(m)
    return result
def cy_Pyrope_Garnet_get_param_value(int index):
    result = Pyrope_Garnet_get_param_value(<int> index)
    return result
def cy_Pyrope_Garnet_set_param_value(int index, double value):
    result = Pyrope_Garnet_set_param_value(<int> index, <double> value)
    return result
def cy_Pyrope_Garnet_dparam_g(double t, double p, int index):
    result = Pyrope_Garnet_dparam_g(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_dgdt(double t, double p, int index):
    result = Pyrope_Garnet_dparam_dgdt(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_dgdp(double t, double p, int index):
    result = Pyrope_Garnet_dparam_dgdp(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d2gdt2(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d2gdt2(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d2gdtdp(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d2gdtdp(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d2gdp2(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d2gdp2(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d3gdt3(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d3gdt3(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d3gdt2dp(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d3gdt2dp(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d3gdtdp2(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d3gdtdp2(<double> t, <double> p, <int> index)
    return result
def cy_Pyrope_Garnet_dparam_d3gdp3(double t, double p, int index):
    result = Pyrope_Garnet_dparam_d3gdp3(<double> t, <double> p, <int> index)
    return result
