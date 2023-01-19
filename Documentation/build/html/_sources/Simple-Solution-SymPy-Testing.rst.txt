
Simple Solution SymPy Testing Theory and Code Generation
========================================================

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()

Simple Solution Properties - Structure of the Equations
-------------------------------------------------------

There are three terms: - Terms describing standard state contributions -
Terms describing the configurational entropy of solution - Terms
describing the excess enthalpy of solution

Assumptions: - There are :math:`c` components in the system - There are
as many endmember species, :math:`s`, as there are components - The
configurational entropy is described as a simple :math:`x_i log(x_i)`
sum - The excess enthalpy is described using an asymmetric regular
solution
:math:`\left[ {{W_{ij}} + \Delta {W_{ij}}\left( {{X_i} - {X_j}} \right)} \right]{X_i}{X_j}`,
where the :math:`W_{ij}` and :math:`\Delta{W_{ij}}` are allowed to be
first order functions of both :math:`T` and :math:`P` - Ternary
interaction terms are permittted, i.e. :math:`W_{ijk}`

Number of solution components
-----------------------------

This notebook illustrates a three component solution

.. code:: ipython3

    c = 3

Define primary compositional variables
--------------------------------------

= :math:`n` is a vector of mole numbers of each component ## Define
derived compositional variables - :math:`n_T` is the total number of
moles in the solution - :math:`X` is a vector of mole fractions of
components in the system

.. code:: ipython3

    component_string = ''
    for i in range(1,c+1):
        component_string += 'n' + str(i) + ' '
    n = sym.Matrix(list(sym.symbols(component_string)))
    nT = (sym.ones(1,c) * n)[0]
    X = n/nT
    nT, X




.. math::

    \left ( n_{1} + n_{2} + n_{3}, \quad \left[\begin{matrix}\frac{n_{1}}{n_{1} + n_{2} + n_{3}}\\\frac{n_{2}}{n_{1} + n_{2} + n_{3}}\\\frac{n_{3}}{n_{1} + n_{2} + n_{3}}\end{matrix}\right]\right )



Define temperature and pressure
-------------------------------

-  :math:`T` is temperature in :math:`K`
-  :math:`P` is pressure in :math:`bars`

.. code:: ipython3

    T,P = sym.symbols('T P')

Define standard state properties
--------------------------------

.. code:: ipython3

    ss_list = []
    for i in range(1,c+1):
        ss_string = 'mu' + str(i)
        ss_list.append(sym.Function(ss_string)(T,P))
    mu = sym.Matrix(ss_list)
    mu




.. math::

    \left[\begin{matrix}\mu_{1}{\left (T,P \right )}\\\mu_{2}{\left (T,P \right )}\\\mu_{3}{\left (T,P \right )}\end{matrix}\right]



.. code:: ipython3

    G_ss = (n.transpose()*mu)[0]
    G_ss




.. math::

    n_{1} \mu_{1}{\left (T,P \right )} + n_{2} \mu_{2}{\left (T,P \right )} + n_{3} \mu_{3}{\left (T,P \right )}



Define configurational entropy and configurational Gibbs free energy
--------------------------------------------------------------------

.. code:: ipython3

    S_config,R = sym.symbols('S_config R')
    S_config = 0
    for i in range(0,c):
        S_config += X[i]*sym.log(X[i])
    S_config *= -R*nT
    S_config




.. math::

    - R \left(n_{1} + n_{2} + n_{3}\right) \left(\frac{n_{1} \log{\left (\frac{n_{1}}{n_{1} + n_{2} + n_{3}} \right )}}{n_{1} + n_{2} + n_{3}} + \frac{n_{2} \log{\left (\frac{n_{2}}{n_{1} + n_{2} + n_{3}} \right )}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3} \log{\left (\frac{n_{3}}{n_{1} + n_{2} + n_{3}} \right )}}{n_{1} + n_{2} + n_{3}}\right)



.. code:: ipython3

    G_config = sym.simplify(-T*S_config)
    G_config




.. math::

    R T \left(n_{1} \log{\left (\frac{n_{1}}{n_{1} + n_{2} + n_{3}} \right )} + n_{2} \log{\left (\frac{n_{2}}{n_{1} + n_{2} + n_{3}} \right )} + n_{3} \log{\left (\frac{n_{3}}{n_{1} + n_{2} + n_{3}} \right )}\right)



Parameterize the excess free energy
-----------------------------------

-  Symmetric terms: :math:`W_{ij} = Wh_{ij} - T Ws_{ij} + P Wv_{ij}`,
   where :math:`Wh_{ij}` is the excess enthalpy along the
   :math:`i`-:math:`j` binary, :math:`Ws_{ij}` is the excess entropy,
   and :math:`Wv_{ij}` is the excess volume
-  Asymetric terms:
   :math:`\Delta W_{ij} = \Delta Wh_{ij} - T \Delta Ws_{ij} + P \Delta Wv_{ij}`
-  Convention:
   :math:`\left[ {{W_{ij}} + \Delta {W_{ij}}\left( {{X_i} - {X_j}} \right)} \right]{X_i}{X_j} = \left( {{W_{ij}} + \Delta {W_{ij}}} \right)X_i^2{X_j} + \left( {{W_{ij}} - \Delta {W_{ij}}} \right){X_i}X_j^2 = {W_{iij}}X_i^2{X_j} + {W_{ijj}}{X_i}X_j^2`
-  Strictly ternary terms: :math:`W_{ijk}X_iX_jX_k`, where
   :math:`i {\ne} j {\ne} k`

.. code:: ipython3

    module = 'asymm_regular'
    params = []
    units = []
    symparam = ()
    G_excess = sym.symbols('G_excess')
    G_excess = 0
    for i in range(1,c):
        for j in range(i+1,c+1):
            param = 'Wh' + str(i) + str(j); params.append(param); units.append('J/m')
            h_term = sym.symbols(param); symparam += (h_term,)
            param = 'Ws' + str(i) + str(j); params.append(param); units.append('J/K-m')
            s_term = sym.symbols(param); symparam += (s_term,)
            param = 'Wv' + str(i) + str(j); params.append(param); units.append('J/bar-m')
            v_term = sym.symbols(param); symparam += (v_term,)
            param = 'dWh' + str(i) + str(j); params.append(param); units.append('J/m')
            dh_term = sym.symbols(param); symparam += (dh_term,)
            param = 'dWs' + str(i) + str(j); params.append(param); units.append('J/K-m')
            ds_term = sym.symbols(param); symparam += (ds_term,)
            param = 'dWv' + str(i) + str(j); params.append(param); units.append('J/bar-m')
            dv_term = sym.symbols(param); symparam += (dv_term,)
            w_term = h_term - T*s_term + P*v_term
            dw_term = dh_term - T*ds_term + P*dv_term
            G_excess += (w_term + dw_term*(n[i-1]-n[j-1])/nT)*n[i-1]*n[j-1]
    G_excess /= nT
    for i in range(1,c-1):
        for j in range(i+1,c):
            for k in range(j+1,c+1):
                param = 'Wh' + str(i) + str(j) + str(k); params.append(param); units.append('J/m')
                h_term = sym.symbols(param); symparam += (h_term,)
                param = 'Ws' + str(i) + str(j) + str(k); params.append(param); units.append('J/K-m')
                s_term = sym.symbols(param); symparam += (s_term,)
                param = 'Wv' + str(i) + str(j) + str(k); params.append(param); units.append('J/bar-m')
                v_term = sym.symbols(param); symparam += (v_term,)
                G_excess += (h_term - T*s_term + P*v_term)*n[i-1]*n[j-1]*n[k-1]/nT/nT
    G_excess




.. math::

    \frac{n_{1} n_{2} n_{3} \left(P Wv_{123} - T Ws_{123} + Wh_{123}\right)}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{n_{1} n_{2} \left(P Wv_{12} - T Ws_{12} + Wh_{12} + \frac{\left(n_{1} - n_{2}\right) \left(P dWv_{12} - T dWs_{12} + dWh_{12}\right)}{n_{1} + n_{2} + n_{3}}\right) + n_{1} n_{3} \left(P Wv_{13} - T Ws_{13} + Wh_{13} + \frac{\left(n_{1} - n_{3}\right) \left(P dWv_{13} - T dWs_{13} + dWh_{13}\right)}{n_{1} + n_{2} + n_{3}}\right) + n_{2} n_{3} \left(P Wv_{23} - T Ws_{23} + Wh_{23} + \frac{\left(n_{2} - n_{3}\right) \left(P dWv_{23} - T dWs_{23} + dWh_{23}\right)}{n_{1} + n_{2} + n_{3}}\right)}{n_{1} + n_{2} + n_{3}}



.. code:: ipython3

    print(params)
    print(units)


.. parsed-literal::

    ['Wh12', 'Ws12', 'Wv12', 'dWh12', 'dWs12', 'dWv12', 'Wh13', 'Ws13', 'Wv13', 'dWh13', 'dWs13', 'dWv13', 'Wh23', 'Ws23', 'Wv23', 'dWh23', 'dWs23', 'dWv23', 'Wh123', 'Ws123', 'Wv123']
    ['J/m', 'J/K-m', 'J/bar-m', 'J/m', 'J/K-m', 'J/bar-m', 'J/m', 'J/K-m', 'J/bar-m', 'J/m', 'J/K-m', 'J/bar-m', 'J/m', 'J/K-m', 'J/bar-m', 'J/m', 'J/K-m', 'J/bar-m', 'J/m', 'J/K-m', 'J/bar-m']


Define the Gibbs free energy of solution
----------------------------------------

.. code:: ipython3

    G = G_ss + G_config + G_excess
    G




.. math::

    R T \left(n_{1} \log{\left (\frac{n_{1}}{n_{1} + n_{2} + n_{3}} \right )} + n_{2} \log{\left (\frac{n_{2}}{n_{1} + n_{2} + n_{3}} \right )} + n_{3} \log{\left (\frac{n_{3}}{n_{1} + n_{2} + n_{3}} \right )}\right) + \frac{n_{1} n_{2} n_{3} \left(P Wv_{123} - T Ws_{123} + Wh_{123}\right)}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + n_{1} \mu_{1}{\left (T,P \right )} + n_{2} \mu_{2}{\left (T,P \right )} + n_{3} \mu_{3}{\left (T,P \right )} + \frac{n_{1} n_{2} \left(P Wv_{12} - T Ws_{12} + Wh_{12} + \frac{\left(n_{1} - n_{2}\right) \left(P dWv_{12} - T dWs_{12} + dWh_{12}\right)}{n_{1} + n_{2} + n_{3}}\right) + n_{1} n_{3} \left(P Wv_{13} - T Ws_{13} + Wh_{13} + \frac{\left(n_{1} - n_{3}\right) \left(P dWv_{13} - T dWs_{13} + dWh_{13}\right)}{n_{1} + n_{2} + n_{3}}\right) + n_{2} n_{3} \left(P Wv_{23} - T Ws_{23} + Wh_{23} + \frac{\left(n_{2} - n_{3}\right) \left(P dWv_{23} - T dWs_{23} + dWh_{23}\right)}{n_{1} + n_{2} + n_{3}}\right)}{n_{1} + n_{2} + n_{3}}



Instantiate C code printers and modify default behavior
-------------------------------------------------------

The C99 code printing class is subclassed and two functions are added: -
The first overrides the pow() function for products with integer
exponent <= 4 - The second implements a version of user function
derivative printing specified to endmember chemical potential
derivatives called externally from a C structure (see implementation
below)

.. code:: ipython3

    from sympy.printing.ccode import C99CodePrinter
    class SubCodePrinter(C99CodePrinter):
        def _print_Pow(self, expr):
            if expr.exp.is_integer and expr.exp > 0 and expr.exp <= 4:
                result = ')*('.join([self._print(expr.base) for i in range(expr.exp)])
                return '((' + result + '))'
            else:
                return super()._print_Pow(expr)
        def _print_Derivative(self, expr):
            number_of_derivatives = len(expr.args) - 1
            function_string_index = int(sym.srepr(expr.args[0].func).split("'")[1][2:]) - 1
            result = ''
            if number_of_derivatives == 1:
                derivative_string = sym.srepr(expr.args[1]).split("'")[1]
                result = '(*endmember[' + str(function_string_index) + '].dmu0d' + derivative_string + ')(T, P)'
            elif number_of_derivatives == 2:
                derivative_string_1 = sym.srepr(expr.args[1]).split("'")[1]
                derivative_string_2 = sym.srepr(expr.args[2]).split("'")[1]
                if derivative_string_1 == 'T' and derivative_string_2 == 'T':
                    result = '(*endmember[' + str(function_string_index) + '].d2mu0dT2)(T, P)'
                elif derivative_string_1 == 'P' and derivative_string_2 == 'T':
                    result = '(*endmember[' + str(function_string_index) + '].d2mu0dTdP)(T, P)'
                elif derivative_string_1 == 'P' and derivative_string_2 == 'P':
                    result = '(*endmember[' + str(function_string_index) + '].d2mu0dP2)(T, P)'
            elif number_of_derivatives == 3:
                derivative_string_1 = sym.srepr(expr.args[1]).split("'")[1]
                derivative_string_2 = sym.srepr(expr.args[2]).split("'")[1]
                derivative_string_3 = sym.srepr(expr.args[3]).split("'")[1]
                if derivative_string_1 == 'T' and derivative_string_2 == 'T' and derivative_string_3 == 'T':
                    result = '(*endmember[' + str(function_string_index) + '].d3mu0dT3)(T, P)'
                elif derivative_string_1 == 'P' and derivative_string_2 == 'T' and derivative_string_3 == 'T':
                    result = '(*endmember[' + str(function_string_index) + '].d3mu0dT2dP)(T, P)'
                elif derivative_string_1 == 'P' and derivative_string_2 == 'P' and derivative_string_3 == 'T':
                    result = '(*endmember[' + str(function_string_index) + '].d3mu0dTdP2)(T, P)'
                elif derivative_string_1 == 'P' and derivative_string_2 == 'P' and derivative_string_3 == 'P':
                    result = '(*endmember[' + str(function_string_index) + '].d3mu0dP3)(T, P)'
            return result

Instantiate the mldified C99 printing class
-------------------------------------------

and add definitions of user functions that describe endmember chemical
potentials

::

    typedef struct _endmembers {
      const char *(*name) (void);
      const char *(*formula) (void);
      const double (*mw) (void);
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
        Akermanite_berman_name,
        Akermanite_berman_formula,
        Akermanite_berman_mw,
        Akermanite_berman_g,
        Akermanite_berman_dgdt,
        Akermanite_berman_dgdp,
        Akermanite_berman_d2gdt2,
        Akermanite_berman_d2gdtdp,
        Akermanite_berman_d2gdp2,
        Akermanite_berman_d3gdt3,
        Akermanite_berman_d3gdt2dp,
        Akermanite_berman_d3gdtdp2,
        Akermanite_berman_d3gdp3
      },
      {
      ...
      }
    };
    static int nc = (sizeof endmember / sizeof(_endmembers));

.. code:: ipython3

    printer = SubCodePrinter(settings={'user_functions':{'mu1':'(*endmember[0].mu0)', 'mu2':'(*endmember[1].mu0)', \
                                       'mu3':'(*endmember[2].mu0)'}})
    printer.doprint(G.diff(n[0]))




.. parsed-literal::

    'R*T*(-n2/(n1 + n2 + n3) - n3/(n1 + n2 + n3) + (-n1/((n1 + n2 + n3)*(n1 + n2 + n3)) + 1.0/(n1 + n2 + n3))*(n1 + n2 + n3) + log(n1/(n1 + n2 + n3))) - 2*n1*n2*n3*(P*Wv123 - T*Ws123 + Wh123)/((n1 + n2 + n3)*(n1 + n2 + n3)*(n1 + n2 + n3)) + n2*n3*(P*Wv123 - T*Ws123 + Wh123)/((n1 + n2 + n3)*(n1 + n2 + n3)) + (*endmember[0].mu0)(T, P) + (n1*n2*(-(n1 - n2)*(P*dWv12 - T*dWs12 + dWh12)/((n1 + n2 + n3)*(n1 + n2 + n3)) + (P*dWv12 - T*dWs12 + dWh12)/(n1 + n2 + n3)) + n1*n3*(-(n1 - n3)*(P*dWv13 - T*dWs13 + dWh13)/((n1 + n2 + n3)*(n1 + n2 + n3)) + (P*dWv13 - T*dWs13 + dWh13)/(n1 + n2 + n3)) - n2*n3*(n2 - n3)*(P*dWv23 - T*dWs23 + dWh23)/((n1 + n2 + n3)*(n1 + n2 + n3)) + n2*(P*Wv12 - T*Ws12 + Wh12 + (n1 - n2)*(P*dWv12 - T*dWs12 + dWh12)/(n1 + n2 + n3)) + n3*(P*Wv13 - T*Ws13 + Wh13 + (n1 - n3)*(P*dWv13 - T*dWs13 + dWh13)/(n1 + n2 + n3)))/(n1 + n2 + n3) - (n1*n2*(P*Wv12 - T*Ws12 + Wh12 + (n1 - n2)*(P*dWv12 - T*dWs12 + dWh12)/(n1 + n2 + n3)) + n1*n3*(P*Wv13 - T*Ws13 + Wh13 + (n1 - n3)*(P*dWv13 - T*dWs13 + dWh13)/(n1 + n2 + n3)) + n2*n3*(P*Wv23 - T*Ws23 + Wh23 + (n2 - n3)*(P*dWv23 - T*dWs23 + dWh23)/(n1 + n2 + n3)))/((n1 + n2 + n3)*(n1 + n2 + n3))'



Use code printers to construct “C” package code
===============================================

.. code:: ipython3

    !mkdir =p working
    %cd working


.. parsed-literal::

    mkdir: =p: File exists
    mkdir: working: File exists
    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working


.. code:: ipython3

    moles_assign_text = ''
    for i in range(0,c):
        moles_assign_text += '    double n' + str(i+1) + ' = n[' + str(i) + '];'
        if i < c-1:
            moles_assign_text += '\n'
    print (moles_assign_text)


.. parsed-literal::

        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];


(1)
~~~

Function template for primary routines in solution_calc.h

.. code:: ipython3

    solution_calc_template = """\
    
    static double {module}_{func}(double T, double P, double n[{number_components}]) {{
    {moles_assign}
        double result = {g_code};
        return result;
    }}
    \
    """
    
    solution_derivative_template = """\
    
    static void {module}_{func}(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {moles_assign}
    {derivative_code}
    }}
    \
    """

| Primary routines in solution_calc.h
| Note that the compositional derivatives augment the standard routines
  inherited from the endmember api

.. code:: ipython3

    f_list = ['g', 'dgdt', 'dgdp', 'd2gdt2', 'd2gdtdp', 'd2gdp2', 'd3gdt3', 'd3gdt2dp', 'd3gdtdp2', 'd3gdp3']
    G_matrix = sym.Matrix([G, G.diff(T), G.diff(P), G.diff(T,T), G.diff(T,P), G.diff(P,P), \
                           G.diff(T,T,T), G.diff(T,T,P), G.diff(T,P,P), G.diff(P,P,P)])
    dn_f_list = ['dgdn', 'd2gdndt', 'd2gdndp', 'd3gdndt2', 'd3gdndtdp', 'd3gdndp2', 'd4gdndt3', 'd4gdndt2dp', \
                 'd4gdndtdp2', 'd4gdndp3']
    d2n_f_list = ['d2gdn2', 'd3gdn2dt', 'd3gdn2dp', 'd4gdn2dt2', 'd4gdn2dtdp', 'd4gdn2dp2', 'd5gdn2dt3', 'd5gdn2dt2dp', \
                 'd5gdn2dtdp2', 'd5gdn2dp3']
    d3n_f_list = ['d3gdn3', 'd4gdn3dt', 'd4gdn3dp', 'd5gdn3dt2', 'd5gdn3dtdp', 'd5gdn3dp2', 'd6gdn3dt3', 'd6gdn3dt2dp', \
                 'd6gdn3dtdp2', 'd6gdn3dp3']

Note that the compositional derivatives are returned in a pre-allocated
one dimensional array that is passed to the function (which returns
void). The higher-oredr compositional derivatives are returned in
compact storage mode, as defined below.

.. code:: ipython3

    solution_calc = '#include <math.h>\n\n'
    for i in range(0,len(f_list)):
        solution_calc += solution_calc_template.format(\
                                  module=module, \
                                  func=f_list[i], \
                                  g_code=printer.doprint(G_matrix[i]), \
                                  number_components=c, \
                                  moles_assign=moles_assign_text)
        
        derivative_code_text = ''
        for j in range(0,c):
            derivative_code_text += '    result[' + str(j) + '] = ' + printer.doprint(G_matrix[i].diff(n[j])) + ';\n'
        solution_calc += solution_derivative_template.format(\
                                  module=module, \
                                  func=dn_f_list[i], \
                                  number_components=c, \
                                  derivative_code=derivative_code_text, \
                                  moles_assign=moles_assign_text)
        
        derivative_code_text = ''
        l = 0
        for j in range(0,c):
            for k in range(j,c):
                derivative_code_text += '    result[' + str(l) + '] = ' \
                                      + printer.doprint(G_matrix[i].diff(n[j]).diff(n[k])) + ';\n'
                l += 1
        solution_calc += solution_derivative_template.format(\
                                  module=module, \
                                  func=d2n_f_list[i], \
                                  number_components=c, \
                                  derivative_code=derivative_code_text, \
                                  moles_assign=moles_assign_text)
        
        derivative_code_text = ''
        m = 0
        for j in range(0,c):
            for k in range(j,c):
                for l in range(k,c):
                    derivative_code_text += '    result[' + str(m) + '] = ' \
                                          + printer.doprint(G_matrix[i].diff(n[j]).diff(n[k]).diff(n[l])) + ';\n'
                    m += 1
        solution_calc += solution_derivative_template.format(\
                                  module=module, \
                                  func=d3n_f_list[i], \
                                  number_components=c, \
                                  derivative_code=derivative_code_text, \
                                  moles_assign=moles_assign_text)

Compact storage scheme for the second order compositional derivatives.
Elements of the upper triangle of the symmetric Hessian matrix are
returned.

.. code:: ipython3

    k = 0
    for i in range(1,c+1):
        for j in range (i,c+1):
            print ((i,j), ' = ', k)
            k += 1


.. parsed-literal::

    (1, 1)  =  0
    (1, 2)  =  1
    (1, 3)  =  2
    (2, 2)  =  3
    (2, 3)  =  4
    (3, 3)  =  5


Compact storage scheme for the third order compositional derivatives.
Elements of the upper wedge of the symmetric third order array are
returned.

.. code:: ipython3

    l = 0
    for i in range (1,c+1):
        for j in range (i,c+1):
            for k in range (j,c+1):
                print ((i,j,k), ' = ', l)
                l += 1


.. parsed-literal::

    (1, 1, 1)  =  0
    (1, 1, 2)  =  1
    (1, 1, 3)  =  2
    (1, 2, 2)  =  3
    (1, 2, 3)  =  4
    (1, 3, 3)  =  5
    (2, 2, 2)  =  6
    (2, 2, 3)  =  7
    (2, 3, 3)  =  8
    (3, 3, 3)  =  9


| Add convenience routiines to soution_calc.h
| These routines mimic those of the endmember API. Convenience routines
  manipulating compositional derivatives are not included, although they
  could be defined if the need arises.

.. code:: ipython3

    convenience_template = """\
    
    static double {module}_s(double T, double P, double n[{number_components}]) {{
        double result = -{module}_dgdt(T, P, n);
        return result;
    }}
    
    static double {module}_v(double T, double P, double n[{number_components}]) {{
        double result = {module}_dgdp(T, P, n);
        return result;
    }}
    
    static double {module}_cv(double T, double P, double n[{number_components}]) {{
        double result = -T*{module}_d2gdt2(T, P, n);
        double dvdt = {module}_d2gdtdp(T, P, n);
        double dvdp = {module}_d2gdp2(T, P, n);
        result += T*dvdt*dvdt/dvdp;
        return result;
    }}
    
    static double {module}_cp(double T, double P, double n[{number_components}]) {{
        double result = -T*{module}_d2gdt2(T, P, n);
        return result;
    }}
    
    static double {module}_dcpdt(double T, double P, double n[{number_components}]) {{
        double result = -T*{module}_d3gdt3(T, P, n) - {module}_d2gdt2(T, P, n);
        return result;
    }}
    
    static double {module}_alpha(double T, double P, double n[{number_components}]) {{
        double result = {module}_d2gdtdp(T, P, n)/{module}_dgdp(T, P, n);
        return result;
    }}
    
    static double {module}_beta(double T, double P, double n[{number_components}]) {{
        double result = -{module}_d2gdp2(T, P, n)/{module}_dgdp(T, P, n);
        return result;
    }}
    
    static double {module}_K(double T, double P, double n[{number_components}]) {{
        double result = -{module}_dgdp(T, P, n)/{module}_d2gdp2(T, P, n);
        return result;
    }}
    
    static double {module}_Kp(double T, double P, double n[{number_components}]) {{
        double result = {module}_dgdp(T, P, n);
        result *= {module}_d3gdp3(T, P, n);
        result /= pow({module}_d2gdp2(T, P, n), 2.0);
        return result - 1.0;
    }}
    
    \
    """

.. code:: ipython3

    solution_calc += convenience_template.format(module=module, number_components=c)

Write contents to solution_calc.h

.. code:: ipython3

    with open('solution_calc.h', 'w') as f:
        f.write(solution_calc)

(2)
~~~

Function template for primary routines in solution_calib.h

.. code:: ipython3

    solution_calib_template = """\
    
    static double {module}_dparam_{func}(double T, double P, double n[{number_components}], int index) {{
    {moles_assign}
        double result = 0.0;
        switch (index) {{
    {switch_code}
        default:
            break;
        }}
            return result;
    }}
    \
    """

.. code:: ipython3

    print (symparam)


.. parsed-literal::

    (Wh12, Ws12, Wv12, dWh12, dWs12, dWv12, Wh13, Ws13, Wv13, dWh13, dWs13, dWv13, Wh23, Ws23, Wv23, dWh23, dWs23, dWv23, Wh123, Ws123, Wv123)


.. code:: ipython3

    G_param_jac = sym.Matrix(G_matrix).jacobian(symparam)

.. code:: ipython3

    solution_calib = '#include <math.h>\n\n'
    
    for j in range(0,len(f_list)):
        G_jac_list = [ printer.doprint(G_param_jac[j,i]) for i in range(0, len(params)) ]
        
        switch_code_text = ''
        for i in range(0, len(params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'
            switch_code_text += '        break;\n'
            
        solution_calib += solution_calib_template.format(\
                                  module=module, \
                                  number_components=c, \
                                  func=f_list[j], \
                                  switch_code=switch_code_text, \
                                  moles_assign=moles_assign_text \
                                 )

Add convenience routines to solution_calib.h

.. code:: ipython3

    extra_template = """\
    
    static int {module}_get_param_number(void) {{
        return {number_params};
    }}
    
    static const char *paramNames[{number_params}] = {names_params};
    static const char *paramUnits[{number_params}] = {units_params};
    
    static const char **{module}_get_param_names(void) {{
        return paramNames;
    }}
    
    static const char **{module}_get_param_units(void) {{
        return paramUnits;
    }}
    
    static void {module}_get_param_values(double *values) {{
    {code_block_one}
    }}
    
    static int {module}_set_param_values(double *values) {{
    {code_block_two}
        return 1;
    }}
    
    static double {module}_get_param_value(int index) {{
        double result = 0.0;
        switch (index) {{
    {code_block_three}
        default:
            break;
        }}
        return result;
    }}
    
    static int {module}_set_param_value(int index, double value) {{
        int result = 1;
        switch (index) {{
    {code_block_four}
        default:
            result = 0;
            break;
        }}
        return result;
    }}
    
    \
    """

.. code:: ipython3

    value_params=[ printer.doprint(symparam[i]) for i in range(0, len(params)) ]
    code_block_one_text = ''
    code_block_two_text = ''
    code_block_three_text = ''
    code_block_four_text = ''
    for i in range(0,len(value_params)):
        code_block_one_text   += '    values[' + str(i) + '] = ' + value_params[i] + ';\n'
        code_block_two_text   += '    ' + value_params[i] + ' = values[' + str(i) + '];\n'
        code_block_three_text += '    case ' + str(i) + ':\n' \
                               + '        result = ' + value_params[i] + ';\n' \
                               + '        break;\n'
        code_block_four_text  += '    case ' + str(i) + ':\n' \
                               + '        ' + value_params[i] + ' = value;\n' \
                               + '        break;\n'
    
    import json
    solution_calib += extra_template.format(module=module, \
                                          number_params=len(params), \
                                          names_params=json.dumps(params).replace('[', '{').replace(']', '}'), \
                                          units_params=json.dumps(units).replace('[', '{').replace(']', '}'), \
                                          code_block_one=code_block_one_text, \
                                          code_block_two=code_block_two_text, \
                                          code_block_three=code_block_three_text, \
                                          code_block_four=code_block_four_text)

.. code:: ipython3

    with open('solution_calib.h', 'w') as f:
        f.write(solution_calib)

Define characteristics of a Feldspar Solution
=============================================

Components 1. albite, :math:`{\rm{NaAlSi_3O_8}}` 2. anorthite,
:math:`{\rm{CaAl_2Si_2O_8}}` 3. sanidine, :math:`{\rm{KAlSi_3O_8}}`

Original calibration from the Elkins-Grove paper using their notation:

::

    whabor   = 18810.0;  /* joules     */
    wsabor   = 10.3;     /* joules/K   */
    wvabor   = 0.4602;   /* joules/bar */
    whorab   = 27320.0;  /* joules     */
    wsorab   = 10.3;     /* joules/K   */
    wvorab   = 0.3264;   /* joules/bar */
    whaban   = 7924.0;   /* joules     */
    whanab   = 0.0;      /* joules     */
    whoran   = 40317.0;  /* joules     */
    whanor   = 38974.0;  /* joules     */
    wvanor   = -0.1037;  /* joules/bar */
    whabanor = 12545.0;  /* joules     */
    wvabanor = -1.095;   /* joules/bar */

| Note that in the notation derveloped in this paper,
  :math:`\left[ {{W_{ij}} + \Delta {W_{ij}}\left( {{X_i} - {X_j}} \right)} \right]{X_i}{X_j}`,
  is related to the convention used in the original paper:
  :math:`{{\tilde W}_{13}}{X_1}{X_3}\left( {{X_3} + \frac{{{X_2}}}{2}} \right) + {{\tilde W}_{31}}{X_1}{X_3}\left( {{X_1} + \frac{{{X_2}}}{2}} \right)`:
| - :math:`{W_{13}} = {{\tilde W}_{13}} + {{\tilde W}_{31}}` -
  :math:`d{W_{13}} = {{\tilde W}_{31}} - {{\tilde W}_{13}}`

.. code:: ipython3

    print (params)
    whabor   = 18810.0
    wsabor   = 10.3
    wvabor   = 0.4602
    whorab   = 27320.0
    wsorab   = 10.3
    wvorab   = 0.3264
    whaban   = 7924.0
    whanab   = 0.0
    whoran   = 40317.0
    whanor   = 38974.0
    wvanor   = -0.1037
    wvoran   = 0.0 
    whabanor = 12545.0
    wvabanor = -1.095
    paramValues = { 'Wh12':whaban+whanab,  'Ws12':0.0,            'Wv12':0.0, \
                    'Wh13':whabor+whorab,  'Ws13':wsabor+wsorab,  'Wv13':wvabor+wvorab, \
                    'Wh23':whanor+whoran,  'Ws23':0.0,            'Wv23':wvanor+wvoran, \
                    'dWh12':whanab-whaban, 'dWs12':0.0,           'dWv12':0.0, \
                    'dWh13':whorab-whabor, 'dWs13':wsorab-wsabor, 'dWv13':wvorab-wvabor, \
                    'dWh23':whoran-whanor, 'dWs23':0.0,           'dWv23':wvoran-wvanor,
                    'Wh123':whabanor, 'Ws123':0.0, 'Wv123':wvabanor}
    print (paramValues)


.. parsed-literal::

    ['Wh12', 'Ws12', 'Wv12', 'dWh12', 'dWs12', 'dWv12', 'Wh13', 'Ws13', 'Wv13', 'dWh13', 'dWs13', 'dWv13', 'Wh23', 'Ws23', 'Wv23', 'dWh23', 'dWs23', 'dWv23', 'Wh123', 'Ws123', 'Wv123']
    {'Wh12': 7924.0, 'Ws12': 0.0, 'Wv12': 0.0, 'Wh13': 46130.0, 'Ws13': 20.6, 'Wv13': 0.7866, 'Wh23': 79291.0, 'Ws23': 0.0, 'Wv23': -0.1037, 'dWh12': -7924.0, 'dWs12': 0.0, 'dWv12': 0.0, 'dWh13': 8510.0, 'dWs13': 0.0, 'dWv13': -0.13379999999999997, 'dWh23': 1343.0, 'dWs23': 0.0, 'dWv23': 0.1037, 'Wh123': 12545.0, 'Ws123': 0.0, 'Wv123': -1.095}


Create templates for instance of an asymmetric regular solution phase
---------------------------------------------------------------------

.. code:: ipython3

    fast_c_template = """\
    
    static const int identifier = {git_identifier};
    {code_block_one}
    
    #include "solution_calc.h"
    
    const int {phase}_{module}_identifier(void) {{
        return identifier;
    }}
    
    const char *{phase}_{module}_name(void) {{
        return "{phase}";
    }}
    
    double {phase}_{module}_g(double T, double P, double n[{number_components}]) {{
        return {module}_g(T, P, n);
    }}
    
    void {phase}_{module}_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_dgdn(T, P, n, result);
    }}
    
    void {phase}_{module}_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d2gdn2(T, P, n, result);
    }}
    
    void {phase}_{module}_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdn3(T, P, n, result);
    }}
    
    double {phase}_{module}_dgdt(double T, double P, double n[{number_components}]) {{
        return {module}_dgdt(T, P, n);
    }}
    
    void {phase}_{module}_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d2gdndt(T, P, n, result);
    }}
    
    void {phase}_{module}_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdn2dt(T, P, n, result);
    }}
    
    void {phase}_{module}_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn3dt(T, P, n, result);
    }}
    
    double {phase}_{module}_dgdp(double T, double P, double n[{number_components}]) {{
        return {module}_dgdp(T, P, n);
    }}
    
    void {phase}_{module}_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d2gdndp(T, P, n, result);
    }}
    
    void {phase}_{module}_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdn2dp(T, P, n, result);
    }}
    
    void {phase}_{module}_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn3dp(T, P, n, result);
    }}
    
    double {phase}_{module}_d2gdt2(double T, double P, double n[{number_components}]) {{
        return {module}_d2gdt2(T, P, n);
    }}
    
    void {phase}_{module}_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdndt2(T, P, n, result);
    }}
    
    void {phase}_{module}_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn2dt2(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn3dt2(T, P, n, result);
    }}
    
    double {phase}_{module}_d2gdtdp(double T, double P, double n[{number_components}]) {{
        return {module}_d2gdtdp(T, P, n);
    }}
    
    void {phase}_{module}_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdndtdp(T, P, n, result);
    }}
    
    void {phase}_{module}_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn2dtdp(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn3dtdp(T, P, n, result);
    }}
    
    double {phase}_{module}_d2gdp2(double T, double P, double n[{number_components}]) {{
        return {module}_d2gdp2(T, P, n);
    }}
    
    void {phase}_{module}_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdndp2(T, P, n, result);
    }}
    
    void {phase}_{module}_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn2dp2(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn3dp2(T, P, n, result);
    }}
    
    double {phase}_{module}_d3gdt3(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdt3(T, P, n);
    }}
    
    void {phase}_{module}_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndt3(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dt3(T, P, n, result);
    }}
    
    void {phase}_{module}_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dt3(T, P, n, result);
    }}
    
    double {phase}_{module}_d3gdt2dp(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdt2dp(T, P, n);
    }}
    
    void {phase}_{module}_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndt2dp(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dt2dp(T, P, n, result);
    }}
    
    void {phase}_{module}_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dt2dp(T, P, n, result);
    }}
    
    double {phase}_{module}_d3gdtdp2(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdtdp2(T, P, n);
    }}
    
    void {phase}_{module}_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndtdp2(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dtdp2(T, P, n, result);
    }}
    
    void {phase}_{module}_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dtdp2(T, P, n, result);
    }}
    
    double {phase}_{module}_d3gdp3(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdp3(T, P, n);
    }}
    
    void {phase}_{module}_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndp3(T, P, n, result);
    }}
    
    void {phase}_{module}_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dp3(T, P, n, result);
    }}
    
    void {phase}_{module}_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dp3(T, P, n, result);
    }}
    
    double {phase}_{module}_s(double T, double P, double n[{number_components}]) {{
        return {module}_s(T, P, n);
    }}
    
    double {phase}_{module}_v(double T, double P, double n[{number_components}]) {{
        return {module}_v(T, P, n);
    }}
    
    double {phase}_{module}_cv(double T, double P, double n[{number_components}]) {{
        return {module}_cv(T, P, n);
    }}
    
    double {phase}_{module}_cp(double T, double P, double n[{number_components}]) {{
        return {module}_cp(T, P, n);
    }}
    
    double {phase}_{module}_dcpdt(double T, double P, double n[{number_components}]) {{
        return {module}_dcpdt(T, P, n);
    }}
    
    double {phase}_{module}_alpha(double T, double P, double n[{number_components}]) {{
        return {module}_alpha(T, P, n);
    }}
    
    double {phase}_{module}_beta(double T, double P, double n[{number_components}]) {{
        return {module}_beta(T, P, n);
    }}
    
    double {phase}_{module}_K(double T, double P, double n[{number_components}]) {{
        return {module}_K(T, P, n);
    }}
    
    double {phase}_{module}_Kp(double T, double P, double n[{number_components}]) {{
        return {module}_Kp(T, P, n);
    }}
    
    \
    """

.. code:: ipython3

    calib_c_template = """\
    
    static int identifier = {git_identifier};
    {code_block_two}
    
    #include "solution_calc.h"
    #include "solution_calib.h"
    
    const int {phase}_{module}_calib_identifier(void) {{
        return identifier;
    }}
    
    const char *{phase}_{module}_calib_name(void) {{
        return "{phase}";
    }}
    
    double {phase}_{module}_calib_g(double T, double P, double n[{number_components}]) {{
        return {module}_g(T, P, n);
    }}
    
    void {phase}_{module}_calib_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_dgdn(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d2gdn2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdn3(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_dgdt(double T, double P, double n[{number_components}]) {{
        return {module}_dgdt(T, P, n);
    }}
    
    void {phase}_{module}_calib_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d2gdndt(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdn2dt(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn3dt(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_dgdp(double T, double P, double n[{number_components}]) {{
        return {module}_dgdp(T, P, n);
    }}
    
    void {phase}_{module}_calib_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d2gdndp(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdn2dp(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn3dp(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d2gdt2(double T, double P, double n[{number_components}]) {{
        return {module}_d2gdt2(T, P, n);
    }}
    
    void {phase}_{module}_calib_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdndt2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn2dt2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn3dt2(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d2gdtdp(double T, double P, double n[{number_components}]) {{
        return {module}_d2gdtdp(T, P, n);
    }}
    
    void {phase}_{module}_calib_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdndtdp(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn2dtdp(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn3dtdp(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d2gdp2(double T, double P, double n[{number_components}]) {{
        return {module}_d2gdp2(T, P, n);
    }}
    
    void {phase}_{module}_calib_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d3gdndp2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdn2dp2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn3dp2(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d3gdt3(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdt3(T, P, n);
    }}
    
    void {phase}_{module}_calib_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndt3(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dt3(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dt3(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d3gdt2dp(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdt2dp(T, P, n);
    }}
    
    void {phase}_{module}_calib_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndt2dp(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dt2dp(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dt2dp(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d3gdtdp2(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdtdp2(T, P, n);
    }}
    
    void {phase}_{module}_calib_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndtdp2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dtdp2(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dtdp2(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_d3gdp3(double T, double P, double n[{number_components}]) {{
        return {module}_d3gdp3(T, P, n);
    }}
    
    void {phase}_{module}_calib_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d4gdndp3(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d5gdn2dp3(T, P, n, result);
    }}
    
    void {phase}_{module}_calib_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
        {module}_d6gdn3dp3(T, P, n, result);
    }}
    
    double {phase}_{module}_calib_s(double T, double P, double n[{number_components}]) {{
        return {module}_s(T, P, n);
    }}
    
    double {phase}_{module}_calib_v(double T, double P, double n[{number_components}]) {{
        return {module}_v(T, P, n);
    }}
    
    double {phase}_{module}_calib_cv(double T, double P, double n[{number_components}]) {{
        return {module}_cv(T, P, n);
    }}
    
    double {phase}_{module}_calib_cp(double T, double P, double n[{number_components}]) {{
        return {module}_cp(T, P, n);
    }}
    
    double {phase}_{module}_calib_dcpdt(double T, double P, double n[{number_components}]) {{
        return {module}_dcpdt(T, P, n);
    }}
    
    double {phase}_{module}_calib_alpha(double T, double P, double n[{number_components}]) {{
        return {module}_alpha(T, P, n);
    }}
    
    double {phase}_{module}_calib_beta(double T, double P, double n[{number_components}]) {{
        return {module}_beta(T, P, n);
    }}
    
    double {phase}_{module}_calib_K(double T, double P, double n[{number_components}]) {{
        return {module}_K(T, P, n);
    }}
    
    double {phase}_{module}_calib_Kp(double T, double P, double n[{number_components}]) {{
        return {module}_Kp(T, P, n);
    }}
    
    int {phase}_{module}_get_param_number(void) {{
        return {module}_get_param_number();
    }}
    
    const char **{phase}_{module}_get_param_names(void) {{
        return {module}_get_param_names();
    }}
    
    const char **{phase}_{module}_get_param_units(void) {{
        return {module}_get_param_units();
    }}
    
    void {phase}_{module}_get_param_values(double *values) {{
        {module}_get_param_values(values);
    }}
    
    int {phase}_{module}_set_param_values(double *values) {{
        return {module}_set_param_values(values);
    }}
    
    double {phase}_{module}_get_param_value(int index) {{
        return {module}_get_param_value(index);
    }}
    
    int {phase}_{module}_set_param_value(int index, double value) {{
        return {module}_set_param_value(index, value);
    }}
    
    double {phase}_{module}_dparam_g(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_g(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_dgdt(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_dgdt(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_dgdp(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_dgdp(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d2gdt2(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d2gdt2(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d2gdtdp(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d2gdtdp(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d2gdp2(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d2gdp2(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d3gdt3(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d3gdt3(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d3gdt2dp(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d3gdt2dp(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d3gdtdp2(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d3gdtdp2(T, P, n, index);
    }}
    
    double {phase}_{module}_dparam_d3gdp3(double T, double P, double n[{number_components}], int index) {{
        return {module}_dparam_d3gdp3(T, P, n, index);
    }}
    
    \
    """

.. code:: ipython3

    fast_h_template = """\
    
    const char *{phase}_{module}_name(void);
    
    double {phase}_{module}_g(double T, double P, double n[{number_components}]);
    double {phase}_{module}_dgdt(double T, double P, double n[{number_components}]);
    double {phase}_{module}_dgdp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d2gdt2(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d2gdtdp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d2gdp2(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d3gdt3(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d3gdt2dp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d3gdtdp2(double T, double P, double n[{number_components}]);
    double {phase}_{module}_d3gdp3(double T, double P, double n[{number_components}]);
    
    void {phase}_{module}_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]);
    
    void {phase}_{module}_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]);
    
    void {phase}_{module}_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]);
    
    double {phase}_{module}_s(double T, double P, double n[{number_components}]);
    double {phase}_{module}_v(double T, double P, double n[{number_components}]);
    double {phase}_{module}_cv(double T, double P, double n[{number_components}]);
    double {phase}_{module}_cp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_dcpdt(double T, double P, double n[{number_components}]);
    double {phase}_{module}_alpha(double T, double P, double n[{number_components}]);
    double {phase}_{module}_beta(double T, double P, double n[{number_components}]);
    double {phase}_{module}_K(double T, double P, double n[{number_components}]);
    double {phase}_{module}_Kp(double T, double P, double n[{number_components}]);
    
    \
    """

.. code:: ipython3

    calib_h_template = """\
    
    const char *{phase}_{module}_calib_name(void);
    
    double {phase}_{module}_calib_g(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_dgdt(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_dgdp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d2gdt2(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d2gdtdp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d2gdp2(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d3gdt3(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d3gdt2dp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d3gdtdp2(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_d3gdp3(double T, double P, double n[{number_components}]);
    
    void {phase}_{module}_calib_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]);
    
    void {phase}_{module}_calib_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]);
    
    void {phase}_{module}_calib_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_calib_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]);
    
    double {phase}_{module}_calib_s(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_v(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_cv(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_cp(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_dcpdt(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_alpha(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_beta(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_K(double T, double P, double n[{number_components}]);
    double {phase}_{module}_calib_Kp(double T, double P, double n[{number_components}]);
    
    
    int {phase}_{module}_get_param_number(void);
    const char *{phase}_{module}_get_param_names(void);
    const char *{phase}_{module}_get_param_units(void);
    void {phase}_{module}_get_param_values(double *values);
    int {phase}_{module}_set_param_values(double *values);
    double {phase}_{module}_get_param_value(int index);
    int {phase}_{module}_set_param_value(int index, double value);
    
    double {phase}_{module}_dparam_g(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_dgdt(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_dgdp(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d2gdt2(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d2gdtdp(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d2gdp2(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdt3(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdt2dp(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdtdp2(double T, double P, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdp3(double T, double P, double n[{number_components}], int index);
    
    \
    """

::

    typedef struct _endmembers {
      const char *(*name) (void);
      const char *(*formula) (void);
      const double (*mw) (void);
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
        Albite_berman_name,
        Albite_berman_formula,
        Albite_berman_mw,
        Albite_berman_g,
        Albite_berman_dgdt,
        Albite_berman_dgdp,
        Albite_berman_d2gdt2,
        Albite_berman_d2gdtdp,
        Albite_berman_d2gdp2,
        Albite_berman_d3gdt3,
        Albite_berman_d3gdt2dp,
        Albite_berman_d3gdtdp2,
        Albite_berman_d3gdp3
      },
      {
        Anorthite_berman_name,
        Anorthite_berman_formula,
        Anorthite_berman_mw,
        Anorthite_berman_g,
        Anorthite_berman_dgdt,
        Anorthite_berman_dgdp,
        Anorthite_berman_d2gdt2,
        Anorthite_berman_d2gdtdp,
        Anorthite_berman_d2gdp2,
        Anorthite_berman_d3gdt3,
        Anorthite_berman_d3gdt2dp,
        Anorthite_berman_d3gdtdp2,
        Anorthite_berman_d3gdp3
      },
      {
        Sanidine_berman_name,
        Sanidine_berman_formula,
        Sanidine_berman_mw,
        Sanidine_berman_g,
        Sanidine_berman_dgdt,
        Sanidine_berman_dgdp,
        Sanidine_berman_d2gdt2,
        Sanidine_berman_d2gdtdp,
        Sanidine_berman_d2gdp2,
        Sanidine_berman_d3gdt3,
        Sanidine_berman_d3gdt2dp,
        Sanidine_berman_d3gdtdp2,
        Sanidine_berman_d3gdp3
      }
    };
    static int nc = (sizeof endmember / sizeof(_endmembers));

.. code:: ipython3

    std_state_h_template = """\
    
    {code_block_four}
    
    typedef struct _endmembers {{
      const char *(*name) (void);
      const char *(*formula) (void);
      const double (*mw) (void);
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
    }} Endmembers;
    
    static Endmembers endmember[] = {{
    {code_block_three}
    }};
    static int nc = (sizeof endmember / sizeof(struct _endmembers));
    
    static const double R=8.3143;
    
    \
    """

Generate both fast computation and calibibration code for the feldspar
solution

.. code:: ipython3

    phase = 'Feldspar'
    endmembers = ['Albite_berman', 'Anorthite_berman', 'Potassium_Feldspar_berman']
    feldspar_phase_c = ''
    feldspar_phase_h = ''
    
    code_block_three_text = ''
    code_block_four_text = ''
    for i in range(0,len(endmembers)):
        code_block_three_text += '  {\n'
        code_block_three_text += '    ' + endmembers[i] + '_name,\n'
        code_block_three_text += '    ' + endmembers[i] + '_formula,\n'
        code_block_three_text += '    ' + endmembers[i] + '_mw,\n'
        code_block_three_text += '    ' + endmembers[i] + '_g,\n'
        code_block_three_text += '    ' + endmembers[i] + '_dgdt,\n'
        code_block_three_text += '    ' + endmembers[i] + '_dgdp,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d2gdt2,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d2gdtdp,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d2gdp2,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d3gdt3,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d3gdt2dp,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d3gdtdp2,\n'
        code_block_three_text += '    ' + endmembers[i] + '_d3gdp3\n'
        code_block_three_text += '  },\n'
        code_block_four_text += '#include "' + endmembers[i] + '.h"\n'
    
    code_block_one_text = ''
    code_block_two_text = ''
    for i in range(0, len(params)):
        symbol = printer.doprint(symparam[i])
        value = str(paramValues[params[i]])
        code_block_one_text += 'static const double ' + symbol + ' = ' + value + ';\n'
        code_block_two_text += 'static double ' + symbol + ' = ' + value + ';\n'
    
    print(code_block_four_text) 
        
    feldspar_phase_c += std_state_h_template.format( \
            code_block_three=code_block_three_text, \
            code_block_four=code_block_four_text, \
            )
    feldspar_phase_c += fast_c_template.format( \
            module=module, \
            phase=phase, \
            param_names=[ printer.doprint(symparam[i]) for i in range(0, len(params)) ], \
            param_values=paramValues, \
            git_identifier=1, \
            number_components=c, \
            code_block_one=code_block_one_text \
            )
    feldspar_phase_h += fast_h_template.format(module=module, phase=phase, number_components=c)
    with open(phase + '_solution.h', 'w') as f:
        f.write(feldspar_phase_h)
    with open(phase + '_solution.c', 'w') as f:
        f.write(feldspar_phase_c)
    feldspar_phase_calib_c = ''
    feldspar_phase_calib_h = ''
    feldspar_phase_calib_c += std_state_h_template.format( \
            code_block_three=code_block_three_text, \
            code_block_four=code_block_four_text, \
            )
    feldspar_phase_calib_c += calib_c_template.format( \
            module=module, \
            phase=phase, \
            param_names=[ printer.doprint(symparam[i]) for i in range(0, len(params)) ], \
            param_values=paramValues, \
            git_identifier=1, \
            number_components=c, \
            code_block_two=code_block_two_text \
            )
    feldspar_phase_calib_h += calib_h_template.format(module=module, phase=phase, number_components=c)
    with open(phase + '_solution_calib.h', 'w') as f:
        f.write(feldspar_phase_calib_h)
    with open(phase + '_solution_calib.c', 'w') as f:
        f.write(feldspar_phase_calib_c)


.. parsed-literal::

    #include "Albite_berman.h"
    #include "Anorthite_berman.h"
    #include "Potassium_Feldspar_berman.h"
    


Load the Cython Jupyter magic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    %load_ext cython

Generate the C-code and python extension

::

    setuptools.extension.Extension
    self, name, sources, include_dirs=None, define_macros=None, undef_macros=None, library_dirs=None, libraries=None, runtime_library_dirs=None, extra_objects=None, extra_compile_args=None, extra_link_args=None, export_symbols=None, swig_opts=None, depends=None, language=None, optional=None, **kw

::

    -O0, -O1, -O2, -O3, -Ofast, -Os, -Oz, -Og, -O, -O4
    Specify which optimization level to use:

    -O0 Means “no optimization”: this level compiles the fastest and generates the most debuggable code.

    -O1 Somewhere between -O0 and -O2.

    -O2 Moderate level of optimization which enables most optimizations.

    -O3 Like -O2, except that it enables optimizations that take longer to perform or that may generate larger code (in an attempt to make the program run faster).

    -Ofast Enables all the optimizations from -O3 along with other aggressive optimizations that may violate strict compliance with language standards.

    -Os Like -O2 with extra optimizations to reduce code size.

    -Oz Like -Os (and thus -O2), but reduces code size further.

    -Og Like -O1. In future versions, this option might disable different optimizations in order to improve debuggability.

    -O Equivalent to -O2.

    -O4 and higher

.. code:: ipython3

    %%writefile feldspar.pyxbld
    import numpy
    
    #            module name specified by `%%cython_pyximport` magic
    #            |        just `modname + ".pyx"`
    #            |        |
    def make_ext(modname, pyxfilename):
        from setuptools.extension import Extension
        return Extension(modname,
                         sources=[pyxfilename, \
                                  'Feldspar_solution.c', 'Feldspar_solution_calib.c', \
                                  'Albite_berman.c', 'Albite_berman_calib.c', \
                                  'Anorthite_berman.c', 'Anorthite_berman_calib.c', \
                                  'Potassium_Feldspar_berman.c', 'Potassium_Feldspar_berman_calib.c', \
                                 ],
                         include_dirs=['.', numpy.get_include()], extra_compile_args=['-O3'])


.. parsed-literal::

    Overwriting feldspar.pyxbld


.. code:: ipython3

    %%cython_pyximport feldspar
    import numpy as np
    cimport numpy as cnp # cimport gives us access to NumPy's C API
    
    # here we just replicate the function signature from the header
    cdef extern from "Feldspar_solution.h":
        double Feldspar_asymm_regular_g(double t, double p, double n[3])
    cdef extern from "Feldspar_solution_calib.h":
        double Feldspar_asymm_regular_calib_g(double t, double p, double n[3])
    
    # here is the "wrapper" signature
    def cy_Feldspar_asymm_regular_g(double t, double p, cnp.ndarray[cnp.double_t,ndim=1] n):
        result = Feldspar_asymm_regular_g(<double> t, <double> p, <double*> n.data)
        return result
    
    def cy_Feldspar_asymm_regular_calib_g(double t, double p, cnp.ndarray[cnp.double_t,ndim=1] n):
        result = Feldspar_asymm_regular_calib_g(<double> t, <double> p, <double*> n.data)
        return result

.. code:: ipython3

    %cd ..


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen


Test and time the generated functions for Feldspar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    t = 2000.00
    p = 1.0
    n = np.array([1.0e-12, 1.0, 1.0])

.. code:: ipython3

    newCodeFast = cy_Feldspar_asymm_regular_g(t, p, n)
    newCodeCalib = cy_Feldspar_asymm_regular_calib_g(t, p, n)
    print (newCodeFast, newCodeCalib, newCodeFast-newCodeCalib)


.. parsed-literal::

    -10194212.10545563 -10194212.10545563 0.0


Instantiate the Objective-C code via Rubicon

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))
    FeldsparRaw = ObjCClass('FeldsparBerman')
    obj = FeldsparRaw.alloc().init()
    print (obj.phaseName)
    print (obj.phaseFormula)


.. parsed-literal::

    Feldspar
    


.. code:: ipython3

    import ctypes
    m = (ctypes.c_double*3)()
    ctypes.cast(m, ctypes.POINTER(ctypes.c_double))
    m[0] = 1.0e-12
    m[1] = 1.0
    m[2] = 1.0
    oldCode = obj.getGibbsFreeEnergyFromMolesOfComponents_andT_andP_(m, t, p)
    print (oldCode, oldCode-newCodeFast)


.. parsed-literal::

    -10210751.991823655 -16539.886368025094


.. code:: ipython3

    AnCorr = 3.7*4184.0 - t*3.7*4184.0/2200.0
    SnCorr = 3400.0
    print (n[1]*AnCorr+n[2]*SnCorr)


.. parsed-literal::

    4807.345454545455


Time the fast code

.. code:: ipython3

    %timeit(cy_Feldspar_asymm_regular_g(t, p, n))


.. parsed-literal::

    666 ns ± 10.5 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)


Time the calibration code

.. code:: ipython3

    %timeit(cy_Feldspar_asymm_regular_calib_g(t, p, n))


.. parsed-literal::

    678 ns ± 9.27 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)


Time the Rubicon wrapped Objective-C code

.. code:: ipython3

    %timeit(obj.getGibbsFreeEnergyFromMolesOfComponents_andT_andP_(m, t, p))


.. parsed-literal::

    52.4 µs ± 2.16 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)

