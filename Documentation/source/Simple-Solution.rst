
Simple Solution SymPy Code Generation
=====================================

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()
    from thermoengine import coder

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

Create a simple solution model
------------------------------

… with the specified number of endmember thermodynamic components

.. code:: ipython3

    model = coder.SimpleSolnModel(nc=c)

Retrieve primary compositional variables
----------------------------------------

-  :math:`n` is a vector of mole numbers of each component
-  :math:`n_T` is the total number of moles in the solution ### and
   construct a derived mole fraction variable
-  :math:`X` is a vector of mole fractions of components in the system

.. code:: ipython3

    n = model.n
    nT = model.nT
    X = n/nT
    n, nT, X




.. math::

    \left ( \left[\begin{matrix}n_{1}\\n_{2}\\n_{3}\end{matrix}\right], \quad n_{1} + n_{2} + n_{3}, \quad \left[\begin{matrix}\frac{n_{1}}{n_{1} + n_{2} + n_{3}}\\\frac{n_{2}}{n_{1} + n_{2} + n_{3}}\\\frac{n_{3}}{n_{1} + n_{2} + n_{3}}\end{matrix}\right]\right )



Retrieve the temperature, pressure, and standard state chemical potentials
--------------------------------------------------------------------------

-  :math:`T` is temperature in :math:`K`
-  :math:`P` is pressure in :math:`bars`
-  :math:`\mu` in Joules

.. code:: ipython3

    T = model.get_symbol_for_t()
    P = model.get_symbol_for_p()
    mu = model.mu
    T,P,mu




.. math::

    \left ( T, \quad P, \quad \left[\begin{matrix}\mu_{1}{\left (T,P \right )}\\\mu_{2}{\left (T,P \right )}\\\mu_{3}{\left (T,P \right )}\end{matrix}\right]\right )



Define the standard state contribution to solution properties
-------------------------------------------------------------

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



Add the Gibbs free energy of solution to the model
--------------------------------------------------

.. code:: ipython3

    model.add_expression_to_model(G, list(zip(params, units, symparam)))

… give the model a unqiue name

.. code:: ipython3

    model.module = "Simple_Solution"

Define Parameters of a Feldspar Solution
========================================

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
                    'Wh123':whabanor, 'Ws123':0.0, 'Wv123':wvabanor, 'T_r':298.15, 'P_r':1.0}
    print (paramValues)


.. parsed-literal::

    ['Wh12', 'Ws12', 'Wv12', 'dWh12', 'dWs12', 'dWv12', 'Wh13', 'Ws13', 'Wv13', 'dWh13', 'dWs13', 'dWv13', 'Wh23', 'Ws23', 'Wv23', 'dWh23', 'dWs23', 'dWv23', 'Wh123', 'Ws123', 'Wv123']
    {'Wh12': 7924.0, 'Ws12': 0.0, 'Wv12': 0.0, 'Wh13': 46130.0, 'Ws13': 20.6, 'Wv13': 0.7866, 'Wh23': 79291.0, 'Ws23': 0.0, 'Wv23': -0.1037, 'dWh12': -7924.0, 'dWs12': 0.0, 'dWv12': 0.0, 'dWh13': 8510.0, 'dWs13': 0.0, 'dWv13': -0.13379999999999997, 'dWh23': 1343.0, 'dWs23': 0.0, 'dWv23': 0.1037, 'Wh123': 12545.0, 'Ws123': 0.0, 'Wv123': -1.095, 'T_r': 298.15, 'P_r': 1.0}


Generate both fast computation and calibibration code for the feldspar
solution

Use code printers to construct “C” package code
===============================================

.. code:: ipython3

    model_working_dir = "working"
    !mkdir -p {model_working_dir}
    %cd {model_working_dir}


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working


Choose model type and create model
----------------------------------

model_type is “fast” or “calib”

.. code:: ipython3

    model_type = "calib"

.. code:: ipython3

    model.create_code_module(phase="Feldspar", params=paramValues, 
                             endmembers=['High_Albite_berman', 'Anorthite_berman', 'Potassium_Feldspar_berman'], 
                             prefix="cy", module_type=model_type, silent=False)


.. parsed-literal::

    Creating generic fast model code file string
    Writing include file to working directory ...
    Creating (once only) generic model calib codetemplate include file string
    Creating (once only) generic model calib codetemplate code file string
    Creating generic calib model code file string
    Writing include file to working directory ...
    Creating code blocks for standard state properties.
    Creating calib code and include files ...
    Writing include file to working directory ...
    Writing code file to working directory ...
    ... done
    Writing pyxbld file to working directory ...
    writing pyx file to working directory ...
    Compiling code and Python bindings ...
    Success! Import the module named  Simple_Solution




.. parsed-literal::

    True



Load the module
---------------

.. code:: ipython3

    import Simple_Solution
    %cd ..


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen


Test and time the generated functions for Feldspar (T in K, P in bars)
----------------------------------------------------------------------

.. code:: ipython3

    t = 2000.00
    p = 1.0
    n = np.array([1.1, 1.2, 1.3])

Available in both “Fast” and “Calib” code versions
--------------------------------------------------

Execute the “fast” or “calibration” code metadata retrieval functions:

.. code:: ipython3

    try:
        print(Simple_Solution.cy_Feldspar_Simple_Solution_identifier())
        print(Simple_Solution.cy_Feldspar_Simple_Solution_name())
    except AttributeError:
        pass
    try:
        print(Simple_Solution.cy_Feldspar_Simple_Solution_calib_identifier())
        print(Simple_Solution.cy_Feldspar_Simple_Solution_calib_name())
    except AttributeError:
        pass


.. parsed-literal::

    Thu Sep 20 12:25:55 2018
    Feldspar


Execute the standard thermodynamic property retrieval functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
    try:
        print(fmt.format('G', Simple_Solution.cy_Feldspar_Simple_Solution_g(t,p,n), 'J'))
        print(fmt.format('dGdT', Simple_Solution.cy_Feldspar_Simple_Solution_dgdt(t,p,n), 'J/K'))
        print(fmt.format('dGdP', Simple_Solution.cy_Feldspar_Simple_Solution_dgdp(t,p,n), 'J/bar'))
        print(fmt.format('d2GdT2', Simple_Solution.cy_Feldspar_Simple_Solution_d2gdt2(t,p,n), 'J/K^2'))
        print(fmt.format('d2GdTdP', Simple_Solution.cy_Feldspar_Simple_Solution_d2gdtdp(t,p,n), 'J/K-bar'))
        print(fmt.format('d2GdP2', Simple_Solution.cy_Feldspar_Simple_Solution_d2gdp2(t,p,n), 'J/bar^2'))
        print(fmt.format('d3GdT3', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdt3(t,p,n), 'J/K^3'))
        print(fmt.format('d3GdT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdt2dp(t,p,n), 'J/K^2-bar'))
        print(fmt.format('d3GdTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdtdp2(t,p,n), 'J/K-bar^2'))
        print(fmt.format('d3GdP3', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdp3(t,p,n), 'J/bar^3'))
        print(fmt.format('S', Simple_Solution.cy_Feldspar_Simple_Solution_s(t,p,n), 'J/K'))
        print(fmt.format('V', Simple_Solution.cy_Feldspar_Simple_Solution_v(t,p,n), 'J/bar'))
        print(fmt.format('Cv', Simple_Solution.cy_Feldspar_Simple_Solution_cv(t,p,n), 'J/K'))
        print(fmt.format('Cp', Simple_Solution.cy_Feldspar_Simple_Solution_cp(t,p,n), 'J/K'))
        print(fmt.format('dCpdT', Simple_Solution.cy_Feldspar_Simple_Solution_dcpdt(t,p,n), 'J/K^2'))
        print(fmt.format('alpha', Simple_Solution.cy_Feldspar_Simple_Solution_alpha(t,p,n), '1/K'))
        print(fmt.format('beta', Simple_Solution.cy_Feldspar_Simple_Solution_beta(t,p,n), '1/bar'))
        print(fmt.format('K', Simple_Solution.cy_Feldspar_Simple_Solution_K(t,p,n), 'bar'))
        print(fmt.format('Kp', Simple_Solution.cy_Feldspar_Simple_Solution_Kp(t,p,n), ''))
    except AttributeError:
        pass
    try:
        print(fmt.format('G', Simple_Solution.cy_Feldspar_Simple_Solution_calib_g(t,p,n), 'J'))
        print(fmt.format('dGdT', Simple_Solution.cy_Feldspar_Simple_Solution_calib_dgdt(t,p,n), 'J/K'))
        print(fmt.format('dGdP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_dgdp(t,p,n), 'J/bar'))
        print(fmt.format('d2GdT2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d2gdt2(t,p,n), 'J/K^2'))
        print(fmt.format('d2GdTdP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d2gdtdp(t,p,n), 'J/K-bar'))
        print(fmt.format('d2GdP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d2gdp2(t,p,n), 'J/bar^2'))
        print(fmt.format('d3GdT3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdt3(t,p,n), 'J/K^3'))
        print(fmt.format('d3GdT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdt2dp(t,p,n), 'J/K^2-bar'))
        print(fmt.format('d3GdTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdtdp2(t,p,n), 'J/K-bar^2'))
        print(fmt.format('d3GdP3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdp3(t,p,n), 'J/bar^3'))
        print(fmt.format('S', Simple_Solution.cy_Feldspar_Simple_Solution_calib_s(t,p,n), 'J/K'))
        print(fmt.format('V', Simple_Solution.cy_Feldspar_Simple_Solution_calib_v(t,p,n), 'J/bar'))
        print(fmt.format('Cv', Simple_Solution.cy_Feldspar_Simple_Solution_calib_cv(t,p,n), 'J/K'))
        print(fmt.format('Cp', Simple_Solution.cy_Feldspar_Simple_Solution_calib_cp(t,p,n), 'J/K'))
        print(fmt.format('dCpdT', Simple_Solution.cy_Feldspar_Simple_Solution_calib_dcpdt(t,p,n), 'J/K^2'))
        print(fmt.format('alpha', Simple_Solution.cy_Feldspar_Simple_Solution_calib_alpha(t,p,n), '1/K'))
        print(fmt.format('beta', Simple_Solution.cy_Feldspar_Simple_Solution_calib_beta(t,p,n), '1/bar'))
        print(fmt.format('K', Simple_Solution.cy_Feldspar_Simple_Solution_calib_K(t,p,n), 'bar'))
        print(fmt.format('Kp', Simple_Solution.cy_Feldspar_Simple_Solution_calib_Kp(t,p,n), ''))
    except AttributeError:
        pass


.. parsed-literal::

    G          -1.821595e+07 J         
    dGdT       -2.816363e+03 J/K       
    dGdP        3.903068e+01 J/bar     
    d2GdT2     -6.171534e-01 J/K^2     
    d2GdTdP     1.195927e-03 J/K-bar   
    d2GdP2     -6.244961e-05 J/bar^2   
    d3GdT3      2.788296e-04 J/K^3     
    d3GdT2dP    3.279670e-07 J/K^2-bar 
    d3GdTdP2    0.000000e+00 J/K-bar^2 
    d3GdP3      0.000000e+00 J/bar^3   
    S           2.816363e+03 J/K       
    V           3.903068e+01 J/bar     
    Cv          1.188502e+03 J/K       
    Cp          1.234307e+03 J/K       
    dCpdT       5.949430e-02 J/K^2     
    alpha       3.064070e-05 1/K       
    beta        1.600013e-06 1/bar     
    K           6.249948e+05 bar       
    Kp         -1.000000e+00           


Execute functions for molar derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First derivative vectors:
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    def printResult(name, result, units):
        print ("{0:<10.10s}".format(name), end=' ')
        [print ("{0:13.6e}".format(x), end=' ') for x in result]
        print ("{0:<10.10s}".format(units))
    def printLabels(n):
        print ("{0:<18.18s}".format(''), end=' ')
        [print ("[{0:3d}]{1:<8.8s}".format(idx, ''), end=' ') for idx in range(len(n))]
        print ()
    printLabels(n)
    try:
        printResult('dGdn', Simple_Solution.cy_Feldspar_Simple_Solution_dgdn(t,p,n), 'J/m')
        printResult('d2GdndT', Simple_Solution.cy_Feldspar_Simple_Solution_d2gdndt(t,p,n), 'J/K-m')
        printResult('d2GdndP', Simple_Solution.cy_Feldspar_Simple_Solution_d2gdndp(t,p,n), 'J/bar-m')
        printResult('d3GdndT2', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdndt2(t,p,n), 'J/K^2-m')
        printResult('d3GdndTdP', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdndtdp(t,p,n), 'J/K-bar-m')
        printResult('d3GdndP2', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdndp2(t,p,n), 'J/bar^2-m')
        printResult('d4GdndT3', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdndt3(t,p,n), 'J/K^3-m')
        printResult('d4GdndT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdndt2dp(t,p,n), 'J/K^2-bar-m')
        printResult('d4GdndTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdndtdp2(t,p,n), 'J/K-bar^2-m')
        printResult('d4GdndP3', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdndp3(t,p,n), 'J/bar^3-m')
    except AttributeError:
        pass
    try:
        printResult('dGdn', Simple_Solution.cy_Feldspar_Simple_Solution_calib_dgdn(t,p,n), 'J/m')
        printResult('d2GdndT', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d2gdndt(t,p,n), 'J/K-m')
        printResult('d2GdndP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d2gdndp(t,p,n), 'J/bar-m')
        printResult('d3GdndT2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdndt2(t,p,n), 'J/K^2-m')
        printResult('d3GdndTdP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdndtdp(t,p,n), 'J/K-bar-m')
        printResult('d3GdndP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdndp2(t,p,n), 'J/bar^2-m')
        printResult('d4GdndT3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdndt3(t,p,n), 'J/K^3-m')
        printResult('d4GdndT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdndt2dp(t,p,n), 'J/K^2-bar-m')
        printResult('d4GdndTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdndtdp2(t,p,n), 'J/K-bar^2-m')
        printResult('d4GdndP3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdndp3(t,p,n), 'J/bar^3-m')
    except AttributeError:
        pass    


.. parsed-literal::

                       [  0]         [  1]         [  2]         
    dGdn       -4.970230e+06 -5.218217e+06 -4.989876e+06 J/m       
    d2GdndT    -7.868998e+02 -7.738057e+02 -7.863130e+02 J/K-m     
    d2GdndP     1.077333e+01  1.024322e+01  1.145242e+01 J/bar-m   
    d3GdndT2   -1.688921e-01 -1.779158e-01 -1.675947e-01 J/K^2-m   
    d3GdndTdP   3.750779e-04  2.540274e-04  3.680836e-04 J/K-bar-m 
    d3GdndP2   -1.960841e-05 -1.281943e-05 -1.961311e-05 J/bar^2-m 
    d4GdndT3    7.680829e-05  7.849093e-05  7.703949e-05 J/K^3-m   
    d4GdndT2dP  6.453120e-08  8.463000e-08  1.195590e-07 J/K^2-bar-
    d4GdndTdP2  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar^2-


The Hessian matrix (molar second derivative matrix) is stored as a compact linear array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A function is provided to map matrix indices to compact storage 1-D
array indices

.. code:: ipython3

    for i in range(1,c+1):
        print ("[ ", end=' ')
        for j in range (1,c+1):
            print ((i,j), end=' ')
        print (']     [', end=' ')
        for j in range (1,c+1):
            print (model.symmetric_index_from_2d_array(elm=(i,j)), end=' ')
        print (']')


.. parsed-literal::

    [  (1, 1) (1, 2) (1, 3) ]     [ 0 1 2 ]
    [  (2, 1) (2, 2) (2, 3) ]     [ 1 3 4 ]
    [  (3, 1) (3, 2) (3, 3) ]     [ 2 4 5 ]


.. code:: ipython3

    def printResult(name, result, units):
        print ("{0:<10.10s}".format(name), end=' ')
        [print ("{0:13.6e}".format(x), end=' ') for x in result]
        print ("{0:<10.10s}".format(units))
    def printLabels(n):
        print ("{0:<18.18s}".format(''), end=' ')
        maxIdx = int(len(n)*(len(n)-1)/2 + len(n))
        [print ("[{0:3d}]{1:<8.8s}".format(idx, ''), end=' ') for idx in range(maxIdx)]
        print ()
    printLabels(n)
    try:
        printResult('d2Gdn2', Simple_Solution.cy_Feldspar_Simple_Solution_d2gdn2(t,p,n), 'J/m^2')
        printResult('d3Gdn2dT', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdn2dt(t,p,n), 'J/K-m^2')
        printResult('d3Gdn2dP', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdn2dp(t,p,n), 'J/bar-m^2')
        printResult('d4Gdn2dT2', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdn2dt2(t,p,n), 'J/K^2-m^2')
        printResult('d4Gdn2dTdP', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdn2dtdp(t,p,n), 'J/K-bar-m^2')
        printResult('d4Gdn2dP2', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdn2dp2(t,p,n), 'J/bar^2-m^2')
        printResult('d5Gdn2dT3', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn2dt3(t,p,n), 'J/K^3-m^2')
        printResult('d5Gdn2dT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn2dt2dp(t,p,n), 'J/K^2-bar-m^2')
        printResult('d5Gdn2dTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn2dtdp2(t,p,n), 'J/K-bar^2-m^2')
        printResult('d5Gdn2dP3', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn2dp3(t,p,n), 'J/bar^3-m^2')
    except AttributeError:
        pass
    try:
        printResult('d2Gdn2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d2gdn2(t,p,n), 'J/m^2')
        printResult('d3Gdn2dT', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdn2dt(t,p,n), 'J/K-m^2')
        printResult('d3Gdn2dP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdn2dp(t,p,n), 'J/bar-m^2')
        printResult('d4Gdn2dT2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdn2dt2(t,p,n), 'J/K^2-m^2')
        printResult('d4Gdn2dTdP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdn2dtdp(t,p,n), 'J/K-bar-m^2')
        printResult('d4Gdn2dP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdn2dp2(t,p,n), 'J/bar^2-m^2')
        printResult('d5Gdn2dT3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn2dt3(t,p,n), 'J/K^3-m^2')
        printResult('d5Gdn2dT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn2dt2dp(t,p,n), 'J/K^2-bar-m^2')
        printResult('d5Gdn2dTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn2dtdp2(t,p,n), 'J/K-bar^2-m^2')
        printResult('d5Gdn2dP3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn2dp3(t,p,n), 'J/bar^3-m^2')
    except AttributeError:
        pass


.. parsed-literal::

                       [  0]         [  1]         [  2]         [  3]         [  4]         [  5]         
    d2Gdn2      1.332988e+04 -6.321989e+03 -5.443443e+03 -2.308700e+03  7.480484e+03 -2.299072e+03 J/m^2     
    d3Gdn2dT    8.118868e+00 -1.505944e+00 -5.479710e+00  3.356281e+00 -1.823845e+00  6.320227e+00 J/K-m^2   
    d3Gdn2dP   -4.982665e-02 -6.268677e-02  1.000257e-01  1.380904e-01 -7.442545e-02 -1.593673e-02 J/bar-m^2 
    d4Gdn2dT2   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^2-m^2 
    d4Gdn2dTdP  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar-m^
    d4Gdn2dP2   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/bar^2-m^
    d5Gdn2dT3   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^3-m^2 
    d5Gdn2dT2d  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^2-bar-
    d5Gdn2dTdP  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar^2-
    d5Gdn2dP3   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/bar^3-m^


The 3-D Tensor (molar third derivative tensor) is stored as a compact linear array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| A function is provided to map matrix indices to compact storage 1-D
  array indices:
| If :math:`n_c` represents the number of components in the solution,
  and
| if :math:`n_d` represents the dimensionality of molar derivative (in
  this case 3), then
| the number of numerically ordered permutations of :math:`n_c` molar
  derivatives taken :math:`n_d` at a time is:

.. code:: ipython3

    n_c,n_d = sym.symbols('n_c n_d')
    q = sym.factorial(n_c+n_d-1)/sym.factorial(n_d)/sym.factorial(n_c-1)
    q




.. math::

    \frac{\left(n_{c} + n_{d} - 1\right)!}{n_{d}! \left(n_{c} - 1\right)!}



Substituting :math:`n_d` equal to 3 and simplifying gives:

.. code:: ipython3

    q = sym.simplify(q.subs(n_d,3))
    q




.. math::

    \frac{n_{c} \left(n_{c} + 1\right) \left(n_{c} + 2\right)}{6}



and, for the number of components in this solution, there will be the
following number of unique terms in the third derivative tensor:

.. code:: ipython3

    q.subs(n_c,c)




.. math::

    10



A function is provided to map matrix indices to compact storage 1-D
array indices

.. code:: ipython3

    for i in range(1,c+1):
        for j in range (1,c+1):
            print ("[", end=' ')
            for k in range (1,c+1):
                print ("{0:1d}{1:1d}{2:1d}".format(i,j,k), end=' ')
            print ('] ', end=' ')
        print ('  ->  ', end=' ')
        for j in range (1,c+1):
            print ("[", end=' ')
            for k in range (1,c+1):
                print (model.symmetric_index_from_3d_array(elm=(i,j,k)), end=' ')
            print ('] ', end=' ')
        print ('')


.. parsed-literal::

    [ 111 112 113 ]  [ 121 122 123 ]  [ 131 132 133 ]    ->   [ 0 1 2 ]  [ 1 3 4 ]  [ 2 4 5 ]  
    [ 211 212 213 ]  [ 221 222 223 ]  [ 231 232 233 ]    ->   [ 1 3 4 ]  [ 3 6 7 ]  [ 4 7 8 ]  
    [ 311 312 313 ]  [ 321 322 323 ]  [ 331 332 333 ]    ->   [ 2 4 5 ]  [ 4 7 8 ]  [ 5 8 9 ]  


.. code:: ipython3

    def printResult(name, result, units):
        print ("{0:<10.10s}".format(name), end=' ')
        [print ("{0:10.3e}".format(x), end=' ') for x in result]
        print ("{0:<14.14s}".format(units))
    def printLabels(n):
        print ("{0:<15.15s}".format(''), end=' ')
        maxIdx = int(len(n)*(len(n)+1)*(len(n)+2)/6)
        [print ("[{0:3d}]{1:<5.5s}".format(idx, ''), end=' ') for idx in range(maxIdx)]
        print ()
    printLabels(n)
    try:
        printResult('d3Gdn3', Simple_Solution.cy_Feldspar_Simple_Solution_d3gdn3(t,p,n), 'J/m^3')
        printResult('d4Gdn3dT', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdn3dt(t,p,n), 'J/K-m^3')
        printResult('d4Gdn3dP', Simple_Solution.cy_Feldspar_Simple_Solution_d4gdn3dp(t,p,n), 'J/bar-m^3')
        printResult('d5Gdn3dT2', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn3dt2(t,p,n), 'J/K^2-m^3')
        printResult('d5Gdn3dTdP', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn3dtdp(t,p,n), 'J/K-bar-m^3')
        printResult('d5Gdn3dP2', Simple_Solution.cy_Feldspar_Simple_Solution_d5gdn3dp2(t,p,n), 'J/bar^2-m^3')
        printResult('d6Gdn3dT3', Simple_Solution.cy_Feldspar_Simple_Solution_d6gdn3dt3(t,p,n), 'J/K^3-m^3')
        printResult('d6Gdn3dT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_d6gdn3dt2dp(t,p,n), 'J/K^2-bar-m^3')
        printResult('d6Gdn3dTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_d6gdn3dtdp2(t,p,n), 'J/K-bar^2-m^3')
        printResult('d6Gdn3dP3', Simple_Solution.cy_Feldspar_Simple_Solution_d6gdn3dp3(t,p,n), 'J/bar^3-m^3')
    except AttributeError:
        pass
    try:
        printResult('d3Gdn3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d3gdn3(t,p,n), 'J/m^3')
        printResult('d4Gdn3dT', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdn3dt(t,p,n), 'J/K-m^3')
        printResult('d4Gdn3dP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d4gdn3dp(t,p,n), 'J/bar-m^3')
        printResult('d5Gdn3dT2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn3dt2(t,p,n), 'J/K^2-m^3')
        printResult('d5Gdn3dTdP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn3dtdp(t,p,n), 'J/K-bar-m^3')
        printResult('d5Gdn3dP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d5gdn3dp2(t,p,n), 'J/bar^2-m^3')
        printResult('d6Gdn3dT3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d6gdn3dt3(t,p,n), 'J/K^3-m^3')
        printResult('d6Gdn3dT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d6gdn3dt2dp(t,p,n), 'J/K^2-bar-m^3')
        printResult('d6Gdn3dTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d6gdn3dtdp2(t,p,n), 'J/K-bar^2-m^3')
        printResult('d6Gdn3dP3', Simple_Solution.cy_Feldspar_Simple_Solution_calib_d6gdn3dp3(t,p,n), 'J/bar^3-m^3')
    except AttributeError:
        pass


.. parsed-literal::

                    [  0]      [  1]      [  2]      [  3]      [  4]      [  5]      [  6]      [  7]      [  8]      [  9]      
    d3Gdn3     -1.450e+04  3.638e+01  1.984e+03  6.227e+03 -9.160e+02  3.354e+03 -1.071e+03 -2.505e+03 -2.667e+03  1.392e+03 J/m^3         
    d4Gdn3dT   -8.621e+00 -6.021e-01  1.606e+00  5.459e-01  1.164e+00  1.782e+00 -4.080e+00  7.225e-01 -2.489e-01 -6.140e+00 J/K-m^3       
    d4Gdn3dP    1.570e-02  7.428e-02 -4.353e-02  1.417e-02 -2.772e-02 -1.453e-02 -1.646e-01  3.375e-02  4.955e-02 -2.119e-02 J/bar-m^3     
    d5Gdn3dT2   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^2-m^3     
    d5Gdn3dTdP  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-bar-m^3   
    d5Gdn3dP2   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar^2-m^3   
    d6Gdn3dT3   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^3-m^3     
    d6Gdn3dT2d  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^2-bar-m^3 
    d6Gdn3dTdP  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-bar^2-m^3 
    d6Gdn3dP3   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar^3-m^3   


Test and time the generated functions for Feldspar
--------------------------------------------------

Time the code

.. code:: ipython3

    try:
        %timeit(Simple_Solution.cy_Feldspar_Simple_Solution_g(t, p, n))
    except AttributeError:
        pass
    try:
        %timeit(Simple_Solution.cy_Feldspar_Simple_Solution_calib_g(t, p, n))
    except AttributeError:
        pass


.. parsed-literal::

    959 ns ± 12 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)


Time the Rubicon wrapped Objective-C code

.. code:: ipython3

    from thermoengine import model as stdmodel
    modelDB = stdmodel.Database()
    FeldsparHC = modelDB.get_phase('Fsp')

.. code:: ipython3

    %timeit(FeldsparHC.gibbs_energy(t,p,mol=n))


.. parsed-literal::

    101 µs ± 1.06 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)


Methods available only in the “Calib” versions of generated code
----------------------------------------------------------------

Execute the parameter value/metadata functions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions are only defined for the “calibration” model code
implementation:

.. code:: ipython3

    Simple_Solution.cy_Feldspar_Simple_Solution_get_param_number()
    try:
        np = Simple_Solution.cy_Feldspar_Simple_Solution_get_param_number()
        names = Simple_Solution.cy_Feldspar_Simple_Solution_get_param_names()
        units = Simple_Solution.cy_Feldspar_Simple_Solution_get_param_units()
        values = Simple_Solution.cy_Feldspar_Simple_Solution_get_param_values()
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], Simple_Solution.cy_Feldspar_Simple_Solution_get_param_value(i), units[i]))
    except AttributeError:
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+00  1.000000e+00 bar       
    Wh12        7.924000e+03  7.924000e+03 J/m       
    Ws12        0.000000e+00  0.000000e+00 J/K-m     
    Wv12        0.000000e+00  0.000000e+00 J/bar-m   
    dWh12      -7.924000e+03 -7.924000e+03 J/m       
    dWs12       0.000000e+00  0.000000e+00 J/K-m     
    dWv12       0.000000e+00  0.000000e+00 J/bar-m   
    Wh13        4.613000e+04  4.613000e+04 J/m       
    Ws13        2.060000e+01  2.060000e+01 J/K-m     
    Wv13        7.866000e-01  7.866000e-01 J/bar-m   
    dWh13       8.510000e+03  8.510000e+03 J/m       
    dWs13       0.000000e+00  0.000000e+00 J/K-m     
    dWv13      -1.338000e-01 -1.338000e-01 J/bar-m   
    Wh23        7.929100e+04  7.929100e+04 J/m       
    Ws23        0.000000e+00  0.000000e+00 J/K-m     
    Wv23       -1.037000e-01 -1.037000e-01 J/bar-m   
    dWh23       1.343000e+03  1.343000e+03 J/m       
    dWs23       0.000000e+00  0.000000e+00 J/K-m     
    dWv23       1.037000e-01  1.037000e-01 J/bar-m   
    Wh123       1.254500e+04  1.254500e+04 J/m       
    Ws123       0.000000e+00  0.000000e+00 J/K-m     
    Wv123      -1.095000e+00 -1.095000e+00 J/bar-m   


Functions that allow modification of the array of parameter values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    try:
        values[1] = 100.0
        Simple_Solution.cy_Feldspar_Simple_Solution_set_param_values(values)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], Simple_Solution.cy_Feldspar_Simple_Solution_get_param_value(i), units[i]))
    except (AttributeError, NameError):
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+02  1.000000e+02 bar       
    Wh12        7.924000e+03  7.924000e+03 J/m       
    Ws12        0.000000e+00  0.000000e+00 J/K-m     
    Wv12        0.000000e+00  0.000000e+00 J/bar-m   
    dWh12      -7.924000e+03 -7.924000e+03 J/m       
    dWs12       0.000000e+00  0.000000e+00 J/K-m     
    dWv12       0.000000e+00  0.000000e+00 J/bar-m   
    Wh13        4.613000e+04  4.613000e+04 J/m       
    Ws13        2.060000e+01  2.060000e+01 J/K-m     
    Wv13        7.866000e-01  7.866000e-01 J/bar-m   
    dWh13       8.510000e+03  8.510000e+03 J/m       
    dWs13       0.000000e+00  0.000000e+00 J/K-m     
    dWv13      -1.338000e-01 -1.338000e-01 J/bar-m   
    Wh23        7.929100e+04  7.929100e+04 J/m       
    Ws23        0.000000e+00  0.000000e+00 J/K-m     
    Wv23       -1.037000e-01 -1.037000e-01 J/bar-m   
    dWh23       1.343000e+03  1.343000e+03 J/m       
    dWs23       0.000000e+00  0.000000e+00 J/K-m     
    dWv23       1.037000e-01  1.037000e-01 J/bar-m   
    Wh123       1.254500e+04  1.254500e+04 J/m       
    Ws123       0.000000e+00  0.000000e+00 J/K-m     
    Wv123      -1.095000e+00 -1.095000e+00 J/bar-m   


Functions that allow modification of a particular parameter value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    try:
        Simple_Solution.cy_Feldspar_Simple_Solution_set_param_value(1, 1.0)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], Simple_Solution.cy_Feldspar_Simple_Solution_get_param_value(i), units[i]))
    except AttributeError:
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+02  1.000000e+00 bar       
    Wh12        7.924000e+03  7.924000e+03 J/m       
    Ws12        0.000000e+00  0.000000e+00 J/K-m     
    Wv12        0.000000e+00  0.000000e+00 J/bar-m   
    dWh12      -7.924000e+03 -7.924000e+03 J/m       
    dWs12       0.000000e+00  0.000000e+00 J/K-m     
    dWv12       0.000000e+00  0.000000e+00 J/bar-m   
    Wh13        4.613000e+04  4.613000e+04 J/m       
    Ws13        2.060000e+01  2.060000e+01 J/K-m     
    Wv13        7.866000e-01  7.866000e-01 J/bar-m   
    dWh13       8.510000e+03  8.510000e+03 J/m       
    dWs13       0.000000e+00  0.000000e+00 J/K-m     
    dWv13      -1.338000e-01 -1.338000e-01 J/bar-m   
    Wh23        7.929100e+04  7.929100e+04 J/m       
    Ws23        0.000000e+00  0.000000e+00 J/K-m     
    Wv23       -1.037000e-01 -1.037000e-01 J/bar-m   
    dWh23       1.343000e+03  1.343000e+03 J/m       
    dWs23       0.000000e+00  0.000000e+00 J/K-m     
    dWv23       1.037000e-01  1.037000e-01 J/bar-m   
    Wh123       1.254500e+04  1.254500e+04 J/m       
    Ws123       0.000000e+00  0.000000e+00 J/K-m     
    Wv123      -1.095000e+00 -1.095000e+00 J/bar-m   


Functions that evaluate parameter derivatives …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    try:
        fmt = "    {0:<10.10s} {1:13.6e}"
        for i in range(0, np):
            print ('Derivative with respect to parameter: ', names[i], ' of')
            print (fmt.format('G', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_g(t, p, n, i)))
            print (fmt.format('dGdT', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_dgdt(t, p, n, i)))
            print (fmt.format('dGdP', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_dgdp(t, p, n, i)))
            print (fmt.format('d2GdT2', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d2gdt2(t, p, n, i)))
            print (fmt.format('d2GdTdP', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d2gdtdp(t, p, n, i)))
            print (fmt.format('d2GdP2', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d2gdp2(t, p, n, i)))
            print (fmt.format('d3GdT3', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d3gdt3(t, p, n, i)))
            print (fmt.format('d3GdT2dP', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d3gdt2dp(t, p, n, i)))
            print (fmt.format('d3GdTdP2', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d3gdtdp2(t, p, n, i)))
            print (fmt.format('d3GdP3', Simple_Solution.cy_Feldspar_Simple_Solution_dparam_d3gdp3(t, p, n, i)))
    except (AttributeError, TypeError):
        pass


.. parsed-literal::

    Derivative with respect to parameter:  T_r  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  P_r  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wh12  of
        G           3.666667e-01
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Ws12  of
        G          -7.333333e+02
        dGdT       -3.666667e-01
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wv12  of
        G           3.666667e-01
        dGdT        0.000000e+00
        dGdP        3.666667e-01
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWh12  of
        G          -1.018519e-02
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWs12  of
        G           2.037037e+01
        dGdT        1.018519e-02
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWv12  of
        G          -1.018519e-02
        dGdT        0.000000e+00
        dGdP       -1.018519e-02
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wh13  of
        G           3.972222e-01
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Ws13  of
        G          -7.944444e+02
        dGdT       -3.972222e-01
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wv13  of
        G           3.972222e-01
        dGdT        0.000000e+00
        dGdP        3.972222e-01
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWh13  of
        G          -2.206790e-02
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWs13  of
        G           4.413580e+01
        dGdT        2.206790e-02
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWv13  of
        G          -2.206790e-02
        dGdT        0.000000e+00
        dGdP       -2.206790e-02
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wh23  of
        G           4.333333e-01
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Ws23  of
        G          -8.666667e+02
        dGdT       -4.333333e-01
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wv23  of
        G           4.333333e-01
        dGdT        0.000000e+00
        dGdP        4.333333e-01
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWh23  of
        G          -1.203704e-02
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWs23  of
        G           2.407407e+01
        dGdT        1.203704e-02
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  dWv23  of
        G          -1.203704e-02
        dGdT        0.000000e+00
        dGdP       -1.203704e-02
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wh123  of
        G           1.324074e-01
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Ws123  of
        G          -2.648148e+02
        dGdT       -1.324074e-01
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Wv123  of
        G           1.324074e-01
        dGdT        0.000000e+00
        dGdP        1.324074e-01
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00


Parameter derivatives of the chemical potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def printResult(name, result, units):
        print ("dmu[*]/d {0:<10.10s}".format(name), end=' ')
        [print ("{0:13.6e}".format(x), end=' ') for x in result]
        print ("{0:<12.12s}".format(units))
    def printLabels(n):
        print ("         {0:<18.18s}".format(''), end=' ')
        [print ("[{0:3d}]{1:<8.8s}".format(idx, ''), end=' ') for idx in range(len(n))]
        print ()
    printLabels(n)
    try:
        for i in range(0, np):
            result = Simple_Solution.cy_Feldspar_Simple_Solution_dparam_dgdn(t,p,n, i)
            printResult(names[i], result, 'J/m^2/p-unit')
    except AttributeError:
        pass    


.. parsed-literal::

                                [  0]         [  1]         [  2]         
    dmu[*]/d T_r         0.000000e+00  0.000000e+00  0.000000e+00 J/m^2/p-unit
    dmu[*]/d P_r         0.000000e+00  0.000000e+00  0.000000e+00 J/m^2/p-unit
    dmu[*]/d Wh12        2.314815e-01  2.037037e-01 -1.018519e-01 J/m^2/p-unit
    dmu[*]/d Ws12       -4.629630e+02 -4.074074e+02  2.037037e+02 J/m^2/p-unit
    dmu[*]/d Wv12        2.314815e-01  2.037037e-01 -1.018519e-01 J/m^2/p-unit
    dmu[*]/d dWh12       9.825103e-02 -1.046811e-01  5.658436e-03 J/m^2/p-unit
    dmu[*]/d dWs12      -1.965021e+02  2.093621e+02 -1.131687e+01 J/m^2/p-unit
    dmu[*]/d dWv12       9.825103e-02 -1.046811e-01  5.658436e-03 J/m^2/p-unit
    dmu[*]/d Wh13        2.507716e-01 -1.103395e-01  1.952160e-01 J/m^2/p-unit
    dmu[*]/d Ws13       -5.015432e+02  2.206790e+02 -3.904321e+02 J/m^2/p-unit
    dmu[*]/d Wv13        2.507716e-01 -1.103395e-01  1.952160e-01 J/m^2/p-unit
    dmu[*]/d dWh13       1.025377e-01  1.225995e-02 -1.150549e-01 J/m^2/p-unit
    dmu[*]/d dWs13      -2.050754e+02 -2.451989e+01  2.301097e+02 J/m^2/p-unit
    dmu[*]/d dWv13       1.025377e-01  1.225995e-02 -1.150549e-01 J/m^2/p-unit
    dmu[*]/d Wh23       -1.203704e-01  2.407407e-01  2.129630e-01 J/m^2/p-unit
    dmu[*]/d Ws23        2.407407e+02 -4.814815e+02 -4.259259e+02 J/m^2/p-unit
    dmu[*]/d Wv23       -1.203704e-01  2.407407e-01  2.129630e-01 J/m^2/p-unit
    dmu[*]/d dWh23       6.687243e-03  1.170267e-01 -1.229424e-01 J/m^2/p-unit
    dmu[*]/d dWs23      -1.337449e+01 -2.340535e+02  2.458848e+02 J/m^2/p-unit
    dmu[*]/d dWv23       6.687243e-03  1.170267e-01 -1.229424e-01 J/m^2/p-unit
    dmu[*]/d Wh123       4.681070e-02  3.677984e-02  2.829218e-02 J/m^2/p-unit
    dmu[*]/d Ws123      -9.362140e+01 -7.355967e+01 -5.658436e+01 J/m^2/p-unit
    dmu[*]/d Wv123       4.681070e-02  3.677984e-02  2.829218e-02 J/m^2/p-unit

