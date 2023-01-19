
Berman Standard State Code Generator + Birch-Murnaghan
======================================================

Required system packages and initialization

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()

.. code:: ipython3

    from thermoengine import coder

.. code:: ipython3

    model = coder.StdStateModel()

.. code:: ipython3

    T = model.get_symbol_for_t()
    P = model.get_symbol_for_p()
    Tr = model.get_symbol_for_tr()
    Pr = model.get_symbol_for_pr()

Define model expressions applicable over all of T,P space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An expression for the Gibbs free energy, :math:`G(T,P)` or the Helmholtz
energy :math:`A(T,V)` is constructed. The expression may have multiple
parts. Often the heat capacity function is postulated, then integrated
to yield expressions for the entahlpy, entropy, and in combination the
energy potential. Then, an equation of state (EOS) is adopted and that
term is integrated in pressure or volume and added to the heat capacity
integrals. This proceedure is follwed here. #### (1) :math:`C_P`
integrals The isobaric heat capacity terms parameterized as: $C_P = k_0
+ k_1 / T^{1/2} + k_2 / T^2 + k_3 / T^3 $, and in addition the reference
condition third law entropy, $ S_{Tr,Pr} $, and enthalpy of formation
from the elements, $ :raw-latex:`\Delta `H_{Tr,Pr} $, constitute
additional parameters:

.. code:: ipython3

    k0,k1,k2,k3 = sym.symbols('k0 k1 k2 k3')
    CpPr = k0+k1/sym.sqrt(T)+k2/T**2+k3/T**3
    STrPr,HTrPr = sym.symbols('S_TrPr H_TrPr')

Specify paramters …

.. code:: ipython3

    params = [('H_TrPr','J',HTrPr), ('S_TrPr','J/K',STrPr), ('k0','J/K-m',k0), ('k1','J-K^(1/2)-m',k1),
              ('k2','J-K/m',k2),  ('k3','J-K^2',k3)]

Define the heat capacity contribution to the Gibbs free energy …

.. code:: ipython3

    GPr = HTrPr + sym.integrate(CpPr,(T,Tr,T)) - T*(STrPr + sym.integrate(CpPr/T,(T,Tr,T)))

… and add this expression to the model

.. code:: ipython3

    model.add_expression_to_model(GPr, params)

(2) :math:`V` (EOS) integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| Next, define a volume-implicit equation of state applicable over the
  whole of temperature and pressure space. We will use the 3rd order
  Birch-Murnaghan expression:
| :math:`P = \frac{{3K}}{2}\left[ {{{\left( {\frac{{{V_{{T_r}.{P_r}}}}}{V}} \right)}^{\frac{7}{3}}} - {{\left( {\frac{{{V_{{T_r}.{P_r}}}}}{V}} \right)}^{\frac{5}{3}}}} \right]\left\{ {\left( {\frac{{3{K_P}}}{4} - 3} \right)\left[ {{{\left( {\frac{{{V_{{T_r}.{P_r}}}}}{V}} \right)}^{\frac{2}{3}}} - 1} \right] + 1} \right\}`
| The parameters in this expression are:

.. code:: ipython3

    VTrPr,K,Kp = sym.symbols('V_TrPr K K_P')
    params = [('V_TrPr', 'J/bar-m', VTrPr), ('K','bar',K), ('K_P','',Kp)]

where *V* is an implicit function of *T* and *P*:

.. code:: ipython3

    V = sym.Function('V')(T,P)

| Define *f*, an implicit function derived from the Birch-Murnaghan
  expression. *f* has a value zero for internally consistent *V*, *T*,
  and *P*:
| :math:`f = 0 = \frac{{3K}}{2}\left[ {{{\left( {\frac{{{V_{{T_r}.{P_r}}}}}{V}} \right)}^{\frac{7}{3}}} - {{\left( {\frac{{{V_{{T_r}.{P_r}}}}}{V}} \right)}^{\frac{5}{3}}}} \right]\left\{ {\left( {\frac{{3{K_P}}}{4} - 3} \right)\left[ {{{\left( {\frac{{{V_{{T_r}.{P_r}}}}}{V}} \right)}^{\frac{2}{3}}} - 1} \right] + 1} \right\} - P`

.. code:: ipython3

    f = (sym.S(3)*K/sym.S(2))*((VTrPr/V)**(sym.S(7)/sym.S(3))-(VTrPr/V)**(sym.S(5)/sym.S(3)))*(1+(sym.S(3)/sym.S(4))*(Kp-4)*((VTrPr/V)**(sym.S(2)/sym.S(3))-1))- P
    f




.. math::

    \frac{3 K \left(\left(\frac{V_{TrPr}}{V{\left (T,P \right )}}\right)^{\frac{7}{3}} - \left(\frac{V_{TrPr}}{V{\left (T,P \right )}}\right)^{\frac{5}{3}}\right) \left(\left(\frac{3 K_{P}}{4} - 3\right) \left(\left(\frac{V_{TrPr}}{V{\left (T,P \right )}}\right)^{\frac{2}{3}} - 1\right) + 1\right)}{2} - P



| Because the EOS is explicit in *P* (and a function of *T* and *V*),
  the natural thermodynamic potential to use is the Helmholtz energy,
  *A*. *A* is obtained by integrating pressure from the reference volume
  to the final volume:
| :math:`{A_{T,P}} - {A_{{T_r},{P_r}}} = \int_{{V_{{T_r},{P_r}}}}^{{V_{T,P}}} {PdV}`
| Note: To perform this integration in SymPy, first define a variable of
  integration, :math:`V_{TP}`, substitute that variable for the function
  :math:`V(T,P)`, and integrate :math:`f` with respect to :math:`V_{TP}`
  over the limits :math:`V_{TrPr}` to :math:`V(T,P)`. This procedure
  generates an expression for the Helmholtz energy that is a function of
  :math:`V(T,P)`.
| Note: The integration is performed on the integrand :math:`f+P`, which
  corresponds to the Birch-Murnaghan expression for :math:`P`.

.. code:: ipython3

    VTP = sym.symbols('V_TP')
    A = sym.integrate((f+P).subs(V,VTP),(VTP,VTrPr,V)).simplify()
    A




.. math::

    \frac{9 K \left(- K_{P} V_{TrPr}^{3} + V_{TrPr}^{\frac{5}{3}} \left(3 K_{P} V_{TrPr}^{\frac{2}{3}} \left(\frac{1}{V{\left (T,P \right )}}\right)^{\frac{2}{3}} - 3 K_{P} - 14 V_{TrPr}^{\frac{2}{3}} \left(\frac{1}{V{\left (T,P \right )}}\right)^{\frac{2}{3}} + 16\right) \left(\frac{1}{V{\left (T,P \right )}}\right)^{\frac{5}{3}} V^{3}{\left (T,P \right )} + 4 V_{TrPr}^{3} - \left(3 K_{P} V_{TrPr}^{\frac{10}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{7}{3}} - 3 K_{P} V_{TrPr}^{\frac{8}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} - K_{P} V_{TrPr} - 14 V_{TrPr}^{\frac{10}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{7}{3}} + 16 V_{TrPr}^{\frac{8}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} + 4 V_{TrPr}\right) V^{2}{\left (T,P \right )}\right)}{16 V^{2}{\left (T,P \right )}}



The Gibbs free energy contribution if given by the identity:
:math:`G = A + PV` and
:math:`{G_{T,P}} - {G_{{T_r},{P_r}}} = {A_{T,P}} + PV - {A_{{T_r},{P_r}}} - {P_r}{V_{{T_r},{P_r}}}`

.. code:: ipython3

    GTr = A+(f+P)*V-A.subs(V,VTrPr)-Pr*VTrPr
    GTr




.. math::

    \frac{3 K \left(\left(\frac{V_{TrPr}}{V{\left (T,P \right )}}\right)^{\frac{7}{3}} - \left(\frac{V_{TrPr}}{V{\left (T,P \right )}}\right)^{\frac{5}{3}}\right) \left(\left(\frac{3 K_{P}}{4} - 3\right) \left(\left(\frac{V_{TrPr}}{V{\left (T,P \right )}}\right)^{\frac{2}{3}} - 1\right) + 1\right) V{\left (T,P \right )}}{2} + \frac{9 K \left(- K_{P} V_{TrPr}^{3} + V_{TrPr}^{\frac{5}{3}} \left(3 K_{P} V_{TrPr}^{\frac{2}{3}} \left(\frac{1}{V{\left (T,P \right )}}\right)^{\frac{2}{3}} - 3 K_{P} - 14 V_{TrPr}^{\frac{2}{3}} \left(\frac{1}{V{\left (T,P \right )}}\right)^{\frac{2}{3}} + 16\right) \left(\frac{1}{V{\left (T,P \right )}}\right)^{\frac{5}{3}} V^{3}{\left (T,P \right )} + 4 V_{TrPr}^{3} - \left(3 K_{P} V_{TrPr}^{\frac{10}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{7}{3}} - 3 K_{P} V_{TrPr}^{\frac{8}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} - K_{P} V_{TrPr} - 14 V_{TrPr}^{\frac{10}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{7}{3}} + 16 V_{TrPr}^{\frac{8}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} + 4 V_{TrPr}\right) V^{2}{\left (T,P \right )}\right)}{16 V^{2}{\left (T,P \right )}} - \frac{9 K \left(- K_{P} V_{TrPr}^{3} + V_{TrPr}^{\frac{14}{3}} \left(3 K_{P} V_{TrPr}^{\frac{2}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{2}{3}} - 3 K_{P} - 14 V_{TrPr}^{\frac{2}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{2}{3}} + 16\right) \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} + 4 V_{TrPr}^{3} - V_{TrPr}^{2} \left(3 K_{P} V_{TrPr}^{\frac{10}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{7}{3}} - 3 K_{P} V_{TrPr}^{\frac{8}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} - K_{P} V_{TrPr} - 14 V_{TrPr}^{\frac{10}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{7}{3}} + 16 V_{TrPr}^{\frac{8}{3}} \left(\frac{1}{V_{TrPr}}\right)^{\frac{5}{3}} + 4 V_{TrPr}\right)\right)}{16 V_{TrPr}^{2}} - P_{r} V_{TrPr}



The *implicit_function* argument of the *add_expression_to_model* method
conveys information about how to compute a value for the implicit
variable contained in the Gibbs free energy expression passed as the
first argument. It is an array of tuples; each tuple has three
components. The first is a sympy expression for the implicit function,
which evaluates to zero for an internally consistent set of *T*, *P*,
*V*. The second component is a symbol for the function definition of the
implicit variable. The third component is is a sympy expression that
initializes *f* in the iterative routine. This expression must be
defined in terms of known parameters and Tr, Pr, T, P..

.. code:: ipython3

    model.add_expression_to_model(GTr, params, implicit_functions=[(f,V,VTrPr/2.0)])

| Note:
| The implicit function will be utilized in code generation not only to
  compute the value of *V*, given a *T* and *P*, but also to evaluate
  derivatives of *V*, Since:

:math:`dP = {\left( {\frac{{\partial P}}{{\partial V}}} \right)_T}dV + {\left( {\frac{{\partial P}}{{\partial T}}} \right)_V}dT`,

:math:`{d^2}P = d{\left( {\frac{{\partial P}}{{\partial V}}} \right)_T}dV + {\left( {\frac{{\partial P}}{{\partial V}}} \right)_T}{d^2}V + d{\left( {\frac{{\partial P}}{{\partial T}}} \right)_V}dT + {\left( {\frac{{\partial P}}{{\partial T}}} \right)_V}{d^2}T`,
or

:math:`{d^2}P = {\left( {\frac{{{\partial ^2}P}}{{\partial {V^2}}}} \right)_T}dVdV + 2\left( {\frac{{{\partial ^2}P}}{{\partial V\partial T}}} \right)dVdT + {\left( {\frac{{{\partial ^2}P}}{{\partial {T^2}}}} \right)_V}dTdT + {\left( {\frac{{\partial P}}{{\partial V}}} \right)_T}{d^2}V + {\left( {\frac{{\partial P}}{{\partial T}}} \right)_V}{d^2}T`

Code Print the Model, compile the code and link a Python module
---------------------------------------------------------------

Name the model class

.. code:: ipython3

    model.set_module_name('berman')
    #model.get_berman_std_state_database()

Make a working sub-directory and move down into the directory. This is
done so that generated files will not clash between alternate model
configurations.

.. code:: ipython3

    model_working_dir = "working"
    !mkdir -p {model_working_dir}
    %cd {model_working_dir}


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working


-  Choose an existing phase from the Berman database
-  Generate an include file and code file for this phase

Note that the call to

::

    model.create_code_module(phase=phase_name, formula=formula, params=param_dict)

generates fast code with unmodifiable model parameters and
“calibration-” related functions. The call to:

::

    model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type='calib')

generates code suitable for model parameter calibration.

.. code:: ipython3

    model_type='fast'

.. code:: ipython3

    param_dict = model.get_berman_std_state_database(42, extend_defs=True)
    #param_dict['K'] = 100.0
    #param_dict['K_P'] = 4.0
    phase_name = param_dict.pop('Phase', None).title()
    formula = param_dict.pop('Formula', None)
    result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type=model_type)


.. parsed-literal::

    Creating (once only) generic fast model code file string
    Creating (once only) generic model fast code template include file string
    Creating (once only) generic model fast code template code file string
    Creating include file ...
    ... done!
    Creating code file ...
    ... done
    Writing include file to working directory ...
    Writing code file to working directory ...
    Writing pyxbld file to working directory ...
    writing pyx file to working directory ...
    Compiling code and Python bindings ...
    Success! Import the module named  berman


.. code:: ipython3

    param_dict




.. parsed-literal::

    {'H_TrPr': -3970790.781,
     'S_TrPr': 214.145,
     'k0': 381.37231,
     'k1': -1941.045,
     'k2': -12037252,
     'k3': 1836425472,
     'V_TrPr': 10.869,
     'v1': -1.8045e-06,
     'v2': 0,
     'v3': 1.51451e-05,
     'v4': 5.5000000000000004e-09,
     'l1': 0.0,
     'l2': 0.0,
     'k_lambda': 0.0,
     'T_lambda_Pr': 0.0,
     'T_lambda_ref': 0.0,
     'H_t': 0.0,
     'd0': 282.98291,
     'd1': -4831.375,
     'd2': 3620706.0,
     'd3': -0.15733000000000003,
     'd4': 3.477e-05,
     'd5': 410629.63,
     'T_D': 1436.15,
     'T_D_ref': 298.15,
     'T_r': 298.15,
     'P_r': 1.0,
     'K': 554170.1302299806,
     'K_P': 4.0}



Import the new module and test the model
----------------------------------------

.. code:: ipython3

    import berman
    %cd ..


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen


Evaluate functions at temperature (K) and pressure (bars)

.. code:: ipython3

    t = 1000.0
    p = 10000.0

Available in both “Fast” and “Calib” code versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execute the “fast” or “calibration” code metadata retrieval functions:

.. code:: ipython3

    try:
        print(berman.cy_Potassium_Feldspar_berman_identifier())
        print(berman.cy_Potassium_Feldspar_berman_name())
        print(berman.cy_Potassium_Feldspar_berman_formula())
        print(berman.cy_Potassium_Feldspar_berman_mw())
        print(berman.cy_Potassium_Feldspar_berman_elements())
    except AttributeError:
        pass
    try:
        print(berman.cy_Potassium_Feldspar_berman_calib_identifier())
        print(berman.cy_Potassium_Feldspar_berman_calib_name())
        print(berman.cy_Potassium_Feldspar_berman_calib_formula())
        print(berman.cy_Potassium_Feldspar_berman_calib_mw())
        print(berman.cy_Potassium_Feldspar_berman_calib_elements())
    except AttributeError:
        pass


.. parsed-literal::

    Thu Sep 20 12:23:25 2018
    Potassium_Feldspar
    KAlSi3O8
    278.33524
    [0. 0. 0. 0. 0. 0. 0. 0. 8. 0. 0. 0. 0. 1. 3. 0. 0. 0. 0. 1. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]


Execute the standard thermodynamic property retrieval functions:

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
    try:
        print(fmt.format('G', berman.cy_Potassium_Feldspar_berman_g(t,p), 'J/m'))
        print(fmt.format('dGdT', berman.cy_Potassium_Feldspar_berman_dgdt(t,p), 'J/K-m'))
        print(fmt.format('dGdP', berman.cy_Potassium_Feldspar_berman_dgdp(t,p), 'J/bar-m'))
        print(fmt.format('d2GdP2', berman.cy_Potassium_Feldspar_berman_d2gdt2(t,p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', berman.cy_Potassium_Feldspar_berman_d2gdtdp(t,p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', berman.cy_Potassium_Feldspar_berman_d2gdp2(t,p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', berman.cy_Potassium_Feldspar_berman_d3gdt3(t,p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', berman.cy_Potassium_Feldspar_berman_d3gdt2dp(t,p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', berman.cy_Potassium_Feldspar_berman_d3gdtdp2(t,p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', berman.cy_Potassium_Feldspar_berman_d3gdp3(t,p), 'J/bar^3-m'))
        print(fmt.format('S', berman.cy_Potassium_Feldspar_berman_s(t,p), 'J/K-m'))
        print(fmt.format('V', berman.cy_Potassium_Feldspar_berman_v(t,p), 'J/bar-m'))
        print(fmt.format('Cv', berman.cy_Potassium_Feldspar_berman_cv(t,p), 'J/K-m'))
        print(fmt.format('Cp', berman.cy_Potassium_Feldspar_berman_cp(t,p), 'J/K-m'))
        print(fmt.format('dCpdT', berman.cy_Potassium_Feldspar_berman_dcpdt(t,p), 'J/K^2-m'))
        print(fmt.format('alpha', berman.cy_Potassium_Feldspar_berman_alpha(t,p), '1/K'))
        print(fmt.format('beta', berman.cy_Potassium_Feldspar_berman_beta(t,p), '1/bar'))
        print(fmt.format('K', berman.cy_Potassium_Feldspar_berman_K(t,p), 'bar'))
        print(fmt.format('Kp', berman.cy_Potassium_Feldspar_berman_Kp(t,p), ''))
    except AttributeError:
        pass
    try:
        print(fmt.format('G', berman.cy_Potassium_Feldspar_berman_calib_g(t,p), 'J/m'))
        print(fmt.format('dGdT', berman.cy_Potassium_Feldspar_berman_calib_dgdt(t,p), 'J/K-m'))
        print(fmt.format('dGdP', berman.cy_Potassium_Feldspar_berman_calib_dgdp(t,p), 'J/bar-m'))
        print(fmt.format('d2GdP2', berman.cy_Potassium_Feldspar_berman_calib_d2gdt2(t,p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', berman.cy_Potassium_Feldspar_berman_calib_d2gdtdp(t,p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', berman.cy_Potassium_Feldspar_berman_calib_d2gdp2(t,p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', berman.cy_Potassium_Feldspar_berman_calib_d3gdt3(t,p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', berman.cy_Potassium_Feldspar_berman_calib_d3gdt2dp(t,p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', berman.cy_Potassium_Feldspar_berman_calib_d3gdtdp2(t,p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', berman.cy_Potassium_Feldspar_berman_calib_d3gdp3(t,p), 'J/bar^3-m'))
        print(fmt.format('S', berman.cy_Potassium_Feldspar_berman_calib_s(t,p), 'J/K-m'))
        print(fmt.format('V', berman.cy_Potassium_Feldspar_berman_calib_v(t,p), 'J/bar-m'))
        print(fmt.format('Cv', berman.cy_Potassium_Feldspar_berman_calib_cv(t,p), 'J/K-m'))
        print(fmt.format('Cp', berman.cy_Potassium_Feldspar_berman_calib_cp(t,p), 'J/K-m'))
        print(fmt.format('dCpdT', berman.cy_Potassium_Feldspar_berman_calib_dcpdt(t,p), 'J/K^2-m'))
        print(fmt.format('alpha', berman.cy_Potassium_Feldspar_berman_calib_alpha(t,p), '1/K'))
        print(fmt.format('beta', berman.cy_Potassium_Feldspar_berman_calib_beta(t,p), '1/bar'))
        print(fmt.format('K', berman.cy_Potassium_Feldspar_berman_calib_K(t,p), 'bar'))
        print(fmt.format('Kp', berman.cy_Potassium_Feldspar_berman_calib_Kp(t,p), ''))
    except AttributeError:
        pass


.. parsed-literal::

    G          -4.206302e+06 J/m       
    dGdT       -5.343985e+02 J/K-m     
    dGdP        1.032149e+01 J/bar-m   
    d2GdP2     -3.097903e-01 J/K^2-m   
    d2GdTdP     0.000000e+00 J/K-bar-m 
    d2GdP2     -5.097131e-05 J/bar^2-m 
    d3GdT3      2.605344e-04 J/K^3-m   
    d3GdT2dP    0.000000e+00 J/K^2-bar-
    d3GdTdP2    0.000000e+00 J/K-bar^2-
    d3GdP3      6.989655e-10 J/bar^3-m 
    S           5.343985e+02 J/K-m     
    V           1.032149e+01 J/bar-m   
    Cv          3.097903e+02 J/K-m     
    Cp          3.097903e+02 J/K-m     
    dCpdT       4.925584e-02 J/K^2-m   
    alpha       0.000000e+00 1/K       
    beta        4.938367e-06 1/bar     
    K           2.024961e+05 bar       
    Kp          1.776812e+00           


Available only in the “Calib” versions of generated code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Execute the parameter value/metadata functions.
| These functions are only defined for the “calibration” model code
  implementation:

.. code:: ipython3

    try:
        np = berman.cy_Potassium_Feldspar_berman_get_param_number()
        names = berman.cy_Potassium_Feldspar_berman_get_param_names()
        units = berman.cy_Potassium_Feldspar_berman_get_param_units()
        values = berman.cy_Potassium_Feldspar_berman_get_param_values()
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], berman.cy_Potassium_Feldspar_berman_get_param_value(i), units[i]))
    except AttributeError:
        pass

Test the functions that allow modification of the array of parameter
values

.. code:: ipython3

    try:
        values[1] = 100.0
        berman.cy_Potassium_Feldspar_berman_set_param_values(values)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], berman.cy_Potassium_Feldspar_berman_get_param_value(i), units[i]))
    except (AttributeError, NameError):
        pass

Test the functions that allow modification of a particular parameter
value

.. code:: ipython3

    try:
        berman.cy_Potassium_Feldspar_berman_set_param_value(1, 1.0)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], berman.cy_Potassium_Feldspar_berman_get_param_value(i), units[i]))
    except AttributeError:
        pass

Evaluate parameter derivatives …

.. code:: ipython3

    try:
        fmt = "    {0:<10.10s} {1:13.6e}"
        for i in range(0, np):
            print ('Derivative with respect to parameter: ', names[i], ' of')
            print (fmt.format('G', berman.cy_Potassium_Feldspar_berman_dparam_g(t, p, i)))
            print (fmt.format('dGdT', berman.cy_Potassium_Feldspar_berman_dparam_dgdt(t, p, i)))
            print (fmt.format('dGdP', berman.cy_Potassium_Feldspar_berman_dparam_dgdp(t, p, i)))
            print (fmt.format('d2GdT2', berman.cy_Potassium_Feldspar_berman_dparam_d2gdt2(t, p, i)))
            print (fmt.format('d2GdTdP', berman.cy_Potassium_Feldspar_berman_dparam_d2gdtdp(t, p, i)))
            print (fmt.format('d2GdP2', berman.cy_Potassium_Feldspar_berman_dparam_d2gdp2(t, p, i)))
            print (fmt.format('d3GdT3', berman.cy_Potassium_Feldspar_berman_dparam_d3gdt3(t, p, i)))
            print (fmt.format('d3GdT2dP', berman.cy_Potassium_Feldspar_berman_dparam_d3gdt2dp(t, p, i)))
            print (fmt.format('d3GdTdP2', berman.cy_Potassium_Feldspar_berman_dparam_d3gdtdp2(t, p, i)))
            print (fmt.format('d3GdP3', berman.cy_Potassium_Feldspar_berman_dparam_d3gdp3(t, p, i)))
    except (AttributeError, TypeError):
        pass

Time execution of the code
--------------------------

.. code:: ipython3

    try:
        %timeit(berman.cy_Potassium_Feldspar_berman_calib_g(t,p))
    except AttributeError:
        pass
    try:
        %timeit(berman.cy_Potassium_Feldspar_berman_g(t,p))
    except AttributeError:
        pass


.. parsed-literal::

    265 ns ± 4.29 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)

