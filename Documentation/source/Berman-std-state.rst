
Berman Standard State Code Generator
====================================

Required system packages and initialization

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()

Required ENKI packages

.. code:: ipython3

    from thermoengine import coder

Types of terms in a standard state properties description
=========================================================

There are three classes of terms: 1. Terms that apply over the whole of
:math:`T`-, :math:`P`-space, :math:`T_r \le T`, :math:`P_r \le P` 2.
Terms that apply over a specified range of :math:`T`-, :math:`P`-space,
:math:`(T_{r_\lambda},P_{r_\lambda}) \le (T,P) \le (T_\lambda,P_\lambda)`
3. Terms that apply to a specific :math:`T_t` and :math:`P_t` and higher
:math:`T`, :math:`P`, :math:`T_t \le T`, :math:`P_t \le P`

Second-order phase transitions (:math:`lambda`-transitions) are an
example of the second type, as are order disorder transformations.
First-order phase transitions are an example of the third type.

Create a model class for the Gibbs free energy
----------------------------------------------

.. code:: ipython3

    model = coder.StdStateModel()

Retrieve sympy symbols for model variables and reference conditions

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

.. code:: ipython3

    CpPr




.. math::

    k_{0} + \frac{k_{2}}{T^{2}} + \frac{k_{3}}{T^{3}} + \frac{k_{1}}{\sqrt{T}}



Specify paramters …

.. code:: ipython3

    params = [('H_TrPr','J',HTrPr), ('S_TrPr','J/K',STrPr), ('k0','J/K-m',k0), ('k1','J-K^(1/2)-m',k1),
              ('k2','J-K/m',k2),  ('k3','J-K^2',k3)]

Define the heat capacity contribution to the Gibbs free energy …

.. code:: ipython3

    GPr = HTrPr + sym.integrate(CpPr,(T,Tr,T)) - T*(STrPr + sym.integrate(CpPr/T,(T,Tr,T)))

.. code:: ipython3

    GPr




.. math::

    H_{TrPr} + 2 \sqrt{T} k_{1} + T k_{0} - T \left(S_{TrPr} + k_{0} \log{\left (T \right )} - k_{0} \log{\left (T_{r} \right )} + \frac{k_{2}}{2 T_{r}^{2}} + \frac{k_{3}}{3 T_{r}^{3}} + \frac{2 k_{1}}{\sqrt{T_{r}}} - \frac{k_{2}}{2 T^{2}} - \frac{k_{3}}{3 T^{3}} - \frac{2 k_{1}}{\sqrt{T}}\right) - 2 \sqrt{T_{r}} k_{1} - T_{r} k_{0} + \frac{k_{2}}{T_{r}} + \frac{k_{3}}{2 T_{r}^{2}} - \frac{k_{2}}{T} - \frac{k_{3}}{2 T^{2}}



… and add this expression to the model

.. code:: ipython3

    model.add_expression_to_model(GPr, params)

(2) :math:`V` (EOS) integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, define a volume-explicit equation of state applicable over the
whole of temperature and pressure space

.. code:: ipython3

    VTrPr,v1,v2,v3,v4 = sym.symbols('V_TrPr v1 v2 v3 v4')
    params = [('V_TrPr', 'J/bar-m', VTrPr), ('v1','1/bar',v1), ('v2','1/bar^2',v2), ('v3','1/K',v3),  ('v4','1/K^2',v4)]

.. code:: ipython3

    GPrToP = sym.integrate(VTrPr*(1+v1*(P-Pr)+v2*(P-Pr)**2+v3*(T-Tr)+v4*(T-Tr)**2),(P,Pr,P))
    model.add_expression_to_model(GPrToP, params)

Define additional lambda heat capacity terms applicable over a restricted range of T,P space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These contributions to the potential function apply over a limited range
of :math:`T,P` or :math:`T,V` space. However, the affects of these
functions propagate beyond the upper limits of the range. Say,
:math:`f(T,P)` is the contribution to the Gibbs free energy that
describes a :math:`\lambda`-like transition over the range
:math:`({T_{\lambda ,ref}},{P_{\lambda ,ref}})` to
:math:`({T_\lambda },{P_\lambda })`. Then, *above the upper limit* of
temperature there is a fixed entropy contribution:

:math:`- {\left. {\frac{{\partial f\left( {T,P} \right)}}{{\partial T}}} \right|_{{T_\lambda },{P_\lambda }}} = {f_S}({T_\lambda },{P_\lambda })`

and *above the upper limit of pressure* there is a fixed volume
contribution:

:math:`{\left. {\frac{{\partial f\left( {T,P} \right)}}{{\partial P}}} \right|_{{T_\lambda },{P_\lambda }}} = {f_V}({T_\lambda },{P_\lambda })`

the consequence of which is that for :math:`T > T_{\lambda}` and
:math:`P > P_{\lambda}`, there is a contribution to the Gibbs free
energy of the form:

:math:`- \left( {T - {T_\lambda }} \right){f_S}({T_\lambda },{P_\lambda }) + \left( {P - {P_\lambda }} \right){f_V}({T_\lambda },{P_\lambda })`

This contribution is linear in :math:`T` and :math:`P`.

The additional energetic contributions applicable above the upper range
limit will be added automatically to the model function, and need not be
explicitly accounted for in building the model expressions.

(1) :math:`\lambda`-transition-like heat capacity integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameters of the Berman (1988) lambda transition model: - :math:`l_1`
and :math:`l_2`, coefficients in Berman (1988)’s :math:`lambda`-heat
capacity model - :math:`k_{\lambda}`,
:math:`\frac{{d{T_\lambda }}}{{dP}}` in Berman (1988)’s
:math:`lambda`-transition model - :math:`T_{\lambda,{P_r}}`, Temperature
of the :math:`\lambda`-transition at reference pressure -
:math:`T_{\lambda,{ref}}`, Temperature of the lower bound of the heat
capacity integral for the :math:`\lambda`-transition at reference
pressure

.. code:: ipython3

    l1,l2 = sym.symbols('l1 l2')
    kl = sym.symbols('k_lambda')
    TlPr, Tlref = sym.symbols('T_lambda_Pr T_lambda_ref')
    params = [('l1','(J/m)^(1/2)-K', l1), ('l2', '(J/m)^(1/2)/K^2', l2), ('k_lambda', 'K/bar', kl), 
              ('T_lambda_Pr', 'K', TlPr), ('T_lambda_ref', 'K', Tlref)]

Calculate the transition temperature, :math:`T_{\lambda}`, at the
pressure, :math:`P`

.. code:: ipython3

    Tl = TlPr + kl*(P-Pr)

Temperature difference between :math:`T_{\lambda}` at :math:`P` and
:math:`P_r`

.. code:: ipython3

    td = TlPr - Tl

Reference temperature for lower limit of heat capacity integral.

.. code:: ipython3

    tr = Tlref - td

Heat capacity due to the :math:`\lambda`-transition at :math:`T` and
:math:`P`. Valid: :math:`T_r \le T \le T_{\lambda}`.

**Note:** The syntax of the arguments to the SymPy Piecewise expression
is of the form:

::

    (e1,c1), (e2,d2), ..., (eDefault, True)

where the sequence is evaluated left to right. *eN* is the value of the
resulting expression if the condition *cN* is *True*. *cN* is always a
logical comparison. The first *cN* that is *True* provides the value of
the expression. If no *cN* evaluate to *True*, then the resulting
expression becomes *eDefault* given by the last tuple in the sequence.
See the example below.

.. code:: ipython3

    Cpl = (T+td)*(l1+l2*(T+td))**2

Compute the Gibbs free energy of the lambda transition

.. code:: ipython3

    Gl = sym.integrate(Cpl,(T,tr,T)) - T*sym.integrate(Cpl/T,(T,tr,T))

… and add this expression to the model

.. code:: ipython3

    model.add_expression_to_model(Gl, params, exp_type='restricted', lower_limits=(tr,Pr), upper_limits=(Tl,None))

(2) First order phase transition terms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Berman terms valid at T :math:`\ge` :math:`T_t`. Parameters (in this
case :math:`T_t` is equivalent to :math:`T_{\lambda}`: -
:math:`{{\Delta}_t}H`, First order enthalpy contribution at
:math:`T_{\lambda}`

.. code:: ipython3

    deltaHt = sym.symbols('H_t')
    params = [('H_t','J/m', deltaHt)]

:math:`{{\Delta}_t}S = {{\Delta}_t}H/T_{\lambda}`, First order enropy
contribution at :math:`T_{\lambda}`:
:math:`{{\Delta}_t}H - T {{\Delta}_t}S = -(T-T_{\lambda}) {{\Delta}_t}H/T_{\lambda}`

.. code:: ipython3

    GaboveTl = -(T-Tl)*deltaHt/Tl

… and add this expression to the model

.. code:: ipython3

    model.add_expression_to_model(GaboveTl , params, exp_type='restricted', lower_limits=(Tl,None), upper_limits=(Tl,None))

(3) Order-disorder contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameters of the Berman (1988) order-disorder model: - :math:`d_0`,
:math:`d_1`, :math:`d_2`, :math:`d_3`, :math:`d_4`, :math:`d_5`,
order-disorder coefficients from the Berman (1988) model -
:math:`T_{D_{ref}}`, :math:`T_D`, minimum, maximum temperature of
ordering interval, :math:`T_{D_{ref}} \le T \le T_D`

.. code:: ipython3

    d0,d1,d2,d3,d4,d5 = sym.symbols('d0 d1 d2 d3 d4 d5')
    TD,TDref = sym.symbols('T_D T_D_ref')
    params = [('d0','J/K-m', d0), ('d1','J/K^(1/2)-m',d1), ('d2','J-K/m',d2), ('d3','J/K^2-m',d3), 
              ('d4','J/K^3-m',d4), ('d5','bar',d5), ('T_D','K',TD), ('T_D_ref','K',TDref)]

.. code:: ipython3

    CpDs = d0 + d1/sym.sqrt(T) + d2/T**2 + d3*T + d4*T**2
    HDs = sym.integrate(CpDs,(T,TDref,T))
    SDs = sym.integrate(CpDs/T,(T,TDref,T))
    VDs = HDs/d5
    GDs = HDs - T*SDs + VDs*(P-Pr)

.. code:: ipython3

    model.add_expression_to_model(GDs , params, exp_type='restricted', lower_limits=(TDref,None), upper_limits=(TD,None))

Code Print the Model, compile the code and link a Python module
---------------------------------------------------------------

Name the model class

.. code:: ipython3

    model.set_module_name('berman')
    model.get_berman_std_state_database()




.. parsed-literal::

    [(0, 'Akermanite', 'MgCa2Si2O7'),
     (1, 'Albite', 'NaAlSi3O8'),
     (2, 'Ca-Al_Pyroxene', 'CaAl2SiO6'),
     (3, 'Calcite', 'CaCO3'),
     (4, 'Chrysotile', 'Mg3Si2O9H4'),
     (5, 'Clinochlore', 'Mg5Al2Si3O18H8'),
     (6, 'Coesite', 'SiO2'),
     (7, 'Cordierite', 'Mg2Al4Si5O18'),
     (8, 'Corundum', 'Al2O3'),
     (9, 'Alpha_Cristobalite', 'SiO2'),
     (10, 'Beta_Cristobalite', 'SiO2'),
     (11, 'Diaspore', 'AlO2H'),
     (12, 'High_Albite', 'NaAlSi3O8'),
     (13, 'Diopside', 'MgCaSi2O6'),
     (14, 'Dolomite', 'MgCaC2O6'),
     (15, 'Clinoenstatite', 'MgSiO3'),
     (16, 'Orthoenstatite', 'MgSiO3'),
     (17, 'Protoenstatite', 'MgSiO3'),
     (18, 'Fayalite', 'Fe2SiO4'),
     (19, 'Ferrosilite', 'SiFeO3'),
     (20, 'Forsterite', 'Mg2SiO4'),
     (21, 'Gehlenite', 'Al2Ca2SiO7'),
     (22, 'Grossular', 'Ca3Al2Si3O12'),
     (23, 'Low_Albite', 'NaAlSi3O8'),
     (24, 'Hematite', 'Fe2O3'),
     (25, 'Ilmenite', 'FeTiO3'),
     (26, 'Jadeite', 'NaAlSi2O6'),
     (27, 'Kaolinite', 'Al2Si2O9H4'),
     (28, 'Kyanite', 'Al2SiO5'),
     (29, 'Lawsonite', 'CaAl2Si2O10H4'),
     (30, 'Lime', 'CaO'),
     (31, 'Magnesite', 'MgCO3'),
     (32, 'Magnetite', 'Fe3O4'),
     (33, 'Margarite', 'CaAl4Si2O12H2'),
     (34, 'Almandine', 'Si3Fe3Al2O12'),
     (35, 'Meionite', 'Ca4Al6Si6O27C'),
     (36, 'Merwinite', 'Ca3MgSi2O8'),
     (37, 'Monticellite', 'CaMgSiO4'),
     (38, 'Muscovite', 'KAl3Si3O12H2'),
     (39, 'Paragonite', 'NaAl3Si3O12H2'),
     (40, 'Periclase', 'MgO'),
     (41, 'Phlogopite', 'KMg3AlSi3O12H2'),
     (42, 'Potassium_Feldspar', 'KAlSi3O8'),
     (43, 'Sanidine', 'KAlSi3O8'),
     (44, 'Microcline', 'KAlSi3O8'),
     (45, 'Andalusite', 'Al2SiO5'),
     (46, 'Prehnite', 'Ca2Al2Si3O12H2'),
     (47, 'Pyrope', 'Mg3Al2Si3O12'),
     (48, 'Pyrophyllite', 'Al2Si4O12H2'),
     (49, 'A-Quartz', 'SiO2'),
     (50, 'B-Quartz', 'SiO2'),
     (51, 'Rutile', 'TiO2'),
     (52, 'Sillimanite', 'Al2SiO5'),
     (53, 'Sphene', 'CaTiSiO5'),
     (54, 'Spinel', 'MgAl2O4'),
     (55, 'Talc', 'Mg3Si4O12H2'),
     (56, 'Anorthite', 'Al2CaSi2O8'),
     (57, 'Tremolite', 'Ca2Mg5Si8O24H2'),
     (58, 'Low_Tridymite', 'SiO2'),
     (59, 'High_Tridymite', 'SiO2'),
     (60, 'Wollastonite', 'CaSiO3'),
     (61, 'Pseudowollastonite', 'CaSiO3'),
     (62, 'Zoisite', 'Ca2Al3Si3O13H'),
     (63, 'Clinozoisite', 'Ca2Al3Si3O13H'),
     (64, 'Carbon_Dioxide', 'CO2'),
     (65, 'Water', 'H2O'),
     (66, 'Oxygen_Gas', 'O2'),
     (67, 'Anthophyllite', 'Mg7Si8O24H2'),
     (68, 'Sulfur_Gas', 'S2'),
     (69, 'Hydrogen_Gas', 'H2'),
     (70, 'Antigorite', 'Mg48Si34O99H62O48'),
     (71, 'Brucite', 'MgO2H2')]



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

| generates code suitable for model parameter calibration.
| model_type is “fast” or “calib”

.. code:: ipython3

    model_type = "fast"

Index 12 in Berman database is High_Albite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    param_dict = model.get_berman_std_state_database(12)
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


Index 56 in Berman database is Anorthite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    param_dict = model.get_berman_std_state_database(56)
    phase_name = param_dict.pop('Phase', None).title()
    formula = param_dict.pop('Formula', None)
    result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type=model_type)


.. parsed-literal::

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


Index 42 in Berman database is Potassium_Feldspar
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    param_dict = model.get_berman_std_state_database(42)
    phase_name = param_dict.pop('Phase', None).title()
    formula = param_dict.pop('Formula', None)
    result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type=model_type)


.. parsed-literal::

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

    Thu Sep 20 12:23:41 2018
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

    G          -4.207188e+06 J/m       
    dGdT       -5.440218e+02 J/K-m     
    dGdP        1.083742e+01 J/bar-m   
    d2GdP2     -3.203108e-01 J/K^2-m   
    d2GdTdP     2.759511e-04 J/K-bar-m 
    d2GdP2     -1.961311e-05 J/bar^2-m 
    d3GdT3      2.898693e-04 J/K^3-m   
    d3GdT2dP    7.416367e-08 J/K^2-bar-
    d3GdTdP2    0.000000e+00 J/K-bar^2-
    d3GdP3      0.000000e+00 J/bar^3-m 
    S           5.440218e+02 J/K-m     
    V           1.083742e+01 J/bar-m   
    Cv          3.164283e+02 J/K-m     
    Cp          3.203108e+02 J/K-m     
    dCpdT       3.044151e-02 J/K^2-m   
    alpha       2.546279e-05 1/K       
    beta        1.809758e-06 1/bar     
    K           5.525601e+05 bar       
    Kp         -1.000000e+00           


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

Test the generated code against the standard Berman code base
-------------------------------------------------------------

.. code:: ipython3

    from thermoengine import model

.. code:: ipython3

    bermanDB = model.Database()

.. code:: ipython3

    abbrv = ""
    for full_name, abbrv in zip(bermanDB.phase_info.phase_name,bermanDB.phase_info.abbrev):
        if phase_name == full_name:
            break
    refModel = bermanDB.get_phase(abbrv)

.. code:: ipython3

    import math
    fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:13.6e} {4:6.2f}%"
    fmts = "{0:<10.10s} {1:13.6e}"
    try:
        x = berman.cy_Potassium_Feldspar_berman_g(t,p)
        y = refModel.gibbs_energy(t,p)
        print(fmt.format('G', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_dgdt(t,p)
        y = -refModel.entropy(t,p)
        print(fmt.format('dGdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_dgdp(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('dGdP', x, y, x-y, 100.0*math.fabs((x-y)/y))) 
        x = berman.cy_Potassium_Feldspar_berman_d2gdt2(t,p)
        print(fmts.format('d2GdT2', x))
        x = berman.cy_Potassium_Feldspar_berman_d2gdtdp(t,p)
        print(fmts.format('d2GdTdP', x))
        x = berman.cy_Potassium_Feldspar_berman_d2gdp2(t,p)
        print(fmts.format('d2GdP2', x))
        x = berman.cy_Potassium_Feldspar_berman_d3gdt3(t,p)
        print(fmts.format('d3GdT3', x))
        x = berman.cy_Potassium_Feldspar_berman_d3gdt2dp(t,p)
        print(fmts.format('d3GdT2dP', x))
        x = berman.cy_Potassium_Feldspar_berman_d3gdtdp2(t,p)
        print(fmts.format('d3GdTdP2', x))
        x = berman.cy_Potassium_Feldspar_berman_d3gdp3(t,p)
        print(fmts.format('d3GdP3', x))
        x = berman.cy_Potassium_Feldspar_berman_s(t,p)
        y = refModel.entropy(t,p)
        print(fmt.format('S', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_v(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('V', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_cv(t,p)
        print(fmts.format('Cv', x))
        x = berman.cy_Potassium_Feldspar_berman_cp(t,p)
        y = refModel.heat_capacity(t,p)
        print(fmt.format('Cp', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_dcpdt(t,p)
        print(fmts.format('dCpdT', x))
        x = berman.cy_Potassium_Feldspar_berman_alpha(t,p)
        print(fmts.format('alpha', x))
        x = berman.cy_Potassium_Feldspar_berman_beta(t,p)
        print(fmts.format('beta', x))
        x = berman.cy_Potassium_Feldspar_berman_K(t,p)
        print(fmts.format('K', x))
        x = berman.cy_Potassium_Feldspar_berman_Kp(t,p)
        print(fmts.format('Kp', x))
    except AttributeError:
        pass
    try:
        x = berman.cy_Potassium_Feldspar_berman_calib_g(t,p)
        y = refModel.gibbs_energy(t,p)
        print(fmt.format('G', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_calib_dgdt(t,p)
        y = -refModel.entropy(t,p)
        print(fmt.format('dGdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_calib_dgdp(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('dGdP', x, y, x-y, 100.0*math.fabs((x-y)/y))) 
        x = berman.cy_Potassium_Feldspar_berman_calib_d2gdt2(t,p)
        print(fmts.format('d2GdT2', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_d2gdtdp(t,p)
        print(fmts.format('d2GdTdP', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_d2gdp2(t,p)
        print(fmts.format('d2GdP2', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_d3gdt3(t,p)
        print(fmts.format('d3GdT3', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_d3gdt2dp(t,p)
        print(fmts.format('d3GdT2dP', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_d3gdtdp2(t,p)
        print(fmts.format('d3GdTdP2', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_d3gdp3(t,p)
        print(fmts.format('d3GdP3', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_s(t,p)
        y = refModel.entropy(t,p)
        print(fmt.format('S', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_calib_v(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('V', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_calib_cv(t,p)
        print(fmts.format('Cv', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_cp(t,p)
        y = refModel.heat_capacity(t,p)
        print(fmt.format('Cp', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = berman.cy_Potassium_Feldspar_berman_calib_dcpdt(t,p)
        print(fmts.format('dCpdT', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_alpha(t,p)
        print(fmts.format('alpha', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_beta(t,p)
        print(fmts.format('beta', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_K(t,p)
        print(fmts.format('K', x))
        x = berman.cy_Potassium_Feldspar_berman_calib_Kp(t,p)
        print(fmts.format('Kp', x))
    except AttributeError:
        pass


.. parsed-literal::

    G          -4.207188e+06 -4.207218e+06  3.024521e+01   0.00%
    dGdT       -5.440218e+02 -5.440935e+02  7.166392e-02   0.01%
    dGdP        1.083742e+01  1.083799e+01 -5.670119e-04   0.01%
    d2GdT2     -3.203108e-01
    d2GdTdP     2.759511e-04
    d2GdP2     -1.961311e-05
    d3GdT3      2.898693e-04
    d3GdT2dP    7.416367e-08
    d3GdTdP2    0.000000e+00
    d3GdP3      0.000000e+00
    S           5.440218e+02  5.440935e+02 -7.166392e-02   0.01%
    V           1.083742e+01  1.083799e+01 -5.670119e-04   0.01%
    Cv          3.164283e+02
    Cp          3.203108e+02  3.203555e+02 -4.466899e-02   0.01%
    dCpdT       3.044151e-02
    alpha       2.546279e-05
    beta        1.809758e-06
    K           5.525601e+05
    Kp         -1.000000e+00


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

    152 ns ± 0.724 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)


.. code:: ipython3

    %timeit(refModel.gibbs_energy(t,p))


.. parsed-literal::

    79.9 µs ± 310 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)

