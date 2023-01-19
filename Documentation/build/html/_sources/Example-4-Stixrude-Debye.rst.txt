Helmholtz energy (Stixrude - Debye integral)
============================================

Required system packages and initialization

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()

.. code:: ipython3

    from thermoengine import coder
    from thermoengine.coder import Debye as db

Declare *T* and *V* to be independent variables
-----------------------------------------------

This specifies that the model expression will be defined in terms of the
Helmholtz free energy. Note, that the defauilt model type, *TP*,
decalres that the model expression will be the Gibbs free energy.

.. code:: ipython3

    model = coder.StdStateModel(model_type='TV')

.. code:: ipython3

    T = model.get_symbol_for_t()
    V = model.get_symbol_for_v()
    Tr = model.get_symbol_for_tr()
    Vr = model.get_symbol_for_vr()

Define model expressions applicable over all of T,P space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An expression for the Gibbs free energy, :math:`G(T,P)` or the Helmholtz
energy :math:`A(T,V)` is constructed. The expression may have multiple
parts. Often the heat capacity function is postulated, then integrated
to yield expressions for the entahlpy, entropy, and in combination the
energy potential. Then, an equation of state (EOS) is adopted and that
term is integrated in pressure or volume and added to the heat capacity
integrals. This proceedure is follwed here. #### (1) Helmholtz free
energy Declare parameters of the Stixrude standard state model:

.. code:: ipython3

    a0,n,v0,k00,k0p,theta0,gamma0,q,refS,R = sym.symbols('a0 n v0 k00 k0p theta0 gamma0 q refS R')

Specify paramters …

.. code:: ipython3

    params = [('a0','J/m',a0), ('n','',n), ('v0','J/bar-m',v0), ('k00','bar',k00), ('k0p','',k0p),  ('theta0','K',theta0), ('gamma0', '',gamma0), ('q', '', q), ('refS', 'J/K-m', refS), ('R', 'J/K-m', R)]

Define the Debye temperature:

.. code:: ipython3

    c1 = sym.S(9)*k00*v0
    c2 = k0p/sym.S(2) - sym.S(2)
    c5 = sym.S(3)*gamma0
    c7 = c5*(-sym.S(2) + sym.S(6)*gamma0 - sym.S(3)*q)
    f = (v0/V)**(sym.S(2)/sym.S(3))/sym.S(2) - sym.S(1)/sym.S(2)
    d0 = theta0*(sym.S(1) + c7*f*f + sym.S(2)*c5*f)**(sym.S(1)/sym.S(2))
    d0




.. math::

    \displaystyle \theta_{0} \sqrt{3 \gamma_{0} \left(\frac{\left(\frac{v_{0}}{V}\right)^{\frac{2}{3}}}{2} - \frac{1}{2}\right)^{2} \left(6 \gamma_{0} - 3 q - 2\right) + 6 \gamma_{0} \left(\frac{\left(\frac{v_{0}}{V}\right)^{\frac{2}{3}}}{2} - \frac{1}{2}\right) + 1}



| Define the Debye Helmholtz free energy …
| db(x) returns the Debye *integral* with upper limit *x*

.. code:: ipython3

    x = d0/T
    A_db = n*R*T*(sym.S(3)*sym.ln(sym.S(1)-sym.exp(-x)) - db(x))

… and from that the quasiharmonic approximation to the Helmholtz energy
…

.. code:: ipython3

    A_quasi = A_db - A_db.subs(T,Tr)

… and finally the Stixrude model expression for the Helmholtz free
energy:

.. code:: ipython3

    A = a0 + c1*f*f*(sym.S(1)/sym.S(2)+c2*f) + A_quasi

… and add this expression to the model

.. code:: ipython3

    model.add_expression_to_model(A, params)

Code Print the Model, compile the code and link a Python module
---------------------------------------------------------------

Name the model class

.. code:: ipython3

    model.set_module_name('stixrude')
    model.set_include_debye_code(include=True)

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

    param_dict = {'a0':-2055.0*1000.0, 'n':7.0, 'v0':43.6/10.0, 'k00':128.0*10000.0, 'k0p':4.2,
                  'theta0':809.0, 'gamma0':0.99, 'q':2.1, 'refS':0.0, 'R':8.314472, 'T_r':300.00, 
                  'V_r':43.6/10.0}
    phase_name = 'Forsterite'
    formula = 'Mg(2)Si(1)O(4)'
    model.set_reference_origin(Vr=param_dict['V_r'])
    result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict)
    #result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type='calib')


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
    Success! Import the module named  stixrude


.. code:: ipython3

    param_dict




.. parsed-literal::

    {'a0': -2055000.0,
     'n': 7.0,
     'v0': 4.36,
     'k00': 1280000.0,
     'k0p': 4.2,
     'theta0': 809.0,
     'gamma0': 0.99,
     'q': 2.1,
     'refS': 0.0,
     'R': 8.314472,
     'T_r': 300.0,
     'V_r': 4.36}



Import the new module and test the model
----------------------------------------

.. code:: ipython3

    import stixrude
    %cd ..


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/Cython/Compiler/Main.py:369: FutureWarning: Cython directive 'language_level' not set, using 2 for now (Py2). This will change in a later release! File: /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working/stixrude.pyx
      tree = Parsing.p_module(s, pxd, full_module_name)


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
        print(stixrude.cy_Forsterite_stixrude_identifier())
        print(stixrude.cy_Forsterite_stixrude_name())
        print(stixrude.cy_Forsterite_stixrude_formula())
        print(stixrude.cy_Forsterite_stixrude_mw())
        print(stixrude.cy_Forsterite_stixrude_elements())
    except AttributeError:
        pass
    try:
        print(stixrude.cy_Forsterite_stixrude_calib_identifier())
        print(stixrude.cy_Forsterite_stixrude_calib_name())
        print(stixrude.cy_Forsterite_stixrude_calib_formula())
        print(stixrude.cy_Forsterite_stixrude_calib_mw())
        print(stixrude.cy_Forsterite_stixrude_calib_elements())
    except AttributeError:
        pass


.. parsed-literal::

    Wed Sep 23 09:56:59 2020
    Forsterite
    Mg2SiO4
    140.6931
    [0. 0. 0. 0. 0. 0. 0. 0. 4. 0. 0. 0. 2. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]


Execute the standard thermodynamic property retrieval functions:

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
    try:
        print(fmt.format('G', stixrude.cy_Forsterite_stixrude_g(t,p), 'J/m'))
        print(fmt.format('dGdT', stixrude.cy_Forsterite_stixrude_dgdt(t,p), 'J/K-m'))
        print(fmt.format('dGdP', stixrude.cy_Forsterite_stixrude_dgdp(t,p), 'J/bar-m'))
        print(fmt.format('d2GdP2', stixrude.cy_Forsterite_stixrude_d2gdt2(t,p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', stixrude.cy_Forsterite_stixrude_d2gdtdp(t,p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', stixrude.cy_Forsterite_stixrude_d2gdp2(t,p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', stixrude.cy_Forsterite_stixrude_d3gdt3(t,p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', stixrude.cy_Forsterite_stixrude_d3gdt2dp(t,p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', stixrude.cy_Forsterite_stixrude_d3gdtdp2(t,p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', stixrude.cy_Forsterite_stixrude_d3gdp3(t,p), 'J/bar^3-m'))
        print(fmt.format('S', stixrude.cy_Forsterite_stixrude_s(t,p), 'J/K-m'))
        print(fmt.format('V', stixrude.cy_Forsterite_stixrude_v(t,p), 'J/bar-m'))
        print(fmt.format('Cv', stixrude.cy_Forsterite_stixrude_cv(t,p), 'J/K-m'))
        print(fmt.format('Cp', stixrude.cy_Forsterite_stixrude_cp(t,p), 'J/K-m'))
        print(fmt.format('dCpdT', stixrude.cy_Forsterite_stixrude_dcpdt(t,p), 'J/K^2-m'))
        print(fmt.format('alpha', stixrude.cy_Forsterite_stixrude_alpha(t,p), '1/K'))
        print(fmt.format('beta', stixrude.cy_Forsterite_stixrude_beta(t,p), '1/bar'))
        print(fmt.format('K', stixrude.cy_Forsterite_stixrude_K(t,p), 'bar'))
        print(fmt.format('Kp', stixrude.cy_Forsterite_stixrude_Kp(t,p), ''))
    except AttributeError:
        pass
    try:
        print(fmt.format('G', stixrude.cy_Forsterite_stixrude_calib_g(t,p), 'J/m'))
        print(fmt.format('dGdT', stixrude.cy_Forsterite_stixrude_calib_dgdt(t,p), 'J/K-m'))
        print(fmt.format('dGdP', stixrude.cy_Forsterite_stixrude_calib_dgdp(t,p), 'J/bar-m'))
        print(fmt.format('d2GdP2', stixrude.cy_Forsterite_stixrude_calib_d2gdt2(t,p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', stixrude.cy_Forsterite_stixrude_calib_d2gdtdp(t,p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', stixrude.cy_Forsterite_stixrude_calib_d2gdp2(t,p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', stixrude.cy_Forsterite_stixrude_calib_d3gdt3(t,p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', stixrude.cy_Forsterite_stixrude_calib_d3gdt2dp(t,p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', stixrude.cy_Forsterite_stixrude_calib_d3gdtdp2(t,p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', stixrude.cy_Forsterite_stixrude_calib_d3gdp3(t,p), 'J/bar^3-m'))
        print(fmt.format('S', stixrude.cy_Forsterite_stixrude_calib_s(t,p), 'J/K-m'))
        print(fmt.format('V', stixrude.cy_Forsterite_stixrude_calib_v(t,p), 'J/bar-m'))
        print(fmt.format('Cv', stixrude.cy_Forsterite_stixrude_calib_cv(t,p), 'J/K-m'))
        print(fmt.format('Cp', stixrude.cy_Forsterite_stixrude_calib_cp(t,p), 'J/K-m'))
        print(fmt.format('dCpdT', stixrude.cy_Forsterite_stixrude_calib_dcpdt(t,p), 'J/K^2-m'))
        print(fmt.format('alpha', stixrude.cy_Forsterite_stixrude_calib_alpha(t,p), '1/K'))
        print(fmt.format('beta', stixrude.cy_Forsterite_stixrude_calib_beta(t,p), '1/bar'))
        print(fmt.format('K', stixrude.cy_Forsterite_stixrude_calib_K(t,p), 'bar'))
        print(fmt.format('Kp', stixrude.cy_Forsterite_stixrude_calib_Kp(t,p), ''))
    except AttributeError:
        pass


.. parsed-literal::

    G          -2.147813e+06 J/m       
    dGdT       -2.653920e+02 J/K-m     
    dGdP        4.410375e+00 J/bar-m   
    d2GdP2     -1.624659e-01 J/K^2-m   
    d2GdTdP     1.347595e-04 J/K-bar-m 
    d2GdP2     -3.718165e-06 J/bar^2-m 
    d3GdT3      1.522700e-04 J/K^3-m   
    d3GdT2dP    3.519210e-08 J/K^2-bar-
    d3GdTdP2   -7.534594e-10 J/K-bar^2-
    d3GdP3      1.715513e-11 J/bar^3-m 
    S           2.653920e+02 J/K-m     
    V           4.410375e+00 J/bar-m   
    Cv          1.575818e+02 J/K-m     
    Cp          1.624659e+02 J/K-m     
    dCpdT       1.019590e-02 J/K^2-m   
    alpha       3.055512e-05 1/K       
    beta        8.430497e-07 1/bar     
    K           1.186170e+06 bar       
    Kp          4.472832e+00           


Available only in the “Calib” versions of generated code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Execute the parameter value/metadata functions.
| These functions are only defined for the “calibration” model code
  implementation:

.. code:: ipython3

    try:
        np = stixrude.cy_Forsterite_stixrude_get_param_number()
        names = stixrude.cy_Forsterite_stixrude_get_param_names()
        units = stixrude.cy_Forsterite_stixrude_get_param_units()
        values = stixrude.cy_Forsterite_stixrude_get_param_values()
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], stixrude.cy_Forsterite_stixrude_get_param_value(i), units[i]))
    except AttributeError:
        pass

Test the functions that allow modification of the array of parameter
values

.. code:: ipython3

    try:
        values[1] = 100.0
        stixrude.cy_Forsterite_stixrude_set_param_values(values)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], stixrude.cy_Forsterite_stixrude_get_param_value(i), units[i]))
    except (AttributeError, NameError):
        pass

Test the functions that allow modification of a particular parameter
value

.. code:: ipython3

    try:
        stixrude.cy_Forsterite_stixrude_set_param_value(1, 1.0)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], stixrude.cy_Forsterite_stixrude_get_param_value(i), units[i]))
    except AttributeError:
        pass

Evaluate parameter derivatives …

.. code:: ipython3

    try:
        fmt = "    {0:<10.10s} {1:13.6e}"
        for i in range(0, np):
            print ('Derivative with respect to parameter: ', names[i], ' of')
            print (fmt.format('G', stixrude.cy_Forsterite_stixrude_dparam_g(t, p, i)))
            print (fmt.format('dGdT', stixrude.cy_Forsterite_stixrude_dparam_dgdt(t, p, i)))
            print (fmt.format('dGdP', stixrude.cy_Forsterite_stixrude_dparam_dgdp(t, p, i)))
            print (fmt.format('d2GdT2', stixrude.cy_Forsterite_stixrude_dparam_d2gdt2(t, p, i)))
            print (fmt.format('d2GdTdP', stixrude.cy_Forsterite_stixrude_dparam_d2gdtdp(t, p, i)))
            print (fmt.format('d2GdP2', stixrude.cy_Forsterite_stixrude_dparam_d2gdp2(t, p, i)))
            print (fmt.format('d3GdT3', stixrude.cy_Forsterite_stixrude_dparam_d3gdt3(t, p, i)))
            print (fmt.format('d3GdT2dP', stixrude.cy_Forsterite_stixrude_dparam_d3gdt2dp(t, p, i)))
            print (fmt.format('d3GdTdP2', stixrude.cy_Forsterite_stixrude_dparam_d3gdtdp2(t, p, i)))
            print (fmt.format('d3GdP3', stixrude.cy_Forsterite_stixrude_dparam_d3gdp3(t, p, i)))
    except (AttributeError, TypeError):
        pass

Time execution of the code
--------------------------

.. code:: ipython3

    try:
        %timeit stixrude.cy_Forsterite_stixrude_calib_g(t,p)
    except AttributeError:
        pass
    try:
        %timeit stixrude.cy_Forsterite_stixrude_g(t,p)
    except AttributeError:
        pass


.. parsed-literal::

    264 ns ± 8.62 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)


Test the generated code against the standard Stixrude code base
---------------------------------------------------------------

.. code:: ipython3

    from thermoengine import model as md
    stixrudeDB = md.Database(database="Stixrude")

.. code:: ipython3

    abbrv = ""
    for full_name, abbrv in zip(stixrudeDB.phase_info.phase_name,stixrudeDB.phase_info.abbrev):
        if phase_name == full_name:
            break
    refModel = stixrudeDB.get_phase(abbrv)

.. code:: ipython3

    import math
    fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:13.6e} {4:6.2f}%"
    fmts = "{0:<10.10s} {1:13.6e}"
    try:
        x = stixrude.cy_Forsterite_stixrude_g(t,p)
        y = refModel.gibbs_energy(t,p)
        print(fmt.format('G', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_dgdt(t,p)
        y = -refModel.entropy(t,p)
        print(fmt.format('dGdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_dgdp(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('dGdP', x, y, x-y, 100.0*math.fabs((x-y)/y))) 
        x = stixrude.cy_Forsterite_stixrude_d2gdt2(t,p)
        print(fmts.format('d2GdT2', x))
        x = stixrude.cy_Forsterite_stixrude_d2gdtdp(t,p)
        print(fmts.format('d2GdTdP', x))
        x = stixrude.cy_Forsterite_stixrude_d2gdp2(t,p)
        print(fmts.format('d2GdP2', x))
        x = stixrude.cy_Forsterite_stixrude_d3gdt3(t,p)
        print(fmts.format('d3GdT3', x))
        x = stixrude.cy_Forsterite_stixrude_d3gdt2dp(t,p)
        print(fmts.format('d3GdT2dP', x))
        x = stixrude.cy_Forsterite_stixrude_d3gdtdp2(t,p)
        print(fmts.format('d3GdTdP2', x))
        x = stixrude.cy_Forsterite_stixrude_d3gdp3(t,p)
        print(fmts.format('d3GdP3', x))
        x = stixrude.cy_Forsterite_stixrude_s(t,p)
        y = refModel.entropy(t,p)
        print(fmt.format('S', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_v(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('V', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_cv(t,p)
        print(fmts.format('Cv', x))
        x = stixrude.cy_Forsterite_stixrude_cp(t,p)
        y = refModel.heat_capacity(t,p)
        print(fmt.format('Cp', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_dcpdt(t,p)
        print(fmts.format('dCpdT', x))
        x = stixrude.cy_Forsterite_stixrude_alpha(t,p)
        print(fmts.format('alpha', x))
        x = stixrude.cy_Forsterite_stixrude_beta(t,p)
        print(fmts.format('beta', x))
        x = stixrude.cy_Forsterite_stixrude_K(t,p)
        print(fmts.format('K', x))
        x = stixrude.cy_Forsterite_stixrude_Kp(t,p)
        print(fmts.format('Kp', x))
    except AttributeError:
        pass
    try:
        x = stixrude.cy_Forsterite_stixrude_calib_g(t,p)
        y = refModel.gibbs_energy(t,p)
        print(fmt.format('G', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_calib_dgdt(t,p)
        y = -refModel.entropy(t,p)
        print(fmt.format('dGdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_calib_dgdp(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('dGdP', x, y, x-y, 100.0*math.fabs((x-y)/y))) 
        x = stixrude.cy_Forsterite_stixrude_calib_d2gdt2(t,p)
        print(fmts.format('d2GdT2', x))
        x = stixrude.cy_Forsterite_stixrude_calib_d2gdtdp(t,p)
        print(fmts.format('d2GdTdP', x))
        x = stixrude.cy_Forsterite_stixrude_calib_d2gdp2(t,p)
        print(fmts.format('d2GdP2', x))
        x = stixrude.cy_Forsterite_stixrude_calib_d3gdt3(t,p)
        print(fmts.format('d3GdT3', x))
        x = stixrude.cy_Forsterite_stixrude_calib_d3gdt2dp(t,p)
        print(fmts.format('d3GdT2dP', x))
        x = stixrude.cy_Forsterite_stixrude_calib_d3gdtdp2(t,p)
        print(fmts.format('d3GdTdP2', x))
        x = stixrude.cy_Forsterite_stixrude_calib_d3gdp3(t,p)
        print(fmts.format('d3GdP3', x))
        x = stixrude.cy_Forsterite_stixrude_calib_s(t,p)
        y = refModel.entropy(t,p)
        print(fmt.format('S', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_calib_v(t,p)
        y = refModel.volume(t,p)
        print(fmt.format('V', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_calib_cv(t,p)
        print(fmts.format('Cv', x))
        x = stixrude.cy_Forsterite_stixrude_calib_cp(t,p)
        y = refModel.heat_capacity(t,p)
        print(fmt.format('Cp', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = stixrude.cy_Forsterite_stixrude_calib_dcpdt(t,p)
        print(fmts.format('dCpdT', x))
        x = stixrude.cy_Forsterite_stixrude_calib_alpha(t,p)
        print(fmts.format('alpha', x))
        x = stixrude.cy_Forsterite_stixrude_calib_beta(t,p)
        print(fmts.format('beta', x))
        x = stixrude.cy_Forsterite_stixrude_calib_K(t,p)
        print(fmts.format('K', x))
        x = stixrude.cy_Forsterite_stixrude_calib_Kp(t,p)
        print(fmts.format('Kp', x))
    except AttributeError:
        pass


.. parsed-literal::

    G          -2.147813e+06 -2.148196e+06  3.829654e+02   0.02%
    dGdT       -2.653920e+02 -2.747216e+02  9.329571e+00   3.40%
    dGdP        4.410375e+00  4.414971e+00 -4.596333e-03   0.10%
    d2GdT2     -1.624659e-01
    d2GdTdP     1.347595e-04
    d2GdP2     -3.718165e-06
    d3GdT3      1.522700e-04
    d3GdT2dP    3.519210e-08
    d3GdTdP2   -7.534594e-10
    d3GdP3      1.715513e-11
    S           2.653920e+02  2.747216e+02 -9.329571e+00   3.40%
    V           4.410375e+00  4.414971e+00 -4.596333e-03   0.10%
    Cv          1.575818e+02
    Cp          1.624659e+02  1.748649e+02 -1.239895e+01   7.09%
    dCpdT       1.019590e-02
    alpha       3.055512e-05
    beta        8.430497e-07
    K           1.186170e+06
    Kp          4.472832e+00


