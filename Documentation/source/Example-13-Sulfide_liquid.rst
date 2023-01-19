Sulfide liquid model of Kress
=============================

This notebook shows how to import the sulfide_liquid module of the
ThermoEngine package in order to calculate the thermodynamic properties
of liquid solutions in the system O-S-Fe-Ni-Cu. The model used is from
Kress (Contributions to Mineralogy and Petrology, 156(6):785-797). DOI:
10.1007/s00410-008-0315-z).

The ThermoEngine sulfide_liquid module emulates a “coder” solution phase
by wrapping the Python package, SulfLiq, which is available on PyPI.org
and as a `GitLab project <https://gitlab.com/ENKI-portal/sulfliq>`__.
The underlying code in the SulfLiq module is written in C++ (by Kress),
which is wrapped using
`PyBind11 <https://pybind11.readthedocs.io/en/stable/>`__ to expose a
Python API.

.. code:: ipython3

    from thermoengine import core, model, phases, sulfide_liquid
    import numpy as np

Test module function calls
--------------------------

This section makes function calls directly to the sulfide_liquid module
API. All methods are tested, and all coded functions are “calibration”
versions. Note that this method of accessing the module functions is not
commonly used. See below for an example of how to load the module into
the ThermoEngine package and access the capabilities of the
sulfide_liquid module using the standard database and phase properties
methods.

Set temperature (K), pressure (bars) and a composition for testing:
moles of (O, S, Fe, Ni, Cu) = (2.0, 5.0, 10.0, 1.0, 1.0)

.. code:: ipython3

    t = 2000.00
    p = 1.0
    n = np.array([2., 5., 10., 1., 1.])

Characteristics
~~~~~~~~~~~~~~~

.. code:: ipython3

    print(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_identifier())
    print(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_name())
    print(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_formula(t,p,n))


.. parsed-literal::

    Version_1_0_0
    Sulfide Liquid
    O2.0S5.0Fe10.0Ni1.0Cu1.0


Composition conversion methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    c = sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_number()
    e = np.zeros(106)
    sum = np.sum(n)
    for index in range(0,c):
        end = sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_elements(index)
        for i in range(0,106):
            e[i] += end[i]*n[index]/sum
    nConv = sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_moles(e)
    for i in range(0,c):
        print ('X[{0:d}] input {1:13.6e}, calc {2:13.6e}, diff {3:13.6e}'.format(
            i, n[i]/sum, nConv[i], nConv[i]-n[i]/sum))
    if not sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_test_moles(nConv):
        print ('Output of intrinsic composition calculation fails tests for permissible values.')


.. parsed-literal::

    X[0] input  1.052632e-01, calc  1.052632e-01, diff  0.000000e+00
    X[1] input  2.631579e-01, calc  2.631579e-01, diff  0.000000e+00
    X[2] input  5.263158e-01, calc  5.263158e-01, diff  0.000000e+00
    X[3] input  5.263158e-02, calc  5.263158e-02, diff  0.000000e+00
    X[4] input  5.263158e-02, calc  5.263158e-02, diff  0.000000e+00


.. code:: ipython3

    print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_moles_to_tot_moles(n))
    print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_moles_to_mole_frac(n))
    e = sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_moles_to_elm(n)
    print (e)
    print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_moles(e))
    print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_tot_moles(e))
    print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_tot_grams(e))


.. parsed-literal::

    19.0
    [0.10526316 0.26315789 0.52631579 0.05263158 0.05263158]
    [ 0.  0.  0.  0.  0.  0.  0.  0.  2.  0.  0.  0.  0.  0.  0.  0.  5.  0.
      0.  0.  0.  0.  0.  0.  0.  0. 10.  0.  1.  1.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
      0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
    [ 2.  5. 10.  1.  1.]
    19.0
    873.0248000000001


.. code:: ipython3

    try:
        sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_g(t,p,n)
    except:
        print ("Exception generated on first call to function.")

Execute the standard thermodynamic property retrieval functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
    print(fmt.format('G', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_g(t,p,n), 'J'))
    print(fmt.format('dGdT', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_dgdt(t,p,n), 'J/K'))
    print(fmt.format('dGdP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_dgdp(t,p,n), 'J/bar'))
    print(fmt.format('d2GdT2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d2gdt2(t,p,n), 'J/K^2'))
    print(fmt.format('d2GdTdP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d2gdtdp(t,p,n), 'J/K-bar'))
    print(fmt.format('d2GdP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d2gdp2(t,p,n), 'J/bar^2'))
    print(fmt.format('d3GdT3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdt3(t,p,n), 'J/K^3'))
    print(fmt.format('d3GdT2dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdt2dp(t,p,n), 'J/K^2-bar'))
    print(fmt.format('d3GdTdP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdtdp2(t,p,n), 'J/K-bar^2'))
    print(fmt.format('d3GdP3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdp3(t,p,n), 'J/bar^3'))
    print(fmt.format('S', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_s(t,p,n), 'J/K'))
    print(fmt.format('V', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_v(t,p,n), 'J/bar'))
    print(fmt.format('Cv', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_cv(t,p,n), 'J/K'))
    print(fmt.format('Cp', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_cp(t,p,n), 'J/K'))
    print(fmt.format('dCpdT', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_dcpdt(t,p,n), 'J/K^2'))
    print(fmt.format('alpha', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_alpha(t,p,n), '1/K'))
    print(fmt.format('beta', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_beta(t,p,n), '1/bar'))
    print(fmt.format('K', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_K(t,p,n), 'bar'))
    print(fmt.format('Kp', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_Kp(t,p,n), ''))


.. parsed-literal::

    G          -3.878956e+06 J         
    dGdT       -2.452130e+03 J/K       
    dGdP        1.884016e+01 J/bar     
    d2GdT2     -3.800446e-01 J/K^2     
    d2GdTdP     2.350344e-03 J/K-bar   
    d2GdP2     -2.574264e-04 J/bar^2   
    d3GdT3      1.880394e-04 J/K^3     
    d3GdT2dP    1.310298e-05 J/K^2-bar 
    d3GdTdP2   -8.854753e-07 J/K-bar^2 
    d3GdP3      6.631852e-08 J/bar^3   
    S           2.452130e+03 J/K       
    V           1.884016e+01 J/bar     
    Cv          7.171711e+02 J/K       
    Cp          7.600892e+02 J/K       
    dCpdT       3.965713e-03 J/K^2     
    alpha       1.247519e-04 1/K       
    beta        1.366371e-05 1/bar     
    K           7.318657e+04 bar       
    Kp          1.785441e+01           


Execute functions that access endmember properties:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<15.15s}"
    print ("number of components", sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_number())
    for index in range(0, c):
        print ("{0:<20.20s}".format(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_name(index)), end=' ')
        print ("{0:<20.20s}".format(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_formula(index)), end=' ')
        print ("mw: {0:10.2f}".format(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_mw(index)))
        print (fmt.format('mu0', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_mu0(index,t,p), 'J/mol'))
        print (fmt.format('dmu0dT', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_dmu0dT(index,t,p), 'J/K-mol'))
        print (fmt.format('dmu0dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_dmu0dP(index,t,p), 'J/bar-mol'))
        print (fmt.format('d2mu0dT2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d2mu0dT2(index,t,p), 'J/K^2-mol'))
        print (fmt.format('d2mu0dTdP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d2mu0dTdP(index,t,p), 'J/K-bar-mol'))
        print (fmt.format('d2mu0dP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d2mu0dP2(index,t,p), 'J/bar^2-mol'))
        print (fmt.format('d3mu0dT3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dT3(index,t,p), 'J/K^3-mol'))
        print (fmt.format('d3mu0dT2dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dT2dP(index,t,p), 'J/K^2-bar-mol'))
        print (fmt.format('d3mu0dTdP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dTdP2(index,t,p), 'J/K-bar^2-mol'))
        print (fmt.format('d3mu0dP3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dP3(index,t,p), 'J/bar^3-mol'))
        print ("Element array:")
        print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_endmember_elements(index))
        print ()


.. parsed-literal::

    number of components 5
    O                    O                    mw:      16.00
    mu0        -1.549132e+05 J/mol          
    dmu0dT     -7.949271e+01 J/K-mol        
    dmu0dP      6.621094e-06 J/bar-mol      
    d2mu0dT2    0.000000e+00 J/K^2-mol      
    d2mu0dTdP   2.208563e-01 J/K-bar-mol    
    d2mu0dP2    4.443341e-03 J/bar^2-mol    
    d3mu0dT3    0.000000e+00 J/K^3-mol      
    d3mu0dT2dP  0.000000e+00 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    0.000000e+00 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    S                    S                    mw:      32.06
    mu0        -1.404959e+05 J/mol          
    dmu0dT     -9.945635e+01 J/K-mol        
    dmu0dP      3.267578e-05 J/bar-mol      
    d2mu0dT2   -3.276800e-02 J/K^2-mol      
    d2mu0dTdP   1.095107e+00 J/K-bar-mol    
    d2mu0dP2    2.192835e-02 J/bar^2-mol    
    d3mu0dT3    0.000000e+00 J/K^3-mol      
    d3mu0dT2dP  0.000000e+00 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    0.000000e+00 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Fe                   Fe                   mw:      55.85
    mu0        -1.276421e+05 J/mol          
    dmu0dT     -1.044423e+02 J/K-mol        
    dmu0dP      1.626953e-05 J/bar-mol      
    d2mu0dT2   -8.192000e-02 J/K^2-mol      
    d2mu0dTdP   5.442765e-01 J/K-bar-mol    
    d2mu0dP2    1.091830e-02 J/bar^2-mol    
    d3mu0dT3    0.000000e+00 J/K^3-mol      
    d3mu0dT2dP  0.000000e+00 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    0.000000e+00 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Ni                   Ni                   mw:      58.71
    mu0        -1.277856e+05 J/mol          
    dmu0dT     -1.014425e+02 J/K-mol        
    dmu0dP      1.543945e-05 J/bar-mol      
    d2mu0dT2    3.276800e-02 J/K^2-mol      
    d2mu0dTdP   5.200282e-01 J/K-bar-mol    
    d2mu0dP2    1.036124e-02 J/bar^2-mol    
    d3mu0dT3    0.000000e+00 J/K^3-mol      
    d3mu0dT2dP  0.000000e+00 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    0.000000e+00 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Cu                   Cu                   mw:      63.55
    mu0        -1.294941e+05 J/mol          
    dmu0dT     -9.668906e+01 J/K-mol        
    dmu0dP      1.678711e-05 J/bar-mol      
    d2mu0dT2    3.276800e-02 J/K^2-mol      
    d2mu0dTdP   5.639373e-01 J/K-bar-mol    
    d2mu0dP2    1.126564e-02 J/bar^2-mol    
    d3mu0dT3    0.000000e+00 J/K^3-mol      
    d3mu0dT2dP  0.000000e+00 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    0.000000e+00 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    


Execute functions that access species properties:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print ("number of species", sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_species_number())
    for index in range(0, sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_species_number()):
        print ("{0:<20.20s}".format(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_species_name(index)), end=' ')
        print ("{0:<20.20s}".format(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_species_formula(index)), end=' ')
        print ("mw: {0:10.2f}".format(sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_species_mw(index)))
        print ("Element array:")
        print (sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_species_elements(index))
        print ()


.. parsed-literal::

    number of species 15
    O                    O                    mw:      16.00
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    S                    S                    mw:      32.06
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Fe                   Fe                   mw:      55.85
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Ni                   Ni                   mw:      58.71
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Cu                   Cu                   mw:      63.55
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    FeO                  FeO                  mw:      71.85
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    FeO1.5               FeO1.5               mw:      79.85
    Element array:
    [0.  0.  0.  0.  0.  0.  0.  0.  1.5 0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0. ]
    
    NiO                  NiO                  mw:      74.71
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    FeS                  FeS                  mw:      87.91
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    NiS                  NiS                  mw:      90.77
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    CuS0.5               CuS0.5               mw:      79.58
    Element array:
    [0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.5 0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0. ]
    
    FeOS                 FeOS                 mw:     103.91
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Ni.25S.25O           Ni.25S.25O           mw:      38.69
    Element array:
    [0.   0.   0.   0.   0.   0.   0.   0.   1.   0.   0.   0.   0.   0.
     0.   0.   0.25 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
     0.25 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
     0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
     0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
     0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
     0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
     0.   0.   0.   0.   0.   0.   0.   0.  ]
    
    CuS                  CuS                  mw:      95.61
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    CuO                  CuO                  mw:      79.55
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    


Execute functions for molar derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Component first derivative vectors:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    printResult('dGdn', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_dgdn(t,p,n), 'J/m')
    printResult('d2GdndT', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d2gdndt(t,p,n), 'J/K-m')
    printResult('d2GdndP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d2gdndp(t,p,n), 'J/bar-m')
    printResult('d3GdndT2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdndt2(t,p,n), 'J/K^2-m')
    printResult('d3GdndTdP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdndtdp(t,p,n), 'J/K-bar-m')
    printResult('d3GdndP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdndp2(t,p,n), 'J/bar^2-m')
    printResult('d4GdndT3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdndt3(t,p,n), 'J/K^3-m')
    printResult('d4GdndT2dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdndt2dp(t,p,n), 'J/K^2-bar-m')
    printResult('d4GdndTdP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdndtdp2(t,p,n), 'J/K-bar^2-m')
    printResult('d4GdndP3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdndp3(t,p,n), 'J/bar^3-m')


.. parsed-literal::

                       [  0]         [  1]         [  2]         [  3]         [  4]         
    dGdn       -4.624664e+05 -2.353277e+05 -1.410271e+05 -1.992151e+05 -1.678986e+05 J/m       
    d2GdndT    -3.279900e+02 -8.095877e+01 -1.117905e+02 -1.290345e+02 -1.173106e+02 J/K-m     
    d2GdndP     6.594476e-01  1.633334e+00  8.133086e-01  7.724063e-01  8.390871e-01 J/bar-m   
    d3GdndT2   -1.573033e-02 -2.062067e-02 -2.301200e-02 -1.945551e-02 -1.642221e-02 J/K^2-m   
    d3GdndTdP  -2.588384e-04  6.446500e-04  7.302876e-05  8.384035e-05  8.892059e-05 J/K-bar-m 
    d3GdndP2   -1.218492e-06 -4.841321e-05 -9.289128e-07 -9.577613e-07 -2.676480e-06 J/bar^2-m 
    d4GdndT3    0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^3-m   
    d4GdndT2dP  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^2-bar-
    d4GdndTdP2  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar^2-
    d4GdndP3    0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/bar^3-m 


The Hessian matrix (molar second derivative matrix) is stored as a compact linear array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    printResult('d2Gdn2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d2gdn2(t,p,n), 'J/m^2')
    printResult('d3Gdn2dT', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdn2dt(t,p,n), 'J/K-m^2')
    printResult('d3Gdn2dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdn2dp(t,p,n), 'J/bar-m^2')
    printResult('d4Gdn2dT2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdn2dt2(t,p,n), 'J/K^2-m^2')
    printResult('d4Gdn2dTdP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdn2dtdp(t,p,n), 'J/K-bar-m^2')
    printResult('d4Gdn2dP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdn2dp2(t,p,n), 'J/bar^2-m^2')
    printResult('d5Gdn2dT3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn2dt3(t,p,n), 'J/K^3-m^2')
    printResult('d5Gdn2dT2dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn2dt2dp(t,p,n), 'J/K^2-bar-m^2')
    printResult('d5Gdn2dTdP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn2dtdp2(t,p,n), 'J/K-bar^2-m^2')
    printResult('d5Gdn2dP3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn2dp3(t,p,n), 'J/bar^3-m^2')


.. parsed-literal::

                       [  0]         [  1]         [  2]         [  3]         [  4]         [  5]         [  6]         [  7]         [  8]         [  9]         [ 10]         [ 11]         [ 12]         [ 13]         [ 14]         
    d2Gdn2      3.631551e+05 -9.909290e+04 -4.091487e+03  1.189443e+05 -1.045312e+04  3.475317e+05 -2.083723e+04 -2.987186e+05  5.169939e+03  2.147133e+04 -1.713196e+04 -2.862025e+04  7.397122e+05  2.155454e+01  1.803476e+05 J/m^2     
    d3Gdn2dT   -1.576708e-09  4.640982e-09 -9.783600e-11 -1.049198e-08 -1.402419e-08  1.943042e-11 -1.780479e-13 -8.672956e-11 -1.443867e-10 -4.148268e-11 -4.472491e-12 -1.000217e-11 -4.992790e-08 -6.762534e-08 -9.142881e-08 J/K-m^2   
    d3Gdn2dP    0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/bar-m^2 
    d4Gdn2dT2   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^2-m^2 
    d4Gdn2dTdP  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar-m^
    d4Gdn2dP2   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/bar^2-m^
    d5Gdn2dT3   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^3-m^2 
    d5Gdn2dT2d  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K^2-bar-
    d5Gdn2dTdP  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar^2-
    d5Gdn2dP3   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 J/bar^3-m^


The 3-D Tensor (molar third derivative tensor) is stored as a compact linear array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    printResult('d3Gdn3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d3gdn3(t,p,n), 'J/m^3')
    printResult('d4Gdn3dT', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdn3dt(t,p,n), 'J/K-m^3')
    printResult('d4Gdn3dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d4gdn3dp(t,p,n), 'J/bar-m^3')
    printResult('d5Gdn3dT2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn3dt2(t,p,n), 'J/K^2-m^3')
    printResult('d5Gdn3dTdP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn3dtdp(t,p,n), 'J/K-bar-m^3')
    printResult('d5Gdn3dP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d5gdn3dp2(t,p,n), 'J/bar^2-m^3')
    printResult('d6Gdn3dT3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d6gdn3dt3(t,p,n), 'J/K^3-m^3')
    printResult('d6Gdn3dT2dP', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d6gdn3dt2dp(t,p,n), 'J/K^2-bar-m^3')
    printResult('d6Gdn3dTdP2', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d6gdn3dtdp2(t,p,n), 'J/K-bar^2-m^3')
    printResult('d6Gdn3dP3', sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_d6gdn3dp3(t,p,n), 'J/bar^3-m^3')


.. parsed-literal::

                    [  0]      [  1]      [  2]      [  3]      [  4]      [  5]      [  6]      [  7]      [  8]      [  9]      [ 10]      [ 11]      [ 12]      [ 13]      [ 14]      [ 15]      [ 16]      [ 17]      [ 18]      [ 19]      [ 20]      [ 21]      [ 22]      [ 23]      [ 24]      [ 25]      [ 26]      [ 27]      [ 28]      [ 29]      [ 30]      [ 31]      [ 32]      [ 33]      [ 34]      
    d3Gdn3      0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/m^3         
    d4Gdn3dT    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-m^3       
    d4Gdn3dP    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar-m^3     
    d5Gdn3dT2   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^2-m^3     
    d5Gdn3dTdP  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-bar-m^3   
    d5Gdn3dP2   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar^2-m^3   
    d6Gdn3dT3   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^3-m^3     
    d6Gdn3dT2d  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^2-bar-m^3 
    d6Gdn3dTdP  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-bar^2-m^3 
    d6Gdn3dP3   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar^3-m^3   


Test and time the Gibbs free energy method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    %timeit sulfide_liquid.cy_SulfLiq_sulfide_liquid_calib_g(t, p, n) 


.. parsed-literal::

    45.9 ms ± 1.5 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)


Import the sulfide_liquid module into the ThermoEngine package
--------------------------------------------------------------

Note that the module is treated as if it was generated using the coder
module. See the Example-8 notebook in this folder for further details.

.. code:: ipython3

    modelDB = model.Database(database="CoderModule", calib="calib", 
                             phase_tuple=('thermoengine.sulfide_liquid', {
                                 'Sulf':['SulfLiq','solution']
                             }))


.. parsed-literal::

    Solution phase code generated by the coder module does not yet provide information on solution species. Species are proxied by components.
    Solution phase code generated by the coder module does not yet provide information on species properties. Species are proxied by components.


.. code:: ipython3

    phase = modelDB.get_phase('Sulf')

Test module methods wrapped by the standard ThermoEngine interface.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print (phase.props['phase_name'])
    print (phase.props['formula'])
    print (phase.props['molwt'])
    print (phase.props['abbrev'])
    print (phase.props['endmember_num'])
    print (phase.props['endmember_name'])


.. parsed-literal::

    Sulfide Liquid
    ['O' 'S' 'Fe' 'Ni' 'Cu']
    [15.9994 32.06   55.847  58.71   63.546 ]
    Sulf
    5
    ['O' 'S' 'Fe' 'Ni' 'Cu']


Test solution composition chosen as above.

.. code:: ipython3

    moles_end = np.array([2., 5., 10., 1., 1.])
    for i in range(0,phase.props['endmember_num']):
        print ("mole number of {0:10.10s} = {1:13.6e}".format(phase.props['endmember_name'][i], moles_end[i]))
    if not phase.test_endmember_comp(moles_end):
        print ("Calculated composition is infeasible!")
    print ('Formula: ', phase.compute_formula(t, p, moles_end))
    print ('Total moles of endmembers: ', phase.covert_endmember_comp(moles_end,output='total_moles'))
    mol_elm = phase.covert_endmember_comp(moles_end,output='moles_elements')
    print ('Mole fractions of endmembers: ', phase.covert_endmember_comp(moles_end,output='mole_fraction'))
    print ('Moles of endmembers: ', phase.convert_elements(mol_elm, output='moles_end'))
    print ('Total moles of endmembers: ', phase.convert_elements(mol_elm, output='total_moles'))
    print ('Total grams of phase: ', phase.convert_elements(mol_elm, output='total_grams'))
    # Check if intrinsic mode fails
    if len(moles_end) == 0:
        print ('Intrinsic mode returned an empty array; estimating composition ...')
        moles_end = np.array([0.20813521, 0.00267478, 0.14968884])
        print (moles_end)


.. parsed-literal::

    mole number of O          =  2.000000e+00
    mole number of S          =  5.000000e+00
    mole number of Fe         =  1.000000e+01
    mole number of Ni         =  1.000000e+00
    mole number of Cu         =  1.000000e+00
    Formula:  O2.0S5.0Fe10.0Ni1.0Cu1.0
    Total moles of endmembers:  19.0
    Mole fractions of endmembers:  [0.10526316 0.26315789 0.52631579 0.05263158 0.05263158]
    Moles of endmembers:  [ 2.  5. 10.  1.  1.]
    Total moles of endmembers:  19.0
    Total grams of phase:  873.0248000000001


Test and output the remainder of the methods for retrieval of
thermodynamic quantities.

.. code:: ipython3

    def test_func(name, func, t, p, units, deriv=None, const=None, endmember=None):
        global moles_end
        try:
            if deriv:
                result = func(t, p, deriv=deriv, mol=moles_end)
                if type(result) is np.ndarray:
                    if len(result.shape) == 2:
                        print ("{0:>10s}".format(name), end=' ')
                        for x in result[0]:
                            print ("{0:15.6e}".format(x), end=' ')
                        print (" {0:<20s}".format(units))
                    elif len(result.shape) == 3:
                        for i in range(0,result.shape[1]):
                            print ("{0:>10s}".format(name), end=' ')
                            for x in result[0][i]:
                                print ("{0:15.6e}".format(x), end=' ')
                            print (" {0:<20s}".format(units))
                    elif len(result.shape) == 4:
                        for i in range(0,result.shape[1]):
                            for j in range(0,result.shape[2]):
                                print ("{0:>10s}".format(name), end=' ')
                                for x in result[0][i][j]:
                                    print ("{0:15.6e}".format(x), end=' ')
                                print (" {0:<20s}".format(units))
                    elif len(result.shape) == 1:
                        print ("{0:>10s}".format(name), end=' ')
                        for x in result:
                            print ("{0:15.6e}".format(x), end=' ')
                        print (" {0:<20s}".format(units))
                    else:
                        print ('A', result.shape)
                else:
                    print ("{0:>10s}{1:15.6e} {2:<20s}".format(name, result, units))
            elif const:
                print ("{0:>10s}{1:15.6e} {2:<20s}".format(name, func(t, p, const=const, mol=moles_end), units))
            else:
                result = func(t, p, mol=moles_end)
                if type(result) is np.ndarray:
                    if len(result.shape) == 2:
                        print ("{0:>10s}".format(name), end=' ')
                        for x in result[0]:
                            print ("{0:15.6e}".format(x), end=' ')
                        print (" {0:<20s}".format(units))
                    elif len(result.shape) == 1:
                        print ("{0:>10s}".format(name), end=' ')
                        for x in result:
                            print ("{0:15.6e}".format(x), end=' ')
                        print (" {0:<20s}".format(units))
                    else:
                        print ('B', len(result.shape))
                else:
                    print ("{0:>10s}{1:15.6e} {2:<20s}".format(name, result, units))
        except AttributeError:
            print ("{0:>10s} is not implemented".format(name))
    
    test_func('G',  phase.gibbs_energy, t, p, 'J/mol')
    test_func('dG/dT', phase.gibbs_energy, t, p, 'J/K-mol', deriv={'dT':1})
    test_func('dG/dP', phase.gibbs_energy, t, p, 'J/bar-mol', deriv={'dP':1})
    test_func('dG/dm', phase.gibbs_energy, t, p, 'J/mol^2', deriv={'dmol':1})
    
    test_func('d2G/dT2', phase.gibbs_energy, t, p, 'J/K^2-mol', deriv={'dT':2})
    test_func('d2G/dTdP', phase.gibbs_energy, t, p, 'J/K-bar-mol', deriv={'dT':1, 'dP':1})
    test_func('d2G/dTdm', phase.gibbs_energy, t, p, 'J/K-mol^2', deriv={'dT':1, 'dmol':1})
    test_func('d2G/dP2', phase.gibbs_energy, t, p, 'J/bar^2-mol', deriv={'dP':2})
    test_func('d2G/dPdm', phase.gibbs_energy, t, p, 'J/bar-mol^2', deriv={'dP':1, 'dmol':1})
    test_func('d2G/dm2', phase.gibbs_energy, t, p, 'J/mol^3', deriv={'dmol':2})
    
    test_func('d3G/dT3', phase.gibbs_energy, t, p, 'J/K^3-mol', deriv={'dT':3})
    test_func('d3G/dT2dP', phase.gibbs_energy, t, p, 'J/K^2-bar-mol', deriv={'dT':2, 'dP':1})
    test_func('d3G/dT2dm', phase.gibbs_energy, t, p, 'J/K^2-mol^2', deriv={'dT':2, 'dmol':1})
    test_func('d3G/dTdP2', phase.gibbs_energy, t, p, 'J/K-bar^2-mol', deriv={'dT':1, 'dP':2})
    test_func('d3G/dTdPdm', phase.gibbs_energy, t, p, 'J/K-bar-mol^2', deriv={'dT':1, 'dP':1, 'dmol':1})
    test_func('d3G/dTdm2', phase.gibbs_energy, t, p, 'J/K-mol^3', deriv={'dT':1, 'dmol':2})
    test_func('d3G/dP3', phase.gibbs_energy, t, p, 'J/bar^3-mol', deriv={'dP':3})
    test_func('d3G/dP2dm', phase.gibbs_energy, t, p, 'J/bar^2-mol^2', deriv={'dP':2, 'dmol':1})
    test_func('d3G/dPdm2', phase.gibbs_energy, t, p, 'J/bar-mol^3', deriv={'dP':1, 'dmol':2})
    test_func('d3G/dm3', phase.gibbs_energy, t, p, 'J/mol^4', deriv={'dmol':3})
    
    test_func('H', phase.enthalpy, t, p, 'J/mol')
    test_func('S', phase.entropy, t, p, 'J/K-mol')
    test_func('dS/dm', phase.entropy, t, p, 'J/K-mol^2', deriv={'dmol':1})
    test_func('d2S/dm2', phase.entropy, t, p, 'J/K-mol^3', deriv={'dmol':2})
    
    test_func('Cv', phase.heat_capacity, t, p, 'J/K-mol', const='V')
    test_func('Cp', phase.heat_capacity, t, p, 'J/K-mol')
    test_func('dCp/dT', phase.heat_capacity, t, p, 'J/-K^2-mol', deriv={'dT':1})
    test_func('dCp/dm', phase.heat_capacity, t, p, 'J/-K-mol^2', deriv={'dmol':1})
    
    test_func('rho', phase.density, t, p, 'gm/cc')
    test_func('alpha', phase.thermal_exp, t, p, '1/K')
    test_func('beta', phase.compressibility, t, p, '1/bar')
    test_func('K', phase.bulk_mod, t, p, '')
    test_func('Kp', phase.bulk_mod, t, p, '1/bar', deriv={'dP':1})
    
    test_func("V", phase.volume, t, p, 'J/bar-mol')
    test_func("dV/dT", phase.volume, t, p, 'J/bar-K-mol', deriv={'dT':1})
    test_func("dv/dP", phase.volume, t, p, 'J/bar^2-mol', deriv={'dP':1})
    test_func("dv/dm", phase.volume, t, p, 'J/bar-mol^2', deriv={'dP':1, 'dmol':1})
    test_func("d2V/dT2", phase.volume, t, p, 'J/bar-K^2-mol', deriv={'dT':2})
    test_func("d2V/dTdP", phase.volume, t, p, 'J/bar^2-K-mol', deriv={'dT':1, 'dP':1})
    test_func("d2V/dP2", phase.volume, t, p, 'J/bar^3-mol', deriv={'dP':2})
    test_func("d2V/dTdm", phase.volume, t, p, 'J/bar-K-mol^2', deriv={'dT':1, 'dmol':1})
    test_func("d2V/dPdm", phase.volume, t, p, 'J/bar^2-mol^2', deriv={'dP':1, 'dmol':1})
    test_func("d2V/dm2", phase.volume, t, p, 'J/bar-mol^3', deriv={'dmol':2})
    
    test_func('mu0', phase.chem_potential, t, p, 'J/mol')
    test_func('activity', phase.activity, t, p, '')
    test_func('da/dm', phase.activity, t, p, '1/mol', deriv={'dmol':1})


.. parsed-literal::

             G  -3.878956e+06 J/mol               
         dG/dT  -2.452130e+03 J/K-mol             
         dG/dP   1.884016e+01 J/bar-mol           
         dG/dm   -4.624664e+05   -2.353277e+05   -1.410271e+05   -1.992151e+05   -1.678986e+05  J/mol^2             
       d2G/dT2  -3.800446e-01 J/K^2-mol           
      d2G/dTdP   2.350344e-03 J/K-bar-mol         
      d2G/dTdm   -3.279900e+02   -8.095877e+01   -1.117905e+02   -1.290345e+02   -1.173106e+02  J/K-mol^2           
       d2G/dP2  -2.574264e-04 J/bar^2-mol         
      d2G/dPdm    6.594476e-01    1.633334e+00    8.133086e-01    7.724063e-01    8.390871e-01  J/bar-mol^2         
       d2G/dm2    3.631551e+05   -9.909290e+04   -4.091487e+03    1.189443e+05   -1.045312e+04    3.475317e+05   -2.083723e+04   -2.987186e+05    5.169939e+03    2.147133e+04   -1.713196e+04   -2.862025e+04    7.397122e+05    2.155454e+01    1.803476e+05  J/mol^3             
       d3G/dT3   1.880394e-04 J/K^3-mol           
     d3G/dT2dP   1.310298e-05 J/K^2-bar-mol       
     d3G/dT2dm   -1.573033e-02   -2.062067e-02   -2.301200e-02   -1.945551e-02   -1.642221e-02  J/K^2-mol^2         
     d3G/dTdP2  -8.854753e-07 J/K-bar^2-mol       
    d3G/dTdPdm   -2.588384e-04    6.446500e-04    7.302876e-05    8.384035e-05    8.892059e-05  J/K-bar-mol^2       
     d3G/dTdm2   -1.576708e-09    4.640982e-09   -9.783600e-11   -1.049198e-08   -1.402419e-08    1.943042e-11   -1.780479e-13   -8.672956e-11   -1.443867e-10   -4.148268e-11   -4.472491e-12   -1.000217e-11   -4.992790e-08   -6.762534e-08   -9.142881e-08  J/K-mol^3           
       d3G/dP3   6.631852e-08 J/bar^3-mol         
     d3G/dP2dm   -1.218492e-06   -4.841321e-05   -9.289128e-07   -9.577613e-07   -2.676480e-06  J/bar^2-mol^2       
     d3G/dPdm2    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  J/bar-mol^3         
       d3G/dm3    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  J/mol^4             
             H   1.025303e+06 J/mol               
             S   2.452130e+03 J/K-mol             
         dS/dm    3.279900e+02    8.095877e+01    1.117905e+02    1.290345e+02    1.173106e+02  J/K-mol^2           
       d2S/dm2    1.576708e-09   -4.640982e-09    9.783600e-11    1.049198e-08    1.402419e-08   -1.942983e-11    1.780292e-13    8.672955e-11    1.443867e-10    4.148261e-11    4.472356e-12    1.000216e-11    4.992790e-08    6.762534e-08    9.142881e-08  J/K-mol^3           
            Cv   7.171711e+02 J/K-mol             
            Cp   7.600892e+02 J/K-mol             
        dCp/dT   3.965713e-03 J/-K^2-mol          
        dCp/dm    7.865164e-06    1.031033e-05    1.150600e-05    9.727757e-06    8.211105e-06  J/-K-mol^2          
           rho is not implemented
         alpha   1.247519e-04 1/K                 
          beta   1.366371e-05 1/bar               
             K   7.318657e+04                     
            Kp   1.785441e+01 1/bar               
             V   1.884016e+01 J/bar-mol           
         dV/dT   2.350344e-03 J/bar-K-mol         
         dv/dP  -2.574264e-04 J/bar^2-mol         
         dv/dm   -1.218492e-06   -4.841321e-05   -9.289128e-07   -9.577613e-07   -2.676480e-06  J/bar-mol^2         
       d2V/dT2   1.310298e-05 J/bar-K^2-mol       
      d2V/dTdP  -8.854753e-07 J/bar^2-K-mol       
       d2V/dP2   6.631852e-08 J/bar^3-mol         
      d2V/dTdm   -2.588384e-04    6.446500e-04    7.302876e-05    8.384035e-05    8.892059e-05  J/bar-K-mol^2       
      d2V/dPdm   -1.218492e-06   -4.841321e-05   -9.289128e-07   -9.577613e-07   -2.676480e-06  J/bar^2-mol^2       
       d2V/dm2    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  J/bar-mol^3         
           mu0   -1.678986e+05  J/mol               
      activity is not implemented
         da/dm is not implemented


