SWIM
====

Standard Water Intergrated Model
================================

The SWIM model is described in a notebook located in the Pure-Phases
folder with the title “Water.ipynb.” In that folder the model is
accessed using server code written in Ovjective-C and bridged to Python
utilizing the Rubicon Python to Objective-C wrappper.

SWIM may also be accessed using Cython wrappers to C/C++ code. That
method of access is illustrated in this notebook. The Cython access
method is faster because it avoids the Rubicon bridge and is more easily
ported to hardware platforms not normally configured to compile and
execute Objective-C. The method of access to SWIM demonstrated here
should be used in preference to the default method, which will soon be
deprecated.

.. code:: ipython3

    from thermoengine import model

The Cython-SWIM module is accessed as …

.. code:: ipython3

    modelDB = model.Database(database="CoderModule", calib=True, phase_tuple=('thermoengine.aqueous', {'SWIM':['SWIM','pure']}))

… the objective-C/Rubicon code is accessed using the default call,
modelDB = model.Database()

The phase is accessed in the standard way once the model database has
been instantiated.

.. code:: ipython3

    SWIM = modelDB.get_phase('SWIM')

Generic Properties of the SWIM phase
------------------------------------

.. code:: ipython3

    print (SWIM.props['phase_name'])
    print (SWIM.props['formula'][0])
    print (SWIM.props['molwt'][0])


.. parsed-literal::

    SWIM
    H2O
    18.01528


Thermodynamic properties of the SWIM phase
------------------------------------------

By default, interpolative smooothing of the integrated model properties
is implemented.

.. code:: ipython3

    t = 1000.0 # K
    p = 1000.0 # bars  

.. code:: ipython3

    import numpy as np
    def test_func(name, func, t, p, units, deriv=None, const=None, endmember=None):
        try:
            if deriv:
                result = func(t, p, deriv=deriv)
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
                print ("{0:>10s}{1:15.6e} {2:<20s}".format(name, func(t, p, const=const), units))
            else:
                result = func(t, p)
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
    
    test_func('G',  SWIM.gibbs_energy, t, p, 'J/mol')
    test_func('dG/dT', SWIM.gibbs_energy, t, p, 'J/K-mol', deriv={'dT':1})
    test_func('dG/dP', SWIM.gibbs_energy, t, p, 'J/bar-mol', deriv={'dP':1})
    
    test_func('d2G/dT2', SWIM.gibbs_energy, t, p, 'J/K^2-mol', deriv={'dT':2})
    test_func('d2G/dTdP', SWIM.gibbs_energy, t, p, 'J/K-bar-mol', deriv={'dT':1, 'dP':1})
    test_func('d2G/dP2', SWIM.gibbs_energy, t, p, 'J/bar^2-mol', deriv={'dP':2})
    
    test_func('d3G/dT3', SWIM.gibbs_energy, t, p, 'J/K^3-mol', deriv={'dT':3})
    test_func('d3G/dT2dP', SWIM.gibbs_energy, t, p, 'J/K^2-bar-mol', deriv={'dT':2, 'dP':1})
    test_func('d3G/dTdP2', SWIM.gibbs_energy, t, p, 'J/K-bar^2-mol', deriv={'dT':1, 'dP':2})
    test_func('d3G/dP3', SWIM.gibbs_energy, t, p, 'J/bar^3-mol', deriv={'dP':3})
    
    test_func('H', SWIM.enthalpy, t, p, 'J/mol')
    test_func('S', SWIM.entropy, t, p, 'J/K-mol')
    
    test_func('Cv', SWIM.heat_capacity, t, p, 'J/K-mol', const='V')
    test_func('Cp', SWIM.heat_capacity, t, p, 'J/K-mol')
    test_func('dCp/dT', SWIM.heat_capacity, t, p, 'J/-K^2-mol', deriv={'dT':1})
    
    test_func('rho', SWIM.density, t, p, 'gm/cc')
    test_func('alpha', SWIM.thermal_exp, t, p, '1/K')
    test_func('beta', SWIM.compressibility, t, p, '1/bar')
    test_func('K', SWIM.bulk_mod, t, p, '')
    test_func('Kp', SWIM.bulk_mod, t, p, '1/bar', deriv={'dP':1})
    
    test_func("V", SWIM.volume, t, p, 'J/bar-mol')
    test_func("dV/dT", SWIM.volume, t, p, 'J/bar-K-mol', deriv={'dT':1})
    test_func("dv/dP", SWIM.volume, t, p, 'J/bar^2-mol', deriv={'dP':1})
    test_func("d2V/dT2", SWIM.volume, t, p, 'J/bar-K^2-mol', deriv={'dT':2})
    test_func("d2V/dTdP", SWIM.volume, t, p, 'J/bar^2-K-mol', deriv={'dT':1, 'dP':1})
    test_func("d2V/dP2", SWIM.volume, t, p, 'J/bar^3-mol', deriv={'dP':2})
    
    test_func('mu0', SWIM.chem_potential, t, p, 'J/mol')


.. parsed-literal::

             G  -3.234830e+05 J/mol               
         dG/dT  -1.671547e+02 J/K-mol             
         dG/dP   6.838333e+00 J/bar-mol           
       d2G/dT2  -7.291385e-02 J/K^2-mol           
      d2G/dTdP   1.431826e-02 J/K-bar-mol         
       d2G/dP2  -7.243395e-03 J/bar^2-mol         
       d3G/dT3   1.971884e-04 J/K^3-mol           
     d3G/dT2dP  -1.818046e-05 J/K^2-bar-mol       
     d3G/dTdP2  -1.299463e-05 J/K-bar^2-mol       
       d3G/dP3   1.744388e-05 J/bar^3-mol         
             H is not implemented
             S   1.671547e+02 J/K-mol             
            Cv   4.461045e+01 J/K-mol             
            Cp   7.291385e+01 J/K-mol             
        dCp/dT  -1.242746e-01 J/-K^2-mol          
           rho is not implemented
         alpha   2.093823e-03 1/K                 
          beta   1.059234e-03 1/bar               
             K   9.440785e+02                     
            Kp  -1.428920e-06 1/bar               
             V   6.838333e+00 J/bar-mol           
         dV/dT   1.431826e-02 J/bar-K-mol         
         dv/dP  -7.243395e-03 J/bar^2-mol         
       d2V/dT2  -1.818046e-05 J/bar-K^2-mol       
      d2V/dTdP  -1.299463e-05 J/bar^2-K-mol       
       d2V/dP2   1.744388e-05 J/bar^3-mol         
           mu0  -3.234830e+05 J/mol               


Smoothing can be turned off using calibration mode
--------------------------------------------------

There is only one calibration parameter for SWIM. It is an integer,
which if non-zero turns off smoothing and forces the use of one of the
four integrated thermodynamic models, even if that model is not
applicable at the T,P conditions specified. See further explanation in
the Pure-Phase notebook mentioned above.

.. code:: ipython3

    try:
        param_props = SWIM.param_props
        supports_calib = param_props['supports_calib']
        print ('This phase supports the Calibration protocol')
        nparam = param_props['param_num']
        print ('... there are', nparam, 'parameters')
        names = param_props['param_names']
        units = param_props['param_units']
        values = param_props['param0']
        t = 1000.0
        p = 1000.0
        for i in range (0, nparam):
            print ("Parameter {0:<15s} has value {1:15.6e}  {2:<20s}".format(names[i], values[i], units[i]))
    except AttributeError:
        print ('This phase does not implement the parameter calibration protocol')


.. parsed-literal::

    This phase supports the Calibration protocol
    ... there are 1 parameters
    Parameter EOS/0-auto/1-DZ2006/2-ZD2005/3-Holten/4-Wagner has value    0.000000e+00  None                


Use the model of Duan and Zhang (2006)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the model adopted in MELTS

.. code:: ipython3

    SWIM.module.cy_SWIM_aqueous_set_param_value(0,1)
    SWIM.gibbs_energy(t, p)




.. parsed-literal::

    -323482.90090449614



Use the model of Zhang and Duan (2005)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the model adopted in DEW

.. code:: ipython3

    SWIM.module.cy_SWIM_aqueous_set_param_value(0,2)
    SWIM.gibbs_energy(1000,1000)




.. parsed-literal::

    -323483.03766402684



Use the model of Holten (2014)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a model applicable to supercoooled water

.. code:: ipython3

    SWIM.module.cy_SWIM_aqueous_set_param_value(0,3)
    SWIM.gibbs_energy(1000,1000)




.. parsed-literal::

    -293932.509702217



Use the model of Wagner (2002)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the standard water model for properties in the liquid and steam
region up to the critical point

.. code:: ipython3

    SWIM.module.cy_SWIM_aqueous_set_param_value(0,4)
    SWIM.gibbs_energy(1000,1000)




.. parsed-literal::

    -323607.85259968095



