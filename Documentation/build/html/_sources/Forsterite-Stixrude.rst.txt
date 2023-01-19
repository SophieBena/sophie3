Forsterite - Stixrude
=====================

.. code:: ipython3

    from thermoengine import phases
    from thermoengine import model

Create a Python reference to the Forsterite phase in teh Stixrude database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    modelDBStix = model.Database(database='Stixrude')
    obj = modelDBStix.get_phase('Fo')

All phases that conform to the Stoichiometric Phase Protocol …
--------------------------------------------------------------

…implement the following functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   (double)getGibbsFreeEnergyFromT:(double)t andP:(double)p;
   (double)getEnthalpyFromT:(double)t andP:(double)p;
   (double)getEntropyFromT:(double)t andP:(double)p;
   (double)getHeatCapacityFromT:(double)t andP:(double)p;
   (double)getDcpDtFromT:(double)t andP:(double)p;
   (double)getVolumeFromT:(double)t andP:(double)p;
   (double)getDvDtFromT:(double)t andP:(double)p;
   (double)getDvDpFromT:(double)t andP:(double)p;
   (double)getD2vDt2FromT:(double)t andP:(double)p;
   (double)getD2vDtDpFromT:(double)t andP:(double)p;
   (double)getD2vDp2FromT:(double)t andP:(double)p;

where *t* (temperature) is in K, and *p* (pressure) is in bars. ### In
Python, these calls are written:

.. code:: ipython3

    def formatted_output(title, value, units, decimals, sci_notation=False):
        format_str = "{0:>10s}{1:15."+str(decimals)
        if sci_notation:
            format_str += "e"
        else:
            format_str += "f"
            
        format_str += "} {2:<20s}"
        
        output = format_str.format(title, value, units)
        return output

.. code:: ipython3

    t = 2100.0 - 273.15
    p = 15000.0
    
    print (formatted_output('G', obj.gibbs_energy(t,p), 'J/mol',2))
    print (formatted_output('H', obj.enthalpy(t,p), 'J/mol',2))
    print (formatted_output('S', obj.entropy(t,p), 'J/K-mol',2))
    print (formatted_output('Cp', obj.heat_capacity(t,p), 'J/K-mol',3))
    print (formatted_output("V", obj.volume(t,p), 'J/bar-mol',3))
    print (formatted_output("dV/dT", obj.volume(t,p, deriv={'dT':1}), 'J/bar-K-mol',6, sci_notation=True))
    print (formatted_output("dv/dP", obj.volume(t,p, deriv={'dP':1}), 'J/bar^2-mol',6, sci_notation=True))
    print (formatted_output("d2V/dT2", obj.volume(t,p, deriv={'dT':2}), 'J/bar-K^2-mol',6, sci_notation=True))
    print (formatted_output("d2V/dTdP", obj.volume(t,p, deriv={'dT':1, 'dP':1}), 'J/bar^2-K-mol',6, sci_notation=True))
    print (formatted_output("d2V/dP2", obj.volume(t,p, deriv={'dP':2}), 'J/bar^3-mol',6, sci_notation=True))


.. parsed-literal::

             G    -2401656.91 J/mol               
             H    -1702496.72 J/mol               
             S         382.71 J/K-mol             
            Cp        187.032 J/K-mol             
             V          4.531 J/bar-mol           
         dV/dT   1.859698e-04 J/bar-K-mol         
         dv/dP  -4.513093e-06 J/bar^2-mol         
       d2V/dT2   6.744500e-08 J/bar-K^2-mol       
      d2V/dTdP  -1.359414e-09 J/bar^2-K-mol       
       d2V/dP2   2.734201e-11 J/bar^3-mol         


.. code:: ipython3

    from scipy.misc import derivative
    t = 2100.0-273.15
    p = 15000.0

.. code:: ipython3

    def test_h(t,p):
        h_est = obj.gibbs_energy(t,p) + t*obj.entropy(t, p)
        h_act = obj.enthalpy(t, p)
        h_err = (h_est-h_act)*100.0/h_act
        print ("H       {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(h_err, h_est, h_act))
    def g(x, doT=True):
        if doT:
            return obj.gibbs_energy(x, p)
        else:
            return obj.gibbs_energy(t, x)
    def test_g_dt(t,p):
        s_est = -derivative(g, t, args=(True,))
        s_act = obj.entropy(t, p)
        s_err = (s_est-s_act)*100.0/s_act
        print ("S       {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(s_err, s_est, s_act))
    def test_g_dp(t,p):
        v_est = derivative(g, p, args=(False,))
        v_act = obj.volume(t, p)
        v_err = (v_est-v_act)*100.0/v_act
        print ("V       {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(v_err, v_est, v_act))
    def s(x, doT=True):
        if doT:
            return obj.entropy(x, p)
        else:
            return obj.entropy(t, x)
    def test_s_dt(t,p):
        cp_est = t*derivative(s, t, args=(True,))
        cp_act = obj.heat_capacity(t, p)
        cp_err = (cp_est-cp_act)*100.0/cp_act
        print ("Cp      {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(cp_err, cp_est, cp_act))
    def cp(x, doT=True):
        if doT:
            return obj.heat_capacity(x, p)
        else:
            return obj.heat_capacity(t, x)
    def test_cp_dt(t,p):
        dcpdt_est = derivative(cp, t, args=(True,))
        dcpdt_act = obj.heat_capacity(t,p, deriv={'dT':1})
        dcpdt_err = (dcpdt_est-dcpdt_act)*100.0/dcpdt_act
        print ("dCpDt   {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(dcpdt_err, dcpdt_est, dcpdt_act))
    def v(x, doT=True):
        if doT:
            return obj.volume(x, p)
        else:
            return obj.volume(t, x)
    def test_v_dt(t,p):
        dvdt_est = derivative(v, t, args=(True,))
        dvdt_act = obj.volume(t,p, deriv={'dT':1})
        dvdt_err = (dvdt_est-dvdt_act)*100.0/dvdt_act
        print ("dVdT    {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(dvdt_err, dvdt_est, dvdt_act))
    def test_v_dp(t,p):
        dvdp_est = derivative(v, p, args=(False,))
        dvdp_act = obj.volume(t,p, deriv={'dP':1})
        dvdp_err = (dvdp_est-dvdp_act)*100.0/dvdp_act
        print ("dVdP    {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(dvdp_err, dvdp_est, dvdp_act))
    def dvdt(x, doT=True):
        if doT:
            return obj.volume(x, p, deriv={'dT':1})
        else:
            return obj.volume(t, x, deriv={'dT':1})
    def dvdp(x, doT=True):
        if doT:
            return obj.volume(x, p, deriv={'dP':1})
        else:
            return obj.volume(t, x, deriv={'dP':1})
    def test_dvdt_dt(t,p):
        d2vdt2_est = derivative(dvdt, t, args=(True,))
        d2vdt2_act = obj.volume(t,p, deriv={'dT':2})
        d2vdt2_err = (d2vdt2_est-d2vdt2_act)*100.0/d2vdt2_act
        print ("d2VdT2  {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(d2vdt2_err, d2vdt2_est, d2vdt2_act))
    def test_dvdt_dp(t,p):
        d2vdtdp_est = derivative(dvdt, p, args=(False,))
        d2vdtdp_act = obj.volume(t,p, deriv={'dT':1, 'dP':1})
        d2vdtdp_err = (d2vdtdp_est-d2vdtdp_act)*100.0/d2vdtdp_act
        print ("d2VdTdP {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(d2vdtdp_err, d2vdtdp_est, d2vdtdp_act))
    def test_dvdp_dt(t,p):
        d2vdtdp_est = derivative(dvdp, t, args=(True,))
        d2vdtdp_act = obj.volume(t,p, deriv={'dT':1, 'dP':1})
        d2vdtdp_err = (d2vdtdp_est-d2vdtdp_act)*100.0/d2vdtdp_act
        print ("d2VdTDp {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(d2vdtdp_err, d2vdtdp_est, d2vdtdp_act))
    def test_dvdp_dp(t,p):
        d2vdp2_est = derivative(dvdp, p, args=(False,))
        d2vdp2_act = obj.volume(t,p, deriv={'dP':2})
        d2vdp2_err = (d2vdp2_est-d2vdp2_act)*100.0/d2vdp2_act
        print ("d2VdP2  {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(d2vdp2_err, d2vdp2_est, d2vdp2_act))

.. code:: ipython3

    test_h(t,p)
    test_g_dt(t,p)
    test_s_dt(t,p)
    test_cp_dt(t,p)
    test_g_dp(t,p)
    test_v_dt(t,p)
    test_v_dp(t,p)
    test_dvdt_dt(t,p)
    test_dvdt_dp(t,p)
    test_dvdp_dt(t,p)
    test_dvdp_dp(t,p)


.. parsed-literal::

    H        -0.000000 % error, est:   -1.702497e+06 act:   -1.702497e+06
    S        -0.000002 % error, est:    3.827135e+02 act:    3.827135e+02
    Cp        0.000009 % error, est:    1.870316e+02 act:    1.870316e+02
    dCpDt    -0.052643 % error, est:    1.544871e-02 act:    1.545684e-02
    V         0.000000 % error, est:    4.530851e+00 act:    4.530851e+00
    dVdT      0.000005 % error, est:    1.859698e-04 act:    1.859698e-04
    dVdP     -0.000000 % error, est:   -4.513093e-06 act:   -4.513093e-06
    d2VdT2   -0.077181 % error, est:    6.744501e-08 act:    6.749710e-08
    d2VdTdP   0.001460 % error, est:   -1.359414e-09 act:   -1.359394e-09
    d2VdTDp  -0.077716 % error, est:   -1.359414e-09 act:   -1.360472e-09
    d2VdP2    0.001319 % error, est:    2.734201e-11 act:    2.734165e-11


.. code:: ipython3

    obj.volume(t,p)




.. parsed-literal::

    4.530851427152421



.. code:: ipython3

    obj.props




.. parsed-literal::

    OrderedDict([('abbrev', 'Fo'),
                 ('phase_name', 'Forsterite'),
                 ('class_name', 'ForsteriteStixrude'),
                 ('identifier', 'Objective-C-base'),
                 ('endmember_name', array(['Forsterite'], dtype='<U10')),
                 ('endmember_ids', [None]),
                 ('formula', array(['Mg2SiO4'], dtype='<U7')),
                 ('atom_num', array([7.])),
                 ('molwt', array([140.6931])),
                 ('elemental_entropy', array([494.47])),
                 ('element_comp',
                  array([[0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 2., 0., 1., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])),
                 ('mol_oxide_comp',
                  array([[1., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.]])),
                 ('endmember_id', array([0]))])



