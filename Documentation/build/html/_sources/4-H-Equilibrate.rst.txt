Enthalpy potential minimization (S, P, constrained)
===================================================

Closed system; crystallization of a rhyolitic liquid using
rhyolite-MELTS

.. code:: ipython3

    import numpy as np
    import scipy.optimize as opt
    import scipy.linalg as lin 
    import sys

.. code:: ipython3

    from thermoengine import core, phases, model, equilibrate

.. code:: ipython3

    np.set_printoptions(linewidth=200, precision=1)

Create phases for equilibrium assemblages
-----------------------------------------

.. code:: ipython3

    modelDB = model.Database(liq_mod='v1.0')

.. code:: ipython3

    Liquid = modelDB.get_phase('Liq')
    Feldspar = modelDB.get_phase('Fsp')
    Quartz = modelDB.get_phase('Qz')

The Berman model database provides the SWIM water model by default.
Instead, override that choice by instantiating the MELTS 1.0.2 water
model directly.

.. code:: ipython3

    Water = phases.PurePhase('WaterMelts', 'H2O', calib=False)

Define elements in system and phases in system
----------------------------------------------

.. code:: ipython3

    elm_sys = ['H','O','Na','Mg','Al','Si','P','K','Ca','Ti','Cr','Mn','Fe','Co','Ni']
    phs_sys = [Liquid, Feldspar, Water, Quartz]

Composition of the system
-------------------------

This is a high-silica rhyolite

.. code:: ipython3

    grm_oxides = {
        'SiO2':  77.5, 
        'TiO2':   0.08, 
        'Al2O3': 12.5, 
        'Fe2O3':  0.207,
        'Cr2O3':  0.0, 
        'FeO':    0.473, 
        'MnO':    0.0,
        'MgO':    0.03, 
        'NiO':    0.0, 
        'CoO':    0.0,
        'CaO':    0.43, 
        'Na2O':   3.98, 
        'K2O':    4.88, 
        'P2O5':   0.0, 
        'H2O':    5.5
    }

Cast this composition as moles of elements for input to the Equilibrate
class

.. code:: ipython3

    mol_oxides = core.chem.format_mol_oxide_comp(grm_oxides, convert_grams_to_moles=True)
    moles_end,oxide_res = Liquid.calc_endmember_comp(
        mol_oxide_comp=mol_oxides, method='intrinsic', output_residual=True)
    if not Liquid.test_endmember_comp(moles_end):
        print ("Calculated composition is infeasible!")
    mol_elm = Liquid.covert_endmember_comp(moles_end,output='moles_elements')

.. code:: ipython3

    blk_cmp = []
    for elm in elm_sys:
        index = core.chem.PERIODIC_ORDER.tolist().index(elm)
        blk_cmp.append(mol_elm[index])
    blk_cmp = np.array(blk_cmp)

Function to constrain the entropy
---------------------------------

Note that the entropy is equivalent to $ -
:raw-latex:`\frac{{\partial G}}{{\partial T}}`$ - Run an equilibration
step at fixed T,P - Calculate the entropy of the system - Define a
function to set the entropy for subsequent equilibration steps

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys)

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Water
    Quad (000) norm:  2.8609503260090e-02 Lin (019) step:  9.4395802431652e-01 func: -1.7280794403704e+06
    Quad (001) norm:  1.7020967377192e-08 Lin (026) step: -3.9941348412705e-01 func: -1.7280794403704e+06
    Quad (002) norm:  2.3819371335756e-08 Lin (037) step:  9.4323348836889e-01 func: -1.7280794403704e+06
    Quad (003) norm:  1.3521403531091e-09 Lin (032) step: -4.7213659552033e-01 func: -1.7280794403704e+06
    Quad (004) norm:  1.9905400629586e-09 Lin (039) step: -1.2872026626417e+00 func: -1.7280794403704e+06
    Quad (005) norm:  4.5527674310389e-09 Lin (039) step: -1.7924442664609e+00 func: -1.7280794403704e+06
    Minimal energy termination of quadratic loop.
    
    Add: Feldspar
    Quad (000) norm:  5.8343755255746e-03 Lin (020) step:  9.7809006932692e-01 func: -1.7280795838328e+06
    Quad (001) norm:  2.2098427302856e-04 Lin (012) step:  1.0471013041669e+00 func: -1.7280795887688e+06
    Quad (002) norm:  9.8888798416875e-05 Lin (025) step:  1.0023917655135e+00 func: -1.7280795888147e+06
    Quad (003) norm:  1.4361892154363e-07 Lin (029) step:  9.9860621867021e-01 func: -1.7280795888147e+06
    Quad (004) norm:  1.6471976415762e-10 Lin (037) step: -1.0717369638414e+00 func: -1.7280795888147e+06
    Quad (005) norm:  3.4122477608691e-10 Lin (034) step: -1.0565415854491e+00 func: -1.7280795888147e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      175.0 MPa
    Liquid          moles:   1.639690 grams: 104.628
               SiO2 form:  SiO2           X:  0.6745  wt%    SiO2   73.73
               TiO2 form:  TiO2           X:  0.0006  wt%    TiO2    0.08
              Al2O3 form:  Al2O3          X:  0.0424  wt%   Al2O3   11.82
              Fe2O3 form:  Fe2O3          X:  0.0008  wt%   Fe2O3    0.20
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.45
            Fe2SiO4 form:  Fe2SiO4        X:  0.0020  wt%     MgO    0.03
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.39
            Mg2SiO4 form:  Mg2SiO4        X:  0.0002  wt%    Na2O    3.76
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.66
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.89
             CaSiO3 form:  CaSiO3         X:  0.0044  
            Na2SiO3 form:  Na2SiO3        X:  0.0387  
            KAlSiO4 form:  KAlSiO4        X:  0.0631  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1732  
    Feldspar        moles:   0.002127 grams:   0.567
             albite form:  NaAlSi3O8      X:  0.7392  wt%    SiO2   63.41
          anorthite form:  CaAl2Si2O8     X:  0.1886  wt%   Al2O3   22.75
           sanidine form:  KAlSi3O8       X:  0.0723  wt%     CaO    3.97
                                                      wt%    Na2O    8.60
                                                      wt%     K2O    1.28
    Water           moles:   0.021370 grams:   0.385
    Quartz          affn:     134.38


.. code:: ipython3

    delta_dGdT = 0.0
    dGdT = state.dGdT(t,p)
    def con(t, p, state):
        return dGdT + delta_dGdT

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[('T',con)])

.. code:: ipython3

    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Water
    Quad (000) norm:  8.4109760121255e-02 Lin (020) step: -7.6247405768250e-01 func: -1.4419214302114e+06
    Quad (001) norm:  5.5409228636560e-07 Lin (040) step:  3.8706352769574e-01 func: -1.4419214302114e+06
    Quad (002) norm:  3.3961756355139e-07 Lin (033) step: -7.4383200643679e-01 func: -1.4419214302114e+06
    Quad (003) norm:  5.9224454039414e-07 Lin (039) step:  5.0200500853460e-01 func: -1.4419214302114e+06
    Quad (004) norm:  2.9492666607495e-07 Lin (045) step: -2.9039028052303e-02 func: -1.4419214302114e+06
    Quad (005) norm:  3.0349128388339e-07 Lin (039) step: -5.0155281000757e-01 func: -1.4419214302114e+06
    Minimal energy termination of quadratic loop.
    
    Add: Feldspar
    Quad (000) norm:  8.2433571075966e-01 Lin (024) step:  9.8350376266276e-01 func: -1.4419216218052e+06
    Quad (001) norm:  1.8528134113517e-02 Lin (010) step:  1.0311516961913e+00 func: -1.4419216272040e+06
    Quad (002) norm:  1.1294075179570e-02 Lin (026) step:  1.0033290064394e+00 func: -1.4419216272453e+06
    Quad (003) norm:  2.5808106760253e-05 Lin (038) step:  9.7213570108658e-01 func: -1.4419216272453e+06
    Quad (004) norm:  7.1851279402016e-07 Lin (039) step:  1.0557287173513e+00 func: -1.4419216272453e+06
    Quad (005) norm:  4.0253001638508e-08 Lin (041) step:  5.5728093004448e-02 func: -1.4419216272453e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      175.0 MPa
    Liquid          moles:   1.639690 grams: 104.628
               SiO2 form:  SiO2           X:  0.6745  wt%    SiO2   73.73
               TiO2 form:  TiO2           X:  0.0006  wt%    TiO2    0.08
              Al2O3 form:  Al2O3          X:  0.0424  wt%   Al2O3   11.82
              Fe2O3 form:  Fe2O3          X:  0.0008  wt%   Fe2O3    0.20
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.45
            Fe2SiO4 form:  Fe2SiO4        X:  0.0020  wt%     MgO    0.03
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.39
            Mg2SiO4 form:  Mg2SiO4        X:  0.0002  wt%    Na2O    3.76
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.66
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.89
             CaSiO3 form:  CaSiO3         X:  0.0044  
            Na2SiO3 form:  Na2SiO3        X:  0.0387  
            KAlSiO4 form:  KAlSiO4        X:  0.0631  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1732  
    Feldspar        moles:   0.002127 grams:   0.567
             albite form:  NaAlSi3O8      X:  0.7392  wt%    SiO2   63.41
          anorthite form:  CaAl2Si2O8     X:  0.1886  wt%   Al2O3   22.75
           sanidine form:  KAlSi3O8       X:  0.0723  wt%     CaO    3.97
                                                      wt%    Na2O    8.60
                                                      wt%     K2O    1.28
    Water           moles:   0.021370 grams:   0.385
    Quartz          affn:     134.38


Pickup runs use previously computed state

.. code:: ipython3

    delta_dGdT = 5.0
    state = equil.execute(t, p, state=state, debug=0)
    state.print_state()


.. parsed-literal::

    Minimal energy termination of quadratic loop.
    Minimal energy termination of quadratic loop.
     
    T =     762.82 °C, P =      175.0 MPa
    Liquid          moles:   1.325625 grams:  84.947
               SiO2 form:  SiO2           X:  0.6770  wt%    SiO2   73.80
               TiO2 form:  TiO2           X:  0.0008  wt%    TiO2    0.09
              Al2O3 form:  Al2O3          X:  0.0399  wt%   Al2O3   11.65
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.24
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.56
            Fe2SiO4 form:  Fe2SiO4        X:  0.0025  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.27
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    3.64
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.90
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.82
             CaSiO3 form:  CaSiO3         X:  0.0031  
            Na2SiO3 form:  Na2SiO3        X:  0.0376  
            KAlSiO4 form:  KAlSiO4        X:  0.0666  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1713  
    Feldspar        moles:   0.025831 grams:   6.875
             albite form:  NaAlSi3O8      X:  0.7556  wt%    SiO2   64.95
          anorthite form:  CaAl2Si2O8     X:  0.1231  wt%   Al2O3   21.51
           sanidine form:  KAlSi3O8       X:  0.1213  wt%     CaO    2.59
                                                      wt%    Na2O    8.80
                                                      wt%     K2O    2.15
    Water           moles:   0.078248 grams:   1.410
    Quartz          moles:   0.107219 grams:   6.442
    Feldspar        moles:   0.021754 grams:   5.906
             albite form:  NaAlSi3O8      X:  0.4246  wt%    SiO2   66.03
          anorthite form:  CaAl2Si2O8     X:  0.0165  wt%   Al2O3   19.09
           sanidine form:  KAlSi3O8       X:  0.5589  wt%     CaO    0.34
                                                      wt%    Na2O    4.85
                                                      wt%     K2O    9.70


