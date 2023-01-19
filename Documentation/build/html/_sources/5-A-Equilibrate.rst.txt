Helmholtz potential minimization (T, V, constrained)
====================================================

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

Function to constrain the volume
--------------------------------

Note that the volume is equivalent to
:math:`\frac{{\partial G}}{{\partial P}}` - Run an equilibration step at
fixed T,P - Calculate the volume of the system - Define a function to
set the volume for subsequent equilibration steps

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

    delta_dGdP = 0.0
    dGdP = state.dGdP(t,p)
    def con(t, p, state):
        return dGdP + delta_dGdP

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[('P',con)])

.. code:: ipython3

    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Water
    Quad (000) norm:  4.0956587294351e+03 Lin (021) step: -2.2696118044255e-01 func: -1.7366366101342e+06
    Quad (001) norm:  4.0702767649697e-04 Lin (034) step: -1.2047661979631e+00 func: -1.7366366101342e+06
    Quad (002) norm:  8.9719750906774e-04 Lin (039) step:  4.5323389607110e-01 func: -1.7366366101342e+06
    Quad (003) norm:  4.9052814275327e-04 Lin (039) step:  9.2151943221000e-01 func: -1.7366366101342e+06
    Quad (004) norm:  3.8643596592639e-05 Lin (051) step:  5.6000724903964e-04 func: -1.7366366101342e+06
    Quad (005) norm:  3.8698640855189e-05 Lin (036) step: -7.8016697776728e-01 func: -1.7366366101342e+06
    Minimal energy termination of quadratic loop.
    
    Add: Feldspar
    Quad (000) norm:  1.2314311837883e+01 Lin (023) step:  9.8213803409058e-01 func: -1.7366367817993e+06
    Quad (001) norm:  5.2977658152546e-01 Lin (028) step:  1.0393811364511e+00 func: -1.7366367867438e+06
    Quad (002) norm:  1.6672951218303e-01 Lin (029) step:  1.0016795265340e+00 func: -1.7366367867830e+06
    Quad (003) norm:  7.8694694907159e-05 Lin (038) step:  1.0204278639913e+00 func: -1.7366367867830e+06
    Quad (004) norm:  1.8885367145348e-07 Lin (031) step:  1.0820505081767e+00 func: -1.7366367867830e+06
    Quad (005) norm:  5.3409589848058e-06 Lin (040) step: -1.6041234936379e+00 func: -1.7366367867830e+06
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

    delta_dGdP = 0.1
    state = equil.execute(t, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  7.5412284472754e+02 Lin (023) step:  7.0294798502808e-01 func: -1.7368007848913e+06
    Quad (001) norm:  5.3333909168720e+00 Lin (021) step:  9.4442344936447e-01 func: -1.7368010144762e+06
    Quad (002) norm:  5.0555604725753e-01 Lin (033) step:  1.0070995028364e+00 func: -1.7368010148712e+06
    Quad (003) norm:  4.3269967709145e-04 Lin (036) step:  9.9782437462778e-01 func: -1.7368010148712e+06
    Quad (004) norm:  5.6091508120282e-06 Lin (034) step: -1.0586635525629e+00 func: -1.7368010148712e+06
    Quad (005) norm:  8.5509193027678e-06 Lin (041) step: -4.8450338014287e-01 func: -1.7368010148712e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      154.5 MPa
    Liquid          moles:   1.587824 grams: 102.153
               SiO2 form:  SiO2           X:  0.6860  wt%    SiO2   74.26
               TiO2 form:  TiO2           X:  0.0006  wt%    TiO2    0.08
              Al2O3 form:  Al2O3          X:  0.0413  wt%   Al2O3   11.68
              Fe2O3 form:  Fe2O3          X:  0.0008  wt%   Fe2O3    0.20
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.46
            Fe2SiO4 form:  Fe2SiO4        X:  0.0021  wt%     MgO    0.03
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.34
            Mg2SiO4 form:  Mg2SiO4        X:  0.0002  wt%    Na2O    3.68
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.74
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.54
             CaSiO3 form:  CaSiO3         X:  0.0039  
            Na2SiO3 form:  Na2SiO3        X:  0.0382  
            KAlSiO4 form:  KAlSiO4        X:  0.0647  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1621  
    Feldspar        moles:   0.009636 grams:   2.565
             albite form:  NaAlSi3O8      X:  0.7539  wt%    SiO2   64.11
          anorthite form:  CaAl2Si2O8     X:  0.1599  wt%   Al2O3   22.22
           sanidine form:  KAlSi3O8       X:  0.0861  wt%     CaO    3.37
                                                      wt%    Na2O    8.78
                                                      wt%     K2O    1.52
    Water           moles:   0.047865 grams:   0.862
    Quartz          affn:      59.08


.. code:: ipython3

    delta_dGdP = 0.2
    state = equil.execute(t, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  2.7393038780015e+02 Lin (015) step:  9.0330055107107e-01 func: -1.7369475154487e+06
    Quad (001) norm:  1.0049268744084e+00 Lin (023) step:  1.0030847801950e+00 func: -1.7369475493666e+06
    Quad (002) norm:  6.3639655636038e-02 Lin (026) step:  1.0030722713724e+00 func: -1.7369475494067e+06
    Quad (003) norm:  2.9042826684625e-05 Lin (034) step:  1.2268736865014e+00 func: -1.7369475494067e+06
    Quad (004) norm:  7.3456492597455e-06 Lin (038) step:  5.7479054798769e-01 func: -1.7369475494067e+06
    Quad (005) norm:  4.5951320971703e-06 Lin (036) step:  1.0556785668262e+00 func: -1.7369475494067e+06
    Minimal energy termination of quadratic loop.
    
    Add: Quartz
    Quad (000) norm:  4.3581885411668e+00 Lin (022) step:  9.9658677128886e-01 func: -1.7369475899600e+06
    Quad (001) norm:  6.0334980393502e-02 Lin (023) step:  1.0016062630753e+00 func: -1.7369475899676e+06
    Quad (002) norm:  1.5289903136215e-05 Lin (031) step:  9.0979700276985e-01 func: -1.7369475899676e+06
    Quad (003) norm:  1.3287186756146e-06 Lin (031) step:  1.2464095147559e+00 func: -1.7369475899676e+06
    Quad (004) norm:  1.7701007889124e-06 Lin (035) step: -4.7741675025314e-01 func: -1.7369475899676e+06
    Quad (005) norm:  2.5408682570821e-06 Lin (036) step: -4.7203595446392e-01 func: -1.7369475899676e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      139.7 MPa
    Liquid          moles:   1.533608 grams:  99.328
               SiO2 form:  SiO2           X:  0.6948  wt%    SiO2   74.64
               TiO2 form:  TiO2           X:  0.0007  wt%    TiO2    0.08
              Al2O3 form:  Al2O3          X:  0.0403  wt%   Al2O3   11.57
              Fe2O3 form:  Fe2O3          X:  0.0008  wt%   Fe2O3    0.21
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.48
            Fe2SiO4 form:  Fe2SiO4        X:  0.0021  wt%     MgO    0.03
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.30
            Mg2SiO4 form:  Mg2SiO4        X:  0.0002  wt%    Na2O    3.59
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.83
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.27
             CaSiO3 form:  CaSiO3         X:  0.0034  
            Na2SiO3 form:  Na2SiO3        X:  0.0376  
            KAlSiO4 form:  KAlSiO4        X:  0.0664  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1536  
    Feldspar        moles:   0.017389 grams:   4.627
             albite form:  NaAlSi3O8      X:  0.7608  wt%    SiO2   64.66
          anorthite form:  CaAl2Si2O8     X:  0.1369  wt%   Al2O3   21.79
           sanidine form:  KAlSi3O8       X:  0.1023  wt%     CaO    2.89
                                                      wt%    Na2O    8.86
                                                      wt%     K2O    1.81
    Water           moles:   0.069672 grams:   1.255
    Quartz          moles:   0.006164 grams:   0.370


.. code:: ipython3

    delta_dGdP = 0.3
    state = equil.execute(t, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  1.9011613487871e+02 Lin (021) step:  9.1968794873612e-01 func: -1.7370822071086e+06
    Quad (001) norm:  3.2315501406743e+00 Lin (033) step:  1.1071739284257e+00 func: -1.7370823269640e+06
    Quad (002) norm:  4.7380397694732e-01 Lin (026) step:  1.0159293218277e+00 func: -1.7370823284221e+06
    Quad (003) norm:  5.4846960932947e-03 Lin (033) step:  1.0026758224404e+00 func: -1.7370823284224e+06
    Quad (004) norm:  1.2018123200807e-05 Lin (036) step:  1.0449521286472e+00 func: -1.7370823284224e+06
    Quad (005) norm:  4.0048377278180e-08 Lin (035) step:  1.1997055567610e+00 func: -1.7370823284224e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      130.3 MPa
    Liquid          moles:   1.436206 grams:  93.694
               SiO2 form:  SiO2           X:  0.6990  wt%    SiO2   74.70
               TiO2 form:  TiO2           X:  0.0007  wt%    TiO2    0.09
              Al2O3 form:  Al2O3          X:  0.0394  wt%   Al2O3   11.57
              Fe2O3 form:  Fe2O3          X:  0.0009  wt%   Fe2O3    0.22
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.50
            Fe2SiO4 form:  Fe2SiO4        X:  0.0023  wt%     MgO    0.03
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.27
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    3.52
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.01
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.09
             CaSiO3 form:  CaSiO3         X:  0.0031  
            Na2SiO3 form:  Na2SiO3        X:  0.0370  
            KAlSiO4 form:  KAlSiO4        X:  0.0694  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1480  
    Feldspar        moles:   0.029239 grams:   7.783
             albite form:  NaAlSi3O8      X:  0.7541  wt%    SiO2   65.23
          anorthite form:  CaAl2Si2O8     X:  0.1104  wt%   Al2O3   21.27
           sanidine form:  KAlSi3O8       X:  0.1355  wt%     CaO    2.33
                                                      wt%    Na2O    8.78
                                                      wt%     K2O    2.40
    Water           moles:   0.092808 grams:   1.672
    Quartz          moles:   0.040472 grams:   2.432


.. code:: ipython3

    delta_dGdP = 0.4
    state = equil.execute(t, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  1.3071901675203e+02 Lin (020) step:  9.5827970688468e-01 func: -1.7372089789564e+06
    Quad (001) norm:  5.4275799328200e+00 Lin (036) step:  1.5873973382285e+00 func: -1.7372092241887e+06
    Quad (002) norm:  2.5782176646416e+00 Lin (028) step:  1.0254507646052e+00 func: -1.7372092338210e+06
    Quad (003) norm:  8.1478184005768e-02 Lin (028) step:  1.0043192189322e+00 func: -1.7372092338364e+06
    Quad (004) norm:  1.5787792364545e-04 Lin (013) step:  7.6393202250021e-01 func: -1.7372092338364e+06
    Quad (005) norm:  3.7768711868793e-05 Lin (036) step:  1.0482418281405e+00 func: -1.7372092338364e+06
    Minimal energy termination of quadratic loop.
    
    Unmixing: Feldspar
    Quad (000) norm:  2.0527253288879e+01 Lin (023) step:  7.7064765119618e-02 func: -1.7372099918662e+06
    Quad (001) norm:  2.8836543431755e+01 Lin (025) step:  9.8904550016763e-01 func: -1.7372108083609e+06
    Quad (002) norm:  1.3528947969788e+00 Lin (019) step:  1.2199012887383e+00 func: -1.7372108656548e+06
    Quad (003) norm:  3.4973655528523e-01 Lin (025) step:  1.0089945715836e+00 func: -1.7372108662996e+06
    Quad (004) norm:  2.4069483516336e-04 Lin (032) step:  9.8822558552838e-01 func: -1.7372108662997e+06
    Quad (005) norm:  1.4108063064080e-06 Lin (036) step: -4.7688092015815e-01 func: -1.7372108662997e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      127.5 MPa
    Liquid          moles:   1.219832 grams:  79.761
               SiO2 form:  SiO2           X:  0.6997  wt%    SiO2   74.65
               TiO2 form:  TiO2           X:  0.0008  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0392  wt%   Al2O3   11.54
              Fe2O3 form:  Fe2O3          X:  0.0011  wt%   Fe2O3    0.26
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.59
            Fe2SiO4 form:  Fe2SiO4        X:  0.0027  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.26
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    3.52
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.01
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.04
             CaSiO3 form:  CaSiO3         X:  0.0030  
            Na2SiO3 form:  Na2SiO3        X:  0.0371  
            KAlSiO4 form:  KAlSiO4        X:  0.0695  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1465  
    Feldspar        moles:   0.033282 grams:   8.860
             albite form:  NaAlSi3O8      X:  0.7515  wt%    SiO2   65.36
          anorthite form:  CaAl2Si2O8     X:  0.1040  wt%   Al2O3   21.14
           sanidine form:  KAlSi3O8       X:  0.1444  wt%     CaO    2.19
                                                      wt%    Na2O    8.75
                                                      wt%     K2O    2.55
    Water           moles:   0.126582 grams:   2.280
    Quartz          moles:   0.120809 grams:   7.259
    Feldspar        moles:   0.027404 grams:   7.420
             albite form:  NaAlSi3O8      X:  0.4701  wt%    SiO2   66.14
          anorthite form:  CaAl2Si2O8     X:  0.0193  wt%   Al2O3   19.19
           sanidine form:  KAlSi3O8       X:  0.5105  wt%     CaO    0.40
                                                      wt%    Na2O    5.38
                                                      wt%     K2O    8.88



