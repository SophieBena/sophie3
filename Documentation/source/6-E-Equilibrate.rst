Internal energy potential minimization (S, V, constrained)
==========================================================

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

Function to constrain the entropy and the volume
------------------------------------------------

Note that the entropy is equivalent to $ -
:raw-latex:`\frac{{\partial G}}{{\partial T}}`$ and that the volume is
equivalent to :math:`\frac{{\partial G}}{{\partial P}}` - Run an
equilibration step at fixed T,P - Calculate the entropy and volume of
the system - Define functions to set the entropy and volume for
subsequent equilibration steps

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
    def conT(t, p, state):
        return dGdT + delta_dGdT
    print (conT(t, p, None), state.dGdT(t,p))
    delta_dGdP = 0.0
    dGdP = state.dGdP(t,p)
    def conP(t, p, state):
        return dGdP + delta_dGdP


.. parsed-literal::

    -272.53139197088126 -272.53139197088126


Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[('T',conT),('P',conP)])

.. code:: ipython3

    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Water
    Quad (000) norm:  2.7440425792652e+03 Lin (018) step: -3.6631021024202e-01 func: -1.4504785905582e+06
    Quad (001) norm:  5.6241198998666e-05 Lin (032) step:  1.0633982368735e+00 func: -1.4504785905583e+06
    Quad (002) norm:  3.8315742232314e-06 Lin (040) step: -7.3760884727322e-01 func: -1.4504785905583e+06
    Quad (003) norm:  6.5187799394262e-06 Lin (030) step:  2.2291077929404e-01 func: -1.4504785905583e+06
    Quad (004) norm:  4.0347687869948e-06 Lin (038) step: -3.8776887618691e-01 func: -1.4504785905583e+06
    Quad (005) norm:  5.8922812452171e-06 Lin (034) step:  4.0330319124792e-01 func: -1.4504785905583e+06
    Minimal energy termination of quadratic loop.
    
    Add: Feldspar
    Quad (000) norm:  1.4255918044039e+01 Lin (024) step:  9.8721565153341e-01 func: -1.4504788192807e+06
    Quad (001) norm:  5.5851657877330e-01 Lin (024) step:  1.0293237248731e+00 func: -1.4504788251705e+06
    Quad (002) norm:  1.6681932555171e-01 Lin (022) step:  1.0038995375858e+00 func: -1.4504788252135e+06
    Quad (003) norm:  3.4194341448514e-04 Lin (007) step:  1.2360679958336e+00 func: -1.4504788252135e+06
    Quad (004) norm:  8.2574899086290e-05 Lin (038) step:  1.6168503100357e+00 func: -1.4504788252135e+06
    Quad (005) norm:  5.4155487274901e-05 Lin (041) step: -3.9575049032453e-01 func: -1.4504788252135e+06
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
    delta_dGdP = 0.1
    state = equil.execute(t, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  1.0665794695480e+03 Lin (017) step:  7.6698918214359e-01 func: -1.4558246702037e+06
    Quad (001) norm:  4.5920341327916e+01 Lin (025) step:  1.0798000483052e+00 func: -1.4558315783073e+06
    Quad (002) norm:  7.8199343002471e-01 Lin (025) step:  1.3227692436755e+00 func: -1.4558318887167e+06
    Quad (003) norm:  1.1841611901740e+00 Lin (028) step:  9.9872636995588e-01 func: -1.4558318924202e+06
    Quad (004) norm:  4.5083625782551e-03 Lin (015) step:  1.0180817909720e+00 func: -1.4558318924203e+06
    Quad (005) norm:  7.8522833105113e-05 Lin (036) step:  1.4160247695405e+00 func: -1.4558318924203e+06
    Minimal energy termination of quadratic loop.
    
    Unmixing: Feldspar
    Add: Quartz
    Quad (000) norm:  9.3008039690192e+01 Lin (017) step:  8.9913103195612e-01 func: -1.4558532895168e+06
    Quad (001) norm:  1.0846094964032e+01 Lin (030) step:  1.1797971119101e+00 func: -1.4558534331925e+06
    Quad (002) norm:  1.3238149403503e+00 Lin (027) step:  1.0154286248998e+00 func: -1.4558534345972e+06
    Quad (003) norm:  2.1525493814175e-02 Lin (034) step:  1.0043899897495e+00 func: -1.4558534345977e+06
    Quad (004) norm:  9.1170103770516e-05 Lin (034) step:  2.6548843415844e-01 func: -1.4558534345977e+06
    Quad (005) norm:  6.6781933106375e-05 Lin (033) step: -9.1130089900093e-01 func: -1.4558534345977e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     764.06 °C, P =      169.5 MPa
    Liquid          moles:   1.298234 grams:  83.372
               SiO2 form:  SiO2           X:  0.6794  wt%    SiO2   73.89
               TiO2 form:  TiO2           X:  0.0008  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0398  wt%   Al2O3   11.63
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.25
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.57
            Fe2SiO4 form:  Fe2SiO4        X:  0.0025  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.27
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    3.63
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.90
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.73
             CaSiO3 form:  CaSiO3         X:  0.0031  
            Na2SiO3 form:  Na2SiO3        X:  0.0376  
            KAlSiO4 form:  KAlSiO4        X:  0.0669  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1686  
    Feldspar        moles:   0.027010 grams:   7.188
             albite form:  NaAlSi3O8      X:  0.7558  wt%    SiO2   65.01
          anorthite form:  CaAl2Si2O8     X:  0.1205  wt%   Al2O3   21.46
           sanidine form:  KAlSi3O8       X:  0.1237  wt%     CaO    2.54
                                                      wt%    Na2O    8.80
                                                      wt%     K2O    2.19
    Water           moles:   0.086351 grams:   1.556
    Quartz          moles:   0.114217 grams:   6.863
    Feldspar        moles:   0.024321 grams:   6.601
             albite form:  NaAlSi3O8      X:  0.4297  wt%    SiO2   66.04
          anorthite form:  CaAl2Si2O8     X:  0.0167  wt%   Al2O3   19.10
           sanidine form:  KAlSi3O8       X:  0.5536  wt%     CaO    0.35
                                                      wt%    Na2O    4.91
                                                      wt%     K2O    9.61


