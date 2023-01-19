Korzhinskii potential minimization (T, P, :math:`\mu`\ SiO2 constrained)
========================================================================

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

Note that quartz saturation will be imposed, so adding quartz as a
system phase is redundant and merely moditors the stauration condition.

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

Function to constrain the chemical potential of SiO2
----------------------------------------------------

.. code:: ipython3

    def muSiO2(t, p, state):
        return Quartz.gibbs_energy(t, p)

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[({'Si':1.0,'O':2.0},muSiO2)])

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    
    Add: Water
    Quad (000) norm:  2.0919339116256e-02 Lin (022) step:  1.1404184587577e+00 func: -6.3463025951454e+05
    Quad (001) norm:  3.3488088202158e-03 Lin (034) step:  1.1687446083504e+00 func: -6.3463005405689e+05
    Quad (002) norm:  5.5404363799077e-04 Lin (024) step:  1.1640376709056e+00 func: -6.3463004850062e+05
    Quad (003) norm:  9.1181655016418e-05 Lin (014) step:  1.1642086859434e+00 func: -6.3463004832265e+05
    Quad (004) norm:  1.4962742159588e-05 Lin (033) step:  1.1664159924080e+00 func: -6.3463004833999e+05
    Quad (005) norm:  2.4903355411261e-06 Lin (035) step:  1.1721712654619e+00 func: -6.3463004834073e+05
    Minimal energy termination of quadratic loop.
    
    Add: Quartz
    Quad (000) norm:  4.2865017842861e-07 Lin (037) step:  1.0227907354696e+00 func: -6.3463989322211e+05
    Quad (001) norm:  9.6509648395378e-09 Lin (038) step: -3.8716905600240e-01 func: -6.3463989322079e+05
    Quad (002) norm:  1.3387516320064e-08 Lin (035) step:  1.6393202493172e+00 func: -6.3463989322079e+05
    Quad (003) norm:  8.5590033366465e-09 Lin (040) step: -1.0509103027852e+00 func: -6.3463989322079e+05
    Quad (004) norm:  1.7554178136231e-08 Lin (038) step:  8.9452579470777e-01 func: -6.3463989322079e+05
    Quad (005) norm:  1.8506933462712e-09 Lin (050) step:  9.0631518900942e-03 func: -6.3463989322079e+05
    Minimal energy termination of quadratic loop.
    
     
    T =     776.85 °C, P =      175.0 MPa
    Liquid          moles:   1.711813 grams: 108.664
               SiO2 form:  SiO2           X:  0.6808  wt%    SiO2   74.34
               TiO2 form:  TiO2           X:  0.0006  wt%    TiO2    0.07
              Al2O3 form:  Al2O3          X:  0.0414  wt%   Al2O3   11.50
              Fe2O3 form:  Fe2O3          X:  0.0008  wt%   Fe2O3    0.19
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.44
            Fe2SiO4 form:  Fe2SiO4        X:  0.0019  wt%     MgO    0.03
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.40
            Mg2SiO4 form:  Mg2SiO4        X:  0.0002  wt%    Na2O    3.66
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.49
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.88
             CaSiO3 form:  CaSiO3         X:  0.0045  
            Na2SiO3 form:  Na2SiO3        X:  0.0375  
            KAlSiO4 form:  KAlSiO4        X:  0.0605  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1719  
    Feldspar        affn:     264.01
             albite form:  NaAlSi3O8      X:  0.7411
          anorthite form:  CaAl2Si2O8     X:  0.1888
           sanidine form:  KAlSi3O8       X:  0.0701
    Water           moles:   0.011074 grams:   0.200
    Quartz          moles:   0.000010 grams:   0.001


Pickup runs use previously computed state

.. code:: ipython3

    state = equil.execute(t-5.0, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Feldspar
    Quad (000) norm:  4.1530329504509e-01 Lin (037) step: -5.5625677045860e-02 func: -6.5238951381008e+05
    Quad (001) norm:  7.0429498604093e-01 Lin (038) step:  4.4370623407854e-03 func: -6.5447193926378e+05
    Quad (002) norm:  5.6770585823860e-01 Lin (037) step: -9.2767951349095e-03 func: -6.5736313225043e+05
    Quad (003) norm:  6.3586523505206e-01 Lin (040) step: -1.0150575486173e-03 func: -6.5737880325226e+05
    Quad (004) norm:  5.7201759044269e-01 Lin (037) step:  3.9933561616004e-03 func: -6.5757808862571e+05
    Quad (005) norm:  9.5799555960644e-01 Lin (038) step: -6.0865053398050e-03 func: -6.5894764733790e+05
    Quad (006) norm:  2.1868223831377e+00 Lin (043) step:  3.5392868187752e-04 func: -6.5902061774685e+05
    Quad (007) norm:  2.1471153681086e+00 Lin (048) step:  3.8726518890751e-05 func: -6.5902800631958e+05
    Quad (008) norm:  2.1369697719763e+00 Lin (058) step:  1.7962899250578e-07 func: -6.5902803994307e+05
    Remove: Water
    Quad (000) norm:  1.8711143615073e+00 Lin (039) step:  3.9933561294381e-03 func: -6.5973127077201e+05
    Quad (001) norm:  7.3595725016869e-01 Lin (037) step: -2.4315330541724e-01 func: -7.7684298032846e+05
    Quad (002) norm:  1.4478457722711e+00 Lin (037) step:  1.4139300658851e-02 func: -7.8112221305756e+05
    Quad (003) norm:  1.0430285089303e+00 Lin (037) step:  1.2725370572772e-02 func: -7.9113032478277e+05
    Quad (004) norm:  6.5323615637559e-01 Lin (037) step: -1.0466952493557e-01 func: -8.2691962540789e+05
    Quad (005) norm:  6.6461155267123e-01 Lin (040) step: -2.6200409720147e-03 func: -8.2796047317863e+05
    Quad (006) norm:  6.8415457564177e-01 Lin (040) step: -2.6200409720147e-03 func: -8.2924105627869e+05
    Quad (007) norm:  6.9378571255856e-01 Lin (040) step: -2.6200409720147e-03 func: -8.3051766289862e+05
    Quad (008) norm:  7.0495145670417e-01 Lin (040) step: -2.6200409720147e-03 func: -8.3175832252417e+05
    Quad (009) norm:  7.2036448957476e-01 Lin (040) step: -2.6200409720147e-03 func: -8.3295174507994e+05
    Quad (010) norm:  7.4171934579853e-01 Lin (040) step: -2.6200409720147e-03 func: -8.3408697267458e+05
    Quad (011) norm:  7.7081352312931e-01 Lin (040) step: -2.6200409647262e-03 func: -8.3515203750756e+05
    Quad (012) norm:  8.0939421471903e-01 Lin (040) step: -2.6200409647262e-03 func: -8.3613426427993e+05
    Quad (013) norm:  8.5880776705343e-01 Lin (040) step: -2.6200409647262e-03 func: -8.3702087385821e+05
    Quad (014) norm:  9.1962373385636e-01 Lin (040) step: -2.6200409566279e-03 func: -8.3779976099576e+05
    Quad (015) norm:  9.9134479158062e-01 Lin (040) step: -2.6200409566279e-03 func: -8.3846000656735e+05
    Quad (016) norm:  1.0724205980903e+00 Lin (040) step: -2.6200409502138e-03 func: -8.3900026475082e+05
    Quad (017) norm:  8.4147166109310e-16
    Add: Water
    Quad (000) norm:  5.9873014135237e-01 Lin (042) step: -4.3029465263668e-05 func: -8.3902618254186e+05
    Quad (001) norm:  5.5438854091657e-01 Lin (060) step: -4.2374165324941e-06 func: -8.3902846091424e+05
    Quad (002) norm:  5.7445840878890e-01 Lin (053) step: -1.3094953623177e-07 func: -8.3902853755442e+05
    Remove: Water
    Quad (000) norm:  8.2723861329952e-16
    Add: Water
    Quad (000) norm:  5.7256978130451e-01 Lin (055) step: -4.7810516978195e-05 func: -8.3905503858508e+05
    Quad (001) norm:  5.7502051822156e-01 Lin (065) step: -2.4640465361673e-07 func: -8.3905517349952e+05
    Remove: Water
    Quad (000) norm:  8.1456766203904e-16
    Add: Water
    Quad (000) norm:  5.9557552355154e-01 Lin (042) step: -4.7810516682650e-05 func: -8.3908321531022e+05
    Quad (001) norm:  5.7081950517278e-01 Lin (055) step: -6.2632755853170e-08 func: -8.3908325171622e+05
    Remove: Water
    Quad (000) norm:  8.0245280288361e-16
    Add: Water
    Quad (000) norm:  7.6453909807180e-16
     
    T =     771.85 °C, P =      175.0 MPa
    Liquid          moles:   1.355664 grams:  81.544
               SiO2 form:  SiO2           X:  0.6603  wt%    SiO2   75.37
               TiO2 form:  TiO2           X:  0.0007  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0186  wt%   Al2O3    9.27
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.25
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.58
            Fe2SiO4 form:  Fe2SiO4        X:  0.0024  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.00
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    2.00
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.65
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    6.74
             CaSiO3 form:  CaSiO3         X:  0.0000  
            Na2SiO3 form:  Na2SiO3        X:  0.0194  
            KAlSiO4 form:  KAlSiO4        X:  0.0722  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.2252  
    Feldspar        moles:   0.089326 grams:  23.638
             albite form:  NaAlSi3O8      X:  0.8499  wt%    SiO2   66.17
          anorthite form:  CaAl2Si2O8     X:  0.0858  wt%   Al2O3   20.92
           sanidine form:  KAlSi3O8       X:  0.0643  wt%     CaO    1.82
                                                      wt%    Na2O    9.95
                                                      wt%     K2O    1.14
    Water           moles:   0.000010 grams:   0.000
    Quartz          moles:   0.000012 grams:   0.001


.. code:: ipython3

    state = equil.execute(t-10.0, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  5.6688245635369e-01 Lin (043) step: -4.7810516655517e-05 func: -8.3831498931258e+05
    Quad (001) norm:  7.5407306663823e-16
     
    T =     766.85 °C, P =      175.0 MPa
    Liquid          moles:   1.342725 grams:  80.765
               SiO2 form:  SiO2           X:  0.6570  wt%    SiO2   75.14
               TiO2 form:  TiO2           X:  0.0007  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0187  wt%   Al2O3    9.35
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.26
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.59
            Fe2SiO4 form:  Fe2SiO4        X:  0.0025  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.00
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    2.01
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.71
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    6.81
             CaSiO3 form:  CaSiO3         X:  0.0000  
            Na2SiO3 form:  Na2SiO3        X:  0.0196  
            KAlSiO4 form:  KAlSiO4        X:  0.0729  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.2274  
    Feldspar        moles:   0.089337 grams:  23.641
             albite form:  NaAlSi3O8      X:  0.8499  wt%    SiO2   66.17
          anorthite form:  CaAl2Si2O8     X:  0.0858  wt%   Al2O3   20.92
           sanidine form:  KAlSi3O8       X:  0.0643  wt%     CaO    1.82
                                                      wt%    Na2O    9.95
                                                      wt%     K2O    1.14
    Water           moles:   0.000001 grams:   0.000
    Quartz          moles:   0.000012 grams:   0.001


.. code:: ipython3

    state = equil.execute(t-15.0, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  7.5572812753852e-16
     
    T =     761.85 °C, P =      175.0 MPa
    Liquid          moles:   1.329912 grams:  79.995
               SiO2 form:  SiO2           X:  0.6537  wt%    SiO2   74.90
               TiO2 form:  TiO2           X:  0.0008  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0189  wt%   Al2O3    9.44
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.26
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.59
            Fe2SiO4 form:  Fe2SiO4        X:  0.0025  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.00
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    2.03
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.76
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    6.88
             CaSiO3 form:  CaSiO3         X:  0.0000  
            Na2SiO3 form:  Na2SiO3        X:  0.0197  
            KAlSiO4 form:  KAlSiO4        X:  0.0736  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.2296  
    Feldspar        moles:   0.089337 grams:  23.641
             albite form:  NaAlSi3O8      X:  0.8499  wt%    SiO2   66.17
          anorthite form:  CaAl2Si2O8     X:  0.0858  wt%   Al2O3   20.92
           sanidine form:  KAlSi3O8       X:  0.0643  wt%     CaO    1.82
                                                      wt%    Na2O    9.95
                                                      wt%     K2O    1.14
    Water           moles:   0.000001 grams:   0.000
    Quartz          moles:   0.000012 grams:   0.001


.. code:: ipython3

    state = equil.execute(t-20.0, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Quad (000) norm:  7.5733316417540e-16
     
    T =     756.85 °C, P =      175.0 MPa
    Liquid          moles:   1.317201 grams:  79.232
               SiO2 form:  SiO2           X:  0.6504  wt%    SiO2   74.66
               TiO2 form:  TiO2           X:  0.0008  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0191  wt%   Al2O3    9.53
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.26
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.60
            Fe2SiO4 form:  Fe2SiO4        X:  0.0025  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.00
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    2.05
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.82
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    6.94
             CaSiO3 form:  CaSiO3         X:  0.0000  
            Na2SiO3 form:  Na2SiO3        X:  0.0199  
            KAlSiO4 form:  KAlSiO4        X:  0.0743  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.2318  
    Feldspar        moles:   0.089337 grams:  23.641
             albite form:  NaAlSi3O8      X:  0.8499  wt%    SiO2   66.17
          anorthite form:  CaAl2Si2O8     X:  0.0858  wt%   Al2O3   20.92
           sanidine form:  KAlSi3O8       X:  0.0643  wt%     CaO    1.82
                                                      wt%    Na2O    9.95
                                                      wt%     K2O    1.14
    Water           moles:   0.000001 grams:   0.000
    Quartz          moles:   0.000012 grams:   0.001


