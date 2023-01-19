Gibbs energy minimization (fixed T, P, bulk composition)
========================================================

Closed system; crystallization of a rhyolitic liquid using
rhyolite-MELTS

.. code:: ipython3

    import numpy as np
    import scipy.optimize as opt
    import scipy.linalg as lin 
    import sys

.. code:: ipython3

    from thermoengine import core, phases, model, equilibrate

Create phases for equilibrium assemblages
-----------------------------------------

.. code:: ipython3

    modelDB = model.Database(liq_mod='v1.0')

.. code:: ipython3

    Liquid = modelDB.get_phase('Liq')
    Feldspar = modelDB.get_phase('Fsp')
    Quartz = modelDB.get_phase('Qz')
    Spinel = modelDB.get_phase('SplS')
    Opx = modelDB.get_phase('Opx')
    RhomOx = modelDB.get_phase('Rhom')

The Berman model database provides the SWIM water model by default.
Instead, override that choice by instantiating the MELTS 1.0.2 water
model directly.

.. code:: ipython3

    Water = phases.PurePhase('WaterMelts', 'H2O', calib=False)

Define elements in system and phases in system
----------------------------------------------

.. code:: ipython3

    elm_sys = ['H','O','Na','Mg','Al','Si','P','K','Ca','Ti','Cr','Mn','Fe','Co','Ni']
    phs_sys = [Liquid, Feldspar, Water, Quartz, Spinel, Opx, RhomOx]

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

Instantiate class instance and run calculation
----------------------------------------------

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
    Spinel          affn:    5163.80
           chromite form:  FeCr2O4        X:  0.0000
          hercynite form:  FeAl2O4        X:  0.0268
          magnetite form:  Fe3O4          X:  0.8177
             spinel form:  MgAl2O4        X:  0.0290
         ulvospinel form:  Fe2TiO4        X:  0.1264
    Orthopyroxene   affn:    9117.94
           diopside form:  CaMgSi2O6      X: -1.3944
     clinoenstatite form:  Mg2Si2O6       X:  0.9774
       hedenbergite form:  CaFeSi2O6      X:  1.3325
    alumino-buffoni form:  CaTi0.5Mg0     X:  0.0568
          buffonite form:  CaTi0.5Mg0     X: -0.0563
           essenite form:  CaFeAlSiO6     X:  0.0788
            jadeite form:  NaAlSi2O6      X:  0.0051
    Ilmenite ss     affn:    8522.12
         geikielite form:  MgTiO3         X:  0.0223
           hematite form:  Fe2O3          X:  0.4727
           ilmenite form:  FeTiO3         X:  0.4815
        pyrophanite form:  MnTiO3         X:  0.0000
           corundum form:  Al2O3          X:  0.0234


Pickup runs use previously computed state

.. code:: ipython3

    state = equil.execute(t-20.0, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Quartz
    Quad (000) norm:  1.7584052453349e-01 Lin (021) step:  7.5129798136470e-01 func: -1.7226695512912e+06
    Quad (001) norm:  5.0942901721106e-02 Lin (017) step:  1.5968504094272e+00 func: -1.7226753142244e+06
    Quad (002) norm:  4.2504295099397e-02 Lin (038) step:  1.9999999656497e+00 func: -1.7226784944406e+06
    Quad (003) norm:  1.1222069968710e+00 Lin (018) step:  4.2609941999829e-01 func: -1.7226940452527e+06
    Quad (004) norm:  2.2723498715356e-01 Lin (019) step:  7.0857168341196e-01 func: -1.7227041627379e+06
    Quad (005) norm:  2.7128481630673e-02 Lin (019) step:  1.0354965272624e+00 func: -1.7227045419501e+06
    Quad (006) norm:  2.3033032581156e-03 Lin (028) step:  1.0060288856278e+00 func: -1.7227045429718e+06
    Quad (007) norm:  9.7478563615133e-07 Lin (036) step:  1.0680437195760e+00 func: -1.7227045429718e+06
    Minimal energy termination of quadratic loop.
    
    Unmixing: Feldspar
    Add: Spinel
    Quad (000) norm:  5.9180648192509e-01 Lin (021) step:  2.5880128423739e-01 func: -1.7227192124839e+06
    Quad (001) norm:  8.6477379612510e-02 Lin (013) step:  1.0392730169330e+00 func: -1.7227251379548e+06
    Quad (002) norm:  1.4062435502107e-02 Lin (021) step:  9.9961786239296e-01 func: -1.7227257460969e+06
    Quad (003) norm:  8.9946821614493e-04 Lin (033) step:  1.1132054779089e+00 func: -1.7227258012406e+06
    Quad (004) norm:  2.4708159733211e-04 Lin (015) step:  9.9995574046988e-01 func: -1.7227258015134e+06
    Quad (005) norm:  6.7982347077598e-07 Lin (033) step:  9.8344720672610e-01 func: -1.7227258015134e+06
    Minimal energy termination of quadratic loop.
    
    Add: Orthopyroxene
    Quad (000) norm:  3.2537713794126e-01 Lin (017) step:  2.3630100340346e-01 func: -1.7227289609526e+06
    Quad (001) norm:  4.8547032241559e-02 Lin (025) step:  1.1168103153307e+00 func: -1.7227310876293e+06
    Quad (002) norm:  2.3452045212057e-03 Lin (016) step:  1.2282298196427e+00 func: -1.7227312111040e+06
    Quad (003) norm:  9.2085985309302e-05 Lin (021) step:  1.0404811807507e+00 func: -1.7227312149335e+06
    Quad (004) norm:  1.1085514763767e-05 Lin (025) step:  1.0010724791354e+00 func: -1.7227312149374e+06
    Quad (005) norm:  1.0596772363073e-08 Lin (045) step:  6.7962862673539e-02 func: -1.7227312149374e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     756.85 °C, P =      175.0 MPa
    Liquid          moles:   0.069015 grams:   4.402
               SiO2 form:  SiO2           X:  0.6693  wt%    SiO2   73.96
               TiO2 form:  TiO2           X:  0.0009  wt%    TiO2    0.12
              Al2O3 form:  Al2O3          X:  0.0334  wt%   Al2O3   10.01
              Fe2O3 form:  Fe2O3          X:  0.0006  wt%   Fe2O3    0.14
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    1.41
            Fe2SiO4 form:  Fe2SiO4        X:  0.0063  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.26
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    4.66
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.31
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    5.08
             CaSiO3 form:  CaSiO3         X:  0.0030  
            Na2SiO3 form:  Na2SiO3        X:  0.0480  
            KAlSiO4 form:  KAlSiO4        X:  0.0584  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1799  
    Feldspar        moles:   0.176184 grams:  47.725
             albite form:  NaAlSi3O8      X:  0.4624  wt%    SiO2   66.20
          anorthite form:  CaAl2Si2O8     X:  0.0154  wt%   Al2O3   19.11
           sanidine form:  KAlSi3O8       X:  0.5222  wt%     CaO    0.32
                                                      wt%    Na2O    5.29
                                                      wt%     K2O    9.08
    Water           moles:   0.292885 grams:   5.276
    Quartz          moles:   0.553322 grams:  33.246
    Spinel          moles:   0.002240 grams:   0.506
           chromite form:  FeCr2O4        X:  0.0000  wt%    TiO2   14.80
          hercynite form:  FeAl2O4        X:  0.0066  wt%   Al2O3    1.31
          magnetite form:  Fe3O4          X:  0.5527  wt%   Fe2O3   39.08
             spinel form:  MgAl2O4        X:  0.0223  wt%     FeO   44.41
         ulvospinel form:  Fe2TiO4        X:  0.4183  wt%     MgO    0.40
    Orthopyroxene   moles:   0.001663 grams:   0.417
           diopside form:  CaMgSi2O6      X: -1.5998  wt%    SiO2   46.82
     clinoenstatite form:  Mg2Si2O6       X:  0.9943  wt%    TiO2    0.01
       hedenbergite form:  CaFeSi2O6      X:  1.5577  wt%   Al2O3    1.46
    alumino-buffoni form:  CaTi0.5Mg0     X:  0.0245  wt%   Fe2O3    0.70
          buffonite form:  CaTi0.5Mg0     X: -0.0240  wt%     FeO   44.64
           essenite form:  CaFeAlSiO6     X:  0.0461  wt%     MgO    6.26
            jadeite form:  NaAlSi2O6      X:  0.0011  wt%     CaO    0.10
                                                      wt%    Na2O    0.01
    Ilmenite ss     affn:     258.13
         geikielite form:  MgTiO3         X:  0.0339
           hematite form:  Fe2O3          X:  0.0567
           ilmenite form:  FeTiO3         X:  0.8836
        pyrophanite form:  MnTiO3         X:  0.0000
           corundum form:  Al2O3          X:  0.0259
    Feldspar        moles:   0.052661 grams:  14.007
             albite form:  NaAlSi3O8      X:  0.7659  wt%    SiO2   65.73
          anorthite form:  CaAl2Si2O8     X:  0.0901  wt%   Al2O3   20.89
           sanidine form:  KAlSi3O8       X:  0.1440  wt%     CaO    1.90
                                                      wt%    Na2O    8.92
                                                      wt%     K2O    2.55


.. code:: ipython3

    state = equil.execute(t-25.0, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Ilmenite ss
    Quad (000) norm:  5.7404029143509e-02 Lin (042) step:  1.4357959490477e-01 func: -1.7214610251210e+06
    Quad (001) norm:  2.9813648107305e-02 Lin (046) step:  1.9395474103949e-02 func: -1.7214610553562e+06
    Quad (002) norm:  2.8228058054619e-02 Lin (051) step:  1.3923968677715e-03 func: -1.7214610573914e+06
    Remove: Ilmenite ss
    Quad (000) norm:  2.8069556006391e-02 Lin (022) step:  5.5356986504905e-01 func: -1.7214615347076e+06
    Quad (001) norm:  1.0076526073353e-03 Lin (026) step:  9.8733594024009e-01 func: -1.7214615415758e+06
    Quad (002) norm:  8.4270612792209e-06 Lin (035) step:  1.0107926159604e+00 func: -1.7214615415778e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     751.85 °C, P =      175.0 MPa
    Liquid          moles:   0.038735 grams:   2.455
               SiO2 form:  SiO2           X:  0.6656  wt%    SiO2   74.42
               TiO2 form:  TiO2           X:  0.0008  wt%    TiO2    0.10
              Al2O3 form:  Al2O3          X:  0.0281  wt%   Al2O3    8.84
              Fe2O3 form:  Fe2O3          X:  0.0006  wt%   Fe2O3    0.16
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    1.45
            Fe2SiO4 form:  Fe2SiO4        X:  0.0064  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.30
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    5.42
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.00
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    5.28
             CaSiO3 form:  CaSiO3         X:  0.0034  
            Na2SiO3 form:  Na2SiO3        X:  0.0554  
            KAlSiO4 form:  KAlSiO4        X:  0.0538  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1856  
    Feldspar        moles:   0.177752 grams:  48.161
             albite form:  NaAlSi3O8      X:  0.4584  wt%    SiO2   66.20
          anorthite form:  CaAl2Si2O8     X:  0.0146  wt%   Al2O3   19.09
           sanidine form:  KAlSi3O8       X:  0.5270  wt%     CaO    0.30
                                                      wt%    Na2O    5.24
                                                      wt%     K2O    9.16
    Water           moles:   0.298109 grams:   5.370
    Quartz          moles:   0.563823 grams:  33.877
    Spinel          moles:   0.002280 grams:   0.515
           chromite form:  FeCr2O4        X:  0.0000  wt%    TiO2   15.03
          hercynite form:  FeAl2O4        X:  0.0039  wt%   Al2O3    1.12
          magnetite form:  Fe3O4          X:  0.5500  wt%   Fe2O3   38.85
             spinel form:  MgAl2O4        X:  0.0210  wt%     FeO   44.63
         ulvospinel form:  Fe2TiO4        X:  0.4252  wt%     MgO    0.37
    Orthopyroxene   moles:   0.001823 grams:   0.458
           diopside form:  CaMgSi2O6      X: -1.6170  wt%    SiO2   46.86
     clinoenstatite form:  Mg2Si2O6       X:  0.9930  wt%    TiO2    0.01
       hedenbergite form:  CaFeSi2O6      X:  1.5842  wt%   Al2O3    1.17
    alumino-buffoni form:  CaTi0.5Mg0     X:  0.0185  wt%   Fe2O3    0.64
          buffonite form:  CaTi0.5Mg0     X: -0.0180  wt%     FeO   45.26
           essenite form:  CaFeAlSiO6     X:  0.0382  wt%     MgO    5.92
            jadeite form:  NaAlSi2O6      X:  0.0012  wt%     CaO    0.13
                                                      wt%    Na2O    0.01
    Ilmenite ss     affn:     164.56
         geikielite form:  MgTiO3         X:  0.0321
           hematite form:  Fe2O3          X:  0.0551
           ilmenite form:  FeTiO3         X:  0.8914
        pyrophanite form:  MnTiO3         X:  0.0000
           corundum form:  Al2O3          X:  0.0215
    Feldspar        moles:   0.055437 grams:  14.742
             albite form:  NaAlSi3O8      X:  0.7693  wt%    SiO2   65.77
          anorthite form:  CaAl2Si2O8     X:  0.0890  wt%   Al2O3   20.88
           sanidine form:  KAlSi3O8       X:  0.1417  wt%     CaO    1.88
                                                      wt%    Na2O    8.96
                                                      wt%     K2O    2.51


