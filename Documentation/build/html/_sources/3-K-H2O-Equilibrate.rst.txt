Korzhinskii potential minimization (T, P, :math:`\mu`\ H2O constrained)
=======================================================================

Open system; crystallization of a rhyolitic liquid using rhyolite-MELTS

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

The Berman model database provides the SWIM water model by default.
Instead, override that choice by instantiating the MELTS 1.0.2 water
model directly.

.. code:: ipython3

    Water = phases.PurePhase('WaterMelts', 'H2O', calib=False)

Define elements in system and phases in system
----------------------------------------------

.. code:: ipython3

    elm_sys = ['H','O','Na','Mg','Al','Si','P','K','Ca','Ti','Cr','Mn','Fe','Co','Ni']
    phs_sys = [Liquid, Feldspar, Quartz, Water]

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

Method to constrain the chemical potential of H2O
-------------------------------------------------

This method is passed to the Equilibrate class and is used to set the
chemical potential of water in the system to that of pure water. This is
equivalent to forcing the system to be saturated with water for all T
and P.

.. code:: ipython3

    def muH2O(t, p, state):
        return Water.gibbs_energy(t, p)

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[({'H':2.0,'O':1.0},muH2O)])

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    
    Add: Feldspar
    Quad (000) norm:  5.3955099899240e-03 Lin (026) step:  9.2108450439213e-01 func: -1.6067853526829e+06
    Quad (001) norm:  5.0689592065047e-04 Lin (028) step:  1.1191143427707e+00 func: -1.6067853659542e+06
    Quad (002) norm:  5.2093020609536e-05 Lin (033) step:  9.6074094554648e-01 func: -1.6067853661133e+06
    Quad (003) norm:  2.1626712406157e-06 Lin (036) step:  9.7782391036993e-01 func: -1.6067853660797e+06
    Quad (004) norm:  4.8163169842462e-08 Lin (035) step:  1.6232198623426e+00 func: -1.6067853660804e+06
    Quad (005) norm:  3.0019986294382e-08 Lin (039) step:  1.3740686088664e+00 func: -1.6067853660804e+06
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
    Quartz          affn:     134.38
    Water           affn:       0.00


.. code:: ipython3

    t = 1030.0
    p = 1750.0
    state = equil.execute(t, p, state=state, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    Add: Quartz
    Quad (000) norm:  1.5738316333786e-01 Lin (024) step:  7.8904074183554e-01 func: -1.6023662509952e+06
    Quad (001) norm:  6.1830006723981e-02 Lin (018) step:  1.5215401431957e+00 func: -1.6023721083608e+06
    Quad (002) norm:  2.8764208858874e-02 Lin (038) step:  1.9999999656497e+00 func: -1.6023752424417e+06
    Quad (003) norm:  2.9818865344059e+00 Lin (016) step: -1.5518893558448e-01 func: -1.6023891660537e+06
    Quad (004) norm:  1.6819900238397e-01 Lin (025) step:  7.7857669133260e-01 func: -1.6024008662770e+06
    Quad (005) norm:  1.4001168241950e-02 Lin (021) step:  9.6028681052471e-01 func: -1.6024013429625e+06
    Quad (006) norm:  4.4326808519568e-03 Lin (025) step:  9.9446260072970e-01 func: -1.6024013517995e+06
    Quad (007) norm:  4.1333324371243e-05 Lin (029) step:  1.0116971344263e+00 func: -1.6024013518454e+06
    Minimal energy termination of quadratic loop.
    
    Unmixing: Feldspar
    Add: Water
    Quad (000) norm:  1.7702750257848e-01 Lin (018) step:  1.5465833493494e-01 func: -1.6024090882168e+06
    Quad (001) norm:  1.3337217702033e-02 Lin (020) step:  1.4524227518911e+00 func: -1.6024112814304e+06
    Quad (002) norm:  1.0349481938031e-02 Lin (011) step:  1.1340244203694e+00 func: -1.6024114495889e+06
    Quad (003) norm:  4.3792689077607e-03 Lin (019) step:  1.0449740135376e+00 func: -1.6024114563752e+06
    Quad (004) norm:  3.2171612864058e-05 Lin (030) step:  9.9362442957517e-01 func: -1.6024114563223e+06
    Quad (005) norm:  1.3238792828513e-07 Lin (037) step:  8.0650453135374e-01 func: -1.6024114563242e+06
    Minimal energy termination of quadratic loop.
    
     
    T =     756.85 °C, P =      175.0 MPa
    Liquid          moles:   0.527063 grams:  33.949
               SiO2 form:  SiO2           X:  0.6693  wt%    SiO2   72.86
               TiO2 form:  TiO2           X:  0.0019  wt%    TiO2    0.24
              Al2O3 form:  Al2O3          X:  0.0403  wt%   Al2O3   11.32
              Fe2O3 form:  Fe2O3          X:  0.0025  wt%   Fe2O3    0.61
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    1.39
            Fe2SiO4 form:  Fe2SiO4        X:  0.0062  wt%     MgO    0.09
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.23
            Mg2SiO4 form:  Mg2SiO4        X:  0.0007  wt%    Na2O    3.82
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.57
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.87
             CaSiO3 form:  CaSiO3         X:  0.0027  
            Na2SiO3 form:  Na2SiO3        X:  0.0397  
            KAlSiO4 form:  KAlSiO4        X:  0.0624  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1743  
    Feldspar        moles:   0.120460 grams:  32.664
             albite form:  NaAlSi3O8      X:  0.4452  wt%    SiO2   66.13
          anorthite form:  CaAl2Si2O8     X:  0.0155  wt%   Al2O3   19.09
           sanidine form:  KAlSi3O8       X:  0.5393  wt%     CaO    0.32
                                                      wt%    Na2O    5.09
                                                      wt%     K2O    9.37
    Quartz          moles:   0.393825 grams:  23.663
    Water           moles:   0.000010 grams:   0.000
    Feldspar        moles:   0.043077 grams:  11.459
             albite form:  NaAlSi3O8      X:  0.7645  wt%    SiO2   65.46
          anorthite form:  CaAl2Si2O8     X:  0.1021  wt%   Al2O3   21.12
           sanidine form:  KAlSi3O8       X:  0.1333  wt%     CaO    2.15
                                                      wt%    Na2O    8.91
                                                      wt%     K2O    2.36


