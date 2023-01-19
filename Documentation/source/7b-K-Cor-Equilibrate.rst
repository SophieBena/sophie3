Khorzhinskii potential minimization (T, P, :math:`\mu`\ Al2O3 constrained)
==========================================================================

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
    Corundum = modelDB.get_phase('Crn')

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

Function to constrain the chemical potential of Al2O3
-----------------------------------------------------

.. code:: ipython3

    def muAl2O3(t, p, state):
        return Corundum.gibbs_energy(t, p)

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[({'Al':2.0,'O':3.0},muAl2O3)])

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    
    Add: Feldspar
    Quad (000) norm:  3.7831240550747e+01 Lin (037) step:  1.0150575567495e-03 func: -1.6258370888059e+06
    Quad (001) norm:  1.6441201922586e+01 Lin (038) step:  1.7190088819108e-03 func: -1.6347976023419e+06
    Quad (002) norm:  9.2516008477522e+01 Lin (037) step: -6.6597925418252e-04 func: -1.6539664291620e+06
    Quad (003) norm:  2.3793134941153e+01 Lin (037) step: -1.0307550195213e-02 func: -1.7278858325724e+06
    Quad (004) norm:  2.4528698592226e+00 Lin (042) step:  3.9933561609542e-03 func: -1.7304635450125e+06
    Quad (005) norm:  1.3938322106641e-01 Lin (047) step:  4.9300693120074e-03 func: -1.7306369722970e+06
    Quad (006) norm:  4.2522083006141e-03 Lin (049) step:  4.9300692965358e-03 func: -1.7306399484230e+06
    Quad (007) norm:  2.1131561990200e-03 Lin (049) step:  4.9300692965358e-03 func: -1.7306400855138e+06
    Quad (008) norm:  2.0124952443650e-03 Lin (049) step:  4.9300692965358e-03 func: -1.7306400954972e+06
    Quad (009) norm:  1.9943095900946e-03 Lin (048) step:  4.9300681410493e-03 func: -1.7306400978945e+06
    Minimal energy termination of quadratic loop.
    
    Add: Water
    Quad (000) norm:  2.1917771332307e-01 Lin (037) step:  5.4778548153093e-03 func: -1.7306419636378e+06
    Quad (001) norm:  2.7758354142296e-03 Lin (048) step:  5.4778536212273e-03 func: -1.7306419648475e+06
    Quad (002) norm:  4.6959561107417e-12
     
    T =     776.85 °C, P =      175.0 MPa
    Liquid          moles:   1.337027 grams:  84.513
               SiO2 form:  SiO2           X:  0.5705  wt%    SiO2   61.17
               TiO2 form:  TiO2           X:  0.0007  wt%    TiO2    0.09
              Al2O3 form:  Al2O3          X:  0.1271  wt%   Al2O3   26.17
              Fe2O3 form:  Fe2O3          X:  0.0010  wt%   Fe2O3    0.24
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.56
            Fe2SiO4 form:  Fe2SiO4        X:  0.0025  wt%     MgO    0.04
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.00
            Mg2SiO4 form:  Mg2SiO4        X:  0.0003  wt%    Na2O    0.00
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    5.24
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    6.49
             CaSiO3 form:  CaSiO3         X:  0.0000  
            Na2SiO3 form:  Na2SiO3        X:  0.0000  
            KAlSiO4 form:  KAlSiO4        X:  0.0703  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.2277  
    Feldspar        moles:   0.145703 grams:  38.484
             albite form:  NaAlSi3O8      X:  0.8815  wt%    SiO2   67.05
          anorthite form:  CaAl2Si2O8     X:  0.0526  wt%   Al2O3   20.32
           sanidine form:  KAlSi3O8       X:  0.0659  wt%     CaO    1.12
                                                      wt%    Na2O   10.34
                                                      wt%     K2O    1.18
    Water           moles:   0.000850 grams:   0.015
    Quartz          affn:    1971.64

