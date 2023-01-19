Korzhinskii potential minimization (T, P, :math:`\mu`\ SiO2, :math:`\mu`\ Al2O3 constrained)
============================================================================================

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
    phs_sys = [Liquid, Feldspar, Water]

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

Function to constrain the chemical potential of Al2O3
-----------------------------------------------------

.. code:: ipython3

    def muAl2O3(t, p, state):
        return Corundum.gibbs_energy(t, p)

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, lagrange_l=[({'Si':1.0,'O':2.0},muSiO2), ({'Al':2.0,'O':3.0},muAl2O3)])

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    
    Add: Feldspar
    Quad (000) norm:  6.0392738026415e-01 Lin (037) step: -1.1629947224609e-01 func: -5.8434811137901e+05
    Quad (001) norm:  7.7699356421318e-01 Lin (037) step: -3.5940204979373e-03 func: -5.8503771603871e+05
    Quad (002) norm:  1.1552906301478e+00 Lin (039) step:  1.9100098532273e-03 func: -5.8673835221826e+05
    Quad (003) norm:  7.7872227963183e-01 Lin (037) step: -2.9111566001470e-03 func: -5.8939279139700e+05
    Quad (004) norm:  7.7328232200283e-01 Lin (037) step: -2.9111566001470e-03 func: -5.9171694455917e+05
    Quad (005) norm:  7.7680107941554e-01 Lin (037) step: -3.2346184481436e-03 func: -5.9325181154855e+05
    Quad (006) norm:  9.7924329200846e-01 Lin (039) step:  1.9100098623831e-03 func: -5.9401976639669e+05
    Quad (007) norm:  7.9165252411658e-01 Lin (037) step: -2.9111566001470e-03 func: -5.9688181278707e+05
    Quad (008) norm:  7.8737701952675e-01 Lin (037) step: -2.9111566001470e-03 func: -5.9961030647662e+05
    Quad (009) norm:  7.7838799160886e-01 Lin (037) step: -2.9111566001470e-03 func: -6.0198562747342e+05
    Quad (010) norm:  7.8164184093681e-01 Lin (037) step: -3.2346184481436e-03 func: -6.0362241373666e+05
    Quad (011) norm:  9.7078321263246e-01 Lin (039) step:  1.9100098669520e-03 func: -6.0426147218625e+05
    Quad (012) norm:  7.9707403813258e-01 Lin (037) step: -2.9111566001470e-03 func: -6.0716080309314e+05
    Quad (013) norm:  7.9459493707071e-01 Lin (037) step: -2.9111566001470e-03 func: -6.0994524233385e+05
    Quad (014) norm:  7.9564981336215e-01 Lin (037) step: -2.9111566001470e-03 func: -6.1244698091678e+05
    Quad (015) norm:  3.5855309006490e-18
     
    T =     776.85 °C, P =      175.0 MPa
    Liquid          moles:   2.671691 grams: 176.084
               SiO2 form:  SiO2           X:  0.6968  wt%    SiO2   68.96
               TiO2 form:  TiO2           X:  0.0004  wt%    TiO2    0.05
              Al2O3 form:  Al2O3          X:  0.1284  wt%   Al2O3   22.86
              Fe2O3 form:  Fe2O3          X:  0.0005  wt%   Fe2O3    0.12
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.27
            Fe2SiO4 form:  Fe2SiO4        X:  0.0012  wt%     MgO    0.02
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.00
            Mg2SiO4 form:  Mg2SiO4        X:  0.0001  wt%    Na2O    1.84
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    2.77
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    3.12
             CaSiO3 form:  CaSiO3         X:  0.0000  
            Na2SiO3 form:  Na2SiO3        X:  0.0195  
            KAlSiO4 form:  KAlSiO4        X:  0.0387  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1143  
    Feldspar        moles:   0.031924 grams:   8.497
             albite form:  NaAlSi3O8      X:  0.7538  wt%    SiO2   62.30
          anorthite form:  CaAl2Si2O8     X:  0.2402  wt%   Al2O3   23.75
           sanidine form:  KAlSi3O8       X:  0.0060  wt%     CaO    5.06
                                                      wt%    Na2O    8.78
                                                      wt%     K2O    0.11
    Water           affn:    4755.55


Pickup runs use previously computed state
