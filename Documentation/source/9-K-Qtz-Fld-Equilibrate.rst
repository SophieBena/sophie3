Korzhinskii potential minimization (T, P, :math:`\mu`\ SiO2, :math:`\mu`\ Albite, :math:`\mu`\ Anorthite, :math:`\mu`\ Sanidine constrained)
============================================================================================================================================

.. code:: ipython3

    import numpy as np
    import scipy.optimize as opt
    import scipy.linalg as lin 
    import scipy as sci
    import sys

.. code:: ipython3

    from thermoengine import core, phases, model, equilibrate

.. code:: ipython3

    #np.set_printoptions(linewidth=200, precision=1)
    np.set_printoptions(linewidth=200, precision=2)

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

Function to constrain the chemical potential of silica
------------------------------------------------------

.. code:: ipython3

    def muSiO2(t, p, state):
        return Quartz.gibbs_energy(t, p)

Function to constrain the chemical potential of feldspar components
-------------------------------------------------------------------

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    moles_fld = np.array([0.4388, 0.0104, 0.5508])
    muFld = Feldspar.gibbs_energy(t,p,mol=moles_fld,deriv={"dmol":1})
    muFld[0][0],muFld[0][1],muFld[0][2]




.. parsed-literal::

    (-4286486.8312951205, -4587946.630682982, -4325071.949624277)



.. code:: ipython3

    def muAb(t, p, state):
        return muFld[0][0]
    def muAn(t, p, state):
        return muFld[0][1]
    def muSn(t, p, state):
        return muFld[0][2]

Instantiate class instance and run calculation
----------------------------------------------

.. code:: ipython3

    equil = equilibrate.Equilibrate(elm_sys, phs_sys, 
                                    lagrange_l=[({'Na':1.0,'Al':1.0,'Si':3.0,'O':8.0},muAb), 
                                                ({'Ca':1.0,'Al':2.0,'Si':2.0,'O':8.0},muAn),
                                                ({ 'K':1.0,'Al':1.0,'Si':3.0,'O':8.0},muSn),
                                                ({'Si':1.0, 'O':2.0},muSiO2)])

.. code:: ipython3

    t = 1050.0
    p = 1750.0
    state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
    state.print_state()


.. parsed-literal::

    
    Add: Quartz
    Quad (000) norm:  0.0000000000000e+00
     
    T =     776.85 °C, P =      175.0 MPa
    Liquid          moles:   1.932513 grams: 124.994
               SiO2 form:  SiO2           X:  0.6898  wt%    SiO2   74.39
               TiO2 form:  TiO2           X:  0.0005  wt%    TiO2    0.06
              Al2O3 form:  Al2O3          X:  0.0400  wt%   Al2O3   11.71
              Fe2O3 form:  Fe2O3          X:  0.0007  wt%   Fe2O3    0.17
            MgCr2O4 form:  MgCr2O4        X:  0.0000  wt%     FeO    0.38
            Fe2SiO4 form:  Fe2SiO4        X:  0.0017  wt%     MgO    0.02
          MnSi0.5O2 form:  MnSi0.5O2      X:  0.0000  wt%     CaO    0.17
            Mg2SiO4 form:  Mg2SiO4        X:  0.0002  wt%    Na2O    3.70
          NiSi0.5O2 form:  NiSi0.5O2      X:  0.0000  wt%     K2O    4.99
          CoSi0.5O2 form:  CoSi0.5O2      X:  0.0000  wt%     H2O    4.40
             CaSiO3 form:  CaSiO3         X:  0.0019  
            Na2SiO3 form:  Na2SiO3        X:  0.0386  
            KAlSiO4 form:  KAlSiO4        X:  0.0686  
          Ca3(PO4)2 form:  Ca3(PO4)2      X:  0.0000  
                H2O form:  H2O            X:  0.1580  
    Feldspar        affn:       0.12
             albite form:  NaAlSi3O8      X:  0.4429
          anorthite form:  CaAl2Si2O8     X:  0.0107
           sanidine form:  KAlSi3O8       X:  0.5464
    Water           affn:     756.37
    Quartz          moles:   0.000010 grams:   0.001


Pickup runs use previously computed state
