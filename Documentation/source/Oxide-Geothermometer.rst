Fe-Ti oxide geothermobarometer
==============================

Constructed based on Fe-Ti exchange for oxide pairs, using the methodology of Ghiorso and Evans (2008)
------------------------------------------------------------------------------------------------------

Required Python code to load the phase library …

.. code:: ipython3

    import numpy as np
    from scipy import optimize as optim
    import thermoengine as thermo

Get access to a thermodynamic database (by default, the Berman (1988)/MELTS database).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    modelDB = thermo.model.Database()

Create a Python reference to the Spinel (“Mag”) and Rhombehedral oxide (“Ilm”) solution phase class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    Mag = modelDB.get_phase('SplS')
    Ilm = modelDB.get_phase('Rhom')

Optional - Obtain some properties of the selected Oxide solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Name, formulas of endmembers, molecular weights of endmembers,
abbreviation, number of endmember components, names of endmembers

.. code:: ipython3

    print (Mag.props['phase_name'])
    print ('Num of endmembers: ', Mag.props['endmember_num'])
    print (Mag.props['formula'])
    print (Mag.props['endmember_name'])
    print ()
    print (Ilm.props['phase_name'])
    print ('Num of endmembers: ', Ilm.props['endmember_num'])
    print (Ilm.props['formula'])
    print (Ilm.props['endmember_name'])


.. parsed-literal::

    Spinel
    Num of endmembers:  5
    ['FeCr2O4' 'FeAl2O4' 'Fe3O4' 'MgAl2O4' 'Fe2TiO4']
    ['chromite' 'hercynite' 'magnetite' 'spinel' 'ulvospinel']
    
    Ilmenite ss
    Num of endmembers:  5
    ['MgTiO3' 'Fe2O3' 'FeTiO3' 'MnTiO3' 'Al2O3']
    ['geikielite' 'hematite' 'ilmenite' 'pyrophanite' 'corundum']


Step 1 - Input compositions of oxide pairs
------------------------------------------

.. code:: ipython3

    case = 1

-  case 0: Brad’s composition
-  case 1: Oxide pair from website geothermometer

.. code:: ipython3

    if case == 0:
        Mag_mol_oxides = thermo.chem.format_mol_oxide_comp(
            {'SiO2':0.0, 'TiO2':4.15, 'Al2O3':2.83, 'Fe2O3':57.89, 'FeO':32.48, 
            'MnO':0.35, 'MgO':1.0, 'CaO':0.0, 'Na2O':0.0, 'Cr2O3':0.07},convert_grams_to_moles=True)
        Ilm_mol_oxides = thermo.chem.format_mol_oxide_comp(
            {'SiO2':0, 'TiO2':37.27, 'Al2O3':0.29, 'Fe2O3':29.62,'FeO':29.19, 
            'MnO':0.39, 'MgO':2.0, 'CaO':0.0, 'Na2O':0.0, 'Cr2O3':0.01},convert_grams_to_moles=True)
    else: # example from website, gives 724 °C, 1.61 delta NNO 0.894 a TiO2
        Mag_mol_oxides = thermo.chem.format_mol_oxide_comp(
            {'SiO2':0.0, 'TiO2':4.35, 'Al2O3':1.94, 'Fe2O3':0.00, 'Cr2O3':0.18, 'FeO':86.34, 
            'MnO':0.44, 'MgO':1.2, 'CaO':0.0, 'Na2O':0.0},convert_grams_to_moles=True)
        Ilm_mol_oxides = thermo.chem.format_mol_oxide_comp(
            {'SiO2':0.0, 'TiO2':28.73, 'Al2O3':0.35, 'Fe2O3':0.00, 'Cr2O3':0.0, 'FeO':65.98, 
            'MnO':0.23, 'MgO':1.02, 'CaO':0.0, 'Na2O':0.0},convert_grams_to_moles=True)

Step 2 - Convert analytical composition to moles of endmember components
------------------------------------------------------------------------

| Note that a method - *test_endmember_comp()* - is called to test the
  validity of the projected composition
| Also note that the “intrinsic” conversion routines now take care of
  FeO-Fe2O3 conversion to balance the cation-anion ratio of the phase.
  Input compositions of Fe2O3 and/or FeO are adjusted to balance charge
  and maintain phase stoichiometry.

.. code:: ipython3

    def validate_endmember_comp(moles_end, phase):
        print(phase.props['phase_name'])
        sum = 0.0
        for i in range(0,phase.props['endmember_num']):
            print ("mole number of {0:10.10s} = {1:13.8f}".format(
                phase.props['endmember_name'][i], moles_end[i]))
            sum += moles_end[i]
        if not phase.test_endmember_comp(moles_end):
            print ("Calculated composition is infeasible!")
            
    Mag_moles_end = Mag.calc_endmember_comp(
        mol_oxide_comp=Mag_mol_oxides, method='intrinsic', normalize=True)
    validate_endmember_comp(Mag_moles_end, Mag)
    print()
    Ilm_moles_end = Ilm.calc_endmember_comp(
        mol_oxide_comp=Ilm_mol_oxides, method='intrinsic', normalize=True)
    validate_endmember_comp(Ilm_moles_end, Ilm)


.. parsed-literal::

    Spinel
    mole number of chromite   =    0.00266617
    mole number of hercynite  =   -0.02419364
    mole number of magnetite  =    0.83193038
    mole number of spinel     =    0.06702845
    mole number of ulvospinel =    0.12256865
    
    Ilmenite ss
    mole number of geikielite =    0.03853892
    mole number of hematite   =    0.44719307
    mole number of ilmenite   =    0.50410315
    mole number of pyrophanit =    0.00493747
    mole number of corundum   =    0.00522739


Implement a Fe-Ti oxide geothermometer.
=======================================

Consider Fe-Ti exchange between oxides
--------------------------------------

Rhom(Ilm) + Spinel(Mag) = Spinel (Ulv) + Rhom(Hm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FeTiO3 + Fe3O4 = Fe2TiO4 + Fe2O3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fix the pressure at 2000 bars

.. code:: ipython3

    P = 2000.0

Corection terms from Ghiorso and Evans (2008) that modify the MELTS models
--------------------------------------------------------------------------

Correction terms for ulvospinel derived in Ghiorso and Evans (2008)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def UlvCorr(t, correctReaction=True):
        tr = 298.15
        h = - 162.0 + 284.5
        s = 0.0
        if correctReaction:
            h += 2039.106175 
            s +=    1.247790
        l1 = - 0.039452*np.sqrt(4.184)
        l2 = 7.54197e-5*np.sqrt(4.184)
        h = h + 0.5*l1*l1*(t*t-tr*tr) + (2.0/3.0)*l1*l2*(t*t*t - tr*tr*tr) + 0.25*l2*l2*(t*t*t*t - tr*tr*tr*tr)
        s = s + l1*l1*(t - tr) + l1*l2*(t*t - tr*tr) + (1.0/3.0)*l2*l2*(t*t*t - tr*tr*tr)
        return h - t*s

Ghiorso and Evans (2008) used the Vinet integral; MELTS uses the Berman integral
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We must substract the latter from computed chemical potentials and add
in the former.

.. code:: ipython3

    def BermanVint(t, p, v0, v1, v2, v3, v4):
        pr = 1.0
        tr = 298.15
        return v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*(t-tr)*(t-tr))*(p-pr))
    def VinetVint(t, p, v0, alpha, K, Kp):
        eta = 3.0*(Kp-1.0)/2.0
        x   = 1.0
        x0  = 1.0
        pr  = 1.0
        tr  = 298.15
        
        iter = 0
        while True:
            fn = x*x*(p/10000.0) - 3.0*K*(1.0-x)*np.exp(eta*(1.0-x)) - x*x*alpha*K*(t-tr)
            dfn = 2.0*x*(p/10000.0) + 3.0*K*(1.0+eta*(1.0-x))*np.exp(eta*(1.0-x)) - 2.0*alpha*K*(t-tr)
            x = x - fn/dfn
            iter += 1
            if ((iter > 500) or (fn*fn < 1.0e-15)):
                break
        # print (iter, x)
        
        iter = 0
        while True:
            fn = x0*x0*(pr/10000.0) - 3.0*K*(1.0-x0)*np.exp(eta*(1.0-x0)) - x0*x0*alpha*K*(t-tr)
            dfn = 2.0*x0*(pr/10000.0) + 3.0*K*(1.0+eta*(1.0-x0))*np.exp(eta*(1.0-x0)) - 2.0*alpha*K*(t-tr)
            x0 = x0 - fn/dfn
            iter += 1
            if ((iter > 500) or (fn*fn < 1.0e-15)):
                break
        # print (iter, x0)
        
        a  = (9.0*v0*K/(eta*eta))*(1.0 - eta*(1.0-x))*np.exp(eta*(1.0-x))
        a += v0*(t-tr)*K*alpha*(x*x*x - 1.0) - 9.0*v0*K/(eta*eta)
        a -= (9.0*v0*K/(eta*eta))*(1.0 - eta*(1.0-x0))*np.exp(eta*(1.0-x0))
        a -= v0*(t-tr)*K*alpha*(x0*x0*x0 - 1.0) - 9.0*v0*K/(eta*eta)
        
        return -a*10000.0 + p*v0*x*x*x - pr*v0*x0*x0*x0

Berman integral for the reaction FeTiO3 + Fe3O4 = Fe2TiO4 + Fe2O3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def rBerVint(T, P):
        vIntBerMag = BermanVint(T, P, 4.452, -0.582E-6, 1.751E-12, 30.291E-6, 138.470E-10)
        vIntBerUlv = BermanVint(T, P, 4.682, 0.0, 0.0, 0.0, 0.0)
        vIntBerHem = BermanVint(T, P, 3.027, -0.479e-6, 0.304e-12, 38.310e-6, 1.650e-10)
        vIntBerIlm = BermanVint(T, P, 3.170, -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10)
        return vIntBerUlv + vIntBerHem - vIntBerMag -vIntBerIlm

Vinet integral for the reaction FeTiO3 + Fe3O4 = Fe2TiO4 + Fe2O3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def rVinetVint(T, P):
        vIntVinetMag = VinetVint(T, P, 4.452, 30.291E-6, 171.821, 9.3387)
        vIntVinetUlv = VinetVint(T, P, 4.682, 30.291E-6, 171.821, 9.3387)
        vIntVinetHem = VinetVint(T, P, 3.027, 38.310E-6, 208.768, 1.64992)
        vIntVinetIlm = VinetVint(T, P, 3.170, 27.248E-6, 171.233, 6.21289)
        return vIntVinetUlv + vIntVinetHem - vIntVinetMag - vIntVinetIlm

This method computes the free energy of the exchange reaction …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def deltaG(T, P, mag_mols, ilm_mols):
        muMag = Mag.chem_potential(T, P, mol=mag_mols, endmember=2)
        muUlv = Mag.chem_potential(T, P, mol=mag_mols, endmember=4) + UlvCorr(T)
        muIlm = Ilm.chem_potential(T, P, mol=ilm_mols, endmember=2)
        muHem = Ilm.chem_potential(T, P, mol=ilm_mols, endmember=1)
        deltaG = muUlv + muHem - muIlm - muMag - rBerVint(T, P) + rVinetVint(T, P)
        return deltaG

This next function is used by the minimizer to zero the free energy of the exchange reaction …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def boundary(P, Tlims, deltaG, mag_mols, ilm_mols):
        Afun = lambda T, P=P: deltaG(T, P, mag_mols, ilm_mols)
        Tbound = optim.brentq(Afun, Tlims[0], Tlims[1])
        return Tbound

Calculate the equilibrium temperature for this oxide pair …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    Teq = boundary(P, [500.,2000.], deltaG, Mag_moles_end, Ilm_moles_end)
    print('Equilibrium Temp = ', Teq-273.15, ' °C')


.. parsed-literal::

    Equilibrium Temp =  724.072806179801  °C


Calculate Oxygen fugacity from the reaction
-------------------------------------------

O2 + 4 Fe3O4 = 6 Fe2O3
~~~~~~~~~~~~~~~~~~~~~~

| Note that the properties of oxygen are defined here for consistency
  instead of using the built-in functions.
| Also note that the chemical potentials of hematite and magnetite are
  adjusted to remove the Berman-type volume integrals and replace them
  with the Vinet-type volume integrals to be consistent with Ghiorso and
  Evans (2008)

.. code:: ipython3

    def muO2(t, p):
        tr = 298.15
        hs = 23.10248*(t-tr) + 2.0*804.8876*(np.sqrt(t)-np.sqrt(tr)) - 1762835.0*(1.0/t-1.0/tr) \
           - 18172.91960*np.log(t/tr) + 0.5*0.002676*(t*t-tr*tr)
        ss = 205.15 + 23.10248*np.log(t/tr)  - 2.0*804.8876*(1.0/np.sqrt(t)-1.0/np.sqrt(tr)) \
           - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr)) + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr)
        return hs - t*ss
    def deltaNNO (T, P, mag_mols, ilm_mols):
        muHem  = Ilm.chem_potential(T, P, mol=ilm_mols, endmember=1)
        muHem -= BermanVint(T, P, 3.027, -0.479e-6, 0.304e-12, 38.310e-6, 1.650e-10)
        muHem += VinetVint(T, P, 3.027, 38.310E-6, 208.768, 1.64992)
        muMag =  Mag.chem_potential(T, P, mol=mag_mols, endmember=2)
        muMag -= BermanVint(T, P, 4.452, -0.582E-6, 1.751E-12, 30.291E-6, 138.470E-10)
        muMag += VinetVint(T, P, 4.452, 30.291E-6, 171.821, 9.3387)
        muOxy  = muO2(T, P)
        logfO2 = (6.0*muHem - 4.0*muMag -  muOxy)/(8.3144598*T)/np.log(10.0)
        return logfO2 - (-25018.7/T + 12.981 + 0.046*(P-1.0)/Teq -0.5117*np.log(T))

Calculate the equilibrium oxygen fugacity for this oxide pair …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print(deltaNNO(Teq, P, Mag_moles_end, Ilm_moles_end))


.. parsed-literal::

    1.613240788544294


Calculate the temperature for Fe-Mg exchange
--------------------------------------------

FeAl2O4 (hercynite) + MgTiO3 (geikielite) = MgAl2O4 (spinel) + FeTiO3 (ilmenite)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The method below is used by the minimizer to evaluate the free energy
change of the Fe-Mg exchange reaction …

.. code:: ipython3

    def deltaGfemg(T, P, mag_mols, ilm_mols):
        muSpn  = Mag.chem_potential(T, P, mol=mag_mols, endmember=3)
        muSpn -= BermanVint(T, P, 3.977, -0.489E-6, 0.0, 21.691E-6, 50.528E-10)
        muSpn += VinetVint(T, P, 3.977, 21.691E-6, 204.499, 4.0)
        
        muHer  = Mag.chem_potential(T, P, mol=mag_mols, endmember=1)
        muHer -= BermanVint(T, P, 0.973948*4.184, 0.0, 0.0, 0.0, 0.0)
        muHer += VinetVint(T, P, 0.973948*4.184, 21.691E-6, 204.499, 4.0)
        
        muIlm  = Ilm.chem_potential(T, P, mol=ilm_mols, endmember=2)
        muIlm -= BermanVint(T, P, 3.170, -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10)
        muIlm += VinetVint(T, P, 3.170, 27.248E-6, 171.233, 6.21289)
        
        muGei  = Ilm.chem_potential(T, P, mol=ilm_mols, endmember=0)
        muGei -= BermanVint(T, P, 3.086, -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10)
        muGei += VinetVint(T, P, 3.086, 27.2476341e-6, 171.240, 6.21527)
        
        deltaG = muSpn + muIlm - muHer - muGei
        return deltaG

Compute the Fe-Mg exchange temperature (if possible) …

.. code:: ipython3

    Tlow  = deltaGfemg(500.0, P, Mag_moles_end, Ilm_moles_end)
    Thigh = deltaGfemg(2000.0, P, Mag_moles_end, Ilm_moles_end)
    if np.sign(Tlow) != np.sign(Thigh):
        Tfemg = boundary(P, [500.,2000.], deltaGfemg, Mag_moles_end, Ilm_moles_end)
        print('Fe-Mg Equilibrium Temp = ', Tfemg-273.15, ' °C')
    else:
        print('No Fe-Mg equilibration temperature found.')


.. parsed-literal::

    No Fe-Mg equilibration temperature found.


Calculate the activity of TiO2 relative to rutile saturation
------------------------------------------------------------

.. code:: ipython3

    Rut = modelDB.get_phase('Rt')
    def aTiO2(T, P, mag_mols, ilm_mols):
        muUlv  = Mag.chem_potential(T, P, mol=mag_mols, endmember=4) + UlvCorr(T, correctReaction=False)
        muUlv -= BermanVint(T, P, 4.682, 0.0, 0.0, 0.0, 0.0)
        muUlv += VinetVint(T, P, 4.682, 30.291E-6, 171.821, 9.3387)
        muIlm  = Ilm.chem_potential(T, P, mol=ilm_mols, endmember=2)
        muIlm -= BermanVint(T, P, 3.170, -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10)
        muIlm += VinetVint(T, P, 3.170, 27.248E-6, 171.233, 6.21289)
        muRut = Rut.chem_potential(T, P)
        return np.exp(-(muRut+muUlv-2.0*muIlm)/(8.3143*T))

.. code:: ipython3

    aTiO2(Teq, P, Mag_moles_end, Ilm_moles_end)




.. parsed-literal::

    0.8938677107462636



