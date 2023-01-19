Feldspar Solid Solution - Berman Compatible
===========================================

Required Python code to load the phase library.

.. code:: ipython3

    from thermoengine import core, phases, model

Listing of solution phases compatible with Berman inherited from the MELTS model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_phase_info returns an ordered list of *pure* and *solution* phases

.. code:: ipython3

    phase_info,info_files = phases.get_phase_info()
    phase_info['solution']




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Abbrev</th>
          <th>Name</th>
          <th>Members</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>AkiS</td>
          <td>Akimotoite Solution</td>
          <td>MgAki/FeAki/AlAki</td>
        </tr>
        <tr>
          <th>1</th>
          <td>Bt</td>
          <td>Biotite</td>
          <td>Ann/Phl</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CfS</td>
          <td>Ca-Ferrite Solution</td>
          <td>MgCf/FeCf/NaCf</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CalS</td>
          <td>Calcite Solution</td>
          <td>Cal/Mgs</td>
        </tr>
        <tr>
          <th>4</th>
          <td>Cam</td>
          <td>Clinoamphibole</td>
          <td>Cum/Gru/Tr</td>
        </tr>
        <tr>
          <th>5</th>
          <td>Fp</td>
          <td>Ferropericlase</td>
          <td>Per/Wus</td>
        </tr>
        <tr>
          <th>6</th>
          <td>Mws</td>
          <td>Magnesiowustite</td>
          <td>Per/Wus</td>
        </tr>
        <tr>
          <th>7</th>
          <td>Fsp</td>
          <td>Feldspar</td>
          <td>Ab/An/Kfs/Or/Mc/Sa</td>
        </tr>
        <tr>
          <th>8</th>
          <td>Grt</td>
          <td>Garnet</td>
          <td>Alm/Grs/Prp/Maj/NaMaj</td>
        </tr>
        <tr>
          <th>9</th>
          <td>Hbl</td>
          <td>Hornblende</td>
          <td>Fprg/Mhs/Prg</td>
        </tr>
        <tr>
          <th>10</th>
          <td>KlsS</td>
          <td>Kalsilite Solution</td>
          <td>Nph/Kls</td>
        </tr>
        <tr>
          <th>11</th>
          <td>LctS</td>
          <td>Leucite Solution</td>
          <td>Anl/Lct/nLct</td>
        </tr>
        <tr>
          <th>12</th>
          <td>Liq</td>
          <td>Liquid</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>13</th>
          <td>Mll</td>
          <td>Melilite</td>
          <td>Ak/Gh/fAk/nMll</td>
        </tr>
        <tr>
          <th>14</th>
          <td>NphS</td>
          <td>Nepheline Solution</td>
          <td>Nph/Kls</td>
        </tr>
        <tr>
          <th>15</th>
          <td>Ol</td>
          <td>Olivine</td>
          <td>CoOl/Fa/Fo/Mtc/NiOl/Tep</td>
        </tr>
        <tr>
          <th>16</th>
          <td>OrOx</td>
          <td>OrthoOxide</td>
          <td>Fpsb/Krt/Psb</td>
        </tr>
        <tr>
          <th>17</th>
          <td>Oam</td>
          <td>Orthoamphibole</td>
          <td>Cum/Gru/Tr</td>
        </tr>
        <tr>
          <th>18</th>
          <td>Opx</td>
          <td>Orthopyroxene</td>
          <td>Abfn/Bfn/CaTs/En/Fs/cEn/Di/Ess/Hd/Jd/oDi</td>
        </tr>
        <tr>
          <th>19</th>
          <td>Cpx</td>
          <td>Clinopyroxene</td>
          <td>Abfn/Bfn/CaTs/En/Fs/cEn/Di/Ess/Hd/Jd</td>
        </tr>
        <tr>
          <th>20</th>
          <td>hpCpx</td>
          <td>High-Pressure Clinopyroxene</td>
          <td>hpcEn/hpcFs</td>
        </tr>
        <tr>
          <th>21</th>
          <td>PrvS</td>
          <td>Perovskite</td>
          <td>Prv/FePrv/AlPrv</td>
        </tr>
        <tr>
          <th>22</th>
          <td>PpvS</td>
          <td>Post-Perovskite</td>
          <td>MgPpv/FePpv/AlPpv</td>
        </tr>
        <tr>
          <th>23</th>
          <td>Rhom</td>
          <td>Rhombohedral</td>
          <td>Crn/Gk/Hem/Ilm/Pph</td>
        </tr>
        <tr>
          <th>24</th>
          <td>Rwd</td>
          <td>Ringwoodite</td>
          <td>MgRwd/FeRwd</td>
        </tr>
        <tr>
          <th>25</th>
          <td>SplS</td>
          <td>Spinel</td>
          <td>Chr/Hc/Mag/Spl/Usp</td>
        </tr>
        <tr>
          <th>26</th>
          <td>Wds</td>
          <td>Wadsleyite</td>
          <td>MgWds/FeWds</td>
        </tr>
        <tr>
          <th>27</th>
          <td>Psu</td>
          <td>PseudoPhase</td>
          <td>NaN</td>
        </tr>
      </tbody>
    </table>
    </div>



Get access to a thermodynamic database (by default, the Berman (1988)/MELTS database).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    modelDB = model.Database()

Create a Python reference to the feldspar solution phase class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The required phase abbreviation method argument is listed above.

.. code:: ipython3

    Feldspar = modelDB.get_phase('Fsp')

Obtain some properties of the selected feldspar solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Name, formulas of endmembers, molecular weights of endmembers,
abbreviation, number of endmember components, names of endmembers

.. code:: ipython3

    print (Feldspar.props['phase_name'])
    print (Feldspar.props['formula'])
    print (Feldspar.props['molwt'])
    print (Feldspar.props['abbrev'])
    print (Feldspar.props['endmember_num'])
    print (Feldspar.props['endmember_name'])


.. parsed-literal::

    Feldspar
    ['NaAlSi3O8' 'CaAl2Si2O8' 'KAlSi3O8']
    [262.22301 278.20928 278.33524]
    Fsp
    3
    ['albite' 'anorthite' 'sanidine']


Specify the composition of the solution
---------------------------------------

Composition input in wt% oxides (treated as grams of oxides)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sandine composition in wt% oxides, taken from Deer, Howie and Zussman,
Table 30, entry 9: Sanidine from a rhyolite in Mitchell Mesa, Texas.

.. code:: ipython3

    mol_oxides = core.chem.format_mol_oxide_comp({'SiO2':67.27,'Al2O3':18.35, 'FeO':0.92, 'CaO':0.15, 
                                                  'Na2O':6.45, 'K2O':7.05, 'H2O':0.16}, convert_grams_to_moles=True)

Convert analytical composition to moles of endmember components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that a method is called - *test_endmember_comp()* - to test the
validity of the projected composition #### (1st method) A default
projection method based on generic least sqaures …

.. code:: ipython3

    moles_end,oxide_res = Feldspar.calc_endmember_comp(mol_oxide_comp=mol_oxides, output_residual=True)
    for i in range(0,Feldspar.props['endmember_num']):
        print ("mole number of {0:10.10s} = {1:13.6e}".format(Feldspar.props['endmember_name'][i], moles_end[i]))
    if not Feldspar.test_endmember_comp(moles_end):
        print ("Calculated composition is infeasible!")


.. parsed-literal::

    mole number of albite     =  2.156163e-01
    mole number of anorthite  = -2.161623e-04
    mole number of sanidine   =  1.571699e-01
    Calculated composition is infeasible!


(2nd method) An intrinsic projection method tailored specifically to the solution phase of interest …
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    moles_end,oxide_res = Feldspar.calc_endmember_comp(mol_oxide_comp=mol_oxides, method='intrinsic', output_residual=True)
    for i in range(0,Feldspar.props['endmember_num']):
        print ("mole number of {0:10.10s} = {1:13.6e}".format(Feldspar.props['endmember_name'][i], moles_end[i]))
    if not Feldspar.test_endmember_comp(moles_end):
        print ("Calculated composition is infeasible!")


.. parsed-literal::

    mole number of albite     =  2.081352e-01
    mole number of anorthite  =  2.674779e-03
    mole number of sanidine   =  1.496888e-01


Some convenience methods to manipulate composition information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| ➜ Moles of endmember components => Molar sum
| ➜ Moles of endmember components => Moles of elements (standard order)
| ➜ Moles of endmember components => Mole fractions of endmember
  components
| ➜ Moles of elements (standard order) => Moles of endmember components
  of the phase
| ➜ Moles of elements (standard order) => Total moles of endmember
  components of the phase
| ➜ Moles of elements (standard order) => Total mass of the phase (g)

.. code:: ipython3

    print ('Total moles of endmembers: ', Feldspar.covert_endmember_comp(moles_end,output='total_moles'))
    mol_elm = Feldspar.covert_endmember_comp(moles_end,output='moles_elements')
    print ('Mole fractions of endmembers: ', Feldspar.covert_endmember_comp(moles_end,output='mole_fraction'))
    print ('Moles of endmembers: ', Feldspar.convert_elements(mol_elm, output='moles_end'))
    print ('Total moles of endmembers: ', Feldspar.convert_elements(mol_elm, output='total_moles'))
    print ('Total grams of phase: ', Feldspar.convert_elements(mol_elm, output='total_grams'))


.. parsed-literal::

    Total moles of endmembers:  0.36049883224427715
    Mole fractions of endmembers:  [0.57735337 0.00741966 0.41522697]
    Moles of endmembers:  [0.20813521 0.00267478 0.14968884]
    Total moles of endmembers:  0.36049883224427715
    Total grams of phase:  96.98566962270826


Compute activities and chemical potentials of endmember components:
-------------------------------------------------------------------

| ➜ Moles of components, T (K), P (bars) => activities of endmember
  components
| ➜ Moles of components, T (K), P (bars) => chemical potentials of
  endmember components (J)

.. code:: ipython3

    t = 1000
    p = 1000

.. code:: ipython3

    names = Feldspar.props['endmember_name']
    forms = Feldspar.props['formula']
    for i in range(0,Feldspar.props['endmember_num']):
        print ("chemical potential of {0:15.15s} ({1:15.15s}) = {2:15.2f} J, activity = {3:10.6f}".format(
            names[i], forms[i], Feldspar.chem_potential(t, p, mol=moles_end, endmember=i),
            Feldspar.activity(t, p, mol=moles_end, endmember=i)))


.. parsed-literal::

    chemical potential of albite          (NaAlSi3O8      ) =     -4265565.12 J, activity =   0.842485
    chemical potential of anorthite       (CaAl2Si2O8     ) =     -4575045.33 J, activity =   0.053515
    chemical potential of sanidine        (KAlSi3O8       ) =     -4305773.90 J, activity =   0.646776


Molar derivatives of activities:
--------------------------------

➜ Moles of components, T (K), P (bars) => d(activities of endmember
components)/d(Moles of components)

.. code:: ipython3

    for i in range(0,Feldspar.props['endmember_num']):
        dadm = Feldspar.activity(t,p,mol=moles_end,deriv={'dmol':1},endmember=i)[0]
        print ("Derivative of a[{0:s}] with respect to:".format(names[i]))
        for j in range(0,Feldspar.props['endmember_num']):
            print ("    {0:15.15s} = {1:13.6e}".format(names[j], dadm[j]))


.. parsed-literal::

    Derivative of a[albite] with respect to:
        albite          =  2.787764e-01
        anorthite       = -7.297958e+00
        sanidine        = -2.572186e-01
    Derivative of a[anorthite] with respect to:
        albite          = -4.635673e-01
        anorthite       =  1.932289e+01
        sanidine        =  2.992890e-01
    Derivative of a[sanidine] with respect to:
        albite          = -1.974669e-01
        anorthite       =  3.617189e+00
        sanidine        =  2.099331e-01


Gibbs free energy and its compositional derivatives:
----------------------------------------------------

| ➜ Moles of components, T (K), P (bars) => Gibbs free energy (J)
| ➜ Moles of components, T (K), P (bars) => d(Gibbs free energy)/d(Moles
  of components) (J)
| ➜ Moles of components, T (K), P (bars) => d^2(Gibbs free
  energy)/d(Moles of components)^2 (J)
| ➜ Moles of components, T (K), P (bars) => d^3(Gibbs free
  energy)/d(Moles of components)^3 (J)

.. code:: ipython3

    print ('Gibbs free energy = {0:12.2f} J'.format(Feldspar.gibbs_energy(t, p, mol=moles_end)))
    print ("")
    dgdm = Feldspar.gibbs_energy(t, p, mol=moles_end, deriv={'dmol':1})[0]
    d2gdm2 = Feldspar.gibbs_energy(t, p, mol=moles_end, deriv={'dmol':2})[0]
    for i in range (0, Feldspar.props['endmember_num']):
        print ('dg/dm[{0:>2d}] = {1:13.6e}     d2gdm2[{0:>2d}][*]:  '.format(i, dgdm[i], i), end=' ')
        for j in range (0, Feldspar.props['endmember_num']):
            print ('{0:13.6e}'.format(d2gdm2[i, j]), end=' ')
        print()
    print ("")
    d3gdm3 = Feldspar.gibbs_energy(t, p, mol=moles_end, deriv={'dmol':3})[0]
    for i in range (0, Feldspar.props['endmember_num']):
        print('d3gdm3[{0:>2d}][*][*]: '.format(i), end=' ')
        for j in range (0, Feldspar.props['endmember_num']):
            for k in range (0, Feldspar.props['endmember_num']):
                print ('{0:13.6e}'.format(d3gdm3[i, j, k]), end=' ')
            print ('  ', end=' ')
        print ()


.. parsed-literal::

    Gibbs free energy =  -1544577.84 J
    
    dg/dm[ 0] = -4.265565e+06     d2gdm2[ 0][*]:    2.751184e+03 -7.202198e+04 -2.538435e+03 
    dg/dm[ 1] = -4.575045e+06     d2gdm2[ 1][*]:   -7.202198e+04  3.002095e+06  4.649893e+04 
    dg/dm[ 2] = -4.305774e+06     d2gdm2[ 2][*]:   -2.538435e+03  4.649893e+04  2.698688e+03 
    
    d3gdm3[ 0][*][*]:  -3.714742e+03  1.866043e+05 -1.613352e+03     5.612260e+04  4.781953e+05  5.992005e+04     7.182150e+02 -1.215085e+05  3.343639e+03    
    d3gdm3[ 1][*][*]:   3.042045e+05  8.211941e+05  2.308674e+04     6.907124e+05 -1.160760e+09  1.827003e+05     2.541830e+04  1.271789e+03 -3.258032e+05    
    d3gdm3[ 2][*][*]:  -3.714742e+03  1.866043e+05 -1.613352e+03     5.612260e+04  4.781953e+05  5.992005e+04     7.182150e+02 -1.215085e+05  3.343639e+03    


Enthalpy, Entropy, and molar derivatives:
-----------------------------------------

| ➜ Moles of components, T (K), P (bars) => enthalpy (J)
| ➜ Moles of components, T (K), P (bars) => entropy (J/K)
| ➜ Moles of components, T (K), P (bars) => d(entropy)/d(Moles of
  components) (J/K)
| ➜ Moles of components, T (K), P (bars) => d^2(entropy)/d(Moles of
  components)^2 (J/K)

.. code:: ipython3

    print ('Entropy = {0:12.2f} J/K'.format(Feldspar.entropy(t, p, mol=moles_end)))
    print ("")
    dsdm = Feldspar.entropy(t, p, mol=moles_end, deriv={'dmol':1})[0]
    d2sdm2 = Feldspar.entropy(t, p, mol=moles_end, deriv={'dmol':2})[0]
    for i in range (0, Feldspar.props['endmember_num']):
        print ('ds/dm[{0:>2d}] = {1:13.6e}     d2sdm2[{0:>2d}][*]:  '.format(i, dsdm[i], i), end=' ')
        for j in range (0, Feldspar.props['endmember_num']):
            print ('{0:13.6e}'.format(d2sdm2[i, j]), end=' ')
        print()


.. parsed-literal::

    Entropy =       199.56 J/K
    
    ds/dm[ 0] =  5.506230e+02     d2sdm2[ 0][*]:   -2.691158e+01  2.489871e+01  3.697437e+01 
    ds/dm[ 1] =  5.697965e+02     d2sdm2[ 1][*]:    2.489871e+01 -3.071644e+03  2.026651e+01 
    ds/dm[ 2] =  5.573604e+02     d2sdm2[ 2][*]:    3.697437e+01  2.026651e+01 -5.177324e+01 


Heat capacity and its derivatives:
----------------------------------

| ➜ Moles of components, T (K), P (bars) => isobaric heat capacity (J/K)
| ➜ Moles of components, T (K), P (bars) => d(isobaric heat capacity)/dT
  (J/K^2)
| ➜ Moles of components, T (K), P (bars) => d(isobaric heat
  capacity)/d(Moles of components) (J/K)

.. code:: ipython3

    print ('{0:<10s}{1:13.6f} {2:<15s}'.format('Cp', Feldspar.heat_capacity(t, p, mol=moles_end), 'J/K'))
    print ('{0:<10s}{1:13.6f} {2:<15s}'.format('dCp/dT', Feldspar.heat_capacity(t, p, mol=moles_end, deriv={'dT':1}), 'J/K^2'))
    print ("")
    dcpdm = Feldspar.heat_capacity(t, p, mol=moles_end, deriv={'dmol':1})[0]
    for i in range (0, Feldspar.props['endmember_num']):
        print ('dCp/dm[{0:>2d}]    = {1:13.6e} J/K-m'.format(i, dcpdm[i]))


.. parsed-literal::

    Cp           117.869612 J/K            
    dCp/dT        -0.069768 J/K^2          
    
    dCp/dm[ 0]    =  3.313154e+02 J/K-m
    dCp/dm[ 1]    =  3.208859e+02 J/K-m
    dCp/dm[ 2]    =  3.210186e+02 J/K-m


Volume and its derivatives:
---------------------------

| ➜ Moles of components, T (K), P (bars) => volume (J/bar)
| ➜ Moles of components, T (K), P (bars) => d(volume)/d(Moles of
  components) (J/bar)
| ➜ Moles of components, T (K), P (bars) => d^2(volume)/d(Moles of
  components)^2 (J/bar)
| ➜ Moles of components, T (K), P (bars) => d(volume)/dT (J/bar-K)
| ➜ Moles of components, T (K), P (bars) => d(volume)/dP (J/bar^2)
| ➜ Moles of components, T (K), P (bars) => d2(volume)/dT^2 (J/bar-K^2)
| ➜ Moles of components, T (K), P (bars) => d2(volume)/dTdP (J/bar^2-K)
| ➜ Moles of components, T (K), P (bars) => d2(volume)/dP^2 (J/bar^3)
| ➜ Moles of components, T (K), P (bars) => d2(volume)/d(Moles of
  components)dT (J/bar-K)
| ➜ Moles of components, T (K), P (bars) => d2(volume)/d(Moles of
  components)dP (J/bar^2)

.. code:: ipython3

    print ('{0:<10s}{1:13.6f} {2:<15s}'.format('Volume', Feldspar.volume(t, p, mol=moles_end), 'J/bar'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('dvdt', Feldspar.volume(t, p, mol=moles_end, deriv={'dT':1}), 'J/bar-K'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('dvdp', Feldspar.volume(t, p, mol=moles_end, deriv={'dP':1}), 'J/bar^2'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('d2vdt2', Feldspar.volume(t, p, mol=moles_end, deriv={'dT':2}), 'J/bar-K^2'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('d2vdtdp', Feldspar.volume(t, p, mol=moles_end, deriv={'dT':1,'dP':1}),  'J/bar^2-K'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('d2vdp2', Feldspar.volume(t, p, mol=moles_end, deriv={'dP':2}), 'J/bar^3'))
    print ("")
    dvdm = Feldspar.volume(t, p, mol=moles_end, deriv={'dmol':1})[0]
    d2vdm2 = Feldspar.volume(t, p, mol=moles_end, deriv={'dmol':2})[0]
    for i in range (0, Feldspar.props['endmember_num']):
        print ('dv/dm[{0:>2d}] = {1:13.6e} J/bar-m     d2vdm2[{0:>2d}][*]:  '.format(i, dvdm[i], i), end=' ')
        for j in range (0, Feldspar.props['endmember_num']):
            print ('{0:13.6e}'.format(d2vdm2[i, j]), end=' ')
        print('J/bar-m^2')


.. parsed-literal::

    Volume         3.843501 J/bar          
    dvdt       1.194789e-04 J/bar-K        
    dvdp      -7.056455e-06 J/bar^2        
    d2vdt2    -2.283345e-07 J/bar-K^2      
    d2vdtdp    7.335540e-10 J/bar^2-K      
    d2vdp2     3.499319e-11 J/bar^3        
    
    dv/dm[ 0] =  1.031408e+01 J/bar-m     d2vdm2[ 0][*]:   -3.493432e-01  4.277645e-01  4.781015e-01 J/bar-m^2
    dv/dm[ 1] =  9.782711e+00 J/bar-m     d2vdm2[ 1][*]:    4.277645e-01  3.692671e+00 -6.607703e-01 J/bar-m^2
    dv/dm[ 2] =  1.116056e+01 J/bar-m     d2vdm2[ 2][*]:    4.781015e-01 -6.607703e-01 -6.529701e-01 J/bar-m^2


Accessing properties of solution species:
-----------------------------------------

| ➜ Moles of components, T (K), P (bars) => formulae as an NSString
  object
| ➜ Retrieves the name of the solution species at the specified index
| ➜ Moles of solution species => moles of endmember components
| ➜ Retrieves an elemental stoichiometry vector for the species at the
  specified index
| ➜ Moles of components, T (K), P (bars) => chemical potentials of
  solution species (J)

.. code:: ipython3

    import numpy as np
    print ('Composition as a formula = ', Feldspar.compute_formula(t, p, moles_end))
    mSpecies = np.zeros(Feldspar.props['species_num'])
    for i in range (0, Feldspar.props['species_num']):
        print ('Species = {0:<15.15s}  elms:'.format(Feldspar.props['species_name'][i]), end=' ')
        elm = Feldspar.props['species_elms'][i]
        for j in range (0, 107):
            if elm[j] > 0.0:
                print ('[{0:>2.2s}] {1:5.2f}'.format(core.chem.PERIODIC_ORDER[j], elm[j]), end=' ')
        print ('   chemical potential = {0:12.2f} J/m'.format(Feldspar.chem_potential(t, p, mol=moles_end, endmember=i, species=True)))
        mSpecies[i] = float(i+1)
    for i in range (0, Feldspar.props['species_num']):
        print ('moles of species   [{0:>2d}] = {1:5.2f}'.format(i, mSpecies[i]))
    mSpToComp = Feldspar.convert_species_to_comp(mSpecies)
    for i in range (0, Feldspar.props['endmember_num']):
        print ('moles of component [{0:>2d}] = {1:5.2f}'.format(i, mSpToComp[i]))


.. parsed-literal::

    Composition as a formula =  K0.42Na0.58Ca0.01Al1.01Si2.99O8
    Species = albite           elms: [ O]  8.00 [Na]  1.00 [Al]  1.00 [Si]  3.00    chemical potential =  -4265565.12 J/m
    Species = anorthite        elms: [ O]  8.00 [Al]  2.00 [Si]  2.00 [Ca]  1.00    chemical potential =  -4575045.33 J/m
    Species = sanidine         elms: [ O]  8.00 [Al]  1.00 [Si]  3.00 [ K]  1.00    chemical potential =  -4305773.90 J/m
    moles of species   [ 0] =  1.00
    moles of species   [ 1] =  2.00
    moles of species   [ 2] =  3.00
    moles of component [ 0] =  1.00
    moles of component [ 1] =  2.00
    moles of component [ 2] =  3.00


Solution model parameters
-------------------------

Accessing, modifying, functional derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    if Feldspar.param_props['supports_calib']:
        print ('This phase supports the Calibration protocol', end=' ')
        np = Feldspar.param_props['param_num']
        print ('and there are', np, 'parameters')
        for i in range (0, np):
            name = Feldspar.param_names[i]
            value = Feldspar.get_param_values(param_names=[name])[0]
            units = Feldspar.param_units(param_names=[name])[0]
            dgdw = Feldspar.gibbs_energy(t, p, mol=moles_end, deriv_param=[name])
            print ('Parameter {0:<10.10s} = {1:13.6e} {2:<10.10s} dggw = {3:13.6e}'.format(
                name, value, units, Feldspar.gibbs_energy(t, p, mol=moles_end, deriv_param=[name])))
            dmudw = Feldspar.gibbs_energy(t, p, mol=moles_end, deriv={'dmol':1}, deriv_param=[name])[0]
            print ('   dmu[*]dw:', end=' ')
            for j in range (0, Feldspar.props['endmember_num']):
                print ('{0:13.6e}'.format(dmudw[j]), end=' ')
            print ()
        Feldspar.set_param_values(param_names=['whabor','whorab'], param_values=['1.0', '2.0'])
        print ()
        print('Reset values of whabor and whorab: ', Feldspar.get_param_values(param_names=['whabor','whorab']))
    else:
        print ('This phase does not implement the parameter calibration protocol')


.. parsed-literal::

    This phase supports the Calibration protocol and there are 13 parameters
    Parameter whabor     =  1.881000e+04 joules     dggw =  3.620592e-02
       dmu[*]dw: -2.691388e-02 -8.100080e-02  2.807416e-01 
    Parameter wsabor     =  1.030000e+01 joules     dggw = -3.620592e+01
       dmu[*]dw:  2.691141e+01  8.100121e+01 -2.807403e+02 
    Parameter wvabor     =  4.602000e-01 joules     dggw =  3.616972e+01
       dmu[*]dw: -2.689048e+01 -8.094307e+01  2.804487e+02 
    Parameter whorab     =  2.732000e+04 joules     dggw =  5.021743e-02
       dmu[*]dw:  2.024044e-01 -1.587322e-01  5.687912e-02 
    Parameter wsorab     =  1.030000e+01 joules     dggw = -5.021743e+01
       dmu[*]dw: -2.024090e+02  1.587318e+02 -5.688107e+01 
    Parameter wvorab     =  3.264000e-01 joules     dggw =  5.016721e+01
       dmu[*]dw:  2.022059e+02 -1.585478e+02  5.687040e+01 
    Parameter whaban     =  7.924000e+03 joules     dggw =  3.320741e-04
       dmu[*]dw: -2.523978e-04  1.265933e-01  2.997224e-04 
    Parameter whanab     =  0.000000e+00 joules     dggw =  1.212219e-03
       dmu[*]dw:  0.000000e+00  4.375000e-01  0.000000e+00 
    Parameter whoran     =  4.031700e+04 joules     dggw =  3.288566e-04
       dmu[*]dw: -2.852395e-04  1.242032e-01  3.720515e-04 
    Parameter whanor     =  3.897400e+04 joules     dggw =  7.817838e-04
       dmu[*]dw: -2.796736e-03  2.879420e-01  3.965785e-03 
    Parameter wvanor     = -1.037000e-01 joules     dggw =  7.810020e-01
       dmu[*]dw: -2.410800e+00  2.874879e+02  4.218901e+00 
    Parameter whabanor   =  1.254500e+04 joules     dggw =  6.412320e-04
       dmu[*]dw: -4.782782e-04  2.361748e-01  7.273814e-04 
    Parameter wvabanor   = -1.095000e+00 joules     dggw =  1.776956e+00
       dmu[*]dw: -4.566210e-01  2.359589e+02  7.420091e-01 
    
    Reset values of whabor and whorab:  [1. 2.]


