PhaseObjC - Solution Phase Example
----------------------------------

DEWFluid- number of endmember species is not equivalent to number of endmember components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Required python code to load the PhaseObjC library. The library
libphaseobjc.dylib (see build instructions in README.md) must be
locatable to teh system in a standard location (by default
/usr/local/lib)

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/rubicon/objc/ctypes_patch.py:24: UserWarning: rubicon.objc.ctypes_patch has only been tested with Python 3.4 through 3.6. The current version is sys.version_info(major=3, minor=7, micro=6, releaselevel='final', serial=0). Most likely things will work properly, but you may experience crashes if Python's internals have changed significantly.
      .format(sys.version_info)




.. parsed-literal::

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7fb233e8df10 at 0x7fb2481c0410>



Create a python reference to the “DEWFluid” solution phase class and instantiate an instance of that class.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In Objective-C the code has the form:

::

   obj = [[DEWFluid alloc] init]

and in python:

.. code:: ipython3

    DEWFluid = ObjCClass('DEWFluid')
    obj = DEWFluid.alloc().init()

Obtain properties of the phase inherited from the PhaseBase class.
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

.. code:: ipython3

    print (obj.phaseName)


.. parsed-literal::

    DEWFluid


Solution Protocol Functions
---------------------------

Solution component and species number …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Retrieves the number of endmember components in the system

::

   (NSUInteger)numberOfSolutionComponents  

–> Retrieves the number of species (dependent endmembers with positive
mole fractions) in the soluton

::

   (NSUInteger)numberOfSolutionSpecies

Note that the number of components (nc) may be the same as the number of
endmember species (ns) if the solution does not involve speciation
(complexing) or if the solid solution is not a reciprocal solution.

.. code:: ipython3

    nc = obj.numberOfSolutionComponents()
    print ('Number of components = ', nc)
    ns = obj.numberOfSolutionSpecies()
    print ('Number of species = ', ns)


.. parsed-literal::

    Number of components =  17
    Number of species =  92


Information about solution components …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Return name, formula and molecular weight of each end-member component
in the solution. Note teh use of the PhaseBase class (for further info
see the Stoichiometric Phase notebook examples)

–> Retrieves superclass instance of PhaseBase object for component at
specified index

::

   (id)componentAtIndex:(NSUInteger)index  

.. code:: ipython3

    PhaseBase = ObjCClass('PhaseBase')
    print ('component name, formula, and molecular weight (g/mol)')
    for i in range(0, nc):
        component = obj.componentAtIndex_(i)
        print (component.phaseName, component.phaseFormula, component.mw)


.. parsed-literal::

    component name, formula, and molecular weight (g/mol)
    Water H2O 18.0152
    CO2,aq CO2 44.0098
    O2,aq O2 31.9988
    HF,aq HF 20.006303
    NaOH,aq NaOH 39.99707
    Mg(OH)2,aq Mg(OH)2 58.319599999999994
    HAlO2,aq HAlO2 59.98824
    SiO2,aq SiO2 60.0843
    H3PO4,aq H3PO4 97.99506
    SO2,aq SO2 64.0588
    HCl,aq HCl 36.4609
    KOH,aq KOH 56.1093
    Ca(OH)2,aq Ca(OH)2 74.0946
    H2CrO4,aq H2CrO4 118.0094
    Mn(OH)2,aq Mn(OH)2 88.9526
    Fe(OH)2,aq Fe(OH)2 89.8616
    Co(OH)2,aq Co(OH)2 92.9478


Create a vector of moles of endmember components …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Allocate a “c”-type pointer to a double precision one-dimensional array,
and initialize the array to hold the moles of each component in the
solution

.. code:: ipython3

    import ctypes
    m = (ctypes.c_double*nc)()
    ctypes.cast(m, ctypes.POINTER(ctypes.c_double))
    m[ 0] = 0.998813
    m[ 1] = 0.0
    m[ 2] = 0.0
    m[ 3] = 0.0
    m[ 4] = 8.45249e-05
    m[ 5] = 0.0
    m[ 6] = 3.49734e-10
    m[ 7] = 0.00110226
    m[ 8] = 0.0
    m[ 9] = 0.0
    m[10] = 0.0
    m[11] = 5.03714e-07
    m[12] = 1.51788e-18
    m[13] = 0.0
    m[14] = 0.0
    m[15] = 0.0
    m[16] = 0.0
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('moles of (', component.phaseName, ') = ', m[i])


.. parsed-literal::

    moles of ( Water ) =  0.998813
    moles of ( CO2,aq ) =  0.0
    moles of ( O2,aq ) =  0.0
    moles of ( HF,aq ) =  0.0
    moles of ( NaOH,aq ) =  8.45249e-05
    moles of ( Mg(OH)2,aq ) =  0.0
    moles of ( HAlO2,aq ) =  3.49734e-10
    moles of ( SiO2,aq ) =  0.00110226
    moles of ( H3PO4,aq ) =  0.0
    moles of ( SO2,aq ) =  0.0
    moles of ( HCl,aq ) =  0.0
    moles of ( KOH,aq ) =  5.03714e-07
    moles of ( Ca(OH)2,aq ) =  1.51788e-18
    moles of ( H2CrO4,aq ) =  0.0
    moles of ( Mn(OH)2,aq ) =  0.0
    moles of ( Fe(OH)2,aq ) =  0.0
    moles of ( Co(OH)2,aq ) =  0.0


Note that moles can be assigned from a vector of element abundances using the following functions …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of elements (standard order) => Moles of end member components
of the phase

::

   (DoubleVector *)convertElementsToMoles:(double *)e

–> Moles of elements (standard order) => Total moles of end member
components of the phase

::

   (double)convertElementsToTotalMoles:(double *)e

–> Moles of elements (standard order) => Total mass of the phase (g)

::

   (double)convertElementsToTotalMass:(double *)e

.. code:: ipython3

    e = (ctypes.c_double*107)()
    ctypes.cast(e, ctypes.POINTER(ctypes.c_double))
    for i in range (0, 107):
        e[i] = 0.0
    e[1]  = 2.0*m[0] + m[4] +     m[6]            + m[11] + 2.0*m[12] # H
    e[8]  =     m[0] + m[4] + 2.0*m[6] + 2.0*m[7] + m[11] + 2.0*m[12] # O
    e[11] =            m[4]                                           # Na
    e[13] =                       m[6]                                # Al
    e[14] =                                  m[7]                     # Si
    e[19] =                                         m[11]             # K
    e[20] =                                                     m[12] # Ca
    mCompute = obj.convertElementsToMoles_(e)
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('assumed moles of', component.phaseName, '= ', m[i], ' computed = ', mCompute.valueAtIndex_(i))
    print ('Computed total number of moles = ', obj.convertElementsToTotalMoles_(e))
    print ('Computed total mass = ', obj.convertElementsToTotalMass_(e))


.. parsed-literal::

    assumed moles of Water =  0.998813  computed =  0.998813
    assumed moles of CO2,aq =  0.0  computed =  0.0
    assumed moles of O2,aq =  0.0  computed =  0.0
    assumed moles of HF,aq =  0.0  computed =  0.0
    assumed moles of NaOH,aq =  8.45249e-05  computed =  8.45249e-05
    assumed moles of Mg(OH)2,aq =  0.0  computed =  0.0
    assumed moles of HAlO2,aq =  3.49734e-10  computed =  3.49734e-10
    assumed moles of SiO2,aq =  0.00110226  computed =  0.00110226
    assumed moles of H3PO4,aq =  0.0  computed =  0.0
    assumed moles of SO2,aq =  0.0  computed =  0.0
    assumed moles of HCl,aq =  0.0  computed =  0.0
    assumed moles of KOH,aq =  5.03714e-07  computed =  5.03714e-07
    assumed moles of Ca(OH)2,aq =  1.51788e-18  computed =  1.51788e-18
    assumed moles of H2CrO4,aq =  0.0  computed =  0.0
    assumed moles of Mn(OH)2,aq =  0.0  computed =  0.0
    assumed moles of Fe(OH)2,aq =  0.0  computed =  0.0
    assumed moles of Co(OH)2,aq =  0.0  computed =  0.0
    Computed total number of moles =  1.000000288963734
    Computed total mass =  18.063453510479906


Test the mole vector and output derived quantities …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of endmember components => validity of input values

::

   (BOOL)testPermissibleValuesOfComponents:(double *)m

–> Moles of endmember components => Moles of elements (standard order)

::

   (DoubleVector *)convertMolesToElements:(double *)m

–> Moles of endmember components => Molar sum

::

   (double)totalMolesFromMolesOfComponents:(double *)m

–> Moles of endmember components => Mole fractions of endmember
components

::

   (DoubleVector *)convertMolesToMoleFractions:(double *)m

.. code:: ipython3

    if (obj.testPermissibleValuesOfComponents_(m) == 1):
        print ('mole input is feasible')
    else:
        print ('mole input is infeasible')
        
    print ('Total moles = ', obj.totalMolesFromMolesOfComponents_(m))
    
    mole_frac_pointer = obj.convertMolesToMoleFractions_(m)
    print ('component name, component formula, mole fraction')
    for i in range (0, nc):
        print (obj.componentAtIndex_(i).phaseName, obj.componentAtIndex_(i).phaseFormula, mole_frac_pointer.valueAtIndex_(i))
    
    moles_pointer = obj.convertMolesToElements_(m)
    ne = moles_pointer.size
    formula = ''
    for i in range(0, ne):
        value = moles_pointer.valueAtIndex_(i)
        if value != 0.0:
            name = PhaseBase.elementNameFromAtomicNumber_(i)
            formula = formula + name + '(' + str(value) + ')'
    print ('Solution formula = ', formula)


.. parsed-literal::

    mole input is feasible
    Total moles =  1.000000288963734
    component name, component formula, mole fraction
    Water H2O 0.9988127113793492
    CO2,aq CO2 0.0
    O2,aq O2 0.0
    HF,aq HF 0.0
    NaOH,aq NaOH 8.452487557537634e-05
    Mg(OH)2,aq Mg(OH)2 0.0
    HAlO2,aq HAlO2 3.4973389893958666e-10
    SiO2,aq SiO2 0.0011022596814869265
    H3PO4,aq H3PO4 0.0
    SO2,aq SO2 0.0
    HCl,aq HCl 0.0
    KOH,aq KOH 5.037138544449638e-07
    Ca(OH)2,aq Ca(OH)2 1.5178795613878541e-18
    H2CrO4,aq H2CrO4 0.0
    Mn(OH)2,aq Mn(OH)2 0.0
    Fe(OH)2,aq Fe(OH)2 0.0
    Co(OH)2,aq Co(OH)2 0.0
    Solution formula =  H(1.997711028963734)O(1.001102549313468)Na(8.45249e-05)Al(3.49734e-10)Si(0.00110226)K(5.03714e-07)Ca(1.51788e-18)


Preliminary Implementation
==========================

Compute activties and chemical potentials of endmember components …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => activities of endmember
components

::

   (DoubleVector *)getActivityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => chemical potentials of
endmember components (J)

::

   (DoubleVector *)getChemicalPotentialFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    t = 1000.0
    p = 1000.0
    activity = obj.getActivityFromMolesOfComponents_andT_andP_(m, t, p)
    potential = obj.getChemicalPotentialFromMolesOfComponents_andT_andP_(m, t, p)
    print ('component, activity, chemical potential')
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print (component.phaseName, activity.valueAtIndex_(i), potential.valueAtIndex_(i))


.. parsed-literal::

    component, activity, chemical potential
    Water nan nan
    CO2,aq nan nan
    O2,aq nan nan
    HF,aq nan nan
    NaOH,aq nan nan
    Mg(OH)2,aq nan nan
    HAlO2,aq nan nan
    SiO2,aq nan nan
    H3PO4,aq nan nan
    SO2,aq nan nan
    HCl,aq nan nan
    KOH,aq nan nan
    Ca(OH)2,aq nan nan
    H2CrO4,aq nan nan
    Mn(OH)2,aq nan nan
    Fe(OH)2,aq nan nan
    Co(OH)2,aq nan nan


Parameter calibration protocol
==============================

.. code:: ipython3

    component = obj.componentAtIndex_(10)
    print ('Does component', component.phaseName, 'implement the CalibrationProtocol?')
    
    try:
        if component.supportsParameterCalibration() == 1:
            print ('This component supports the Calibration protocol')
        np = component.getNumberOfFreeParameters()
        print ('... there are', np, 'parameters')
        names = component.getArrayOfNamesOfFreeParameters()
        for i in range (0, np):
            name = names.objectAtIndex_(i)
            value = component.getValueForParameterName_(name)
            units = component.getUnitsForParameterName_(name)
            print ('Parameter ', name, 'has value ', value, units)
            dmudw = component.getChemicalPotentialDerivativesForParameter_usingMolesOfComponents_andT_andP_(name, m, t, p)
            print ('   dmu0/dw = ', dmudw.valueAtIndex_(0))
    except AttributeError:
        print ('This phase does not implement the parameter calibration protocol')


.. parsed-literal::

    Does component HCl,aq implement the CalibrationProtocol?
    This component supports the Calibration protocol
    ... there are 12 parameters
    Parameter  deltaGibbsFreeEnergyOfFormationInTheReferenceState has value  -127235.44 joules/mol
       dmu0/dw =  1.0000000076752398
    Parameter  deltaEnthalpyOfFormationInTheReferenceState has value  -175953.93600000002 joules/mol
       dmu0/dw =  -0.0
    Parameter  entropyInTheReferenceState has value  13.388800000000002 joules/K/mol
       dmu0/dw =  -701.8500291260896
    Parameter  volumeInTheReferenceState has value  1.5 joules/bar/mol
       dmu0/dw =  0.0
    Parameter  heatCapacityInTheReferenceState has value  163.17600000000002 joules/K/mol
       dmu0/dw =  0.0
    Parameter  a1HKF has value  1.8148518400000002 joules/bar/mol
       dmu0/dw =  999.0000693391625
    Parameter  a2HKF has value  -64.14072 joules/mol
       dmu0/dw =  0.3250299673993947
    Parameter  a3HKF has value  20.7802544 joules-K/bar/mol
       dmu0/dw =  1.2940469605881244
    Parameter  a4HKF has value  -116009.768 joules-K/mol
       dmu0/dw =  0.0004210313581320377
    Parameter  c1HKF has value  113.7227936 joules/K/mol
       dmu0/dw =  -508.3085779471905
    Parameter  c2HKF has value  205421.848 joules-K/mol
       dmu0/dw =  -0.02102468523565826
    Parameter  omegaHKF has value  -83680.00000000001 joules/mol
       dmu0/dw =  0.34932168118128587


Higher derivatives NOT IMPLEMENTED
==================================

Execution of the code block following this cell will crash the python
kernel

Gibbs free energy and its compositional derivatives …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => Gibbs free energy (J)

::

   (double)getGibbsFreeEnergyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(Gibbs free energy)/d(Moles
of components) (J)

::

   (DoubleVector *)getDgDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d^2(Gibbs free
energy)/d(Moles of components)^2 (J)

::

   (DoubleMatrix *)getD2gDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d^3(Gibbs free
energy)/d(Moles of components)^3 (J)

::

   (DoubleTensor *)getD3gDm3FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

print (‘Gibbs free energy (J) =’,
obj.getGibbsFreeEnergyFromMolesOfComponents_andT_andP_(m, t, p)) dgdm =
obj.getDgDmFromMolesOfComponents_andT_andP_(m, t, p) for i in range (0,
nc): print (‘dg/dm (’, i, ‘) =’, dgdm.valueAtIndex_(i)) d2gdm2 =
obj.getD2gDm2FromMolesOfComponents_andT_andP_(m, t, p) for i in range
(0, nc): for j in range (0, nc): print (‘d2g/dm2 (’, i, ‘) (’, j, ‘) =’,
d2gdm2.valueAtRowIndex_andColIndex_(i, j)) d3gdm3 =
obj.getD3gDm3FromMolesOfComponents_andT_andP_(m, t, p) for i in range
(0, nc): for j in range (0, nc): for k in range (0, nc): print (‘d3g/dm3
(’, i, ‘) (’, j, ‘) (’, k, ‘) =’,
d3gdm3.valueAtFirstIndex_andSecondIndex_andThirdIndex_(i, j, k))

Higher derivatives NOT IMPLEMENTED
==================================

Execution of the code block following this cell will crash the python
kernel

Molar derivatives of activities …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => d(activities of endmember
components)/d(Moles of components)

::

   (DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

dadm = obj.getDaDmFromMolesOfComponents_andT_andP_(m, t, p) for i in
range (0, nc): for j in range (0, nc): print (‘da/dm (’, i, ‘) (’, j, ‘)
=’, dadm.valueAtRowIndex_andColIndex_(i, j))

Higher derivatives NOT IMPLEMENTED
==================================

Execution of the code block following this cell will crash the python
kernel

Enthalpy, Entroipy and molar derivatives …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => enthalpy (J)

::

   (double)getEnthalpyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => entropy (J/K)

::

   (double)getEntropyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(entropy)/d(Moles of
components) (J/K)

::

   (DoubleVector *)getDsDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d^2(entropy)/d(Moles of
components)^2 (J/K)

::

   (DoubleMatrix *)getD2sDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

print (‘Enthalpy (J) =’,
obj.getEnthalpyFromMolesOfComponents_andT_andP_(m, t, p)) print
(‘Entropy (J/K) =’, obj.getEntropyFromMolesOfComponents_andT_andP_(m, t,
p)) dsdm = obj.getDsDmFromMolesOfComponents_andT_andP_(m, t, p) for i in
range (0, nc): print (‘ds/dm (’, i, ‘) =’, dsdm.valueAtIndex_(i)) d2sdm2
= obj.getD2sDm2FromMolesOfComponents_andT_andP_(m, t, p) for i in range
(0, nc): for j in range (0, nc): print (‘d2s/dm2 (’, i, ‘) (’, j, ‘) =’,
d2sdm2.valueAtRowIndex_andColIndex_(i, j))

Higher derivatives NOT IMPLEMENTED
==================================

Execution of the code block following this cell will crash the python
kernel

Heat capcaity and its derivatives …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => isobaric heat capacity (J/K)

::

   (double)getHeatCapacityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(isobaric heat capacity)/dT
(J/K^2)

::

   (double)getDcpDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(isobaric heat
capacity)/d(Moles of components) (J/K)

::

   (DoubleVector *)getDCpDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

print (‘Heat capacity (J/K) =’,
obj.getHeatCapacityFromMolesOfComponents_andT_andP_(m, t, p)) print
(‘dcpdt (J/K^2) =’, obj.getDcpDtFromMolesOfComponents_andT_andP_(m, t,
p)) dcpdm = obj.getDCpDmFromMolesOfComponents_andT_andP_(m, t, p) for i
in range (0, nc): print (‘dcp/dm (’, i, ‘) =’, dcpdm.valueAtIndex_(i))

Higher derivatives NOT IMPLEMENTED
==================================

Execution of the code block following this cell will crash the python
kernel

Volume and its derivatives …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => volume (J/bar)

::

   (double)getVolumeFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(volume)/d(Moles of
components) (J/bar)

::

   (DoubleVector *)getDvDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d^2(volume)/d(Moles of
components)^2 (J/bar)

::

   (DoubleMatrix *)getD2vDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(volume)/dT (J/bar-K)

::

   (double)getDvDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d(volume)/dP (J/bar^2)

::

   (double)getDvDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d2(volume)/dT^2 (J/bar-K^2)

::

   (double)getD2vDt2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d2(volume)/dTdP (J/bar^2-K)

::

   (double)getD2vDtDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d2(volume)/dP^2 (J/bar^3)

::

   (double)getD2vDp2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d2(volume)/d(Moles of
components)dT (J/bar-K)

::

   (DoubleVector *)getD2vDmDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Moles of components, T (K), P (bars) => d2(volume)/d(Moles of
components)dP (J/bar^2)

::

   (DoubleVector *)getD2vDmDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

print (‘Volume (J/bar) =’,
obj.getVolumeFromMolesOfComponents_andT_andP_(m, t, p)) dvdm =
obj.getDvDmFromMolesOfComponents_andT_andP_(m, t, p) for i in range (0,
nc): print (‘dv/dm (’, i, ‘) =’, dvdm.valueAtIndex_(i)) d2vdm2 =
obj.getD2vDm2FromMolesOfComponents_andT_andP_(m, t, p) for i in range
(0, nc): for j in range (0, nc): print (‘d2v/dm2 (’, i, ‘) (’, j, ‘) =’,
d2vdm2.valueAtRowIndex_andColIndex_(i, j)) print (‘dvdt (J/bar-K) =’,
obj.getDvDtFromMolesOfComponents_andT_andP_(m, t, p)) print (‘dvdp
(J/bar^2) =’, obj.getDvDpFromMolesOfComponents_andT_andP_(m, t, p))
print (‘d2vdt2 (J/bar-K^2) =’,
obj.getD2vDt2FromMolesOfComponents_andT_andP_(m, t, p)) print (‘d2vdtdp
(J/bar^2-K) =’, obj.getD2vDtDpFromMolesOfComponents_andT_andP_(m, t, p))
print (‘d2vdp2 (J/bar^3) =’,
obj.getD2vDp2FromMolesOfComponents_andT_andP_(m, t, p)) d2vdmdt =
obj.getD2vDmDtFromMolesOfComponents_andT_andP_(m, t, p) for i in range
(0, nc): print (‘d2vdmdt (’, i, ‘) =’, d2vdmdt.valueAtIndex_(i)) d2vdmdp
= obj.getD2vDmDpFromMolesOfComponents_andT_andP_(m, t, p) for i in range
(0, nc): print (‘d2vdmdp (’, i, ‘) =’, d2vdmdp.valueAtIndex_(i))

Higher derivatives NOT IMPLEMENTED
==================================

Execution of the code block following this cell will crash the python
kernel

Accessing properties of solution species …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

–> Moles of components, T (K), P (bars) => formulae as an NSString
object

::

   (NSString *)getFormulaFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

–> Retrieves the name of the solution species at the specified index

::

   (NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index

–> Moles of solution species => moles of endmember components

::

   (DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies

–> Retrieves an elemental stoichiometry vector for the species at the
specified index

::

   (DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index

–> Moles of components, T (K), P (bars) => chemical potentials of
solution species (J)

::

   (DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

Note that the first nc species are identical to the solution components

print (‘formula =’, obj.getFormulaFromMolesOfComponents_andT_andP_(m, t,
p)) muSpecies =
obj.chemicalPotentialsOfSpeciesFromMolesOfComponents_andT_andP_(m, t, p)
for i in range (0, ns): print (‘species =’,
obj.nameOfSolutionSpeciesAtIndex_(i)) elm =
obj.elementalCompositionOfSpeciesAtIndex_(i) for j in range (0, 107): if
elm.valueAtIndex_(j) > 0.0: print (’ element (‘, j,’) = ‘,
elm.valueAtIndex_(j)) print (’ chemical potential = ’,
muSpecies.valueAtIndex_(i)) mSpecies = (ctypes.c_double*3)()
ctypes.cast(mSpecies, ctypes.POINTER(ctypes.c_double)) mSpecies[0] = 1
mSpecies[1] = 2 mSpecies[2] = 3 mSpToComp =
obj.convertMolesOfSpeciesToMolesOfComponents_(mSpecies) for i in range
(0, nc): print (‘moles of component (’, i, ‘) =’,
mSpToComp.valueAtIndex_(i))
