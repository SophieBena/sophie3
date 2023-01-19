
PhaseObjC - Solution Phase Example 1
====================================

Feldspar - number of endmember species equivalent to number of endmember components
-----------------------------------------------------------------------------------

Required Python code to load the phase library.

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))




.. parsed-literal::

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7f9bae5f7b10 at 0x106aa92b0>



Create a Python reference to the ``Feldspar`` solution phase class, and instantiate an instance of that class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    Feldspar = ObjCClass('FeldsparBerman')
    obj = Feldspar.alloc().init()

Obtain properties of the phase inherited from the PhaseBase class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print (obj.phaseName)


.. parsed-literal::

    Feldspar


Solution Protocol Functions
---------------------------

Solution component and species number:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Retrieves the number of endmember components in the system

::

    (NSUInteger)numberOfSolutionComponents  

➜ Retrieves the number of species (dependent endmembers with positive
mole fractions) in the solution

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

    Number of components =  3
    Number of species =  3


Information about solution components:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Return name, formula, and molecular weight of each endmember component
in the solution. Note the use of the ``PhaseBase`` class (For further
information, see the Stoichiometric Phase notebook examples).

➜ Retrieves superclass instance of ``PhaseBase`` object for component at
specified index

::

    (id)componentAtIndex:(NSUInteger)index  

.. code:: ipython3

    PhaseBase = ObjCClass('PhaseBase')
    print ("{0:>20s} {1:>20s} {2:>15s}".format('component name', 'formula', 'MW (g/mol)'))
    for i in range(0, nc):
        component = obj.componentAtIndex_(i)
        print ("{0:>20s} {1:>20s} {2:15.3f}".format(component.phaseName, component.phaseFormula, component.mw))


.. parsed-literal::

          component name              formula      MW (g/mol)
                  albite            NaAlSi3O8         262.223
               anorthite           CaAl2Si2O8         278.209
                sanidine             KAlSi3O8         278.335


Create a vector of moles of endmember components.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Allocate a “c”-type pointer to a double precision one-dimensional array,
and initialize the array to hold the moles of each component in the
solution.

.. code:: ipython3

    import ctypes
    m = (ctypes.c_double*nc)()
    ctypes.cast(m, ctypes.POINTER(ctypes.c_double))
    m[0] = 1.0
    m[1] = 2.0
    m[2] = 3.0
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('moles of (', component.phaseName, ') = ', m[i])


.. parsed-literal::

    moles of ( albite ) =  1.0
    moles of ( anorthite ) =  2.0
    moles of ( sanidine ) =  3.0


Note that moles can be assigned from a vector of element abundances using the following functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of elements (standard order) => Moles of endmember components of
the phase

::

    (DoubleVector *)convertElementsToMoles:(double *)e

➜ Moles of elements (standard order) => Total moles of endmember
components of the phase

::

    (double)convertElementsToTotalMoles:(double *)e

➜ Moles of elements (standard order) => Total mass of the phase (g)

::

    (double)convertElementsToTotalMass:(double *)e

.. code:: ipython3

    e = (ctypes.c_double*107)()
    ctypes.cast(e, ctypes.POINTER(ctypes.c_double))
    for i in range (0, 107):
        e[i] = 0.0
    e[8]  = m[0]*8.0 + m[1]*8.0 + m[2]*8.0 # O
    e[11] = m[0]                           # Na
    e[13] = m[0] + 2.0*m[1] + m[2]         # Al
    e[14] = m[0]*3.0 + m[1]*2.0 + m[2]*3.0 # Si
    e[19] = m[2]                           # K
    e[20] = m[1]                           # Ca
    mCompute = obj.convertElementsToMoles_(e)
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('assumed moles of {0:<10s} = {1:5.1f}  computed = {2:5.1f}'.format(component.phaseName, m[i], 
                                                                                  mCompute.valueAtIndex_(i)))
    print ('Computed total number of moles = ', obj.convertElementsToTotalMoles_(e))
    print ('Computed total mass = ', obj.convertElementsToTotalMass_(e))


.. parsed-literal::

    assumed moles of albite     =   1.0  computed =   1.0
    assumed moles of anorthite  =   2.0  computed =   2.0
    assumed moles of sanidine   =   3.0  computed =   3.0
    Computed total number of moles =  6.0
    Computed total mass =  1653.6472899999999


Test the mole vector and output derived quantities:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of endmember components => validity of input values

::

    (BOOL)testPermissibleValuesOfComponents:(double *)m

➜ Moles of endmember components => Moles of elements (standard order)

::

    (DoubleVector *)convertMolesToElements:(double *)m

➜ Moles of endmember components => Molar sum

::

    (double)totalMolesFromMolesOfComponents:(double *)m

➜ Moles of endmember components => Mole fractions of endmember
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
    print ('{0:<20s}{1:<20s}{2:>20s}'.format('component name', 'component formula', 'mole fraction'))
    for i in range (0, nc):
        print ('{0:<20s}{1:<20s}{2:20.13e}'.format(obj.componentAtIndex_(i).phaseName, 
                                                 obj.componentAtIndex_(i).phaseFormula, 
                                                 mole_frac_pointer.valueAtIndex_(i)))
    
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
    Total moles =  6.0
    component name      component formula          mole fraction
    albite              NaAlSi3O8            1.6666666666667e-01
    anorthite           CaAl2Si2O8           3.3333333333333e-01
    sanidine            KAlSi3O8             5.0000000000000e-01
    Solution formula =  O(48.0)Na(1.0)Al(8.0)Si(16.0)K(3.0)Ca(2.0)


Compute activities and chemical potentials of endmember components:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => activities of endmember
components

::

    (DoubleVector *)getActivityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => chemical potentials of
endmember components (J)

::

    (DoubleVector *)getChemicalPotentialFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    t = 1000.0
    p = 1000.0
    activity = obj.getActivityFromMolesOfComponents_andT_andP_(m, t, p)
    potential = obj.getChemicalPotentialFromMolesOfComponents_andT_andP_(m, t, p)
    print ('{0:<10s} {1:>20s} {2:>20s}'.format('component', 'activity', 'chemical potential'))
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('{0:<10s} {1:20.13e} {2:20.13e}'.format(component.phaseName, 
                                                   activity.valueAtIndex_(i), 
                                                   potential.valueAtIndex_(i)))


.. parsed-literal::

    component              activity   chemical potential
    albite      1.9039980487749e-01 -4.2779303912811e+06
    anorthite   1.6804137737887e+00 -4.5463872785053e+06
    sanidine    1.1510158436500e+00 -4.3009815409060e+06


Gibbs free energy and its compositional derivatives:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => Gibbs free energy (J)

::

    (double)getGibbsFreeEnergyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(Gibbs free energy)/d(Moles
of components) (J)

::

    (DoubleVector *)getDgDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d^2(Gibbs free energy)/d(Moles
of components)^2 (J)

::

    (DoubleMatrix *)getD2gDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d^3(Gibbs free energy)/d(Moles
of components)^3 (J)

::

    (DoubleTensor *)getD3gDm3FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    print ('Gibbs free energy (J) = ', obj.getGibbsFreeEnergyFromMolesOfComponents_andT_andP_(m, t, p))
    print ("")
    dgdm = obj.getDgDmFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('dg/dm (', i, ') = ', dgdm.valueAtIndex_(i))
    print ("")
    d2gdm2 = obj.getD2gDm2FromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        for j in range (0, nc):
            print ('d2g/dm2 (', i, ') (', j, ') = ', d2gdm2.valueAtRowIndex_andColIndex_(i, j))
    print ("")
    d3gdm3 = obj.getD3gDm3FromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        for j in range (0, nc):
            for k in range (0, nc):
                print ('d3g/dm3 (', i, ') (', j, ') (', k, ') = ', d3gdm3.valueAtFirstIndex_andSecondIndex_andThirdIndex_(i, j, k))


.. parsed-literal::

    Gibbs free energy (J) =  -26273649.571009792
    
    dg/dm ( 0 ) =  -4277930.391281143
    dg/dm ( 1 ) =  -4546387.278505273
    dg/dm ( 2 ) =  -4300981.540906034
    
    d2g/dm2 ( 0 ) ( 0 ) =  6370.3737597222225
    d2g/dm2 ( 0 ) ( 1 ) =  -2388.3560152777777
    d2g/dm2 ( 0 ) ( 2 ) =  -531.2205763888887
    d2g/dm2 ( 1 ) ( 0 ) =  -2388.356015277779
    d2g/dm2 ( 1 ) ( 1 ) =  -1746.5741152777791
    d2g/dm2 ( 1 ) ( 2 ) =  1960.5014152777783
    d2g/dm2 ( 2 ) ( 0 ) =  -531.2205763888888
    d2g/dm2 ( 2 ) ( 1 ) =  1960.501415277778
    d2g/dm2 ( 2 ) ( 2 ) =  -1129.9274180555558
    
    d3g/dm3 ( 0 ) ( 0 ) ( 0 ) =  -73.9575655092593
    d3g/dm3 ( 0 ) ( 0 ) ( 1 ) =  -268.6991104166666
    d3g/dm3 ( 0 ) ( 0 ) ( 2 ) =  -2.5631199845679262
    d3g/dm3 ( 0 ) ( 1 ) ( 0 ) =  -680.9961331018519
    d3g/dm3 ( 0 ) ( 1 ) ( 1 ) =  224.4508567129629
    d3g/dm3 ( 0 ) ( 1 ) ( 2 ) =  -334.5324746141975
    d3g/dm3 ( 0 ) ( 2 ) ( 0 ) =  272.3015618055556
    d3g/dm3 ( 0 ) ( 2 ) ( 1 ) =  -471.9648155092592
    d3g/dm3 ( 0 ) ( 2 ) ( 2 ) =  567.2570045524692
    d3g/dm3 ( 1 ) ( 0 ) ( 0 ) =  -161.60744212962982
    d3g/dm3 ( 1 ) ( 0 ) ( 1 ) =  1156.44336712963
    d3g/dm3 ( 1 ) ( 0 ) ( 2 ) =  -304.39633317901234
    d3g/dm3 ( 1 ) ( 1 ) ( 0 ) =  744.1463444444445
    d3g/dm3 ( 1 ) ( 1 ) ( 1 ) =  786.7901300925932
    d3g/dm3 ( 1 ) ( 1 ) ( 2 ) =  51.21902746913571
    d3g/dm3 ( 1 ) ( 2 ) ( 0 ) =  -29.531651388888903
    d3g/dm3 ( 1 ) ( 2 ) ( 1 ) =  -86.213313425926
    d3g/dm3 ( 1 ) ( 2 ) ( 2 ) =  -619.442536882716
    d3g/dm3 ( 2 ) ( 0 ) ( 0 ) =  -73.9575655092593
    d3g/dm3 ( 2 ) ( 0 ) ( 1 ) =  -268.6991104166666
    d3g/dm3 ( 2 ) ( 0 ) ( 2 ) =  -2.5631199845679262
    d3g/dm3 ( 2 ) ( 1 ) ( 0 ) =  -680.9961331018519
    d3g/dm3 ( 2 ) ( 1 ) ( 1 ) =  224.4508567129629
    d3g/dm3 ( 2 ) ( 1 ) ( 2 ) =  -334.5324746141975
    d3g/dm3 ( 2 ) ( 2 ) ( 0 ) =  272.3015618055556
    d3g/dm3 ( 2 ) ( 2 ) ( 1 ) =  -471.9648155092592
    d3g/dm3 ( 2 ) ( 2 ) ( 2 ) =  567.2570045524692


Molar derivatives of activities:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => d(activities of endmember
components)/d(Moles of components)

::

    (DoubleMatrix *)getDaDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    dadm = obj.getDaDmFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        for j in range (0, nc):
            print ('da/dm (', i, ') (', j, ') = ', dadm.valueAtRowIndex_andColIndex_(i, j))


.. parsed-literal::

    da/dm ( 0 ) ( 0 ) =  0.14588334806872377
    da/dm ( 0 ) ( 1 ) =  -0.054694023464015876
    da/dm ( 0 ) ( 2 ) =  -0.012165100380230668
    da/dm ( 1 ) ( 0 ) =  -0.4827136794178558
    da/dm ( 1 ) ( 1 ) =  -0.3530023213325952
    da/dm ( 1 ) ( 2 ) =  0.3962394406943487
    da/dm ( 2 ) ( 0 ) =  -0.0735411640061708
    da/dm ( 2 ) ( 1 ) =  0.27140807891019186
    da/dm ( 2 ) ( 2 ) =  -0.15642499793807094


Enthalpy, Entropy, and molar derivatives:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => enthalpy (J)

::

    (double)getEnthalpyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => entropy (J/K)

::

    (double)getEntropyFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(entropy)/d(Moles of
components) (J/K)

::

    (DoubleVector *)getDsDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d^2(entropy)/d(Moles of
components)^2 (J/K)

::

    (DoubleMatrix *)getD2sDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    print ('Enthalpy (J) = ', obj.getEnthalpyFromMolesOfComponents_andT_andP_(m, t, p))
    print ('Entropy (J/K) = ', obj.getEntropyFromMolesOfComponents_andT_andP_(m, t, p))
    print ("")
    dsdm = obj.getDsDmFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('ds/dm (', i, ') = ', dsdm.valueAtIndex_(i))
    print ("")
    d2sdm2 = obj.getD2sDm2FromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        for j in range (0, nc):
            print ('d2s/dm2 (', i, ') (', j, ') = ', d2sdm2.valueAtRowIndex_andColIndex_(i, j))


.. parsed-literal::

    Enthalpy (J) =  -22971079.86850238
    Entropy (J/K) =  3302.569702507414
    
    ds/dm ( 0 ) =  563.4372036636106
    ds/dm ( 1 ) =  539.7714091733154
    ds/dm ( 2 ) =  553.1965601657242
    
    d2s/dm2 ( 0 ) ( 0 ) =  -8.35913888888889
    d2s/dm2 ( 0 ) ( 1 ) =  0.8134944444444441
    d2s/dm2 ( 0 ) ( 2 ) =  2.2440499999999997
    d2s/dm2 ( 1 ) ( 0 ) =  0.8134944444444445
    d2s/dm2 ( 1 ) ( 1 ) =  -2.485322222222223
    d2s/dm2 ( 1 ) ( 2 ) =  1.3857166666666672
    d2s/dm2 ( 2 ) ( 0 ) =  2.24405
    d2s/dm2 ( 2 ) ( 1 ) =  1.385716666666667
    d2s/dm2 ( 2 ) ( 2 ) =  -1.6718277777777777


Heat capacity and its derivatives:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => isobaric heat capacity (J/K)

::

    (double)getHeatCapacityFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(isobaric heat capacity)/dT
(J/K^2)

::

    (double)getDcpDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(isobaric heat
capacity)/d(Moles of components) (J/K)

::

    (DoubleVector *)getDCpDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    print ('Heat capacity (J/K) = ', obj.getHeatCapacityFromMolesOfComponents_andT_andP_(m, t, p))
    print ("")
    print ('dcpdt (J/K^2) = ', obj.getDcpDtFromMolesOfComponents_andT_andP_(m, t, p))
    print ("")
    dcpdm = obj.getDCpDmFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('dcp/dm (', i, ') = ', dcpdm.valueAtIndex_(i))


.. parsed-literal::

    Heat capacity (J/K) =  1936.142955434807
    
    dcpdt (J/K^2) =  -0.14642280679830652
    
    dcp/dm ( 0 ) =  331.3154455676114
    dcp/dm ( 1 ) =  320.8858547164275
    dcp/dm ( 2 ) =  321.0186001447802


Volume and its derivatives:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => volume (J/bar)

::

    (double)getVolumeFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(volume)/d(Moles of
components) (J/bar)

::

    (DoubleVector *)getDvDmFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d^2(volume)/d(Moles of
components)^2 (J/bar)

::

    (DoubleMatrix *)getD2vDm2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(volume)/dT (J/bar-K)

::

    (double)getDvDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d(volume)/dP (J/bar^2)

::

    (double)getDvDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d2(volume)/dT^2 (J/bar-K^2)

::

    (double)getD2vDt2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d2(volume)/dTdP (J/bar^2-K)

::

    (double)getD2vDtDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d2(volume)/dP^2 (J/bar^3)

::

    (double)getD2vDp2FromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d2(volume)/d(Moles of
components)dT (J/bar-K)

::

    (DoubleVector *)getD2vDmDtFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Moles of components, T (K), P (bars) => d2(volume)/d(Moles of
components)dP (J/bar^2)

::

    (DoubleVector *)getD2vDmDpFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

.. code:: ipython3

    print ('{0:<10s}{1:13.6f} {2:<15s}'.format('Volume', obj.getVolumeFromMolesOfComponents_andT_andP_(m, t, p), 'J/bar'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('dvdt', obj.getDvDtFromMolesOfComponents_andT_andP_(m, t, p), 'J/bar-K'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('dvdp', obj.getDvDpFromMolesOfComponents_andT_andP_(m, t, p), 'J/bar^2'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('d2vdt2', obj.getD2vDt2FromMolesOfComponents_andT_andP_(m, t, p), 'J/bar-K^2'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('d2vdtdp', obj.getD2vDtDpFromMolesOfComponents_andT_andP_(m, t, p),  'J/bar^2-K'))
    print ('{0:<10s}{1:13.6e} {2:<15s}'.format('d2vdp2', obj.getD2vDp2FromMolesOfComponents_andT_andP_(m, t, p), 'J/bar^3'))
    print ("")
    dvdm = obj.getDvDmFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('dv/dm (', i, ') = ', dvdm.valueAtIndex_(i))
    print ("")
    d2vdm2 = obj.getD2vDm2FromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        for j in range (0, nc):
            print ('d2v/dm2 (', i, ') (', j, ') = ', d2vdm2.valueAtRowIndex_andColIndex_(i, j))
    print ("")
    d2vdmdt = obj.getD2vDmDtFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('d2vdmdt (', i, ') = ', d2vdmdt.valueAtIndex_(i))
    print ("")
    d2vdmdp = obj.getD2vDmDpFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('d2vdmdp (', i, ') = ', d2vdmdp.valueAtIndex_(i))


.. parsed-literal::

    Volume        63.586940 J/bar          
    dvdt       1.539713e-03 J/bar-K        
    dvdp      -1.037353e-04 J/bar^2        
    d2vdt2    -7.606511e-07 J/bar-K^2      
    d2vdtdp    3.524411e-09 J/bar^2-K      
    d2vdp2     5.487518e-10 J/bar^3        
    
    dv/dm ( 0 ) =  10.315177941459686
    dv/dm ( 1 ) =  10.08327841108531
    dv/dm ( 2 ) =  11.035068569246032
    
    d2v/dm2 ( 0 ) ( 0 ) =  0.01929583333333334
    d2v/dm2 ( 0 ) ( 1 ) =  -0.05370694444444446
    d2v/dm2 ( 0 ) ( 2 ) =  0.02937268518518519
    d2v/dm2 ( 1 ) ( 0 ) =  -0.05370694444444446
    d2v/dm2 ( 1 ) ( 1 ) =  0.05328194444444446
    d2v/dm2 ( 1 ) ( 2 ) =  -0.017618981481481485
    d2v/dm2 ( 2 ) ( 0 ) =  0.02937268518518519
    d2v/dm2 ( 2 ) ( 1 ) =  -0.017618981481481482
    d2v/dm2 ( 2 ) ( 2 ) =  0.0019550925925925933
    
    d2vdmdt ( 0 ) =  0.0003735005077224484
    d2vdmdt ( 1 ) =  0.00016937520208375
    d2vdmdt ( 2 ) =  0.0002758206068054751
    
    d2vdmdp ( 0 ) =  -1.9709725291649575e-05
    d2vdmdp ( 1 ) =  -1.2751467596399999e-05
    d2vdmdp ( 2 ) =  -1.9507531468656002e-05


Accessing properties of solution species:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Moles of components, T (K), P (bars) => formulae as an NSString object

::

    (NSString *)getFormulaFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

➜ Retrieves the name of the solution species at the specified index

::

    (NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index

➜ Moles of solution species => moles of endmember components

::

    (DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies

➜ Retrieves an elemental stoichiometry vector for the species at the
specified index

::

    (DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index

➜ Moles of components, T (K), P (bars) => chemical potentials of
solution species (J)

::

    (DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p

Note that the first nc species are identical to the solution components

.. code:: ipython3

    print ('formula = ', obj.getFormulaFromMolesOfComponents_andT_andP_(m, t, p))
    muSpecies = obj.chemicalPotentialsOfSpeciesFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, ns):
        print ('species = ', obj.nameOfSolutionSpeciesAtIndex_(i))
        elm = obj.elementalCompositionOfSpeciesAtIndex_(i)
        for j in range (0, 107):
            if elm.valueAtIndex_(j) > 0.0:
                print ('   element (', j, ') = ', elm.valueAtIndex_(j))
        print ('   chemical potential = ', muSpecies.valueAtIndex_(i))
    mSpecies = (ctypes.c_double*3)()
    ctypes.cast(mSpecies, ctypes.POINTER(ctypes.c_double))
    mSpecies[0] = 1
    mSpecies[1] = 2
    mSpecies[2] = 3
    mSpToComp = obj.convertMolesOfSpeciesToMolesOfComponents_(mSpecies)
    for i in range (0, nc):
        print ('moles of component (', i, ') = ', mSpToComp.valueAtIndex_(i))


.. parsed-literal::

    formula =  K0.50Na0.17Ca0.33Al1.33Si2.67O8
    species =  albite
       element ( 8 ) =  8.0
       element ( 11 ) =  1.0
       element ( 13 ) =  1.0
       element ( 14 ) =  3.0
       chemical potential =  -4277930.391281143
    species =  anorthite
       element ( 8 ) =  8.0
       element ( 13 ) =  2.0
       element ( 14 ) =  2.0
       element ( 20 ) =  1.0
       chemical potential =  -4546387.278505273
    species =  sanidine
       element ( 8 ) =  8.0
       element ( 13 ) =  1.0
       element ( 14 ) =  3.0
       element ( 19 ) =  1.0
       chemical potential =  -4300981.540906034
    moles of component ( 0 ) =  1.0
    moles of component ( 1 ) =  2.0
    moles of component ( 2 ) =  3.0


Parameter calibration protocol
------------------------------

.. code:: ipython3

    try:
        if obj.supportsParameterCalibration() == 1:
            print ('This phase supports the Calibration protocol')
        np = obj.getNumberOfFreeParameters()
        print ('... there are', np, 'parameters')
        names = obj.getArrayOfNamesOfFreeParameters()
        dgdw = obj.getDgDwFromMolesOfComponents_andT_andP_(m, t, p)
        for i in range (0, np):
            name = names.objectAtIndex_(i)
            value = obj.getValueForParameterName_(name)
            units = obj.getUnitsForParameterName_(name)
            print ('Parameter ', name, 'has value ', value, units, 'dgdw = ', dgdw.valueAtIndex_(i))
            dmudw = obj.getChemicalPotentialDerivativesForParameter_usingMolesOfComponents_andT_andP_(name, m, t, p)
            for j in range (0, nc):
                print ('   dmu (', j, ')dw = ', dmudw.valueAtIndex_(j))
    except AttributeError:
        print ('This phase does not implement the parameter calibration protocol')


.. parsed-literal::

    This phase supports the Calibration protocol
    ... there are 13 parameters
    Parameter  whabor has value  18810.0 joules dgdw =  0.3333333333333333
       dmu ( 0 )dw =  0.2222222222222222
       dmu ( 1 )dw =  -0.06944444444444445
       dmu ( 2 )dw =  0.08333333333333333
    Parameter  wsabor has value  10.3 joules/K dgdw =  -333.3333333333333
       dmu ( 0 )dw =  -222.22087481517642
       dmu ( 1 )dw =  69.44781585548834
       dmu ( 2 )dw =  -83.33131106539287
    Parameter  wvabor has value  0.4602 joules/bar dgdw =  333.0
       dmu ( 0 )dw =  222.05019529967933
       dmu ( 1 )dw =  -69.39917418846247
       dmu ( 2 )dw =  83.2518469227544
    Parameter  whorab has value  27320.0 joules dgdw =  0.16666666666666666
       dmu ( 0 )dw =  0.19444546120058565
       dmu ( 1 )dw =  -0.013888634699853587
       dmu ( 2 )dw =  0.0
    Parameter  wsorab has value  10.3 joules/K dgdw =  -166.66666666666666
       dmu ( 0 )dw =  -194.4417484729811
       dmu ( 1 )dw =  13.889563171097668
       dmu ( 2 )dw =  0.0
    Parameter  wvorab has value  0.3264 joules/bar dgdw =  166.5
       dmu ( 0 )dw =  194.35508542584833
       dmu ( 1 )dw =  -13.786764680454265
       dmu ( 2 )dw =  0.0
    Parameter  whaban has value  7924.0 joules dgdw =  0.1944444444444444
       dmu ( 0 )dw =  0.12962992175668855
       dmu ( 1 )dw =  0.08796851337708228
       dmu ( 2 )dw =  -0.03703937405350833
    Parameter  whanab has value  0.0 joules dgdw =  0.13888888888888887
       dmu ( 0 )dw =  0.1875
       dmu ( 1 )dw =  0.0625
       dmu ( 2 )dw =  0.0
    Parameter  whoran has value  40317.0 joules dgdw =  0.41666666666666663
       dmu ( 0 )dw =  -0.055555038817372324
       dmu ( 1 )dw =  0.23611162784929435
       dmu ( 2 )dw =  0.0
    Parameter  whanor has value  38974.0 joules dgdw =  0.5833333333333334
       dmu ( 0 )dw =  -0.1111109329296454
       dmu ( 1 )dw =  0.09722346949248217
       dmu ( 2 )dw =  0.1666672012110638
    Parameter  wvanor has value  -0.1037 joules/bar dgdw =  582.75
       dmu ( 0 )dw =  -111.49951798892664
       dmu ( 1 )dw =  97.03471565522804
       dmu ( 2 )dw =  166.3452268375338
    Parameter  whabanor has value  12545.0 joules dgdw =  0.16666666666666666
       dmu ( 0 )dw =  0.11111498605021922
       dmu ( 1 )dw =  0.02777999202869669
       dmu ( 2 )dw =  0.0
    Parameter  wvabanor has value  -1.095 joules/bar dgdw =  27.75
       dmu ( 0 )dw =  110.9589042303866
       dmu ( 1 )dw =  27.73972605759665
       dmu ( 2 )dw =  -0.0

