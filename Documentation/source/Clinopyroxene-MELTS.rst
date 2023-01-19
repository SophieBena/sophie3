Clinopyroxene - number of endmember species equivalent to number of endmember components
========================================================================================

Required Python code to load the phase library.

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/rubicon/objc/ctypes_patch.py:24: UserWarning: rubicon.objc.ctypes_patch has only been tested with Python 3.4 through 3.6. The current version is sys.version_info(major=3, minor=7, micro=6, releaselevel='final', serial=0). Most likely things will work properly, but you may experience crashes if Python's internals have changed significantly.
      .format(sys.version_info)




.. parsed-literal::

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7fc0e7cd3c80 at 0x7fc09854ab10>



Create a Python reference to the ``Clinopyroxene`` solution phase class, and instantiate an instance of that class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    Cpx = ObjCClass('CpxBerman')
    obj = Cpx.alloc().init()

Obtain properties of the phase inherited from the ``PhaseBase`` class.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print (obj.phaseName)


.. parsed-literal::

    Clinopyroxene


Solution Protocol Functions
---------------------------

Solution component and species number:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

➜ Retrieves the number of endmember components in the system

::

   (NSUInteger)numberOfSolutionComponents  

➜ Retrieves the number of species (dependent endmembers with positive
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

    Number of components =  7
    Number of species =  14


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
                diopside            CaMgSi2O6         216.552
          clinoenstatite             Mg2Si2O6         200.777
            hedenbergite            CaFeSi2O6         248.094
       alumino-buffonite   CaTi0.5Mg0.5AlSiO6         227.246
               buffonite   CaTi0.5Mg0.5FeSiO6         256.111
                essenite           CaFeAlSiO6         246.990
                 jadeite            NaAlSi2O6         202.139


Create a vector of moles of endmember components.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Allocate a “c”-type pointer to a double precision one-dimensional array,
and initialize the array to hold the moles of each component in the
solution.

.. code:: ipython3

    import ctypes
    m = (ctypes.c_double*nc)()
    ctypes.cast(m, ctypes.POINTER(ctypes.c_double))
    m[0] =  1.0
    m[1] =  2.0
    m[2] =  3.0
    m[3] =  1.5
    m[4] = -1.3
    m[5] =  1.4
    m[6] =  0.5
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('moles of (', component.phaseName.ljust(20), ') = ', m[i])


.. parsed-literal::

    moles of ( diopside             ) =  1.0
    moles of ( clinoenstatite       ) =  2.0
    moles of ( hedenbergite         ) =  3.0
    moles of ( alumino-buffonite    ) =  1.5
    moles of ( buffonite            ) =  -1.3
    moles of ( essenite             ) =  1.4
    moles of ( jadeite              ) =  0.5


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
    e[8]  = m[0]*6.0 + m[1]*6.0 + m[2]*6.0 + m[3]*6.0 + m[4]*6.0 + m[5]*6.0 + m[6]*6.0 # O
    e[11] = m[6]                                                                       # Na
    e[12] = m[0] + m[1]*2.0 + m[3]*0.5 + m[4]*0.5                                      # Mg
    e[13] = m[3] + m[5] + m[6]                                                         # Al
    e[14] = m[0]*2.0 + m[1]*2.0 + m[2]*2.0 + m[3] + m[4] + m[5] + 2.0*m[6]             # Si
    e[20] = m[0] + m[2] + m[3] + m[4] + m[5]                                           # Ca
    e[22] = m[3]*0.5 + m[4]*0.5                                                        # Ti
    e[26] = m[2] + m[4] + m[5]                                                         # Fe
    mCompute = obj.convertElementsToMoles_(e)
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('assumed moles of {0:<10s} = {1:5.1f}  computed = {2:5.1f}'.format(component.phaseName, m[i], 
                                                                                  mCompute.valueAtIndex_(i)))
    print ('Computed total number of moles = ', obj.convertElementsToTotalMoles_(e))
    print ('Computed total mass = ', obj.convertElementsToTotalMass_(e))


.. parsed-literal::

    assumed moles of diopside   =   1.0  computed =   1.0
    assumed moles of clinoenstatite =   2.0  computed =   2.0
    assumed moles of hedenbergite =   3.0  computed =   3.0
    assumed moles of alumino-buffonite =   1.5  computed =   1.5
    assumed moles of buffonite  =  -1.3  computed =  -1.3
    assumed moles of essenite   =   1.4  computed =   1.4
    assumed moles of jadeite    =   0.5  computed =   0.5
    Computed total number of moles =  8.1
    Computed total mass =  1817.1704610000002


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
    Total moles =  8.1
    component name      component formula          mole fraction
    diopside            CaMgSi2O6            1.2345679012346e-01
    clinoenstatite      Mg2Si2O6             2.4691358024691e-01
    hedenbergite        CaFeSi2O6            3.7037037037037e-01
    alumino-buffonite   CaTi0.5Mg0.5AlSiO6   1.8518518518519e-01
    buffonite           CaTi0.5Mg0.5FeSiO6  -1.6049382716049e-01
    essenite            CaFeAlSiO6           1.7283950617284e-01
    jadeite             NaAlSi2O6            6.1728395061728e-02
    Solution formula =  O(48.6)Na(0.5)Mg(5.1)Al(3.4)Si(14.6)Ca(5.6)Ti(0.09999999999999998)Fe(3.0999999999999996)


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
    print ('{0:<20s} {1:>20s} {2:>20s}'.format('component', 'activity', 'chemical potential'))
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        print ('{0:<20s} {1:20.13e} {2:20.13e}'.format(component.phaseName, 
                                                   activity.valueAtIndex_(i), 
                                                   potential.valueAtIndex_(i)))


.. parsed-literal::

    component                        activity   chemical potential
    diopside              5.4339210428822e-01 -3.4450276296895e+06
    clinoenstatite        3.6620618253059e-01 -3.3265693414830e+06
    hedenbergite          2.0710734522041e-01 -3.1282135794144e+06
    alumino-buffonite     8.2583290211142e-01 -3.5283913273026e+06
    buffonite             1.4172341259723e-02 -3.1491507685069e+06
    essenite              4.1697879900223e-02 -3.1743250543351e+06
    jadeite               4.0521374174147e-02 -3.2794771577169e+06


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
    dgdm = obj.getDgDmFromMolesOfComponents_andT_andP_(m, t, p)
    print ("")
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

    Gibbs free energy (J) =  -26765291.69772121
    
    dg/dm ( 0 ) =  -3445027.629689464
    dg/dm ( 1 ) =  -3326569.3414830053
    dg/dm ( 2 ) =  -3128213.579414352
    dg/dm ( 3 ) =  -3528391.327302648
    dg/dm ( 4 ) =  -3149150.7685068524
    dg/dm ( 5 ) =  -3174325.0543351285
    dg/dm ( 6 ) =  -3279477.157716875
    
    d2g/dm2 ( 0 ) ( 0 ) =  577.451325640708
    d2g/dm2 ( 0 ) ( 1 ) =  944.5339793407175
    d2g/dm2 ( 0 ) ( 2 ) =  -216.5516977318639
    d2g/dm2 ( 0 ) ( 3 ) =  -955.4593837687955
    d2g/dm2 ( 0 ) ( 4 ) =  -612.1951461999807
    d2g/dm2 ( 0 ) ( 5 ) =  -490.2278652876815
    d2g/dm2 ( 0 ) ( 6 ) =  -986.4195882611632
    d2g/dm2 ( 1 ) ( 0 ) =  944.5339793407175
    d2g/dm2 ( 1 ) ( 1 ) =  4090.4566487442453
    d2g/dm2 ( 1 ) ( 2 ) =  -2992.4628738391184
    d2g/dm2 ( 1 ) ( 3 ) =  1455.0547636848783
    d2g/dm2 ( 1 ) ( 4 ) =  3207.2896103602493
    d2g/dm2 ( 1 ) ( 5 ) =  1302.2196754686713
    d2g/dm2 ( 1 ) ( 6 ) =  31.456293946018263
    d2g/dm2 ( 2 ) ( 0 ) =  -216.55169773186395
    d2g/dm2 ( 2 ) ( 1 ) =  -2992.4628738391175
    d2g/dm2 ( 2 ) ( 2 ) =  2989.9683051394177
    d2g/dm2 ( 2 ) ( 3 ) =  -1563.9134780782463
    d2g/dm2 ( 2 ) ( 4 ) =  -2427.296116307686
    d2g/dm2 ( 2 ) ( 5 ) =  -2182.6605201358548
    d2g/dm2 ( 2 ) ( 6 ) =  -1044.6349518011739
    d2g/dm2 ( 3 ) ( 0 ) =  -955.4593837687956
    d2g/dm2 ( 3 ) ( 1 ) =  1455.0547636848776
    d2g/dm2 ( 3 ) ( 2 ) =  -1563.9134780782456
    d2g/dm2 ( 3 ) ( 3 ) =  15820.252805335742
    d2g/dm2 ( 3 ) ( 4 ) =  12789.958058397251
    d2g/dm2 ( 3 ) ( 5 ) =  -2860.7542813212517
    d2g/dm2 ( 3 ) ( 6 ) =  -722.5748952073358
    d2g/dm2 ( 4 ) ( 0 ) =  -612.1951461999807
    d2g/dm2 ( 4 ) ( 1 ) =  3207.289610360248
    d2g/dm2 ( 4 ) ( 2 ) =  -2427.2961163076857
    d2g/dm2 ( 4 ) ( 3 ) =  12789.95805839725
    d2g/dm2 ( 4 ) ( 4 ) =  95035.04688101422
    d2g/dm2 ( 4 ) ( 5 ) =  76900.64205519359
    d2g/dm2 ( 4 ) ( 6 ) =  -3641.5414902917464
    d2g/dm2 ( 5 ) ( 0 ) =  -490.2278652876821
    d2g/dm2 ( 5 ) ( 1 ) =  1302.2196754686709
    d2g/dm2 ( 5 ) ( 2 ) =  -2182.660520135853
    d2g/dm2 ( 5 ) ( 3 ) =  -2860.754281321254
    d2g/dm2 ( 5 ) ( 4 ) =  76900.64205519359
    d2g/dm2 ( 5 ) ( 5 ) =  79661.21511506289
    d2g/dm2 ( 5 ) ( 6 ) =  -5659.929985193169
    d2g/dm2 ( 6 ) ( 0 ) =  -986.4195882611608
    d2g/dm2 ( 6 ) ( 1 ) =  31.456293946020423
    d2g/dm2 ( 6 ) ( 2 ) =  -1044.6349518011712
    d2g/dm2 ( 6 ) ( 3 ) =  -722.574895207335
    d2g/dm2 ( 6 ) ( 4 ) =  -3641.5414902917446
    d2g/dm2 ( 6 ) ( 5 ) =  -5659.929985193171
    d2g/dm2 ( 6 ) ( 6 ) =  16662.344480949632
    
    d3g/dm3 ( 0 ) ( 0 ) ( 0 ) =  55.78920067831597
    d3g/dm3 ( 0 ) ( 0 ) ( 1 ) =  -8.041065459832481
    d3g/dm3 ( 0 ) ( 0 ) ( 2 ) =  -14.779774876809526
    d3g/dm3 ( 0 ) ( 0 ) ( 3 ) =  105.32514451958218
    d3g/dm3 ( 0 ) ( 0 ) ( 4 ) =  157.86949583597743
    d3g/dm3 ( 0 ) ( 0 ) ( 5 ) =  77.02793387052434
    d3g/dm3 ( 0 ) ( 0 ) ( 6 ) =  69.34693506409621
    d3g/dm3 ( 0 ) ( 1 ) ( 0 ) =  121.98560008470508
    d3g/dm3 ( 0 ) ( 1 ) ( 1 ) =  442.5969532343207
    d3g/dm3 ( 0 ) ( 1 ) ( 2 ) =  -60.455339751644146
    d3g/dm3 ( 0 ) ( 1 ) ( 3 ) =  254.61064105072796
    d3g/dm3 ( 0 ) ( 1 ) ( 4 ) =  464.15322653668466
    d3g/dm3 ( 0 ) ( 1 ) ( 5 ) =  124.92431338399336
    d3g/dm3 ( 0 ) ( 1 ) ( 6 ) =  -4.90652726410849
    d3g/dm3 ( 0 ) ( 2 ) ( 0 ) =  -110.58815255726262
    d3g/dm3 ( 0 ) ( 2 ) ( 1 ) =  101.73354662920016
    d3g/dm3 ( 0 ) ( 2 ) ( 2 ) =  -256.46501627844964
    d3g/dm3 ( 0 ) ( 2 ) ( 3 ) =  -95.31491702534586
    d3g/dm3 ( 0 ) ( 2 ) ( 4 ) =  -155.1506035906118
    d3g/dm3 ( 0 ) ( 2 ) ( 5 ) =  -43.412918013079604
    d3g/dm3 ( 0 ) ( 2 ) ( 6 ) =  -13.678765434634236
    d3g/dm3 ( 0 ) ( 3 ) ( 0 ) =  167.61538971859648
    d3g/dm3 ( 0 ) ( 3 ) ( 1 ) =  -65.40127235080325
    d3g/dm3 ( 0 ) ( 3 ) ( 2 ) =  62.78370585412148
    d3g/dm3 ( 0 ) ( 3 ) ( 3 ) =  248.94344348687147
    d3g/dm3 ( 0 ) ( 3 ) ( 4 ) =  324.75016727121556
    d3g/dm3 ( 0 ) ( 3 ) ( 5 ) =  224.06232499161675
    d3g/dm3 ( 0 ) ( 3 ) ( 6 ) =  210.19502607655667
    d3g/dm3 ( 0 ) ( 4 ) ( 0 ) =  263.10961480748244
    d3g/dm3 ( 0 ) ( 4 ) ( 1 ) =  13.144198129056718
    d3g/dm3 ( 0 ) ( 4 ) ( 2 ) =  45.89789306134633
    d3g/dm3 ( 0 ) ( 4 ) ( 3 ) =  367.7000410437063
    d3g/dm3 ( 0 ) ( 4 ) ( 4 ) =  484.158601968835
    d3g/dm3 ( 0 ) ( 4 ) ( 5 ) =  329.22390408240983
    d3g/dm3 ( 0 ) ( 4 ) ( 6 ) =  265.8672419394313
    d3g/dm3 ( 0 ) ( 5 ) ( 0 ) =  120.47755700895874
    d3g/dm3 ( 0 ) ( 5 ) ( 1 ) =  -137.62370273276946
    d3g/dm3 ( 0 ) ( 5 ) ( 2 ) =  95.84508280580798
    d3g/dm3 ( 0 ) ( 5 ) ( 3 ) =  205.22170293103696
    d3g/dm3 ( 0 ) ( 5 ) ( 4 ) =  267.4334082493392
    d3g/dm3 ( 0 ) ( 5 ) ( 5 ) =  194.52065983676627
    d3g/dm3 ( 0 ) ( 5 ) ( 6 ) =  154.63238457639204
    d3g/dm3 ( 0 ) ( 6 ) ( 0 ) =  89.18516790991765
    d3g/dm3 ( 0 ) ( 6 ) ( 1 ) =  -195.43980298840142
    d3g/dm3 ( 0 ) ( 6 ) ( 2 ) =  101.96784509164031
    d3g/dm3 ( 0 ) ( 6 ) ( 3 ) =  167.74301372336387
    d3g/dm3 ( 0 ) ( 6 ) ( 4 ) =  180.4653558137479
    d3g/dm3 ( 0 ) ( 6 ) ( 5 ) =  131.0209942837789
    d3g/dm3 ( 0 ) ( 6 ) ( 6 ) =  93.35866218052266
    d3g/dm3 ( 1 ) ( 0 ) ( 0 ) =  -124.8382159330601
    d3g/dm3 ( 1 ) ( 0 ) ( 1 ) =  -113.69723106171182
    d3g/dm3 ( 1 ) ( 0 ) ( 2 ) =  90.39715951007256
    d3g/dm3 ( 1 ) ( 0 ) ( 3 ) =  -275.92284464892185
    d3g/dm3 ( 1 ) ( 0 ) ( 4 ) =  -436.183116163927
    d3g/dm3 ( 1 ) ( 0 ) ( 5 ) =  -189.61635216929537
    d3g/dm3 ( 1 ) ( 0 ) ( 6 ) =  -166.1893058096897
    d3g/dm3 ( 1 ) ( 1 ) ( 0 ) =  16.329434482825718
    d3g/dm3 ( 1 ) ( 1 ) ( 1 ) =  -2534.4965273252656
    d3g/dm3 ( 1 ) ( 1 ) ( 2 ) =  572.7743009836909
    d3g/dm3 ( 1 ) ( 1 ) ( 3 ) =  -388.17694046354376
    d3g/dm3 ( 1 ) ( 1 ) ( 4 ) =  -1027.2818261957111
    d3g/dm3 ( 1 ) ( 1 ) ( 5 ) =  7.366358919996799
    d3g/dm3 ( 1 ) ( 1 ) ( 6 ) =  403.3504228967072
    d3g/dm3 ( 1 ) ( 2 ) ( 0 ) =  -5.411218170380522
    d3g/dm3 ( 1 ) ( 2 ) ( 1 ) =  734.9631873645354
    d3g/dm3 ( 1 ) ( 2 ) ( 2 ) =  439.51321617923963
    d3g/dm3 ( 1 ) ( 2 ) ( 3 ) =  -51.99458654272644
    d3g/dm3 ( 1 ) ( 2 ) ( 4 ) =  130.5042574813345
    d3g/dm3 ( 1 ) ( 2 ) ( 5 ) =  -210.2956835301387
    d3g/dm3 ( 1 ) ( 2 ) ( 6 ) =  -300.9848488943971
    d3g/dm3 ( 1 ) ( 3 ) ( 0 ) =  -213.63259944990764
    d3g/dm3 ( 1 ) ( 3 ) ( 1 ) =  -708.1888538650752
    d3g/dm3 ( 1 ) ( 3 ) ( 2 ) =  106.1040363367409
    d3g/dm3 ( 1 ) ( 3 ) ( 3 ) =  -461.68316344314644
    d3g/dm3 ( 1 ) ( 3 ) ( 4 ) =  -692.8936709853957
    d3g/dm3 ( 1 ) ( 3 ) ( 5 ) =  -385.7957520326195
    d3g/dm3 ( 1 ) ( 3 ) ( 6 ) =  -343.50049034168626
    d3g/dm3 ( 1 ) ( 4 ) ( 0 ) =  -330.9429971924221
    d3g/dm3 ( 1 ) ( 4 ) ( 1 ) =  -1478.290854603339
    d3g/dm3 ( 1 ) ( 4 ) ( 2 ) =  331.5527541332925
    d3g/dm3 ( 1 ) ( 4 ) ( 3 ) =  -649.9437972129052
    d3g/dm3 ( 1 ) ( 4 ) ( 4 ) =  -1005.1424080345486
    d3g/dm3 ( 1 ) ( 4 ) ( 5 ) =  -532.5915794809509
    d3g/dm3 ( 1 ) ( 4 ) ( 6 ) =  -339.35375994486697
    d3g/dm3 ( 1 ) ( 5 ) ( 0 ) =  -146.1667290308609
    d3g/dm3 ( 1 ) ( 5 ) ( 1 ) =  -255.18165719676597
    d3g/dm3 ( 1 ) ( 5 ) ( 2 ) =  -71.03768271125107
    d3g/dm3 ( 1 ) ( 5 ) ( 3 ) =  -404.63637409319915
    d3g/dm3 ( 1 ) ( 5 ) ( 4 ) =  -594.3820753140214
    d3g/dm3 ( 1 ) ( 5 ) ( 5 ) =  -371.9981926556737
    d3g/dm3 ( 1 ) ( 5 ) ( 6 ) =  -250.33895311153236
    d3g/dm3 ( 1 ) ( 6 ) ( 0 ) =  -146.35107296386838
    d3g/dm3 ( 1 ) ( 6 ) ( 1 ) =  212.81714717241428
    d3g/dm3 ( 1 ) ( 6 ) ( 2 ) =  -185.33823836812266
    d3g/dm3 ( 1 ) ( 6 ) ( 3 ) =  -385.9525026948792
    d3g/dm3 ( 1 ) ( 6 ) ( 4 ) =  -424.7556460705507
    d3g/dm3 ( 1 ) ( 6 ) ( 5 ) =  -273.9503434041451
    d3g/dm3 ( 1 ) ( 6 ) ( 6 ) =  -159.08023048921365
    d3g/dm3 ( 2 ) ( 0 ) ( 0 ) =  55.78920067831597
    d3g/dm3 ( 2 ) ( 0 ) ( 1 ) =  -8.041065459832481
    d3g/dm3 ( 2 ) ( 0 ) ( 2 ) =  -14.779774876809526
    d3g/dm3 ( 2 ) ( 0 ) ( 3 ) =  105.32514451958218
    d3g/dm3 ( 2 ) ( 0 ) ( 4 ) =  157.86949583597743
    d3g/dm3 ( 2 ) ( 0 ) ( 5 ) =  77.02793387052434
    d3g/dm3 ( 2 ) ( 0 ) ( 6 ) =  69.34693506409621
    d3g/dm3 ( 2 ) ( 1 ) ( 0 ) =  121.98560008470508
    d3g/dm3 ( 2 ) ( 1 ) ( 1 ) =  442.5969532343207
    d3g/dm3 ( 2 ) ( 1 ) ( 2 ) =  -60.455339751644146
    d3g/dm3 ( 2 ) ( 1 ) ( 3 ) =  254.61064105072796
    d3g/dm3 ( 2 ) ( 1 ) ( 4 ) =  464.15322653668466
    d3g/dm3 ( 2 ) ( 1 ) ( 5 ) =  124.92431338399336
    d3g/dm3 ( 2 ) ( 1 ) ( 6 ) =  -4.90652726410849
    d3g/dm3 ( 2 ) ( 2 ) ( 0 ) =  -110.58815255726262
    d3g/dm3 ( 2 ) ( 2 ) ( 1 ) =  101.73354662920016
    d3g/dm3 ( 2 ) ( 2 ) ( 2 ) =  -256.46501627844964
    d3g/dm3 ( 2 ) ( 2 ) ( 3 ) =  -95.31491702534586
    d3g/dm3 ( 2 ) ( 2 ) ( 4 ) =  -155.1506035906118
    d3g/dm3 ( 2 ) ( 2 ) ( 5 ) =  -43.412918013079604
    d3g/dm3 ( 2 ) ( 2 ) ( 6 ) =  -13.678765434634236
    d3g/dm3 ( 2 ) ( 3 ) ( 0 ) =  167.61538971859648
    d3g/dm3 ( 2 ) ( 3 ) ( 1 ) =  -65.40127235080325
    d3g/dm3 ( 2 ) ( 3 ) ( 2 ) =  62.78370585412148
    d3g/dm3 ( 2 ) ( 3 ) ( 3 ) =  248.94344348687147
    d3g/dm3 ( 2 ) ( 3 ) ( 4 ) =  324.75016727121556
    d3g/dm3 ( 2 ) ( 3 ) ( 5 ) =  224.06232499161675
    d3g/dm3 ( 2 ) ( 3 ) ( 6 ) =  210.19502607655667
    d3g/dm3 ( 2 ) ( 4 ) ( 0 ) =  263.10961480748244
    d3g/dm3 ( 2 ) ( 4 ) ( 1 ) =  13.144198129056718
    d3g/dm3 ( 2 ) ( 4 ) ( 2 ) =  45.89789306134633
    d3g/dm3 ( 2 ) ( 4 ) ( 3 ) =  367.7000410437063
    d3g/dm3 ( 2 ) ( 4 ) ( 4 ) =  484.158601968835
    d3g/dm3 ( 2 ) ( 4 ) ( 5 ) =  329.22390408240983
    d3g/dm3 ( 2 ) ( 4 ) ( 6 ) =  265.8672419394313
    d3g/dm3 ( 2 ) ( 5 ) ( 0 ) =  120.47755700895874
    d3g/dm3 ( 2 ) ( 5 ) ( 1 ) =  -137.62370273276946
    d3g/dm3 ( 2 ) ( 5 ) ( 2 ) =  95.84508280580798
    d3g/dm3 ( 2 ) ( 5 ) ( 3 ) =  205.22170293103696
    d3g/dm3 ( 2 ) ( 5 ) ( 4 ) =  267.4334082493392
    d3g/dm3 ( 2 ) ( 5 ) ( 5 ) =  194.52065983676627
    d3g/dm3 ( 2 ) ( 5 ) ( 6 ) =  154.63238457639204
    d3g/dm3 ( 2 ) ( 6 ) ( 0 ) =  89.18516790991765
    d3g/dm3 ( 2 ) ( 6 ) ( 1 ) =  -195.43980298840142
    d3g/dm3 ( 2 ) ( 6 ) ( 2 ) =  101.96784509164031
    d3g/dm3 ( 2 ) ( 6 ) ( 3 ) =  167.74301372336387
    d3g/dm3 ( 2 ) ( 6 ) ( 4 ) =  180.4653558137479
    d3g/dm3 ( 2 ) ( 6 ) ( 5 ) =  131.0209942837789
    d3g/dm3 ( 2 ) ( 6 ) ( 6 ) =  93.35866218052266
    d3g/dm3 ( 3 ) ( 0 ) ( 0 ) =  55.78920067831597
    d3g/dm3 ( 3 ) ( 0 ) ( 1 ) =  -8.041065459832481
    d3g/dm3 ( 3 ) ( 0 ) ( 2 ) =  -14.779774876809526
    d3g/dm3 ( 3 ) ( 0 ) ( 3 ) =  105.32514451958218
    d3g/dm3 ( 3 ) ( 0 ) ( 4 ) =  157.86949583597743
    d3g/dm3 ( 3 ) ( 0 ) ( 5 ) =  77.02793387052434
    d3g/dm3 ( 3 ) ( 0 ) ( 6 ) =  69.34693506409621
    d3g/dm3 ( 3 ) ( 1 ) ( 0 ) =  121.98560008470508
    d3g/dm3 ( 3 ) ( 1 ) ( 1 ) =  442.5969532343207
    d3g/dm3 ( 3 ) ( 1 ) ( 2 ) =  -60.455339751644146
    d3g/dm3 ( 3 ) ( 1 ) ( 3 ) =  254.61064105072796
    d3g/dm3 ( 3 ) ( 1 ) ( 4 ) =  464.15322653668466
    d3g/dm3 ( 3 ) ( 1 ) ( 5 ) =  124.92431338399336
    d3g/dm3 ( 3 ) ( 1 ) ( 6 ) =  -4.90652726410849
    d3g/dm3 ( 3 ) ( 2 ) ( 0 ) =  -110.58815255726262
    d3g/dm3 ( 3 ) ( 2 ) ( 1 ) =  101.73354662920016
    d3g/dm3 ( 3 ) ( 2 ) ( 2 ) =  -256.46501627844964
    d3g/dm3 ( 3 ) ( 2 ) ( 3 ) =  -95.31491702534586
    d3g/dm3 ( 3 ) ( 2 ) ( 4 ) =  -155.1506035906118
    d3g/dm3 ( 3 ) ( 2 ) ( 5 ) =  -43.412918013079604
    d3g/dm3 ( 3 ) ( 2 ) ( 6 ) =  -13.678765434634236
    d3g/dm3 ( 3 ) ( 3 ) ( 0 ) =  167.61538971859648
    d3g/dm3 ( 3 ) ( 3 ) ( 1 ) =  -65.40127235080325
    d3g/dm3 ( 3 ) ( 3 ) ( 2 ) =  62.78370585412148
    d3g/dm3 ( 3 ) ( 3 ) ( 3 ) =  248.94344348687147
    d3g/dm3 ( 3 ) ( 3 ) ( 4 ) =  324.75016727121556
    d3g/dm3 ( 3 ) ( 3 ) ( 5 ) =  224.06232499161675
    d3g/dm3 ( 3 ) ( 3 ) ( 6 ) =  210.19502607655667
    d3g/dm3 ( 3 ) ( 4 ) ( 0 ) =  263.10961480748244
    d3g/dm3 ( 3 ) ( 4 ) ( 1 ) =  13.144198129056718
    d3g/dm3 ( 3 ) ( 4 ) ( 2 ) =  45.89789306134633
    d3g/dm3 ( 3 ) ( 4 ) ( 3 ) =  367.7000410437063
    d3g/dm3 ( 3 ) ( 4 ) ( 4 ) =  484.158601968835
    d3g/dm3 ( 3 ) ( 4 ) ( 5 ) =  329.22390408240983
    d3g/dm3 ( 3 ) ( 4 ) ( 6 ) =  265.8672419394313
    d3g/dm3 ( 3 ) ( 5 ) ( 0 ) =  120.47755700895874
    d3g/dm3 ( 3 ) ( 5 ) ( 1 ) =  -137.62370273276946
    d3g/dm3 ( 3 ) ( 5 ) ( 2 ) =  95.84508280580798
    d3g/dm3 ( 3 ) ( 5 ) ( 3 ) =  205.22170293103696
    d3g/dm3 ( 3 ) ( 5 ) ( 4 ) =  267.4334082493392
    d3g/dm3 ( 3 ) ( 5 ) ( 5 ) =  194.52065983676627
    d3g/dm3 ( 3 ) ( 5 ) ( 6 ) =  154.63238457639204
    d3g/dm3 ( 3 ) ( 6 ) ( 0 ) =  89.18516790991765
    d3g/dm3 ( 3 ) ( 6 ) ( 1 ) =  -195.43980298840142
    d3g/dm3 ( 3 ) ( 6 ) ( 2 ) =  101.96784509164031
    d3g/dm3 ( 3 ) ( 6 ) ( 3 ) =  167.74301372336387
    d3g/dm3 ( 3 ) ( 6 ) ( 4 ) =  180.4653558137479
    d3g/dm3 ( 3 ) ( 6 ) ( 5 ) =  131.0209942837789
    d3g/dm3 ( 3 ) ( 6 ) ( 6 ) =  93.35866218052266
    d3g/dm3 ( 4 ) ( 0 ) ( 0 ) =  55.78920067831597
    d3g/dm3 ( 4 ) ( 0 ) ( 1 ) =  -8.041065459832481
    d3g/dm3 ( 4 ) ( 0 ) ( 2 ) =  -14.779774876809526
    d3g/dm3 ( 4 ) ( 0 ) ( 3 ) =  105.32514451958218
    d3g/dm3 ( 4 ) ( 0 ) ( 4 ) =  157.86949583597743
    d3g/dm3 ( 4 ) ( 0 ) ( 5 ) =  77.02793387052434
    d3g/dm3 ( 4 ) ( 0 ) ( 6 ) =  69.34693506409621
    d3g/dm3 ( 4 ) ( 1 ) ( 0 ) =  121.98560008470508
    d3g/dm3 ( 4 ) ( 1 ) ( 1 ) =  442.5969532343207
    d3g/dm3 ( 4 ) ( 1 ) ( 2 ) =  -60.455339751644146
    d3g/dm3 ( 4 ) ( 1 ) ( 3 ) =  254.61064105072796
    d3g/dm3 ( 4 ) ( 1 ) ( 4 ) =  464.15322653668466
    d3g/dm3 ( 4 ) ( 1 ) ( 5 ) =  124.92431338399336
    d3g/dm3 ( 4 ) ( 1 ) ( 6 ) =  -4.90652726410849
    d3g/dm3 ( 4 ) ( 2 ) ( 0 ) =  -110.58815255726262
    d3g/dm3 ( 4 ) ( 2 ) ( 1 ) =  101.73354662920016
    d3g/dm3 ( 4 ) ( 2 ) ( 2 ) =  -256.46501627844964
    d3g/dm3 ( 4 ) ( 2 ) ( 3 ) =  -95.31491702534586
    d3g/dm3 ( 4 ) ( 2 ) ( 4 ) =  -155.1506035906118
    d3g/dm3 ( 4 ) ( 2 ) ( 5 ) =  -43.412918013079604
    d3g/dm3 ( 4 ) ( 2 ) ( 6 ) =  -13.678765434634236
    d3g/dm3 ( 4 ) ( 3 ) ( 0 ) =  167.61538971859648
    d3g/dm3 ( 4 ) ( 3 ) ( 1 ) =  -65.40127235080325
    d3g/dm3 ( 4 ) ( 3 ) ( 2 ) =  62.78370585412148
    d3g/dm3 ( 4 ) ( 3 ) ( 3 ) =  248.94344348687147
    d3g/dm3 ( 4 ) ( 3 ) ( 4 ) =  324.75016727121556
    d3g/dm3 ( 4 ) ( 3 ) ( 5 ) =  224.06232499161675
    d3g/dm3 ( 4 ) ( 3 ) ( 6 ) =  210.19502607655667
    d3g/dm3 ( 4 ) ( 4 ) ( 0 ) =  263.10961480748244
    d3g/dm3 ( 4 ) ( 4 ) ( 1 ) =  13.144198129056718
    d3g/dm3 ( 4 ) ( 4 ) ( 2 ) =  45.89789306134633
    d3g/dm3 ( 4 ) ( 4 ) ( 3 ) =  367.7000410437063
    d3g/dm3 ( 4 ) ( 4 ) ( 4 ) =  484.158601968835
    d3g/dm3 ( 4 ) ( 4 ) ( 5 ) =  329.22390408240983
    d3g/dm3 ( 4 ) ( 4 ) ( 6 ) =  265.8672419394313
    d3g/dm3 ( 4 ) ( 5 ) ( 0 ) =  120.47755700895874
    d3g/dm3 ( 4 ) ( 5 ) ( 1 ) =  -137.62370273276946
    d3g/dm3 ( 4 ) ( 5 ) ( 2 ) =  95.84508280580798
    d3g/dm3 ( 4 ) ( 5 ) ( 3 ) =  205.22170293103696
    d3g/dm3 ( 4 ) ( 5 ) ( 4 ) =  267.4334082493392
    d3g/dm3 ( 4 ) ( 5 ) ( 5 ) =  194.52065983676627
    d3g/dm3 ( 4 ) ( 5 ) ( 6 ) =  154.63238457639204
    d3g/dm3 ( 4 ) ( 6 ) ( 0 ) =  89.18516790991765
    d3g/dm3 ( 4 ) ( 6 ) ( 1 ) =  -195.43980298840142
    d3g/dm3 ( 4 ) ( 6 ) ( 2 ) =  101.96784509164031
    d3g/dm3 ( 4 ) ( 6 ) ( 3 ) =  167.74301372336387
    d3g/dm3 ( 4 ) ( 6 ) ( 4 ) =  180.4653558137479
    d3g/dm3 ( 4 ) ( 6 ) ( 5 ) =  131.0209942837789
    d3g/dm3 ( 4 ) ( 6 ) ( 6 ) =  93.35866218052266
    d3g/dm3 ( 5 ) ( 0 ) ( 0 ) =  55.78920067831597
    d3g/dm3 ( 5 ) ( 0 ) ( 1 ) =  -8.041065459832481
    d3g/dm3 ( 5 ) ( 0 ) ( 2 ) =  -14.779774876809526
    d3g/dm3 ( 5 ) ( 0 ) ( 3 ) =  105.32514451958218
    d3g/dm3 ( 5 ) ( 0 ) ( 4 ) =  157.86949583597743
    d3g/dm3 ( 5 ) ( 0 ) ( 5 ) =  77.02793387052434
    d3g/dm3 ( 5 ) ( 0 ) ( 6 ) =  69.34693506409621
    d3g/dm3 ( 5 ) ( 1 ) ( 0 ) =  121.98560008470508
    d3g/dm3 ( 5 ) ( 1 ) ( 1 ) =  442.5969532343207
    d3g/dm3 ( 5 ) ( 1 ) ( 2 ) =  -60.455339751644146
    d3g/dm3 ( 5 ) ( 1 ) ( 3 ) =  254.61064105072796
    d3g/dm3 ( 5 ) ( 1 ) ( 4 ) =  464.15322653668466
    d3g/dm3 ( 5 ) ( 1 ) ( 5 ) =  124.92431338399336
    d3g/dm3 ( 5 ) ( 1 ) ( 6 ) =  -4.90652726410849
    d3g/dm3 ( 5 ) ( 2 ) ( 0 ) =  -110.58815255726262
    d3g/dm3 ( 5 ) ( 2 ) ( 1 ) =  101.73354662920016
    d3g/dm3 ( 5 ) ( 2 ) ( 2 ) =  -256.46501627844964
    d3g/dm3 ( 5 ) ( 2 ) ( 3 ) =  -95.31491702534586
    d3g/dm3 ( 5 ) ( 2 ) ( 4 ) =  -155.1506035906118
    d3g/dm3 ( 5 ) ( 2 ) ( 5 ) =  -43.412918013079604
    d3g/dm3 ( 5 ) ( 2 ) ( 6 ) =  -13.678765434634236
    d3g/dm3 ( 5 ) ( 3 ) ( 0 ) =  167.61538971859648
    d3g/dm3 ( 5 ) ( 3 ) ( 1 ) =  -65.40127235080325
    d3g/dm3 ( 5 ) ( 3 ) ( 2 ) =  62.78370585412148
    d3g/dm3 ( 5 ) ( 3 ) ( 3 ) =  248.94344348687147
    d3g/dm3 ( 5 ) ( 3 ) ( 4 ) =  324.75016727121556
    d3g/dm3 ( 5 ) ( 3 ) ( 5 ) =  224.06232499161675
    d3g/dm3 ( 5 ) ( 3 ) ( 6 ) =  210.19502607655667
    d3g/dm3 ( 5 ) ( 4 ) ( 0 ) =  263.10961480748244
    d3g/dm3 ( 5 ) ( 4 ) ( 1 ) =  13.144198129056718
    d3g/dm3 ( 5 ) ( 4 ) ( 2 ) =  45.89789306134633
    d3g/dm3 ( 5 ) ( 4 ) ( 3 ) =  367.7000410437063
    d3g/dm3 ( 5 ) ( 4 ) ( 4 ) =  484.158601968835
    d3g/dm3 ( 5 ) ( 4 ) ( 5 ) =  329.22390408240983
    d3g/dm3 ( 5 ) ( 4 ) ( 6 ) =  265.8672419394313
    d3g/dm3 ( 5 ) ( 5 ) ( 0 ) =  120.47755700895874
    d3g/dm3 ( 5 ) ( 5 ) ( 1 ) =  -137.62370273276946
    d3g/dm3 ( 5 ) ( 5 ) ( 2 ) =  95.84508280580798
    d3g/dm3 ( 5 ) ( 5 ) ( 3 ) =  205.22170293103696
    d3g/dm3 ( 5 ) ( 5 ) ( 4 ) =  267.4334082493392
    d3g/dm3 ( 5 ) ( 5 ) ( 5 ) =  194.52065983676627
    d3g/dm3 ( 5 ) ( 5 ) ( 6 ) =  154.63238457639204
    d3g/dm3 ( 5 ) ( 6 ) ( 0 ) =  89.18516790991765
    d3g/dm3 ( 5 ) ( 6 ) ( 1 ) =  -195.43980298840142
    d3g/dm3 ( 5 ) ( 6 ) ( 2 ) =  101.96784509164031
    d3g/dm3 ( 5 ) ( 6 ) ( 3 ) =  167.74301372336387
    d3g/dm3 ( 5 ) ( 6 ) ( 4 ) =  180.4653558137479
    d3g/dm3 ( 5 ) ( 6 ) ( 5 ) =  131.0209942837789
    d3g/dm3 ( 5 ) ( 6 ) ( 6 ) =  93.35866218052266
    d3g/dm3 ( 6 ) ( 0 ) ( 0 ) =  55.78920067831597
    d3g/dm3 ( 6 ) ( 0 ) ( 1 ) =  -8.041065459832481
    d3g/dm3 ( 6 ) ( 0 ) ( 2 ) =  -14.779774876809526
    d3g/dm3 ( 6 ) ( 0 ) ( 3 ) =  105.32514451958218
    d3g/dm3 ( 6 ) ( 0 ) ( 4 ) =  157.86949583597743
    d3g/dm3 ( 6 ) ( 0 ) ( 5 ) =  77.02793387052434
    d3g/dm3 ( 6 ) ( 0 ) ( 6 ) =  69.34693506409621
    d3g/dm3 ( 6 ) ( 1 ) ( 0 ) =  121.98560008470508
    d3g/dm3 ( 6 ) ( 1 ) ( 1 ) =  442.5969532343207
    d3g/dm3 ( 6 ) ( 1 ) ( 2 ) =  -60.455339751644146
    d3g/dm3 ( 6 ) ( 1 ) ( 3 ) =  254.61064105072796
    d3g/dm3 ( 6 ) ( 1 ) ( 4 ) =  464.15322653668466
    d3g/dm3 ( 6 ) ( 1 ) ( 5 ) =  124.92431338399336
    d3g/dm3 ( 6 ) ( 1 ) ( 6 ) =  -4.90652726410849
    d3g/dm3 ( 6 ) ( 2 ) ( 0 ) =  -110.58815255726262
    d3g/dm3 ( 6 ) ( 2 ) ( 1 ) =  101.73354662920016
    d3g/dm3 ( 6 ) ( 2 ) ( 2 ) =  -256.46501627844964
    d3g/dm3 ( 6 ) ( 2 ) ( 3 ) =  -95.31491702534586
    d3g/dm3 ( 6 ) ( 2 ) ( 4 ) =  -155.1506035906118
    d3g/dm3 ( 6 ) ( 2 ) ( 5 ) =  -43.412918013079604
    d3g/dm3 ( 6 ) ( 2 ) ( 6 ) =  -13.678765434634236
    d3g/dm3 ( 6 ) ( 3 ) ( 0 ) =  167.61538971859648
    d3g/dm3 ( 6 ) ( 3 ) ( 1 ) =  -65.40127235080325
    d3g/dm3 ( 6 ) ( 3 ) ( 2 ) =  62.78370585412148
    d3g/dm3 ( 6 ) ( 3 ) ( 3 ) =  248.94344348687147
    d3g/dm3 ( 6 ) ( 3 ) ( 4 ) =  324.75016727121556
    d3g/dm3 ( 6 ) ( 3 ) ( 5 ) =  224.06232499161675
    d3g/dm3 ( 6 ) ( 3 ) ( 6 ) =  210.19502607655667
    d3g/dm3 ( 6 ) ( 4 ) ( 0 ) =  263.10961480748244
    d3g/dm3 ( 6 ) ( 4 ) ( 1 ) =  13.144198129056718
    d3g/dm3 ( 6 ) ( 4 ) ( 2 ) =  45.89789306134633
    d3g/dm3 ( 6 ) ( 4 ) ( 3 ) =  367.7000410437063
    d3g/dm3 ( 6 ) ( 4 ) ( 4 ) =  484.158601968835
    d3g/dm3 ( 6 ) ( 4 ) ( 5 ) =  329.22390408240983
    d3g/dm3 ( 6 ) ( 4 ) ( 6 ) =  265.8672419394313
    d3g/dm3 ( 6 ) ( 5 ) ( 0 ) =  120.47755700895874
    d3g/dm3 ( 6 ) ( 5 ) ( 1 ) =  -137.62370273276946
    d3g/dm3 ( 6 ) ( 5 ) ( 2 ) =  95.84508280580798
    d3g/dm3 ( 6 ) ( 5 ) ( 3 ) =  205.22170293103696
    d3g/dm3 ( 6 ) ( 5 ) ( 4 ) =  267.4334082493392
    d3g/dm3 ( 6 ) ( 5 ) ( 5 ) =  194.52065983676627
    d3g/dm3 ( 6 ) ( 5 ) ( 6 ) =  154.63238457639204
    d3g/dm3 ( 6 ) ( 6 ) ( 0 ) =  89.18516790991765
    d3g/dm3 ( 6 ) ( 6 ) ( 1 ) =  -195.43980298840142
    d3g/dm3 ( 6 ) ( 6 ) ( 2 ) =  101.96784509164031
    d3g/dm3 ( 6 ) ( 6 ) ( 3 ) =  167.74301372336387
    d3g/dm3 ( 6 ) ( 6 ) ( 4 ) =  180.4653558137479
    d3g/dm3 ( 6 ) ( 6 ) ( 5 ) =  131.0209942837789
    d3g/dm3 ( 6 ) ( 6 ) ( 6 ) =  93.35866218052266


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

    da/dm ( 0 ) ( 0 ) =  0.03774009729789971
    da/dm ( 0 ) ( 1 ) =  0.06173127101568136
    da/dm ( 0 ) ( 2 ) =  -0.01415302343164224
    da/dm ( 0 ) ( 3 ) =  -0.06244531531314135
    da/dm ( 0 ) ( 4 ) =  -0.04001082577350318
    da/dm ( 0 ) ( 5 ) =  -0.03203949235646955
    da/dm ( 0 ) ( 6 ) =  -0.06446876054224095
    da/dm ( 1 ) ( 0 ) =  0.04160232164400971
    da/dm ( 1 ) ( 1 ) =  0.18016555983588442
    da/dm ( 1 ) ( 2 ) =  -0.13180404909531004
    da/dm ( 1 ) ( 3 ) =  0.06408838391469956
    da/dm ( 1 ) ( 4 ) =  0.14126616606088868
    da/dm ( 1 ) ( 5 ) =  0.05735671026660135
    da/dm ( 1 ) ( 6 ) =  0.0013855032080308596
    da/dm ( 2 ) ( 0 ) =  -0.005394254142888497
    da/dm ( 2 ) ( 1 ) =  -0.07454157794059005
    da/dm ( 2 ) ( 2 ) =  0.07447943879467711
    da/dm ( 2 ) ( 3 ) =  -0.03895673341101451
    da/dm ( 2 ) ( 4 ) =  -0.0604634009732974
    da/dm ( 2 ) ( 5 ) =  -0.05436958322922292
    da/dm ( 2 ) ( 6 ) =  -0.026021621975631005
    da/dm ( 3 ) ( 0 ) =  -0.09490273333261592
    da/dm ( 3 ) ( 1 ) =  0.1445259490546316
    da/dm ( 3 ) ( 2 ) =  -0.15533853797102762
    da/dm ( 3 ) ( 3 ) =  1.5713752554474445
    da/dm ( 3 ) ( 4 ) =  1.2703857427864635
    da/dm ( 3 ) ( 5 ) =  -0.284149598928495
    da/dm ( 3 ) ( 6 ) =  -0.07177105982487131
    da/dm ( 4 ) ( 0 ) =  -0.0010435320507429205
    da/dm ( 4 ) ( 1 ) =  0.0054670631173746815
    da/dm ( 4 ) ( 2 ) =  -0.004137506332308416
    da/dm ( 4 ) ( 3 ) =  0.021801432507985698
    da/dm ( 4 ) ( 4 ) =  0.16199428887957767
    da/dm ( 4 ) ( 5 ) =  0.1310828503058563
    da/dm ( 4 ) ( 6 ) =  -0.00620727766761529
    da/dm ( 5 ) ( 0 ) =  -0.002458590939767461
    da/dm ( 5 ) ( 1 ) =  0.006530892514270619
    da/dm ( 5 ) ( 2 ) =  -0.010946479707441815
    da/dm ( 5 ) ( 3 ) =  -0.014347255745713168
    da/dm ( 5 ) ( 4 ) =  0.3856721235302447
    da/dm ( 5 ) ( 5 ) =  0.3995169503835242
    da/dm ( 5 ) ( 6 ) =  -0.028385682590988632
    da/dm ( 6 ) ( 0 ) =  -0.004807509619407412
    da/dm ( 6 ) ( 1 ) =  0.00015330842730217542
    da/dm ( 6 ) ( 2 ) =  -0.005091233628486803
    da/dm ( 6 ) ( 3 ) =  -0.0035216106825038686
    da/dm ( 6 ) ( 4 ) =  -0.017747767737367574
    da/dm ( 6 ) ( 5 ) =  -0.027584780526260674
    da/dm ( 6 ) ( 6 ) =  0.08120720870441273


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

    Enthalpy (J) =  -23301155.288087707
    Entropy (J/K) =  3464.1364096335037
    
    ds/dm ( 0 ) =  410.0507455804627
    ds/dm ( 1 ) =  412.5313344657526
    ds/dm ( 2 ) =  450.39944723645357
    ds/dm ( 3 ) =  432.12104337722974
    ds/dm ( 4 ) =  510.7162366847344
    ds/dm ( 5 ) =  492.1167130045531
    ds/dm ( 6 ) =  409.2215956602216
    
    d2s/dm2 ( 0 ) ( 0 ) =  -1.139460183279871
    d2s/dm2 ( 0 ) ( 1 ) =  -0.7132633761330777
    d2s/dm2 ( 0 ) ( 2 ) =  0.34822106132360464
    d2s/dm2 ( 0 ) ( 3 ) =  -0.09353835618269568
    d2s/dm2 ( 0 ) ( 4 ) =  -0.1181719409433798
    d2s/dm2 ( 0 ) ( 5 ) =  0.6935404961682898
    d2s/dm2 ( 0 ) ( 6 ) =  1.0741021359745244
    d2s/dm2 ( 1 ) ( 0 ) =  -0.7132633761330776
    d2s/dm2 ( 1 ) ( 1 ) =  -3.9505233145892458
    d2s/dm2 ( 1 ) ( 2 ) =  1.8295535654490198
    d2s/dm2 ( 1 ) ( 3 ) =  1.3090103211277588
    d2s/dm2 ( 1 ) ( 4 ) =  2.090402979359636
    d2s/dm2 ( 1 ) ( 5 ) =  2.431871612273436
    d2s/dm2 ( 1 ) ( 6 ) =  0.9500748865151925
    d2s/dm2 ( 2 ) ( 0 ) =  0.34822106132360486
    d2s/dm2 ( 2 ) ( 1 ) =  1.8295535654490198
    d2s/dm2 ( 2 ) ( 2 ) =  -1.9390778040162355
    d2s/dm2 ( 2 ) ( 3 ) =  0.4640610956454114
    d2s/dm2 ( 2 ) ( 4 ) =  0.17422867581520995
    d2s/dm2 ( 2 ) ( 5 ) =  0.5798359198817742
    d2s/dm2 ( 2 ) ( 6 ) =  1.0570811341684845
    d2s/dm2 ( 3 ) ( 0 ) =  -0.09353835618269367
    d2s/dm2 ( 3 ) ( 1 ) =  1.3090103211277602
    d2s/dm2 ( 3 ) ( 2 ) =  0.4640610956454124
    d2s/dm2 ( 3 ) ( 3 ) =  -21.625910274128053
    d2s/dm2 ( 3 ) ( 4 ) =  -20.415841356800787
    d2s/dm2 ( 3 ) ( 5 ) =  0.5399382856524635
    d2s/dm2 ( 3 ) ( 6 ) =  2.4513849488571307
    d2s/dm2 ( 4 ) ( 0 ) =  -0.11817194094338022
    d2s/dm2 ( 4 ) ( 1 ) =  2.0904029793596353
    d2s/dm2 ( 4 ) ( 2 ) =  0.1742286758152095
    d2s/dm2 ( 4 ) ( 3 ) =  -20.415841356800783
    d2s/dm2 ( 4 ) ( 4 ) =  -105.09799765973003
    d2s/dm2 ( 4 ) ( 5 ) =  -80.79579774539434
    d2s/dm2 ( 4 ) ( 6 ) =  5.050323751765411
    d2s/dm2 ( 5 ) ( 0 ) =  0.6935404961682854
    d2s/dm2 ( 5 ) ( 1 ) =  2.4318716122734316
    d2s/dm2 ( 5 ) ( 2 ) =  0.5798359198817701
    d2s/dm2 ( 5 ) ( 3 ) =  0.5399382856524583
    d2s/dm2 ( 5 ) ( 4 ) =  -80.79579774539437
    d2s/dm2 ( 5 ) ( 5 ) =  -82.01275177725545
    d2s/dm2 ( 5 ) ( 6 ) =  3.353233020611556
    d2s/dm2 ( 6 ) ( 0 ) =  1.0741021359745242
    d2s/dm2 ( 6 ) ( 1 ) =  0.950074886515192
    d2s/dm2 ( 6 ) ( 2 ) =  1.0570811341684845
    d2s/dm2 ( 6 ) ( 3 ) =  2.4513849488571324
    d2s/dm2 ( 6 ) ( 4 ) =  5.050323751765407
    d2s/dm2 ( 6 ) ( 5 ) =  3.3532330206115497
    d2s/dm2 ( 6 ) ( 6 ) =  -15.903356172714398


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
    print ('dcpdt (J/K^2) = ', obj.getDcpDtFromMolesOfComponents_andT_andP_(m, t, p))
    print ("")
    dcpdm = obj.getDCpDmFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, nc):
        print ('dcp/dm (', i, ') = ', dcpdm.valueAtIndex_(i))


.. parsed-literal::

    Heat capacity (J/K) =  2025.399710653537
    dcpdt (J/K^2) =  0.2917833486953511
    
    dcp/dm ( 0 ) =  249.43464059836083
    dcp/dm ( 1 ) =  257.2017126736729
    dcp/dm ( 2 ) =  250.41000411248376
    dcp/dm ( 3 ) =  248.12793496076208
    dcp/dm ( 4 ) =  253.36431868719188
    dcp/dm ( 5 ) =  247.0084370131884
    dcp/dm ( 6 ) =  243.40306480824367


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

    Volume        53.938741 J/bar          
    dvdt       1.963440e-03 J/bar-K        
    dvdp      -4.679171e-05 J/bar^2        
    d2vdt2     7.538216e-07 J/bar-K^2      
    d2vdtdp   -1.306170e-11 J/bar^2-K      
    d2vdp2     1.516196e-10 J/bar^3        
    
    dv/dm ( 0 ) =  6.7755194155226075
    dv/dm ( 1 ) =  6.43106968036745
    dv/dm ( 2 ) =  6.957252912273878
    dv/dm ( 3 ) =  6.4562808603554265
    dv/dm ( 4 ) =  6.82793556118577
    dv/dm ( 5 ) =  6.827933180819789
    dv/dm ( 6 ) =  6.124224794563601
    
    d2v/dm2 ( 0 ) ( 0 ) =  -0.0012385807484713918
    d2v/dm2 ( 0 ) ( 1 ) =  0.0004735999736928697
    d2v/dm2 ( 0 ) ( 2 ) =  -0.00046187456982968943
    d2v/dm2 ( 0 ) ( 3 ) =  0.0003977407458464053
    d2v/dm2 ( 0 ) ( 4 ) =  0.00019906311283538868
    d2v/dm2 ( 0 ) ( 5 ) =  0.0006964427331411663
    d2v/dm2 ( 0 ) ( 6 ) =  0.0007283112241869619
    d2v/dm2 ( 1 ) ( 0 ) =  0.0004735999736928695
    d2v/dm2 ( 1 ) ( 1 ) =  0.012007219676400182
    d2v/dm2 ( 1 ) ( 2 ) =  -0.005254205916742371
    d2v/dm2 ( 1 ) ( 3 ) =  -0.004077118146881666
    d2v/dm2 ( 1 ) ( 4 ) =  -0.0041732315599436465
    d2v/dm2 ( 1 ) ( 5 ) =  -0.0042341585381709575
    d2v/dm2 ( 1 ) ( 6 ) =  -0.00421424686086206
    d2v/dm2 ( 2 ) ( 0 ) =  -0.0004618745698296892
    d2v/dm2 ( 2 ) ( 1 ) =  -0.005254205916742371
    d2v/dm2 ( 2 ) ( 2 ) =  0.0029631129403313054
    d2v/dm2 ( 2 ) ( 3 ) =  0.0005128950067821089
    d2v/dm2 ( 2 ) ( 4 ) =  0.00019067601668224237
    d2v/dm2 ( 2 ) ( 5 ) =  0.0007970701832875802
    d2v/dm2 ( 2 ) ( 6 ) =  0.0008871712744633012
    d2v/dm2 ( 3 ) ( 0 ) =  0.0003977407458464052
    d2v/dm2 ( 3 ) ( 1 ) =  -0.004077118146881666
    d2v/dm2 ( 3 ) ( 2 ) =  0.0005128950067821085
    d2v/dm2 ( 3 ) ( 3 ) =  0.0026970640393349414
    d2v/dm2 ( 3 ) ( 4 ) =  0.0026257011833396363
    d2v/dm2 ( 3 ) ( 5 ) =  0.0029320856175640408
    d2v/dm2 ( 3 ) ( 6 ) =  0.0029614122846401123
    d2v/dm2 ( 4 ) ( 0 ) =  0.00019906311283538858
    d2v/dm2 ( 4 ) ( 1 ) =  -0.004173231559943646
    d2v/dm2 ( 4 ) ( 2 ) =  0.00019067601668224204
    d2v/dm2 ( 4 ) ( 3 ) =  0.0026257011833396363
    d2v/dm2 ( 4 ) ( 4 ) =  0.0017935901374497658
    d2v/dm2 ( 4 ) ( 5 ) =  0.003121134146018646
    d2v/dm2 ( 4 ) ( 6 ) =  0.0031977991125086195
    d2v/dm2 ( 5 ) ( 0 ) =  0.0006964427331411663
    d2v/dm2 ( 5 ) ( 1 ) =  -0.004234158538170956
    d2v/dm2 ( 5 ) ( 2 ) =  0.0007970701832875798
    d2v/dm2 ( 5 ) ( 3 ) =  0.00293208561756404
    d2v/dm2 ( 5 ) ( 4 ) =  0.003121134146018647
    d2v/dm2 ( 5 ) ( 5 ) =  0.002703269445002389
    d2v/dm2 ( 5 ) ( 6 ) =  0.0025108650676256786
    d2v/dm2 ( 6 ) ( 0 ) =  0.0007283112241869621
    d2v/dm2 ( 6 ) ( 1 ) =  -0.00421424686086206
    d2v/dm2 ( 6 ) ( 2 ) =  0.0008871712744633012
    d2v/dm2 ( 6 ) ( 3 ) =  0.002961412284640112
    d2v/dm2 ( 6 ) ( 4 ) =  0.00319779911250862
    d2v/dm2 ( 6 ) ( 5 ) =  0.002510865067625679
    d2v/dm2 ( 6 ) ( 6 ) =  0.0024769559975446774
    
    d2vdmdt ( 0 ) =  0.00025991217886138084
    d2vdmdt ( 1 ) =  0.0002211892459485272
    d2vdmdt ( 2 ) =  0.0002906842482211079
    d2vdmdt ( 3 ) =  0.0001878798178765582
    d2vdmdt ( 4 ) =  0.00019644509560098943
    d2vdmdt ( 5 ) =  0.00020073471801652963
    d2vdmdt ( 6 ) =  0.00016325409543757537
    
    d2vdmdp ( 0 ) =  -5.751021867420583e-06
    d2vdmdp ( 1 ) =  -4.734136173691597e-06
    d2vdmdp ( 2 ) =  -6.720436349014606e-06
    d2vdmdp ( 3 ) =  -5.502950451885919e-06
    d2vdmdp ( 4 ) =  -5.8222304398726725e-06
    d2vdmdp ( 5 ) =  -5.8177438204639245e-06
    d2vdmdp ( 6 ) =  -5.161483887143959e-06


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

Note that the first nc species are identical to the solution components.

.. code:: ipython3

    import numpy as np
    print ('formula = ', obj.getFormulaFromMolesOfComponents_andT_andP_(m, t, p))
    muSpecies = obj.chemicalPotentialsOfSpeciesFromMolesOfComponents_andT_andP_(m, t, p)
    for i in range (0, ns):
        print ('species = ', obj.nameOfSolutionSpeciesAtIndex_(i))
        elm = obj.elementalCompositionOfSpeciesAtIndex_(i)
        for j in range (0, 107):
            if elm.valueAtIndex_(j) > 0.0:
                print ('   element (', j, ') = ', elm.valueAtIndex_(j))
        print ('   chemical potential = ', muSpecies.valueAtIndex_(i))
    mSpecies = (ctypes.c_double*ns)()
    ctypes.cast(mSpecies, ctypes.POINTER(ctypes.c_double))
    for i in range (0, ns):
        mSpecies[i] = np.random.rand()
        name = obj.nameOfSolutionSpeciesAtIndex_(i)
        print ('Mole fraction of species (', name.ljust(20), ') = ', mSpecies[i])
    mSpToComp = obj.convertMolesOfSpeciesToMolesOfComponents_(mSpecies)
    for i in range (0, nc):
        name = obj.nameOfSolutionSpeciesAtIndex_(i)
        print ('moles of component (', name.ljust(20), ') = ', mSpToComp.valueAtIndex_(i))


.. parsed-literal::

    formula =  Na0.06Ca0.69Fe''0.37Mg0.63Fe'''0.01Ti0.01Al0.42Si1.80O6
    species =  diopside
       element ( 8 ) =  6.0
       element ( 12 ) =  1.0
       element ( 14 ) =  2.0
       element ( 20 ) =  1.0
       chemical potential =  -3445027.629689464
    species =  clinoenstatite
       element ( 8 ) =  6.0
       element ( 12 ) =  2.0
       element ( 14 ) =  2.0
       chemical potential =  -3326569.3414830053
    species =  hedenbergite
       element ( 8 ) =  6.0
       element ( 14 ) =  2.0
       element ( 20 ) =  1.0
       element ( 26 ) =  1.0
       chemical potential =  -3128213.579414352
    species =  alumino-buffonite
       element ( 8 ) =  6.0
       element ( 12 ) =  0.5
       element ( 13 ) =  1.0
       element ( 14 ) =  1.0
       element ( 20 ) =  1.0
       element ( 22 ) =  0.5
       chemical potential =  -3528391.327302648
    species =  buffonite
       element ( 8 ) =  6.0
       element ( 12 ) =  0.5
       element ( 14 ) =  1.0
       element ( 20 ) =  1.0
       element ( 22 ) =  0.5
       element ( 26 ) =  1.0
       chemical potential =  -3149150.7685068524
    species =  essenite
       element ( 8 ) =  6.0
       element ( 13 ) =  1.0
       element ( 14 ) =  1.0
       element ( 20 ) =  1.0
       element ( 26 ) =  1.0
       chemical potential =  -3174325.0543351285
    species =  jadeite
       element ( 8 ) =  6.0
       element ( 11 ) =  1.0
       element ( 13 ) =  1.0
       element ( 14 ) =  2.0
       chemical potential =  -3279477.157716875
    species =  fe-aluminobuffonite
       element ( 8 ) =  6.0
       element ( 13 ) =  1.0
       element ( 14 ) =  1.0
       element ( 20 ) =  1.0
       element ( 22 ) =  0.5
       element ( 26 ) =  0.5
       chemical potential =  -3369984.302165092
    species =  fe-buffonite
       element ( 8 ) =  6.0
       element ( 14 ) =  1.0
       element ( 20 ) =  1.0
       element ( 22 ) =  0.5
       element ( 26 ) =  1.5
       chemical potential =  -2990743.743369296
    species =  ca-tschermaks
       element ( 8 ) =  6.0
       element ( 13 ) =  2.0
       element ( 14 ) =  1.0
       element ( 20 ) =  1.0
       chemical potential =  -3553565.6131309243
    species =  mg-tschermaks
       element ( 8 ) =  6.0
       element ( 12 ) =  1.0
       element ( 13 ) =  2.0
       element ( 14 ) =  1.0
       chemical potential =  -3435107.3249244657
    species =  fe-tschermaks
       element ( 8 ) =  6.0
       element ( 13 ) =  2.0
       element ( 14 ) =  1.0
       element ( 26 ) =  1.0
       chemical potential =  -3118293.2746493537
    species =  acmite
       element ( 8 ) =  6.0
       element ( 11 ) =  1.0
       element ( 14 ) =  2.0
       element ( 26 ) =  1.0
       chemical potential =  -2900236.5989210787
    species =  ferrosilite
       element ( 8 ) =  6.0
       element ( 14 ) =  2.0
       element ( 26 ) =  2.0
       chemical potential =  -2692941.2409327812
    Mole fraction of species ( diopside             ) =  0.5890059304031954
    Mole fraction of species ( clinoenstatite       ) =  0.4255077604633264
    Mole fraction of species ( hedenbergite         ) =  0.6594832970644128
    Mole fraction of species ( alumino-buffonite    ) =  0.6859388889199427
    Mole fraction of species ( buffonite            ) =  0.6773027641140398
    Mole fraction of species ( essenite             ) =  0.7261516233724139
    Mole fraction of species ( jadeite              ) =  0.040846849153800546
    Mole fraction of species ( fe-aluminobuffonite  ) =  0.6130994635571386
    Mole fraction of species ( fe-buffonite         ) =  0.4338562556842901
    Mole fraction of species ( ca-tschermaks        ) =  0.791706498940424
    Mole fraction of species ( mg-tschermaks        ) =  0.2720723724391767
    Mole fraction of species ( fe-tschermaks        ) =  0.7544402417014177
    Mole fraction of species ( acmite               ) =  0.34512250102027053
    Mole fraction of species ( ferrosilite          ) =  0.5341198179795654
    moles of component ( diopside             ) =  -2.7836644210186616
    moles of component ( clinoenstatite       ) =  1.9861401925834863
    moles of component ( hedenbergite         ) =  3.005641034345676
    moles of component ( alumino-buffonite    ) =  2.7721349645378295
    moles of component ( buffonite            ) =  -0.36193759226241795
    moles of component ( essenite             ) =  2.5443707364534323
    moles of component ( jadeite              ) =  0.3859693501740711

