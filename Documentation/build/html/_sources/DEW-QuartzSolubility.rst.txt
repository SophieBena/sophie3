PhaseObjC - DEW Quartz Solubility
---------------------------------

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

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7fe0a04dcfb0 at 0x7fe0785a8650>



Create a python reference to the DEWFluid class and instantiate an instance of that class.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In Objective-C the code has the form:

::

   obj = [[CpxBerman alloc] init]

and in python:

.. code:: ipython3

    DEWFluid = ObjCClass('DEWFluid')
    obj = DEWFluid.alloc().init()

Create a python reference to the Quartz class and instantiate the phase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Quartz = ObjCClass('QuartzBerman')
    qtz = Quartz.alloc().init()

Obtain properties of the phase inherited from the PhaseBase class.
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

.. code:: ipython3

    print (obj.phaseName)
    print (qtz.phaseName)


.. parsed-literal::

    DEWFluid
    Quartz


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
    print ('Component name'.ljust(20), '\t', 'Formula'.ljust(20), '\t', 'Molecular weight (g/mol)')
    for i in range(0, nc):
        component = obj.componentAtIndex_(i)
        print (component.phaseName.ljust(20), '\t', component.phaseFormula.ljust(20), '\t', component.mw)


.. parsed-literal::

    Component name       	 Formula              	 Molecular weight (g/mol)
    Water                	 H2O                  	 18.0152
    CO2,aq               	 CO2                  	 44.0098
    O2,aq                	 O2                   	 31.9988
    HF,aq                	 HF                   	 20.006303
    NaOH,aq              	 NaOH                 	 39.99707
    Mg(OH)2,aq           	 Mg(OH)2              	 58.319599999999994
    HAlO2,aq             	 HAlO2                	 59.98824
    SiO2,aq              	 SiO2                 	 60.0843
    H3PO4,aq             	 H3PO4                	 97.99506
    SO2,aq               	 SO2                  	 64.0588
    HCl,aq               	 HCl                  	 36.4609
    KOH,aq               	 KOH                  	 56.1093
    Ca(OH)2,aq           	 Ca(OH)2              	 74.0946
    H2CrO4,aq            	 H2CrO4               	 118.0094
    Mn(OH)2,aq           	 Mn(OH)2              	 88.9526
    Fe(OH)2,aq           	 Fe(OH)2              	 89.8616
    Co(OH)2,aq           	 Co(OH)2              	 92.9478


Create a vector of moles of endmember components …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Allocate a “c”-type pointer to a double precision one-dimensional array,
and initialize the array to hold the moles of each component in the
solution

.. code:: ipython3

    import ctypes
    m = (ctypes.c_double*nc)()
    ctypes.cast(m, ctypes.POINTER(ctypes.c_double))
    m[ 0] =  55.55
    m[ 1] =   0.0
    m[ 2] =   0.0
    m[ 3] =   0.0
    m[ 4] =   0.0
    m[ 5] =   0.0
    m[ 6] =   0.0
    m[ 7] =   0.04223 # 0.03 at 400 C and 2 kb, 0.3 at 600 C and 10 kb
    m[ 8] =   0.0
    m[ 9] =   0.0
    m[10] =   0.0
    m[11] =   0.0
    m[12] =   0.0
    m[13] =   0.0
    m[14] =   0.0
    m[15] =   0.0
    m[16] =   0.0
    for i in range (0, nc):
        component = obj.componentAtIndex_(i)
        if m[i] > 0.0:
            print ('moles of (', component.phaseName.ljust(20), ') = ', m[i])


.. parsed-literal::

    moles of ( Water                ) =  55.55
    moles of ( SiO2,aq              ) =  0.04223


Test the mole vector and output derived quantities …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    if (obj.testPermissibleValuesOfComponents_(m) == 1):
        print ('mole input is feasible')
    else:
        print ('mole input is infeasible')
        
    print ('Total moles = ', obj.totalMolesFromMolesOfComponents_(m))
    
    mole_frac_pointer = obj.convertMolesToMoleFractions_(m)
    print ('component name'.ljust(20), 'component formula'.ljust(20), 'mole fraction')
    for i in range (0, nc):
        if m[i] > 0.0:
            print (obj.componentAtIndex_(i).phaseName.ljust(20), obj.componentAtIndex_(i).phaseFormula.ljust(20), mole_frac_pointer.valueAtIndex_(i))
    
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
    Total moles =  55.59223
    component name       component formula    mole fraction
    Water                H2O                  0.9992403614677806
    SiO2,aq              SiO2                 0.0007596385322193406
    Solution formula =  H(111.1)O(55.63446)Si(0.04223)


Change convention for Berman properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| Use Helgeson (SUPCRT) convention of Gibbs free energy of formation
  rather than enthalpy of formation at Tr, Pr
| Do not implement quartz enthalpy correction that is used in
  rhyolite-MELTS

.. code:: ipython3

    Quartz.enableGibbsFreeEnergyReferenceStateUsed()
    Quartz.disableQuartzCorrectionUsed()

Compute activties and chemical potentials of endmember components/species …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    BermanGTrPr = -856288.0
    t = 673.15
    p = 1000.0
    
    activity = obj.getActivityFromMolesOfComponents_andT_andP_(m, t, p)
    potential = obj.getChemicalPotentialFromMolesOfComponents_andT_andP_(m, t, p)
    print ('Gibbs free energy (J) = ', obj.getGibbsFreeEnergyFromMolesOfComponents_andT_andP_(m, t, p))
    print ('component, activity, chemical potential, mu0')
    
    mu0_quartz = qtz.getGibbsFreeEnergyFromT_andP_(t, p)
    
    for i in range (0, nc):
        if m[i] > 0.0:
            component = obj.componentAtIndex_(i)
            g = component.getGibbsFreeEnergyFromT_andP_(t, p)
            print ("{0:>20s} {1:10.6f} {2:15.2f} {3:15.2f}".format(component.phaseName, activity.valueAtIndex_(i), potential.valueAtIndex_(i), g))
    print ("{0:>20s} {1:15.2f}".format('Quartz', mu0_quartz))
    import numpy as np
    muSpecies = obj.chemicalPotentialsOfSpeciesFromMolesOfComponents_andT_andP_(m, t, p)
    print ('species, activity, mu, mu0')
    for i in range (0, ns):
        if muSpecies.valueAtIndex_(i) != 0.0:
            pure = obj.componentAtIndex_(i)
            g = pure.getGibbsFreeEnergyFromT_andP_(t, p)
            print("{0:>20s} {1:10.6e} {2:15.2f} {3:15.2f}".format(obj.nameOfSolutionSpeciesAtIndex_(i), np.exp((muSpecies.valueAtIndex_(i)-g)/8.3143/t), muSpecies.valueAtIndex_(i), g))
    print ('deltaG quartz -> SiO2 (aq)', muSpecies.valueAtIndex_(7)-mu0_quartz)
    print ("at log (molality) = {0:10.6e}".format(np.log10(m[7])))


.. parsed-literal::

    Gibbs free energy (J) =  -15269944.699001541
    component, activity, chemical potential, mu0
                   Water   0.999394      -274219.43      -274216.04
                 SiO2,aq   0.025098      -877459.59      -856835.60
                  Quartz      -878750.13
    species, activity, mu, mu0
                   Water 9.993938e-01      -274219.43      -274216.04
                 SiO2,aq 2.509769e-02      -877459.59      -856835.60
                      H+ 7.652831e-06       -65932.40            0.00
                     OH- 4.452470e-06      -208287.04      -139323.34
                  HSiO3- 3.200361e-06     -1085746.63     -1014934.90
                Si2O4,aq 8.535953e-03     -1754919.18     -1728259.14
    deltaG quartz -> SiO2 (aq) 1290.539609355852
    at log (molality) = -1.374379e+00


.. code:: ipython3

    print ("{0:>10s} {1:10.6f} {2:10.6f} {3:10.6f}".format('G TP TPr TrPr', qtz.getGibbsFreeEnergyFromT_andP_(t, p), qtz.getGibbsFreeEnergyFromT_andP_(t, 1.0), qtz.getGibbsFreeEnergyFromT_andP_(298.15, 1.0)))
    print ("{0:>10s} {1:10.6f} {2:10.6f} {3:10.6f}".format('H TP TPr TrPr', qtz.getEnthalpyFromT_andP_(t, p), qtz.getEnthalpyFromT_andP_(t, 1.0), qtz.getEnthalpyFromT_andP_(298.15, 1.0)))
    print ("{0:>10s} {1:10.6f} {2:10.6f} {3:10.6f}".format('S TP TPr TrPr', qtz.getEntropyFromT_andP_(t, p), qtz.getEntropyFromT_andP_(t, 1.0), qtz.getEntropyFromT_andP_(298.15, 1.0)))
    print ("{0:>10s} {1:10.6f} {2:10.6f} {3:10.6f}".format('C TP TPr TrPr', qtz.getHeatCapacityFromT_andP_(t, p), qtz.getHeatCapacityFromT_andP_(t, 1.0), qtz.getHeatCapacityFromT_andP_(298.15, 1.0)))
    
    HTP = qtz.getEnthalpyFromT_andP_(t, p)
    HTrPr = qtz.getEnthalpyFromT_andP_(298.15, 1.0)
    STP = qtz.getEntropyFromT_andP_(t,p)
    STrPr = qtz.getEntropyFromT_andP_(298.15, 1.0)
    GTP = qtz.getGibbsFreeEnergyFromT_andP_(t, p)
    GTrPr = qtz.getGibbsFreeEnergyFromT_andP_(298.15, 1.0)
    
    print ('G(T,P)-G(Tr,Pr) = ', GTP - GTrPr)
    
    BermanG = -856288
    BermanH = -910700
    
    print ('Berman deltaG-deltaH formation at 298 K and 1 bar', BermanG - BermanH)
    print ('BermanG + integral from PhaseObjC', BermanG+GTP-GTrPr)
    
    print ('Dimitri at 400C and 1KB deltaG quartz = ', -210026.5*4.184)
    print ('Dimitri at 400C and 1KB deltaG SiO2,aq = ', -204755.13*4.184)
    print ('Dimitri at 400C and 1KB deltaG Si2O4 = ', -413054.54*4.184)


.. parsed-literal::

    G TP TPr TrPr -878750.130116 -881047.101142 -856287.786356
    H TP TPr TrPr -820003.153422 -822172.288621 -843929.046342
    S TP TPr TrPr  87.271747  87.461654  41.451417
    C TP TPr TrPr  67.531020  68.196114  44.844355
    G(T,P)-G(Tr,Pr) =  -22462.343760243733
    Berman deltaG-deltaH formation at 298 K and 1 bar 54412
    BermanG + integral from PhaseObjC -878750.3437602438
    Dimitri at 400C and 1KB deltaG quartz =  -878750.876
    Dimitri at 400C and 1KB deltaG SiO2,aq =  -856695.46392
    Dimitri at 400C and 1KB deltaG Si2O4 =  -1728220.19536


