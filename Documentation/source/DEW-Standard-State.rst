PhaseObjC - DEW Standard State properties
=========================================

Required python code to load the PhaseObjC library. The library
libphaseobjc.dylib (see build instructions in README.md) must be
locatable to the system in a standard location (by default
/usr/local/lib)

.. code:: ipython3

    import numpy as np
    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/rubicon/objc/ctypes_patch.py:24: UserWarning: rubicon.objc.ctypes_patch has only been tested with Python 3.4 through 3.6. The current version is sys.version_info(major=3, minor=7, micro=6, releaselevel='final', serial=0). Most likely things will work properly, but you may experience crashes if Python's internals have changed significantly.
      .format(sys.version_info)




.. parsed-literal::

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7fe86bdc1980 at 0x7fe830975610>



.. code:: ipython3

    DEWFluid = ObjCClass('DEWFluid')
    obj = DEWFluid.alloc().init()
    print (obj.phaseName)
    ns = obj.numberOfSolutionSpecies()
    print ('Number of species = ', ns)


.. parsed-literal::

    DEWFluid
    Number of species =  92


Compute standard state chemical potentials of endmember components/species …
----------------------------------------------------------------------------

| The standard state Gibbs free energies of all the species known to my
  implementation of DEW are printed and stopred in a python dictionary,
  called species
| Indoividual species can be accessed by specifying
  e.g. species[‘CO2,aq’]

.. code:: ipython3

    t = 399.85 + 273.15 # K
    p = 23599.0 # bars
    print("{0:>20s} {1:>15s}".format('species', 'mu0'))
    species = {}
    for i in range (0, ns):
        pure = obj.componentAtIndex_(i)
        g = pure.getGibbsFreeEnergyFromT_andP_(t, p)
        print("{0:>20s} {1:15.2f}".format(obj.nameOfSolutionSpeciesAtIndex_(i), g))
        species[obj.nameOfSolutionSpeciesAtIndex_(i)] = g


.. parsed-literal::

                 species             mu0
                   Water      -235328.77
                  CO2,aq      -379028.11
                   O2,aq        11393.60
                   HF,aq      -254324.32
                 NaOH,aq      -420012.64
              Mg(OH)2,aq      -707343.06
                HAlO2,aq      -847967.29
                 SiO2,aq      -814208.46
                H3PO4,aq     -1103926.66
                  SO2,aq      -294911.74
                  HCl,aq      -110085.29
                  KOH,aq      -434338.03
              Ca(OH)2,aq      -822271.04
               H2CrO4,aq      -642448.25
              Mn(OH)2,aq      -502620.19
              Fe(OH)2,aq      -354307.95
              Co(OH)2,aq      -319921.11
                      H+            0.00
                     OH-      -145885.52
                   H2,aq        30308.72
                   CO3-2      -469979.24
                   HCO3-      -558973.70
                   CO,aq      -109680.58
                      F-      -274324.32
                 NaCl,aq      -377240.15
                     Na+      -280827.69
                  NaCO3-      -749930.05
               NaHCO3,aq      -845356.91
              NaHSiO3,aq     -1267235.39
                MgCO3,aq      -903280.40
              Mg(HSiO3)+     -1448610.09
               Mg(HCO3)+     -1009274.60
                    Mg+2      -435572.03
                   MgCl+      -551634.66
                   MgOH+      -606765.02
                MgSO4,aq     -1137741.77
                    Al+3      -423414.42
                   AlO2-      -793847.25
                  HSiO3-      -973875.71
                Si2O4,aq     -1632620.84
                  H2PO4-     -1092847.49
                  HPO4-2     -1043229.94
                   PO4-3      -955730.19
                 H3P2O7-     -2007817.47
                H2P2O7-2     -1956585.87
                  H2S,aq       -14246.36
                     HS-        42970.36
                    S2-2       161146.51
                  S2O3-2      -438902.52
                  S2O4-2      -511494.08
                  S2O5-2      -700442.16
                  S2O6-2      -874015.55
                  S2O8-2      -990899.32
                    S3-2       159583.09
                  S3O6-2      -863971.93
                    S4-2       159205.88
                  S4O6-2      -932750.30
                    S5-2       160124.78
                  S5O6-2      -860623.27
                   SO3-2      -413501.98
                   HSO3-      -519924.47
                   SO4-2      -678674.84
                   HSO4-      -729088.61
                   HSO5-      -621918.51
                     Cl-      -107280.74
                      K+      -291322.11
                  KCl,aq      -386524.80
                   KSO4-      -964999.30
                CaCO3,aq     -1062603.26
               Ca(HCO3)+     -1141421.73
                 Ca(OH)+      -716294.10
                    Ca+2      -550500.01
                   CaCl+      -667477.48
                CaCl2,aq      -757382.71
                CaSO4,aq     -1284398.52
                    Cr+2      -158023.64
                    Cr+3      -177774.85
                 Cr2O7-2     -1240458.10
                  CrO4-2      -662448.25
                  HCrO4-      -748576.33
                    Mn+2      -230849.16
                   MnCl+      -371282.95
                   MnO4-      -449920.93
                  MnO4-2      -415762.76
                MnSO4,aq      -912339.54
                    Fe+2       -82536.92
                    Fe+3         2636.34
                   FeCl+      -198234.95
                  FeCl+2      -150792.01
                FeCl2,aq      -320420.79
                    Co+2       -48150.08
                    Co+3       147941.51


.. code:: ipython3

    print (species['CO2,aq'])


.. parsed-literal::

    -379028.1103565087


