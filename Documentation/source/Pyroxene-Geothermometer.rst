Pyroxene Geothermometer
=======================

Using opx and cpx solution models from Sack and Ghiorso (1994a, b, c):
----------------------------------------------------------------------

| Sack RO, Ghiorso MS (1994a) Thermodynamics of multicomponent
  pyroxenes: I. Formulation of a general model. *Contrib Mineral Petrol*
  116, 277-286
| Sack RO, Ghiorso MS (1994b) Thermodynamics of multicomponent
  pyroxenes: II. Phase relations in the quadrilateral. *Contrib Mineral
  Petrol* 116, 287-300
| Sack RO, Ghiorso MS (1994c) Thermodynamics of multicomponent
  pyroxenes: III. Calibration of Fe2+(Mg)-1, TiAl(MgSi)-1,
  TiFe3+(MgSi)-1, AlFe3+(MgSi)-1, NaAl(CaMg)-1, Al2(MgSi)-1 and Ca(Mg)-1
  exchange reactions between pyroxenes and silicate melts. *Contrib
  Mineral Petrol* 118, 271-296

Initialize some required packages, and load the phase library.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/rubicon/objc/ctypes_patch.py:24: UserWarning: rubicon.objc.ctypes_patch has only been tested with Python 3.4 through 3.6. The current version is sys.version_info(major=3, minor=7, micro=6, releaselevel='final', serial=0). Most likely things will work properly, but you may experience crashes if Python's internals have changed significantly.
      .format(sys.version_info)




.. parsed-literal::

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7f9754c053b0 at 0x7f977852c090>



Define some conversion functions that …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

… take dictionaries of oxide names and oxides values and return
molecular weights and arrays of molar concentrations.

.. code:: ipython3

    def oxide_mw (formulas=["H2O"]):
        result = {}
        PhaseBase = ObjCClass('PhaseBase')
        for formula in formulas:
            obj = PhaseBase.alloc().init()
            obj.setPhaseFormula_(formula)
            result[formula]= obj.mw
        return result
    
    import ctypes
    def oxides_wts_to_element_moles (oxides={"H2O" : 100.0}):
        e = (ctypes.c_double*107)()
        ctypes.cast(e, ctypes.POINTER(ctypes.c_double))
        for i in range (0, 107):
            e[i] = 0.0
        PhaseBase = ObjCClass('PhaseBase')
        for formula, value in oxides.items():
            obj = PhaseBase.alloc().init()
            obj.setPhaseFormula_(formula)
            moles = value/obj.mw
            elements = obj.formulaAsElementArray
            for i in range (0, 107):
                coeff = elements.valueAtIndex_(i)
                if coeff != 0.0:
                    e[i] += coeff*moles
        return e
    
    def element_moles_to_pyx_moles(e):
        m = (ctypes.c_double*nc)()
        ctypes.cast(m, ctypes.POINTER(ctypes.c_double))
        p = 2000.0
        Na = 11
        Mg = 12
        Al = 13
        Si = 14
        Ca = 20
        Ti = 22
        Cr = 24
        Mn = 25
        Fe = 26
        sumcat  = e[Na] +     e[Mg] +     e[Al] +     e[Si] +     e[Ca] +     e[Ti] +     e[Cr] +     e[Mn] + e[Fe]
        sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn]
        if e[Na]+e[Ca] > 0.25*sumcat:
            corrSi = 4.0*(e[Na]+e[Ca]) - sumcat
        else: 
            corrSi = 0.0
        sumcat += corrSi;
    
        # catch low-P oxidized samples and acmites
        if (p < 1000.0) or (e[Na] > e[Al]): 
            fe3 = 3.0*sumcat - sumchg - 2.0*e[Fe]
            fe2 = e[Fe] - fe3
            if fe3 < 0.01*e[Fe]:
                fe3 = 0.01*e[Fe]
                fe2 = 0.99*e[Fe]
            if fe2 < 0.01*e[Fe]:
                fe2 = 0.01*e[Fe]
                fe3 = 0.99*e[Fe]
        else:
            fe2 = e[Fe]
            fe3 = 0.0
    
        m[0] = -fe3/2.0 - fe2 - e[Mn] - e[Al]/2.0 - e[Cr]/2.0 + e[Ca] + e[Na]/2.0 - e[Ti]
        m[1] =  fe3/4.0 + fe2/2.0 + e[Mn]/2.0 + e[Al]/4.0 + e[Cr]/4.0 - e[Ca]/2.0 + e[Mg]/2.0 - e[Na]/4.0
        m[2] =  fe2 + e[Mn]
        m[3] = -fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 + e[Ti]
        m[4] =  fe3/2.0 - e[Al]/2.0 - e[Cr]/2.0 + e[Na]/2.0 + e[Ti]
        m[5] =  fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 - e[Ti]
        m[6] =  e[Na]
        return m

Test the oxide formula to molecular weight method.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print (oxide_mw(["Al2O3", "SiO2"]))


.. parsed-literal::

    {'Al2O3': 101.96127999999999, 'SiO2': 60.0843}


Implement a two-pyroxene geothermometer.
========================================

Reference pyroxene compositions from the Bishop Tuff (Hildreth, 1977).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+--------+-------+-------+--------+-----+-----+-----+-----+-------+
| Phase  | SiO2  | TiO2  | Al2O3  | FeO | MnO | MgO | CaO | Na2O  |
+========+=======+=======+========+=====+=====+=====+=====+=======+
| cpx    | 5     | 0     | 0.150  | 12  | 0.  | 12  | 20  | 0     |
|        | 1.946 | .7111 | 769231 | .80 | 556 | .69 | .63 | .3811 |
|        | 15385 | 53846 |        | 846 | 923 | 576 | 307 | 53846 |
|        |       |       |        | 154 | 077 | 923 | 692 |       |
+--------+-------+-------+--------+-----+-----+-----+-----+-------+
| :m     | 0     | 0     | 0.025  | 0.  | 0.  | 0.  | 0.  | 0     |
| ath:`\ | .3252 | .1598 | 443754 | 305 | 038 | 250 | 200 | .0158 |
| sigma` | 45469 | 08058 |        | 754 | 860 | 984 | 434 | 30837 |
|        |       |       |        | 049 | 698 | 829 | 912 |       |
+--------+-------+-------+--------+-----+-----+-----+-----+-------+
| opx    | 5     | 0     | 0.128  | 28  | 1.  | 18  | 0   | 0     |
|        | 0.929 | .4255 | 888889 | .49 | 103 | .33 | .97 | .0251 |
|        | 25926 | 55556 |        | 518 | 703 | 037 | 962 | 85185 |
|        |       |       |        | 519 | 704 | 037 | 963 |       |
+--------+-------+-------+--------+-----+-----+-----+-----+-------+
| :m     | 0     | 0     | 0.023  | 0.  | 0.  | 0.  | 0.  | 0     |
| ath:`\ | .4603 | .1076 | 912233 | 493 | 045 | 257 | 022 | .0070 |
| sigma` | 25353 | 43762 |        | 233 | 161 | 241 | 951 | 00203 |
|        |       |       |        | 993 | 943 | 005 | 702 |       |
+--------+-------+-------+--------+-----+-----+-----+-----+-------+

Values in wt%. Averages and standard deviations computed from analyzed
pyroxenes found in the late eruptive units.

Instantiate a clinopyroxene with the specified composition.
-----------------------------------------------------------

| As an illustration of use, compute and print properties at 800 °C and
  200 MPa. Properties are output as a Python dictionary.
| Output the number of components, their names, and their formulas.

.. code:: ipython3

    CpxBerman = ObjCClass('CpxBerman')
    cpx = CpxBerman.alloc().init()
    nc = cpx.numberOfSolutionComponents()
    e = oxides_wts_to_element_moles ({'SiO2':51.94615385, 'TiO2':0.711153846, 'Al2O3':0.150769231, 'FeO':12.80846154, 
                                      'MnO':0.556923077, 'MgO':12.69576923, 'CaO':20.63307692, 'Na2O':0.381153846})
    mCpx = element_moles_to_pyx_moles(e)
    
    if (cpx.testPermissibleValuesOfComponents_(mCpx) == 1):
        print ('Cpx composition is feasible')
    else:
        print ('Cpx composition is infeasible')
        
    t = 1073.15 # K
    p = 2000.0  # bars
    potential = cpx.getChemicalPotentialFromMolesOfComponents_andT_andP_(mCpx, t, p)
    
    for i in range (0, nc):
        component = cpx.componentAtIndex_(i)
        print("{0:>20s}{1:15.2f}".format(component.phaseName, potential.valueAtIndex_(i)))


.. parsed-literal::

    Cpx composition is feasible
                diopside    -3466683.65
          clinoenstatite    -3352376.72
            hedenbergite    -3150901.60
       alumino-buffonite    -3575864.14
               buffonite    -3151406.39
                essenite    -3221797.55
                 jadeite    -3336129.22


Instantiate an orthopyroxene with the specified composition.
------------------------------------------------------------

| As an illustration of use, compute and print properties at 800 °C and
  200 MPa. Properties are output as a Python dictionary.
| Output the number of components, their names, and their formulas.

.. code:: ipython3

    OpxBerman = ObjCClass('OpxBerman')
    opx = OpxBerman.alloc().init()
    nc = opx.numberOfSolutionComponents()
    e = oxides_wts_to_element_moles ({'SiO2':50.92925926, 'TiO2':0.425555556, 'Al2O3':0.128888889, 'FeO':28.49518519, 
                                      'MnO':1.103703704, 'MgO':18.33037037, 'CaO':0.97962963, 'Na2O':0.025185185})
    mOpx = element_moles_to_pyx_moles(e)
    
    if (opx.testPermissibleValuesOfComponents_(mOpx) == 1):
        print ('Opx composition is feasible')
    else:
        print ('Opx composition is infeasible')
        
    t = 1073.15 # K
    p = 2000.0  # bars
    potential = opx.getChemicalPotentialFromMolesOfComponents_andT_andP_(mOpx, t, p)
    
    for i in range (0, nc):
        component = opx.componentAtIndex_(i)
        print("{0:>20s}{1:15.2f}".format(component.phaseName, potential.valueAtIndex_(i)))


.. parsed-literal::

    Opx composition is feasible
                diopside    -3470369.82
          clinoenstatite    -3356457.02
            hedenbergite    -3153990.55
       alumino-buffonite    -3556103.19
               buffonite    -3144518.00
                essenite    -3487657.71
                 jadeite    -3582560.72


Define an Fe-Mg exchange reaction between opx and cpx.
------------------------------------------------------

CaMgSi2O6 [cpx] + CaFeSi2O6 [opx] = CaMgSi2O6 [opx] +CaFeSi2O6 [cpx]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the ``get_properties`` function for the class instance returns
a Python dictionary. The chemical potential of the endmember components
are retrieved from this dictionary by using the name of the component as
a key. Otherwise, the other thermodynamic properties are extensive (mass
dependent) quantities and pertain to the phase as a whole.

.. code:: ipython3

    def deltaG(t, p):
        cpxPotentials = cpx.getChemicalPotentialFromMolesOfComponents_andT_andP_(mCpx, t, p)
        opxPotentials = opx.getChemicalPotentialFromMolesOfComponents_andT_andP_(mOpx, t, p)
        return opxPotentials.valueAtIndex_(0) + cpxPotentials.valueAtIndex_(2) - cpxPotentials.valueAtIndex_(0) - opxPotentials.valueAtIndex_(2)

The Gibbs free energy computed by the ``deltaG`` function defined above
must be zero at equilibrium. In order to find this zero, we . . . ## . .
. import a minimizer routine from SciPy called *BrentQ.* We will use
BrentQ to find the temperature that zeroes the Gibbs free energy of a
reaction within a specified range of values.

.. code:: ipython3

    from scipy.optimize import brentq

Solve for the temperature that zeroes the exchange free energy.
---------------------------------------------------------------

Upper and lower bounds on T are specified by Tmin and Tmax (both in K).
The pressure is specified in bars.

.. code:: ipython3

    Tmin = 500.0
    Tmax = 1500.0
    p = 2000.0
    print ('Equilibrium T (°C) = ', brentq(deltaG, Tmin, Tmax, args=(p)) - 273.15)


.. parsed-literal::

    Equilibrium T (°C) =  695.2941439766115


