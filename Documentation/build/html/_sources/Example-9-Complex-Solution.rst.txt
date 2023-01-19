Complex Solution SymPy Code Generation
======================================

Solutions with Cation Ordering
------------------------------

Generation of configurational entropy (:math:`{S^{config}}`) using
combinatorics coupled with :math:`n^{th}`-order Taylor expansion of the
non-configurational Gibbs free energy (:math:`{G^*}`). The final
expression for the Gibbs free energy of solution is given by
:math:`G = -{T}{S^{config}} + {G^*}`.

This notebook illustrates construction of this problem using the coder
module of the thermoengine package.

Generally, the Taylor expansion of :math:`{G^*}` is taken to order two
(equivalent to regular solution theory) and cation-ordering between
symmetrically non-equivalent crystallographic sites is assumed to be
non-convergent, i.e. the random ordering state is not acheived at finite
temperature. Alternately, cation-ordering may be modeled as convergent,
inducing a symmetry breaking phase transition at finite temperature,
which necessitates Taylor expansion of :math:`{G^*}` to at least
:math:`4^{th}` order (in ordering parameter) with retention of only even
powers of the ordering variable(s) in the expansion.

This notebook illustrates non-convergent ordering in a reciprocal
solution model for orthpyroxene in the compositional space of the
pyroxene quadrilateral: Mg2Si2O6-Fe2Si2O6-CaMgSi2O6-CaFeSi2O6

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()
    from thermoengine import coder
    from thermoengine import core

Complex Solution Properties - General structure of the model
------------------------------------------------------------

There are three terms: - Terms describing standard state contributions -
Terms describing the configurational entropy of solution - Terms
describing the excess enthalpy of solution

Assumptions: - There are :math:`c` components in the system - There may
be more endmember species, :math:`w`, than there are components, thereby
allowing for reciprocal solutions - Cation ordering is permitted, which
may be either convergent or non-convergent. There may be zero or more
cation ordering variables, :math:`s`. - The configurational entropy
formulation assumes random mixing on symmetrically distinct
crystallographic sites - The excess enthalpy is described using a Taylor
series expansion in compositional and ordering variables. The order of
the expansion is :math:`\nu`.

Number of solution components and number of solution species
------------------------------------------------------------

| Note, that the example illustrated in this notebook - orthopyroxene in
  the system Mg2Si2O6-Fe2Si2O6-CaMgSi2O6-CaFeSi2O6, requires three
  endmember thermodynamic components but clearly has four endmember
  species. This is an example of a recipocal solution. One of the
  species endmembers is compositionally redundent, but *not*
  energetically redundant. Hence the Gibbs free energy change of the
  reaction:
| Mg2Si2O6 + 2 CaFeSi2O6 = Fe2Si2O6 + 2 CaMgSi2O6
| is not zero, even though the concentration of the species Fe2Si2O6 may
  be expressed as:
| 2 CaFeSi2O6 - 2 CaMgSi2O6 + Mg2Si2O6
| Additionally, there is one variable that denotes the degree of cation
  ordering of Fe++ and Mg over the M1 and M2 crystallographic sites in
  the pyroxene structure.

.. code:: ipython3

    nc = 3
    nw = 4
    ns = 1

Create a complex solution model
-------------------------------

| A *complex* solution is one that includes ordering parameters as well
  as endmember thermodynamic components.
| Instantiate the class with the specified number of endmember
  thermodynamic components and species

.. code:: ipython3

    model = coder.ComplexSolnModel(nc=nc, ns=ns, nw=nw)

Retrieve primary compositional variables
----------------------------------------

-  :math:`n` is a vector of mole numbers of each component
-  :math:`n_T` is the total number of moles in the solution
-  :math:`s` is a vector of ordering parameters

and construct a derived mole fraction variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  :math:`X` is a vector of mole fractions of thermodynamic components
   in the system

and a reduced set of independent composition variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  :math:`r` is a vector of independent mole fractions in the system. By
   convention, :math:`r_{i-1}=X_i`, where :math:`i` ranges from the
   second index of :math:`X` up to :math:`c`. Hence the length of the
   vector :math:`r` is :math:`c-1`

.. code:: ipython3

    n = model.n
    nT = model.nT
    s = model.s
    X = n/nT
    r = X[1:]
    n, nT, X, r, s




.. math::

    \displaystyle \left( \left[\begin{matrix}n_{1}\\n_{2}\\n_{3}\end{matrix}\right], \  n_{1} + n_{2} + n_{3}, \  \left[\begin{matrix}\frac{n_{1}}{n_{1} + n_{2} + n_{3}}\\\frac{n_{2}}{n_{1} + n_{2} + n_{3}}\\\frac{n_{3}}{n_{1} + n_{2} + n_{3}}\end{matrix}\right], \  \left[ \frac{n_{2}}{n_{1} + n_{2} + n_{3}}, \  \frac{n_{3}}{n_{1} + n_{2} + n_{3}}\right], \  \left[\begin{matrix}s_{1}\end{matrix}\right]\right)



Retrieve the temperature, pressure, and standard state chemical potentials
--------------------------------------------------------------------------

-  :math:`T` is temperature in :math:`K`
-  :math:`P` is pressure in :math:`bars`
-  :math:`\mu^o` in Joules

.. code:: ipython3

    T = model.get_symbol_for_t()
    P = model.get_symbol_for_p()
    mu = model.mu
    T,P,mu




.. math::

    \displaystyle \left( T, \  P, \  \left[\begin{matrix}\mu_{1}{\left(T,P \right)}\\\mu_{2}{\left(T,P \right)}\\\mu_{3}{\left(T,P \right)}\end{matrix}\right]\right)



Define the standard state contribution to solution properties
-------------------------------------------------------------

.. code:: ipython3

    G_ss = (n.transpose()*mu)[0]
    G_ss




.. math::

    \displaystyle n_{1} \mu_{1}{\left(T,P \right)} + n_{2} \mu_{2}{\left(T,P \right)} + n_{3} \mu_{3}{\left(T,P \right)}



Define configurational entropy and configurational Gibbs free energy
--------------------------------------------------------------------

| Configurational enropy is calculated by counting site configurations,
  that is the number of ways of mixing Fe++, Mg and Ca on the M2 site,
  :math:`\Omega^{M2}`, times the number of ways of mixing Fe++ and Mg on
  the M1 site, :math:`\Omega^{M1}`; configurations (:math:`\Omega`)
  equal :math:`\Omega^{M1}\Omega^{M2}`. The assumption is made that the
  mixing on each site is random, i.e.
| If there are two cations on site M1, and their mole fractions on that
  site are denoted :math:`X` and :math:`Y`, and if there is one such
  sites in the formula unit, then the number of configurations,
  :math:`\Omega`, associted with **random** mixing of cations on that
  site is:
| :math:`\Omega = \left[ {\frac{{\left( {X + Y} \right)!}}{{X!Y!}}} \right]`
| and the molar configurational entropy conribution associated with
  these configurations is given by Boltzmann’s law:
  :math:`{{\hat S}^{conf}} = R\log \Omega`:
| :math:`{{\hat S}^{conf}} = cR\log \left[ {\frac{{\left( {X + Y} \right)!}}{{X!Y!}}} \right]`
| Using Stirlings approximation for large factorials,
  :math:`\log X! = X\log X - X`, the configurational entropy can be
  written:
| :math:`{{\hat S}^{conf}} = cR\left[ - {X\log X - Y\log Y + \left( {X + Y} \right)\log \left( {X + Y} \right)} \right]`

| Consequently, to utilize this appropach we must define site mole
  fractions in terms of our chosen set of independent compositional
  variables and ordering parameters. #### There are 5 site mole
  fractions:
| :math:`X_{Ca}^{M2}`, :math:`X_{Mg}^{M2}`, :math:`X_{{Fe}^{2+}}^{M2}`,
  :math:`X_{Mg}^{M1}`, :math:`X_{{Fe}^{2+}}^{M1}`
| #### The requirement of filled sites requires: 1. :math:`X_{Ca}^{M2}`
  + :math:`X_{Mg}^{M2}` + :math:`X_{{Fe}^{2+}}^{M2}` = 1 2.
  :math:`X_{Mg}^{M1}` + :math:`X_{{Fe}^{2+}}^{M1}` =1

Asuuming the endmembers are ordered as:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  :math:`n_1`, :math:`X_1`, CaMgSi2O6
-  :math:`n_2`, :math:`X_2`, CaFeSi2O6
-  :math:`n_3`, :math:`X_3`, Mg2Si2O6

There are two independent compositional variables:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  CaMgSi2O6 = :math:`1-r_1-r_2`
-  CaFeSi2O6 = :math:`r_1`
-  Mg2Si2O6 = :math:`r_2`

The requirement of mass balance requires:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

3. :math:`r_1` = :math:`X_{{Fe}^{2+}}^{M2}` + :math:`X_{{Fe}^{2+}}^{M1}`
4. :math:`r_2` = 1 - :math:`X_{Ca}^{M2}`

There is one ordering parameter:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

5. :math:`s_1` = :math:`X_{{Fe}^{2+}}^{M2}` - :math:`X_{Mg}^{M2}`

Relations 1-5 may be solved simultaneously to give the following site mole fraction definitions:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  :math:`X_{Ca}^{M2}` = :math:`1-r_2`
-  :math:`X_{Mg}^{M2}` = :math:`\frac{r_2-s_1}{2}`
-  :math:`X_{{Fe}^{2+}}^{M2}` = :math:`\frac{r_2+s_1}{2}`
-  :math:`X_{Mg}^{M1}` = :math:`1-r_1+\frac{r_2+s_1}{2}`
-  :math:`X_{{Fe}^{2+}}^{M1}` = :math:`r_1-\frac{r_2+s_1}{2}`

While this system is fairly easy to solve by inspection, for more
complex situations, assemble the relations in a list of equations that
evaluate to zero, and automatiocally solve that system of equations
using the sympy routine linsolve, i.e. 

::

   system = [xCaM2 + xMgM2 + xFeM2 - 1, xMgM1 + xFeM1 - 1, xFeM2 + xFeM1 - r[0], 1 - xCaM2 - r[1],
             xFeM2 - xMgM2 - s[0]]
   ans = sym.linsolve(system, xCaM2, xMgM2, xFeM2, xMgM1, xFeM1)

.. code:: ipython3

    xCaM2 = 1 - r[1]
    xMgM2 = (r[1]-s[0])/2
    xFeM2 = (r[1]+s[0])/2
    xMgM1 = 1 - r[0] + (r[1]+s[0])/2
    xFeM1 = r[0] - (r[1]+s[0])/2

The following functions implement random mixing configurational entropy on the M1 and M2 sites:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def Sconf_M1_random(X, Y):
        A = X*sym.log(X) - X
        B = Y*sym.log(Y) - Y
        ApB = (X+Y)*sym.log(X+Y) - (X+Y)
        return ApB - A - B
    def Sconf_M2_random(X, Y, Z):
        A = X*sym.log(X) - X
        B = Y*sym.log(Y) - Y
        C = Z*sym.log(Z) - Z
        ApBpC = (X+Y+Z)*sym.log(X+Y+Z) - (X+Y+Z)
        return ApBpC - A - B - C

Configurational entropy
-----------------------

:math:`R` is the gas constant

.. code:: ipython3

    R = sym.symbols('R')
    S_config = Sconf_M1_random(xMgM1, xFeM1) + Sconf_M2_random(xCaM2, xMgM2, xFeM2)
    S_config *= R*nT
    S_config




.. math::

    \displaystyle R \left(n_{1} + n_{2} + n_{3}\right) \left(- \left(- \frac{n_{3}}{n_{1} + n_{2} + n_{3}} + 1\right) \log{\left(- \frac{n_{3}}{n_{1} + n_{2} + n_{3}} + 1 \right)} - \left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2}\right) \log{\left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2} \right)} - \left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2}\right) \log{\left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} \right)} - \left(\frac{n_{2}}{n_{1} + n_{2} + n_{3}} - \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2}\right) \log{\left(\frac{n_{2}}{n_{1} + n_{2} + n_{3}} - \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2} \right)} - \left(- \frac{n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} + 1\right) \log{\left(- \frac{n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} + 1 \right)}\right)



Configurational Gibbs free energy
---------------------------------

Note that this quantity is extensive, with units of J, *not J/mole*

.. code:: ipython3

    G_config = -T*S_config
    G_config




.. math::

    \displaystyle - R T \left(n_{1} + n_{2} + n_{3}\right) \left(- \left(- \frac{n_{3}}{n_{1} + n_{2} + n_{3}} + 1\right) \log{\left(- \frac{n_{3}}{n_{1} + n_{2} + n_{3}} + 1 \right)} - \left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2}\right) \log{\left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2} \right)} - \left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2}\right) \log{\left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} \right)} - \left(\frac{n_{2}}{n_{1} + n_{2} + n_{3}} - \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2}\right) \log{\left(\frac{n_{2}}{n_{1} + n_{2} + n_{3}} - \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2} \right)} - \left(- \frac{n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} + 1\right) \log{\left(- \frac{n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} + 1 \right)}\right)



:math:`\hat G^*` - Non-configurational molar Gibbs free energy
--------------------------------------------------------------

:math:`\hat G^*` includes all standard state and excess Gibbs free
energy contributions. It is generally modeled as a Taylor expansion in
composition (:math:`r`) and ordering (:math:`s`) variables of order 2, 3
or 4. Here, we choose a model of order two. #### Taylor expansion of
:math:`\hat G^*` For a second order expansion, the number of Taylor
expansion coefficients is: - 1 for :math:`G_{0}` - nc-1 for
:math:`G_{r_i}`, :math:`i=1...nc-1` - ns for :math:`G_{s_i}`,
:math:`i=1...ns` - (nc-1)(nc-2)/2 for :math:`G_{{r_i},{r_{i+1}}}`,
:math:`i= 1...nc-2` - ns(ns-1)/2 for :math:`G_{{s_i},{s_{i+1}}}`,
:math:`i= 1...ns-1` - ns(nc-1) for :math:`G_{{r_i},{s_j}}`,
:math:`i= 1...nc-1`, :math:`j=1...ns` - nc-1 for
:math:`G_{{r_i},{r_i}}`, :math:`i= 1...nc-1` - ns for
:math:`G_{{s_i},{s_i}}`, :math:`i= 1...ns`

.. code:: ipython3

    (count, taylor, taylor_coeff, taylor_terms) = model.taylor_expansion()
    print ('Number of Taylor expansion terms = {0:.0f}'.format(count))
    taylor


.. parsed-literal::

    Number of Taylor expansion terms = 10




.. math::

    \displaystyle G_{0} + \frac{Gr_{1} n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{Gr_{2} n_{3}}{n_{1} + n_{2} + n_{3}} + \frac{Grr_{11} n_{2}^{2}}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{Grr_{12} n_{2} n_{3}}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{Grr_{22} n_{3}^{2}}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{Grs_{11} n_{2} s_{1}}{n_{1} + n_{2} + n_{3}} + \frac{Grs_{21} n_{3} s_{1}}{n_{1} + n_{2} + n_{3}} + Gs_{1} s_{1} + Gss_{11} s_{1}^{2}



Identify Taylor terms of :math:`\hat{G}^*` corresponding to component endmembers:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. diopside, CaMgSi2O6
2. hedenbergite, CaFeSi2O6
3. enstatite, Mg2Si2O6

.. code:: ipython3

    eqn1 = model.eval_endmember([1,0,0],[0],taylor) - mu[0]
    eqn2 = model.eval_endmember([0,1,0],[0],taylor) - mu[1]
    eqn3 = model.eval_endmember([0,0,1],[-1],taylor) - mu[2]
    params = []
    units = []
    symparams = []
    eqn1, eqn2, eqn3




.. math::

    \displaystyle \left( G_{0} - \mu_{1}{\left(T,P \right)}, \  G_{0} + Gr_{1} + Grr_{11} - \mu_{2}{\left(T,P \right)}, \  G_{0} + Gr_{2} + Grr_{22} - Grs_{21} - Gs_{1} + Gss_{11} - \mu_{3}{\left(T,P \right)}\right)



Identify Taylor terms of :math:`\hat{G}^*` corresponding to dependent species endmembers:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ferrosilite, Fe2Si2O6

.. code:: ipython3

    gFs = model.eval_endmember([-2,2,1],[1],taylor)
    gFs




.. math::

    \displaystyle G_{0} + 2 Gr_{1} + Gr_{2} + 4 Grr_{11} + 2 Grr_{12} + Grr_{22} + 2 Grs_{11} + Grs_{21} + Gs_{1} + Gss_{11}



Identify the free energy of the reciprocal reaction between endmember species:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Mg2Si2O6 + 2 CaFeSi2O6 = Fe2Si2O6 + 2 CaMgSi2O6
| is defined as the “reciprocal energy,” :math:`F`, denoting the
  non-co-planarity of the non-configurational Gibbs free energy of the
  endmember species. In general, all reciprocal solutions have non-zero
  :math:`F`. In the paper on pyroxene thermodynamics by Sack and Ghiorso
  (Contributions to Mineralogy and Petrology 116: 277-286, 1994)
  :math:`F` is notated as :math:`\Delta \bar G_{27}^o`

:math:`F` is conveniently defined in terms of expressions 1-4:

.. code:: ipython3

    Fh,Fs,Fv = sym.symbols('Fh Fs Fv')
    params.append('Fh')
    units.append('J/mol')
    symparams.append(Fh)
    params.append('Fs')
    units.append('J/K-mol')
    symparams.append(Fs)
    params.append('Fv')
    units.append('J/bar-mol')
    symparams.append(Fv)
    F = Fh - T*Fs + P*Fv
    eqn4  =   model.eval_endmember([-2,2,1],[ 1],taylor)
    eqn4 += 2*model.eval_endmember([ 1,0,0],[ 0],taylor)
    eqn4 -=   model.eval_endmember([ 0,0,1],[-1],taylor)
    eqn4 -= 2*model.eval_endmember([ 0,1,0],[ 0],taylor)
    eqn4 -= F
    eqn4




.. math::

    \displaystyle - Fh + Fs T - Fv P + 2 Grr_{11} + 2 Grr_{12} + 2 Grs_{11} + 2 Grs_{21} + 2 Gs_{1}



Identify the free energy of the ordering reaction:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| MgFeSi2O6 = FeMgSi2O6
| which will be notated as :math:`Gex`, (in Sack and Ghiorso, 1994,
  :math:`\Delta \bar G_{EX}^o`)
| Note that both compositions, MgFeSi2O6 and FeMgSi2O6, are equivalent
  and defined by CaFeSi2O6 - CaMgSi2O6 + Mg2Si2O6. They differ only by
  the sign of the ordering parameter.

.. code:: ipython3

    Hex,Vex = sym.symbols('Hex Vex')
    params.append('Hex')
    units.append('J/mol')
    symparams.append(Hex)
    params.append('Vex')
    units.append('J/bar-mol')
    symparams.append(Vex)
    Gex = Hex + P*Vex
    eqn5  = model.eval_endmember([-1,1,1],[ 1],taylor)
    eqn5 -= model.eval_endmember([-1,1,1],[-1],taylor)
    eqn5 -= Gex
    eqn5




.. math::

    \displaystyle 2 Grs_{11} + 2 Grs_{21} + 2 Gs_{1} - Hex - P Vex



Identify the free energy of the reciprocal ordering reaction:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Mg2Si2O6 + Fe2Si2O6 = MgFeSi2O6 + FeMgSi2O6
| which will be notated as :math:`Gx`, (in Sack and Ghiorso, 1994,
  :math:`\Delta \bar G_{X}^o`)

.. code:: ipython3

    Hx,Vx = sym.symbols('Hx Vx')
    params.append('Hx')
    units.append('J/mol')
    symparams.append(Hx)
    params.append('Vx')
    units.append('J/bar-mol')
    symparams.append(Vx)
    Gx = Hx + P*Vx
    eqn6  = model.eval_endmember([-1,1,1],[ 1],taylor)
    eqn6 += model.eval_endmember([-1,1,1],[-1],taylor)
    eqn6 -= model.eval_endmember([ 0,0,1],[-1],taylor)
    eqn6 -= model.eval_endmember([-2,2,1],[ 1],taylor)
    eqn6 -= Gx
    eqn6




.. math::

    \displaystyle - 2 Grr_{11} - 2 Grs_{11} - Hx - P Vx



Identify regular solution interaction parameters:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Ca-Mg interaction on the M2 site, the join Mg2Si2O6 - CaMgSi2O6,
   denoted WM2CaMg
-  Ca-Fe interaction on the M2 site, the join Fe2Si2O6 - CaFeSi2O6,
   denoted WM2CaFe
-  Fe-Mg interaction on the M1 site, the joins Mg2Si2O6 - MgFeSi2O6 or
   Fe2Si2O6 - FeMgSi2O6 or CaMgSi2O6 - CaFeSi2O6, which are assumed to
   be energetically equivalent, denoted WM1FeMg (in Sack and Ghiorso,
   1994, :math:`W_{12}`
-  Fe-Mg interaction on the M2 site, the joins FeMgSi2O6 - Mg2Si2O6 or
   Fe2Si2O6 - MgFeSi2O6, which are assumed to be energetically
   equivalent, denoted WM2FeMg

| Along the A-B join, described using a regular solution parameter,
  :math:`W`, :math:`\hat G^*` is given by
| :math:`{\hat G}^*(X_A,X_B)={X_A}{\hat G}^*(A)+{X_B}{\hat G}^*(B)+W{X_A}{X_B}`,
  so
| :math:`W = \frac{{\hat G}^*(X_A,X_B) - {X_A}{\hat G}^*(A) - {X_B}{\hat G}^*(B)}{{X_A}{X_B}}`
| Taking the midpoint of the join provides a way to define the
  parameter:
| :math:`W = \frac{{\hat G}^*(\frac{1}{2},\frac{1}{2}) - {\frac{1}{2}}{\hat G}^*(A) - {\frac{1}{2}}{\hat G}^*(B)}{{\frac{1}{2}}{\frac{1}{2}}} = 4{\hat G}^*(\frac{1}{2},\frac{1}{2}) - 2{\hat G}^*(A) - 2{\hat G}^*(B)`

-  Ca-Mg interaction on the M2 site, the join Mg2Si2O6 - CaMgSi2O6,
   denoted WM2CaMg

.. code:: ipython3

    WhM2CaMg,WvM2CaMg = sym.symbols('WhM2CaMg WvM2CaMg')
    params.append('WhM2CaMg')
    units.append('J/mol')
    symparams.append(WhM2CaMg)
    params.append('WvM2CaMg')
    units.append('J/bar-mol')
    symparams.append(WvM2CaMg)
    WM2CaMg = WhM2CaMg + P*WvM2CaMg
    eqn7 = model.eval_regular_param([1,0,0],[0],[0,0,1],[-1],taylor) - WM2CaMg
    eqn7




.. math::

    \displaystyle - Grr_{22} + Grs_{21} - Gss_{11} - P WvM2CaMg - WhM2CaMg



-  Ca-Fe interaction on the M2 site, the join Fe2Si2O6 - CaFeSi2O6,
   denoted WM2CaFe

.. code:: ipython3

    WhM2CaFe,WvM2CaFe = sym.symbols('WhM2CFe WvM2CaFe')
    params.append('WhM2CaFe')
    units.append('J/mol')
    symparams.append(WhM2CaFe)
    params.append('WvM2CaFe')
    units.append('J/bar-mol')
    symparams.append(WvM2CaFe)
    WM2CaFe = WhM2CaFe + P*WvM2CaFe
    eqn8 = model.eval_regular_param([0,1,0],[0],[-2,2,1],[1],taylor) - WM2CaFe
    eqn8




.. math::

    \displaystyle - Grr_{11} - Grr_{12} - Grr_{22} - Grs_{11} - Grs_{21} - Gss_{11} - P WvM2CaFe - WhM2CFe



-  Fe-Mg interaction on the M1 site, the joins Mg2Si2O6 - MgFeSi2O6 or
   Fe2Si2O6 - FeMgSi2O6 or CaMgSi2O6 - CaFeSi2O6, which are assumed to
   be energetically equivalent, denoted WM1FeMg (in Sack and Ghiorso,
   1994, :math:`W_{12}`

.. code:: ipython3

    WhM1FeMg,WvM1FeMg = sym.symbols('WhM1FeMg WvM1FeMg')
    params.append('WhM1FeMg')
    units.append('J/mol')
    symparams.append(WhM1FeMg)
    params.append('WvM1FeMg')
    units.append('J/bar-mol')
    symparams.append(WvM1FeMg)
    WM1FeMg = WhM1FeMg + P*WvM1FeMg
    eqn9 = model.eval_regular_param([1,0,0],[0],[0,1,0],[0],taylor) - WM1FeMg
    eqn9




.. math::

    \displaystyle - Grr_{11} - P WvM1FeMg - WhM1FeMg



-  Fe-Mg interaction on the M2 site, the joins FeMgSi2O6 - Mg2Si2O6 or
   Fe2Si2O6 - MgFeSi2O6, which are assumed to be energetically
   equivalent, denoted WM2FeMg

.. code:: ipython3

    WhM2FeMg,WvM2FeMg = sym.symbols('WhM2FeMg WvM2FeMg')
    params.append('WhM2FeMg')
    units.append('J/mol')
    symparams.append(WhM2FeMg)
    params.append('WvM2FeMg')
    units.append('J/bar-mol')
    symparams.append(WvM2FeMg)
    WM2FeMg = WhM2FeMg + P*WvM2FeMg
    eqn10 = model.eval_regular_param([-1,1,1],[1],[0,0,1],[-1],taylor) - WM2FeMg
    eqn10




.. math::

    \displaystyle - Grr_{11} - 2 Grs_{11} - 4 Gss_{11} - P WvM2FeMg - WhM2FeMg



Solve for the Taylor coefficients in terms of the preferred parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    system = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10]
    system




.. math::

    \displaystyle \left[ G_{0} - \mu_{1}{\left(T,P \right)}, \  G_{0} + Gr_{1} + Grr_{11} - \mu_{2}{\left(T,P \right)}, \  G_{0} + Gr_{2} + Grr_{22} - Grs_{21} - Gs_{1} + Gss_{11} - \mu_{3}{\left(T,P \right)}, \  - Fh + Fs T - Fv P + 2 Grr_{11} + 2 Grr_{12} + 2 Grs_{11} + 2 Grs_{21} + 2 Gs_{1}, \  2 Grs_{11} + 2 Grs_{21} + 2 Gs_{1} - Hex - P Vex, \  - 2 Grr_{11} - 2 Grs_{11} - Hx - P Vx, \  - Grr_{22} + Grs_{21} - Gss_{11} - P WvM2CaMg - WhM2CaMg, \  - Grr_{11} - Grr_{12} - Grr_{22} - Grs_{11} - Grs_{21} - Gss_{11} - P WvM2CaFe - WhM2CFe, \  - Grr_{11} - P WvM1FeMg - WhM1FeMg, \  - Grr_{11} - 2 Grs_{11} - 4 Gss_{11} - P WvM2FeMg - WhM2FeMg\right]



.. code:: ipython3

    taylor_soln = sym.linsolve(system, taylor_coeff).args[0]
    taylor_soln




.. math::

    \displaystyle \left( \mu_{1}{\left(T,P \right)}, \  P WvM1FeMg + WhM1FeMg - \mu_{1}{\left(T,P \right)} + \mu_{2}{\left(T,P \right)}, \  \frac{Fh}{4} - \frac{Fs T}{4} + \frac{Fv P}{4} + \frac{Hex}{4} + \frac{Hx}{4} + \frac{P Vex}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{2} + \frac{P WvM2CaFe}{2} + \frac{P WvM2CaMg}{2} - \frac{WhM1FeMg}{2} + \frac{WhM2CFe}{2} + \frac{WhM2CaMg}{2} - \mu_{1}{\left(T,P \right)} + \mu_{3}{\left(T,P \right)}, \  \frac{Fh}{4} - \frac{Fs T}{4} + \frac{Fv P}{4} + \frac{Hex}{4} + \frac{Hx}{4} + \frac{P Vex}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{2} + \frac{P WvM2CaFe}{2} - \frac{P WvM2CaMg}{2} - \frac{WhM1FeMg}{2} + \frac{WhM2CFe}{2} - \frac{WhM2CaMg}{2}, \  - P WvM1FeMg - WhM1FeMg, \  \frac{Fh}{2} - \frac{Fs T}{2} + \frac{Fv P}{2} - \frac{Hex}{2} - \frac{P Vex}{2} + P WvM1FeMg + WhM1FeMg, \  - \frac{Hx}{2} - \frac{P Vx}{2} + P WvM1FeMg + WhM1FeMg, \  - \frac{Fh}{4} + \frac{Fs T}{4} - \frac{Fv P}{4} + \frac{Hex}{4} + \frac{P Vex}{4} - \frac{P WvM1FeMg}{4} - \frac{P WvM2CaFe}{2} - \frac{P WvM2CaMg}{2} + \frac{P WvM2FeMg}{4} - \frac{WhM1FeMg}{4} - \frac{WhM2CFe}{2} - \frac{WhM2CaMg}{2} + \frac{WhM2FeMg}{4}, \  - \frac{Fh}{4} + \frac{Fs T}{4} - \frac{Fv P}{4} + \frac{Hex}{4} + \frac{Hx}{4} + \frac{P Vex}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{2} - \frac{P WvM2CaFe}{2} + \frac{P WvM2CaMg}{2} - \frac{WhM1FeMg}{2} - \frac{WhM2CFe}{2} + \frac{WhM2CaMg}{2}, \  \frac{Hx}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{4} - \frac{P WvM2FeMg}{4} - \frac{WhM1FeMg}{4} - \frac{WhM2FeMg}{4}\right)



Substitute terms into :math:`\hat G^*`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    sub_list = []
    for a,b in zip(taylor_coeff,taylor_soln):
        sub_list.append((a,b))
    G_star_molar = taylor.subs(sub_list)

.. code:: ipython3

    print(params)
    print(units)
    print(symparams)


.. parsed-literal::

    ['Fh', 'Fs', 'Fv', 'Hex', 'Vex', 'Hx', 'Vx', 'WhM2CaMg', 'WvM2CaMg', 'WhM2CaFe', 'WvM2CaFe', 'WhM1FeMg', 'WvM1FeMg', 'WhM2FeMg', 'WvM2FeMg']
    ['J/mol', 'J/K-mol', 'J/bar-mol', 'J/mol', 'J/bar-mol', 'J/mol', 'J/bar-mol', 'J/mol', 'J/bar-mol', 'J/mol', 'J/bar-mol', 'J/mol', 'J/bar-mol', 'J/mol', 'J/bar-mol']
    [Fh, Fs, Fv, Hex, Vex, Hx, Vx, WhM2CaMg, WvM2CaMg, WhM2CFe, WvM2CaFe, WhM1FeMg, WvM1FeMg, WhM2FeMg, WvM2FeMg]


Define the Gibbs free energy of solution
----------------------------------------

.. code:: ipython3

    G = G_config + nT*G_star_molar
    G




.. math::

    \displaystyle - R T \left(n_{1} + n_{2} + n_{3}\right) \left(- \left(- \frac{n_{3}}{n_{1} + n_{2} + n_{3}} + 1\right) \log{\left(- \frac{n_{3}}{n_{1} + n_{2} + n_{3}} + 1 \right)} - \left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2}\right) \log{\left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2} \right)} - \left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2}\right) \log{\left(\frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} \right)} - \left(\frac{n_{2}}{n_{1} + n_{2} + n_{3}} - \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2}\right) \log{\left(\frac{n_{2}}{n_{1} + n_{2} + n_{3}} - \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} - \frac{s_{1}}{2} \right)} - \left(- \frac{n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} + 1\right) \log{\left(- \frac{n_{2}}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}}{2 \left(n_{1} + n_{2} + n_{3}\right)} + \frac{s_{1}}{2} + 1 \right)}\right) + \left(n_{1} + n_{2} + n_{3}\right) \left(\frac{n_{2}^{2} \left(- P WvM1FeMg - WhM1FeMg\right)}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{n_{2} n_{3} \left(\frac{Fh}{2} - \frac{Fs T}{2} + \frac{Fv P}{2} - \frac{Hex}{2} - \frac{P Vex}{2} + P WvM1FeMg + WhM1FeMg\right)}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{n_{2} s_{1} \left(- \frac{Hx}{2} - \frac{P Vx}{2} + P WvM1FeMg + WhM1FeMg\right)}{n_{1} + n_{2} + n_{3}} + \frac{n_{2} \left(P WvM1FeMg + WhM1FeMg - \mu_{1}{\left(T,P \right)} + \mu_{2}{\left(T,P \right)}\right)}{n_{1} + n_{2} + n_{3}} + \frac{n_{3}^{2} \left(- \frac{Fh}{4} + \frac{Fs T}{4} - \frac{Fv P}{4} + \frac{Hex}{4} + \frac{P Vex}{4} - \frac{P WvM1FeMg}{4} - \frac{P WvM2CaFe}{2} - \frac{P WvM2CaMg}{2} + \frac{P WvM2FeMg}{4} - \frac{WhM1FeMg}{4} - \frac{WhM2CFe}{2} - \frac{WhM2CaMg}{2} + \frac{WhM2FeMg}{4}\right)}{\left(n_{1} + n_{2} + n_{3}\right)^{2}} + \frac{n_{3} s_{1} \left(- \frac{Fh}{4} + \frac{Fs T}{4} - \frac{Fv P}{4} + \frac{Hex}{4} + \frac{Hx}{4} + \frac{P Vex}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{2} - \frac{P WvM2CaFe}{2} + \frac{P WvM2CaMg}{2} - \frac{WhM1FeMg}{2} - \frac{WhM2CFe}{2} + \frac{WhM2CaMg}{2}\right)}{n_{1} + n_{2} + n_{3}} + \frac{n_{3} \left(\frac{Fh}{4} - \frac{Fs T}{4} + \frac{Fv P}{4} + \frac{Hex}{4} + \frac{Hx}{4} + \frac{P Vex}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{2} + \frac{P WvM2CaFe}{2} + \frac{P WvM2CaMg}{2} - \frac{WhM1FeMg}{2} + \frac{WhM2CFe}{2} + \frac{WhM2CaMg}{2} - \mu_{1}{\left(T,P \right)} + \mu_{3}{\left(T,P \right)}\right)}{n_{1} + n_{2} + n_{3}} + s_{1}^{2} \left(\frac{Hx}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{4} - \frac{P WvM2FeMg}{4} - \frac{WhM1FeMg}{4} - \frac{WhM2FeMg}{4}\right) + s_{1} \left(\frac{Fh}{4} - \frac{Fs T}{4} + \frac{Fv P}{4} + \frac{Hex}{4} + \frac{Hx}{4} + \frac{P Vex}{4} + \frac{P Vx}{4} - \frac{P WvM1FeMg}{2} + \frac{P WvM2CaFe}{2} - \frac{P WvM2CaMg}{2} - \frac{WhM1FeMg}{2} + \frac{WhM2CFe}{2} - \frac{WhM2CaMg}{2}\right) + \mu_{1}{\left(T,P \right)}\right)



Find the condition of homogeneous equilibrium:
----------------------------------------------

:math:`\frac{{\partial \hat G*}}{{\partial {s_1}}} = 0`

.. code:: ipython3

    dgds = (nT*G_star_molar+G_config).diff(s[0]).simplify()
    dgds




.. math::

    \displaystyle - \frac{R T \left(n_{1} + n_{2} + n_{3}\right) \left(\log{\left(\frac{n_{3} - s_{1} \left(n_{1} + n_{2} + n_{3}\right)}{n_{1} + n_{2} + n_{3}} \right)} - \log{\left(\frac{n_{3} + s_{1} \left(n_{1} + n_{2} + n_{3}\right)}{n_{1} + n_{2} + n_{3}} \right)} + \log{\left(- \frac{- 2 n_{2} + n_{3} + s_{1} \left(n_{1} + n_{2} + n_{3}\right)}{n_{1} + n_{2} + n_{3}} \right)} - \log{\left(\frac{- 2 n_{2} + n_{3} + \left(s_{1} + 2\right) \left(n_{1} + n_{2} + n_{3}\right)}{n_{1} + n_{2} + n_{3}} \right)}\right)}{2} - \frac{n_{2} \left(Hx + P Vx - 2 P WvM1FeMg - 2 WhM1FeMg\right)}{2} + \frac{n_{3} \left(- Fh + Fs T - Fv P + Hex + Hx + P Vex + P Vx - 2 P WvM1FeMg - 2 P WvM2CaFe + 2 P WvM2CaMg - 2 WhM1FeMg - 2 WhM2CFe + 2 WhM2CaMg\right)}{4} + \frac{\left(n_{1} + n_{2} + n_{3}\right) \left(Fh - Fs T + Fv P + Hex + Hx + P Vex + P Vx - 2 P WvM1FeMg + 2 P WvM2CaFe - 2 P WvM2CaMg - 2 WhM1FeMg + 2 WhM2CFe - 2 WhM2CaMg - 2 s_{1} \left(- Hx - P Vx + P WvM1FeMg + P WvM2FeMg + WhM1FeMg + WhM2FeMg\right)\right)}{4}



Identify bounds on the ordering parameter
-----------------------------------------

The code generated to implement this model must compute numerical values
of the ordering parameter as a functrion of compositrion, temperature
and pressure. This task requires an iterative procedure. To construct
this procedure the model must have information on the permissble domain
of the ordering parameter.

Values of the ordering parameter are bounded by the composition of the
solution. We contrain all site mole fractions to have values in the
range 0 to 1, and solve this system of inequality constraints to obtain
a logical expression that embodies the feasible domain for the numerical
procedure.

.. code:: ipython3

    out = sym.reduce_inequalities(inequalities=[
        0 <= (r[1]-s[0])/2, (r[1]-s[0])/2 <= 1,
        0 <= (r[1]+s[0])/2, (r[1]+s[0])/2 <= 1, 
        0 <= 1-r[0]+(r[1]+s[0])/2,  1-r[0]+(r[1]+s[0])/2 <= 1, 
        0 <= r[0]-(r[1]+s[0])/2, r[0]-(r[1]+s[0])/2 <= 1], symbols=[s[0]])
    out




.. math::

    \displaystyle s_{1} \geq - \frac{2 n_{3}}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \geq \frac{2 \left(- 2 n_{1} - 3 n_{3}\right)}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \geq - \frac{2 \left(2 n_{1} + 3 n_{3}\right)}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \geq \frac{2 \left(- 2 n_{1} - 2 n_{2} - n_{3}\right)}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \leq \frac{2 n_{3}}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \leq - \frac{2 \left(- 2 n_{2} + n_{3}\right)}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \leq \frac{2 \left(2 n_{2} - n_{3}\right)}{2 n_{1} + 2 n_{2} + 2 n_{3}} \wedge s_{1} \leq - \frac{2 \left(- 2 n_{1} - 2 n_{2} - n_{3}\right)}{2 n_{1} + 2 n_{2} + 2 n_{3}}



Add the Gibbs free energy of solution to the model
--------------------------------------------------

.. code:: ipython3

    model.add_expression_to_model(G, list(zip(params, units, symparams)), ordering_functions=([dgds],s,[0],out))

… give the model a unqiue name

.. code:: ipython3

    model.module = "Complex_Solution"

| … assign a formula string for code generation
| … assign a conversion string to map element concentrations to moles of
  end members … assign a test string to evaluate the feasibility of
  input compositions

.. code:: ipython3

    model.formula_string = 'Ca[Ca]Mg[Mg]Fe[Fe]Si[Si]O6'
    model.conversion_string = ['[0]=[Ca]-[Fe]', '[1]=[Fe]', '[2]=([Mg]-[Ca])/2.0']
    model.test_string = ['[0]+[1] > 0.0', '[1] > 0.0', '[0]+2.0*[2] > 0.0']

Define Parameters of an Orthopyroxene Solution
==============================================

Components 1. diopside, CaMgSi2O6 2. hedenbergite, CaFeSi2O6 3.
enstatite, Mg2Si2O6

Original calibration from Sack and Ghiorso (Contributions to Mineralogy
and Petrology, 116:287-300, 1994):

::

   F       = -13807 + 2.319*T - 0.05878*P;  /* joules     */
   Gex     =  -7824           - 0.1213*P;   /* joules/K   */
   Gx      =  -1883           + 0.02824;    /* joules/bar */
   WM2CaMg =  31631           + 0.03347*P;  /* joules     */
   WM2CaFe =  17238           + 0.04602*P;  /* joules/K   */
   WM1FeMg =   8368           + 0.01412*P;  /* joules/bar */
   WM2FeMg =   8368           + 0.01412*P;  /* joules     */

Asymmetry along the Ca-Mg and Ca-Fe joins (considered by Sack and
Ghiorso, 1994) is not considered in this example in order to simplify
the presentation.

.. code:: ipython3

    print (params)
    paramValues = {'Fh':-13807.0, 'Fs':-2.319, 'Fv':-0.05878, \
                   'Hex':-7824.0, 'Vex':-0.1213, \
                   'Hx':-1883.0, 'Vx':0.02824, \
                   'WhM2CaMg':31631.0, 'WvM2CaMg':0.03347, \
                   'WhM2CaFe':17238.0, 'WvM2CaFe':0.04602, \
                   'WhM1FeMg':8368.0, 'WvM1FeMg':0.01412, \
                   'WhM2FeMg':8368.0, 'WvM2FeMg':0.01412, \
                   'T_r':298.15, 'P_r':1.0}
    print (paramValues)


.. parsed-literal::

    ['Fh', 'Fs', 'Fv', 'Hex', 'Vex', 'Hx', 'Vx', 'WhM2CaMg', 'WvM2CaMg', 'WhM2CaFe', 'WvM2CaFe', 'WhM1FeMg', 'WvM1FeMg', 'WhM2FeMg', 'WvM2FeMg']
    {'Fh': -13807.0, 'Fs': -2.319, 'Fv': -0.05878, 'Hex': -7824.0, 'Vex': -0.1213, 'Hx': -1883.0, 'Vx': 0.02824, 'WhM2CaMg': 31631.0, 'WvM2CaMg': 0.03347, 'WhM2CaFe': 17238.0, 'WvM2CaFe': 0.04602, 'WhM1FeMg': 8368.0, 'WvM1FeMg': 0.01412, 'WhM2FeMg': 8368.0, 'WvM2FeMg': 0.01412, 'T_r': 298.15, 'P_r': 1.0}


Generate both fast computation and calibibration code for the feldspar
solution

Use code printers to construct “C” package code
===============================================

.. code:: ipython3

    model_working_dir = "working"
    !mkdir -p {model_working_dir}
    %cd {model_working_dir}


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working


Choose model type and create model
----------------------------------

model_type is “fast” or “calib”

.. code:: ipython3

    model_type = "calib"

.. code:: ipython3

    model.create_code_module(phase="Orthopyroxene", params=paramValues, 
                             endmembers=['Diopside_berman', 'Hedenbergite_berman', 'Enstatite_berman'], 
                             prefix="cy", module_type=model_type, silent=False, 
                             add_code_to_access_order_paramater=True)


.. parsed-literal::

    Creating generic fast model code file string
    .Writing include file to working directory ...
    Creating (once only) generic model calib codetemplate include file string
    Creating (once only) generic model calib codetemplate code file string
    Creating generic calib model code file string
    Writing include file to working directory ...
    Creating code blocks for standard state properties.
    Creating calib code and include files ...
    Writing include file to working directory ...
    Writing code file to working directory ...
    ... done
    Writing pyxbld file to working directory ...
    writing pyx file to working directory ...
    Compiling code and Python bindings ...
    Success! Import the module named  Complex_Solution
    Generating additional methods for order parameter access.
    Done!


Load the module
---------------

.. code:: ipython3

    import Complex_Solution
    %cd ..


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/Cython/Compiler/Main.py:369: FutureWarning: Cython directive 'language_level' not set, using 2 for now (Py2). This will change in a later release! File: /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working/Complex_Solution.pyx
      tree = Parsing.p_module(s, pxd, full_module_name)


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen


Test and time the generated functions for Orthopyroxene (T in K, P in bars)
---------------------------------------------------------------------------

.. code:: ipython3

    t = 2000.00
    p = 1.0
    n = np.array([1.1, 1.2, 1.3])

Available in both “Fast” and “Calib” code versions
--------------------------------------------------

Execute the “fast” or “calibration” code metadata retrieval functions:

.. code:: ipython3

    try:
        print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_identifier())
        print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_name())
        print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_formula(t,p,n))
    except AttributeError:
        pass
    try:
        print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_identifier())
        print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_name())
        print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_formula(t,p,n))
    except AttributeError:
        pass


.. parsed-literal::

    Wed Sep 23 09:58:22 2020
    Orthopyroxene
    Ca0.639Mg1.028Fe0.333Si2.000O6


Test intrinsic element conversion routine …

.. code:: ipython3

    try:
        e = np.zeros(106)
        sum = np.sum(n)
        for index in range(0,nc):
            end = Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_elements(index)
            for i in range(0,106):
                e[i] += end[i]*n[index]/sum
        nConv = Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_elm_to_moles(e)
        for i in range(0,nc):
            print ('X[{0:d}] input {1:13.6e}, calc {2:13.6e}, diff {3:13.6e}'.format(
            i, n[i]/sum, nConv[i], nConv[i]-n[i]/sum))
        if not Complex_Solution.cy_Orthopyroxene_Complex_Solution_test_moles(nConv):
            print ('Output of intrinsic composition calculation fails tests for permissible values.')
    except AttributeError:
        pass
    try:
        e = np.zeros(106)
        sum = np.sum(n)
        for index in range(0,nc):
            end = Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_elements(index)
            for i in range(0,106):
                e[i] += end[i]*n[index]/sum
        nConv = Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_elm_to_moles(e)
        for i in range(0,nc):
            print ('X[{0:d}] input {1:13.6e}, calc {2:13.6e}, diff {3:13.6e}'.format(
            i, n[i]/sum, nConv[i], nConv[i]-n[i]/sum))
        if not Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_test_moles(nConv):
            print ('Output of intrinsic composition calculation fails tests for permissible values.')
    except AttributeError:
        pass


.. parsed-literal::

    X[0] input  3.055556e-01, calc  3.055556e-01, diff  5.551115e-17
    X[1] input  3.333333e-01, calc  3.333333e-01, diff  0.000000e+00
    X[2] input  3.611111e-01, calc  1.944444e-01, diff -1.666667e-01


Test various conversion routines …

.. code:: ipython3

    try:
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_moles_to_tot_moles(n))
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_moles_to_mole_frac(n))
        e = Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_moles_to_elm(n)
        print (e)
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_elm_to_moles(e))
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_elm_to_tot_moles(e))
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_conv_elm_to_tot_grams(e))
    except AttributeError:
        pass
    try:
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_moles_to_tot_moles(n))
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_moles_to_mole_frac(n))
        e = Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_moles_to_elm(n)
        print (e)
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_elm_to_moles(e))
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_elm_to_tot_moles(e))
        print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_conv_elm_to_tot_grams(e))
    except AttributeError:
        pass


.. parsed-literal::

    3.5999999999999996
    [0.30555556 0.33333333 0.36111111]
    [ 0.   0.   0.   0.   0.   0.   0.   0.  21.6  0.   0.   0.   3.7  0.
      7.2  0.   0.   0.   0.   0.   2.3  0.   0.   0.   0.   0.   1.2  0.
      0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
      0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
      0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
      0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
      0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
      0.   0.   0.   0.   0.   0.   0.   0. ]
    [1.1 1.2 0.7]
    3.0
    676.4650999999999


Execute a method that retrieves the ordering parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This method is normally hidden from the Phases module implementation.

.. code:: ipython3

    nn = np.array([-0.99, 1, 1])
    print(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_formula(t,p,nn))
    print("Kd = (XFeM1 XMgM2)/(XFeM2 XMgM1), s = XFeM2 - XMgM2")
    for tc in [600.0, 900.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0]:
        ss = Complex_Solution.cy_Orthopyroxene_Complex_Solution_order_params(tc+273.15,p,nn)
        vCaM2 = xCaM2.subs([(model.n[0],nn[0]),(model.n[1],nn[1]),(model.n[2],nn[2]),(s[0],ss[0])])
        vMgM2 = xMgM2.subs([(model.n[0],nn[0]),(model.n[1],nn[1]),(model.n[2],nn[2]),(s[0],ss[0])])
        vFeM2 = xFeM2.subs([(model.n[0],nn[0]),(model.n[1],nn[1]),(model.n[2],nn[2]),(s[0],ss[0])])
        vMgM1 = xMgM1.subs([(model.n[0],nn[0]),(model.n[1],nn[1]),(model.n[2],nn[2]),(s[0],ss[0])])
        vFeM1 = xFeM1.subs([(model.n[0],nn[0]),(model.n[1],nn[1]),(model.n[2],nn[2]),(s[0],ss[0])])
        K = vFeM1*vMgM2/(vMgM1*vFeM2)
        print ("T {0:8.2f}  s {1:8.4f}  RT ln(Kd) {2:10.3f} kJ".format(tc, ss[0], 8.3143*(tc+273.15)*(np.log(float(K))/1000.0))) 


.. parsed-literal::

    Ca0.010Mg1.000Fe0.990Si2.000O6
    Kd = (XFeM1 XMgM2)/(XFeM2 XMgM1), s = XFeM2 - XMgM2
    T   600.00  s   0.5516  RT ln(Kd)    -18.350 kJ
    T   900.00  s   0.3513  RT ln(Kd)    -14.616 kJ
    T  1000.00  s   0.3089  RT ln(Kd)    -13.825 kJ
    T  1200.00  s   0.2467  RT ln(Kd)    -12.665 kJ
    T  1400.00  s   0.2042  RT ln(Kd)    -11.871 kJ
    T  1600.00  s   0.1736  RT ln(Kd)    -11.299 kJ
    T  1800.00  s   0.1507  RT ln(Kd)    -10.870 kJ
    T  2000.00  s   0.1329  RT ln(Kd)    -10.536 kJ


Execute the standard thermodynamic property retrieval functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
    try:
        print(fmt.format('G', Complex_Solution.cy_Orthopyroxene_Complex_Solution_g(t,p,n), 'J'))
        print(fmt.format('dGdT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dgdt(t,p,n), 'J/K'))
        print(fmt.format('dGdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dgdp(t,p,n), 'J/bar'))
        print(fmt.format('d2GdT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d2gdt2(t,p,n), 'J/K^2'))
        print(fmt.format('d2GdTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d2gdtdp(t,p,n), 'J/K-bar'))
        print(fmt.format('d2GdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d2gdp2(t,p,n), 'J/bar^2'))
        print(fmt.format('d3GdT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdt3(t,p,n), 'J/K^3'))
        print(fmt.format('d3GdT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdt2dp(t,p,n), 'J/K^2-bar'))
        print(fmt.format('d3GdTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdtdp2(t,p,n), 'J/K-bar^2'))
        print(fmt.format('d3GdP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdp3(t,p,n), 'J/bar^3'))
        print(fmt.format('S', Complex_Solution.cy_Orthopyroxene_Complex_Solution_s(t,p,n), 'J/K'))
        print(fmt.format('V', Complex_Solution.cy_Orthopyroxene_Complex_Solution_v(t,p,n), 'J/bar'))
        print(fmt.format('Cv', Complex_Solution.cy_Orthopyroxene_Complex_Solution_cv(t,p,n), 'J/K'))
        print(fmt.format('Cp', Complex_Solution.cy_Orthopyroxene_Complex_Solution_cp(t,p,n), 'J/K'))
        print(fmt.format('dCpdT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dcpdt(t,p,n), 'J/K^2'))
        print(fmt.format('alpha', Complex_Solution.cy_Orthopyroxene_Complex_Solution_alpha(t,p,n), '1/K'))
        print(fmt.format('beta', Complex_Solution.cy_Orthopyroxene_Complex_Solution_beta(t,p,n), '1/bar'))
        print(fmt.format('K', Complex_Solution.cy_Orthopyroxene_Complex_Solution_K(t,p,n), 'bar'))
        print(fmt.format('Kp', Complex_Solution.cy_Orthopyroxene_Complex_Solution_Kp(t,p,n), ''))
    except AttributeError:
        pass
    try:
        print(fmt.format('G', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_g(t,p,n), 'J'))
        print(fmt.format('dGdT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_dgdt(t,p,n), 'J/K'))
        print(fmt.format('dGdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_dgdp(t,p,n), 'J/bar'))
        print(fmt.format('d2GdT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d2gdt2(t,p,n), 'J/K^2'))
        print(fmt.format('d2GdTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d2gdtdp(t,p,n), 'J/K-bar'))
        print(fmt.format('d2GdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d2gdp2(t,p,n), 'J/bar^2'))
        print(fmt.format('d3GdT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdt3(t,p,n), 'J/K^3'))
        print(fmt.format('d3GdT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdt2dp(t,p,n), 'J/K^2-bar'))
        print(fmt.format('d3GdTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdtdp2(t,p,n), 'J/K-bar^2'))
        print(fmt.format('d3GdP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdp3(t,p,n), 'J/bar^3'))
        print(fmt.format('S', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_s(t,p,n), 'J/K'))
        print(fmt.format('V', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_v(t,p,n), 'J/bar'))
        print(fmt.format('Cv', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_cv(t,p,n), 'J/K'))
        print(fmt.format('Cp', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_cp(t,p,n), 'J/K'))
        print(fmt.format('dCpdT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_dcpdt(t,p,n), 'J/K^2'))
        print(fmt.format('alpha', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_alpha(t,p,n), '1/K'))
        print(fmt.format('beta', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_beta(t,p,n), '1/bar'))
        print(fmt.format('K', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_K(t,p,n), 'bar'))
        print(fmt.format('Kp', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_Kp(t,p,n), ''))
    except AttributeError:
        pass


.. parsed-literal::

    G          -1.376489e+07 J         
    dGdT       -2.178302e+03 J/K       
    dGdP        2.533452e+01 J/bar     
    d2GdT2     -4.926731e-01 J/K^2     
    d2GdTdP     1.323450e-03 J/K-bar   
    d2GdP2     -2.071512e-05 J/bar^2   
    d3GdT3      2.271906e-04 J/K^3     
    d3GdT2dP    3.633866e-07 J/K^2-bar 
    d3GdTdP2    5.963507e-11 J/K-bar^2 
    d3GdP3      5.644345e-11 J/bar^3   
    S           2.178302e+03 J/K       
    V           2.533452e+01 J/bar     
    Cv          8.162407e+02 J/K       
    Cp          9.853461e+02 J/K       
    dCpdT       3.829190e-02 J/K^2     
    alpha       5.223899e-05 1/K       
    beta        8.176639e-07 1/bar     
    K           1.222996e+06 bar       
    Kp          2.332355e+00           


Execute functions that access endmember properties:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<15.15s}"
    try:
        print ("number of components", Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_number())
        for index in range(0, nc):
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_name(index)), end=' ')
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_formula(index)))
            print ("mw: {0:10.2f}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_mw(index)))
            print (fmt.format('mu0', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_mu0(index,t,p), 'J/mol'))
            print (fmt.format('dmu0dT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_dmu0dT(index,t,p), 'J/K-mol'))
            print (fmt.format('dmu0dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_dmu0dP(index,t,p), 'J/bar-mol'))
            print (fmt.format('d2mu0dT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d2mu0dT2(index,t,p), 'J/K^2-mol'))
            print (fmt.format('d2mu0dTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d2mu0dTdP(index,t,p), 'J/K-bar-mol'))
            print (fmt.format('d2mu0dP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d2mu0dP2(index,t,p), 'J/bar^2-mol'))
            print (fmt.format('d3mu0dT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d3mu0dT3(index,t,p), 'J/K^3-mol'))
            print (fmt.format('d3mu0dT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d3mu0dT2dP(index,t,p), 'J/K^2-bar-mol'))
            print (fmt.format('d3mu0dTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d3mu0dTdP2(index,t,p), 'J/K-bar^2-mol'))
            print (fmt.format('d3mu0dP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_d3mu0dP3(index,t,p), 'J/bar^3-mol'))
            print ("Element array:")
            print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_endmember_elements(index))
            print ()
    except AttributeError:
        pass
    try:
        print ("number of components", Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_number())
        for index in range(0, nc):
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_name(index)), end=' ')
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_formula(index)), end=' ')
            print ("mw: {0:10.2f}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_mw(index)))
            print (fmt.format('mu0', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_mu0(index,t,p), 'J/mol'))
            print (fmt.format('dmu0dT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_dmu0dT(index,t,p), 'J/K-mol'))
            print (fmt.format('dmu0dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_dmu0dP(index,t,p), 'J/bar-mol'))
            print (fmt.format('d2mu0dT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d2mu0dT2(index,t,p), 'J/K^2-mol'))
            print (fmt.format('d2mu0dTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d2mu0dTdP(index,t,p), 'J/K-bar-mol'))
            print (fmt.format('d2mu0dP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d2mu0dP2(index,t,p), 'J/bar^2-mol'))
            print (fmt.format('d3mu0dT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d3mu0dT3(index,t,p), 'J/K^3-mol'))
            print (fmt.format('d3mu0dT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d3mu0dT2dP(index,t,p), 'J/K^2-bar-mol'))
            print (fmt.format('d3mu0dTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d3mu0dTdP2(index,t,p), 'J/K-bar^2-mol'))
            print (fmt.format('d3mu0dP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_d3mu0dP3(index,t,p), 'J/bar^3-mol'))
            print ("Element array:")
            print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_endmember_elements(index))
            print ()
    except AttributeError:
        pass


.. parsed-literal::

    number of components 3
    Diopside             MgCaSi2O6            mw:     216.55
    mu0        -3.947955e+06 J/mol          
    dmu0dT     -5.818145e+02 J/K-mol        
    dmu0dP      7.092442e+00 J/bar-mol      
    d2mu0dT2   -1.339235e-01 J/K^2-mol      
    d2mu0dTdP   3.712074e-04 J/K-bar-mol    
    d2mu0dP2   -5.772640e-06 J/bar^2-mol    
    d3mu0dT3    6.166661e-05 J/K^3-mol      
    d3mu0dT2dP  1.100006e-07 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    2.260068e-11 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 6. 0. 0. 0. 1. 0. 2. 0. 0. 0. 0. 0. 1. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Hedenbergite         CaFeSi2O6            mw:     248.09
    mu0        -3.660551e+06 J/mol          
    dmu0dT     -6.198971e+02 J/K-mol        
    dmu0dP      7.316411e+00 J/bar-mol      
    d2mu0dT2   -1.352710e-01 J/K^2-mol      
    d2mu0dTdP   4.063486e-04 J/K-bar-mol    
    d2mu0dP2   -6.738479e-06 J/bar^2-mol    
    d3mu0dT3    6.238454e-05 J/K^3-mol      
    d3mu0dT2dP  1.136165e-07 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    2.014415e-11 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 6. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 1. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Enstatite            Mg2Si2O6             mw:     200.78
    mu0        -3.821935e+06 J/mol          
    dmu0dT     -5.805518e+02 J/K-mol        
    dmu0dP      6.730275e+00 J/bar-mol      
    d2mu0dT2   -1.392010e-01 J/K^2-mol      
    d2mu0dTdP   3.168470e-04 J/K-bar-mol    
    d2mu0dP2   -4.739597e-06 J/bar^2-mol    
    d3mu0dT3    6.237364e-05 J/K^3-mol      
    d3mu0dT2dP  9.450086e-08 J/K^2-bar-mol  
    d3mu0dTdP2  0.000000e+00 J/K-bar^2-mol  
    d3mu0dP3    5.657143e-12 J/bar^3-mol    
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 6. 0. 0. 0. 2. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    


Execute functions that access species properties:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<15.15s}"
    try:
        print ("number of species", Complex_Solution.cy_Orthopyroxene_Complex_Solution_species_number())
        for index in range(0, nc):
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_species_name(index)), end=' ')
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_species_formula(index)))
            print ("mw: {0:10.2f}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_species_mw(index)))
            print ("Element array:")
            print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_species_elements(index))
            print ()
    except AttributeError:
        pass
    try:
        print ("number of species", Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_species_number())
        for index in range(0, nc):
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_species_name(index)), end=' ')
            print ("{0:<20.20s}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_species_formula(index)), end=' ')
            print ("mw: {0:10.2f}".format(Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_species_mw(index)))
            print ("Element array:")
            print (Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_species_elements(index))
            print ()
    except AttributeError:
        pass


.. parsed-literal::

    number of species 3
    Diopside             MgCaSi2O6            mw:     216.55
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 6. 0. 0. 0. 1. 0. 2. 0. 0. 0. 0. 0. 1. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Hedenbergite         CaFeSi2O6            mw:     248.09
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 6. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 1. 0. 0. 0.
     0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    
    Enstatite            Mg2Si2O6             mw:     200.78
    Element array:
    [0. 0. 0. 0. 0. 0. 0. 0. 6. 0. 0. 0. 2. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    


Execute functions for molar derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First derivative vectors:
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    def printResult(name, result, units):
        print ("{0:<10.10s}".format(name), end=' ')
        [print ("{0:13.6e}".format(x), end=' ') for x in result]
        print ("{0:<10.10s}".format(units))
    def printLabels(n):
        print ("{0:<18.18s}".format(''), end=' ')
        [print ("[{0:3d}]{1:<8.8s}".format(idx, ''), end=' ') for idx in range(len(n))]
        print ()
    printLabels(n)
    try:
        printResult('dGdn', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dgdn(t,p,n), 'J/m')
        printResult('d2GdndT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d2gdndt(t,p,n), 'J/K-m')
        printResult('d2GdndP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d2gdndp(t,p,n), 'J/bar-m')
        printResult('d3GdndT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdndt2(t,p,n), 'J/K^2-m')
        printResult('d3GdndTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdndtdp(t,p,n), 'J/K-bar-m')
        printResult('d3GdndP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdndp2(t,p,n), 'J/bar^2-m')
        printResult('d4GdndT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdndt3(t,p,n), 'J/K^3-m')
        printResult('d4GdndT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdndt2dp(t,p,n), 'J/K^2-bar-m')
        printResult('d4GdndTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdndtdp2(t,p,n), 'J/K-bar^2-m')
        printResult('d4GdndP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdndp3(t,p,n), 'J/bar^3-m')
    except AttributeError:
        pass
    try:
        printResult('dGdn', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_dgdn(t,p,n), 'J/m')
        printResult('d2GdndT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d2gdndt(t,p,n), 'J/K-m')
        printResult('d2GdndP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d2gdndp(t,p,n), 'J/bar-m')
        printResult('d3GdndT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdndt2(t,p,n), 'J/K^2-m')
        printResult('d3GdndTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdndtdp(t,p,n), 'J/K-bar-m')
        printResult('d3GdndP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdndp2(t,p,n), 'J/bar^2-m')
        printResult('d4GdndT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdndt3(t,p,n), 'J/K^3-m')
        printResult('d4GdndT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdndt2dp(t,p,n), 'J/K^2-bar-m')
        printResult('d4GdndTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdndtdp2(t,p,n), 'J/K-bar^2-m')
        printResult('d4GdndP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdndp3(t,p,n), 'J/bar^3-m')
    except AttributeError:
        pass    


.. parsed-literal::

                       [  0]         [  1]         [  2]         
    dGdn       -3.955206e+06 -3.688954e+06 -3.836474e+06 J/m       
    d2GdndT    -5.876726e+02 -6.339219e+02 -5.931968e+02 J/K-m     
    d2GdndP     7.099392e+00  7.309661e+00  6.733534e+00 J/bar-m   
    d3GdndT2   -1.343084e-01 -1.363411e-01 -1.394804e-01 J/K^2-m   
    d3GdndTdP            nan           nan           nan J/K-bar-m 
    d3GdndP2   -5.761793e-06 -6.779023e-06 -4.801784e-06 J/bar^2-m 
    d4GdndT3    0.000000e+00  0.000000e+00  0.000000e+00 J/K^3-m   
    d4GdndT2dP  0.000000e+00  0.000000e+00  0.000000e+00 J/K^2-bar-
    d4GdndTdP2  0.000000e+00  0.000000e+00  0.000000e+00 J/K-bar^2-
    d4GdndP3    0.000000e+00  0.000000e+00  0.000000e+00 J/bar^3-m 


The Hessian matrix (molar second derivative matrix) is stored as a compact linear array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A function is provided to map matrix indices to compact storage 1-D
array indices

.. code:: ipython3

    for i in range(1,nc+1):
        print ("[ ", end=' ')
        for j in range (1,nc+1):
            print ((i,j), end=' ')
        print (']     [', end=' ')
        for j in range (1,nc+1):
            print (model.symmetric_index_from_2d_array(elm=(i,j)), end=' ')
        print (']')


.. parsed-literal::

    [  (1, 1) (1, 2) (1, 3) ]     [ 0 1 2 ]
    [  (2, 1) (2, 2) (2, 3) ]     [ 1 3 4 ]
    [  (3, 1) (3, 2) (3, 3) ]     [ 2 4 5 ]


.. code:: ipython3

    def printResult(name, result, units):
        print ("{0:<10.10s}".format(name), end=' ')
        [print ("{0:13.6e}".format(x), end=' ') for x in result]
        print ("{0:<10.10s}".format(units))
    def printLabels(n):
        print ("{0:<18.18s}".format(''), end=' ')
        maxIdx = int(len(n)*(len(n)-1)/2 + len(n))
        [print ("[{0:3d}]{1:<8.8s}".format(idx, ''), end=' ') for idx in range(maxIdx)]
        print ()
    printLabels(n)
    try:
        printResult('d2Gdn2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d2gdn2(t,p,n), 'J/m^2')
        printResult('d3Gdn2dT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdn2dt(t,p,n), 'J/K-m^2')
        printResult('d3Gdn2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdn2dp(t,p,n), 'J/bar-m^2')
        printResult('d4Gdn2dT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdn2dt2(t,p,n), 'J/K^2-m^2')
        printResult('d4Gdn2dTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdn2dtdp(t,p,n), 'J/K-bar-m^2')
        printResult('d4Gdn2dP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdn2dp2(t,p,n), 'J/bar^2-m^2')
        printResult('d5Gdn2dT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn2dt3(t,p,n), 'J/K^3-m^2')
        printResult('d5Gdn2dT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn2dt2dp(t,p,n), 'J/K^2-bar-m^2')
        printResult('d5Gdn2dTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn2dtdp2(t,p,n), 'J/K-bar^2-m^2')
        printResult('d5Gdn2dP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn2dp3(t,p,n), 'J/bar^3-m^2')
    except AttributeError:
        pass
    try:
        printResult('d2Gdn2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d2gdn2(t,p,n), 'J/m^2')
        printResult('d3Gdn2dT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdn2dt(t,p,n), 'J/K-m^2')
        printResult('d3Gdn2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdn2dp(t,p,n), 'J/bar-m^2')
        printResult('d4Gdn2dT2', Complex_Solution.cy_Orthopyroxene_ComplexSolution_calib_d4gdn2dt2(t,p,n), 'J/K^2-m^2')
        printResult('d4Gdn2dTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdn2dtdp(t,p,n), 'J/K-bar-m^2')
        printResult('d4Gdn2dP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdn2dp2(t,p,n), 'J/bar^2-m^2')
        printResult('d5Gdn2dT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn2dt3(t,p,n), 'J/K^3-m^2')
        printResult('d5Gdn2dT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn2dt2dp(t,p,n), 'J/K^2-bar-m^2')
        printResult('d5Gdn2dTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn2dtdp2(t,p,n), 'J/K-bar^2-m^2')
        printResult('d5Gdn2dP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn2dp3(t,p,n), 'J/bar^3-m^2')
    except AttributeError:
        pass


.. parsed-literal::

                       [  0]         [  1]         [  2]         [  3]         [  4]         [  5]         
    d2Gdn2      1.597539e+03 -1.949313e+03  4.476024e+02  1.234486e+04 -9.745836e+03  8.617416e+03 J/m^2     
    d3Gdn2dT    1.735760e+00 -4.137925e-01 -1.086758e+00  5.012103e+00 -4.276425e+00  4.867033e+00 J/K-m^2   
    d3Gdn2dP   -2.723683e-03 -3.119739e-04  2.592631e-03  1.279253e-02 -1.154451e-02  8.462708e-03 J/bar-m^2 


The 3-D Tensor (molar third derivative tensor) is stored as a compact linear array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| A function is provided to map matrix indices to compact storage 1-D
  array indices:
| If :math:`n_c` represents the number of components in the solution,
  and
| if :math:`n_d` represents the dimensionality of molar derivative (in
  this case 3), then
| the number of numerically ordered permutations of :math:`n_c` molar
  derivatives taken :math:`n_d` at a time is:

.. code:: ipython3

    n_c,n_d = sym.symbols('n_c n_d')
    q = sym.factorial(n_c+n_d-1)/sym.factorial(n_d)/sym.factorial(n_c-1)
    q




.. math::

    \displaystyle \frac{\left(n_{c} + n_{d} - 1\right)!}{n_{d}! \left(n_{c} - 1\right)!}



Substituting :math:`n_d` equal to 3 and simplifying gives:

.. code:: ipython3

    q = sym.simplify(q.subs(n_d,3))
    q




.. math::

    \displaystyle \frac{n_{c} \left(n_{c} + 1\right) \left(n_{c} + 2\right)}{6}



and, for the number of components in this solution, there will be the
following number of unique terms in the third derivative tensor:

.. code:: ipython3

    q.subs(n_c,nc)




.. math::

    \displaystyle 10



A function is provided to map matrix indices to compact storage 1-D
array indices

.. code:: ipython3

    for i in range(1,nc+1):
        for j in range (1,nc+1):
            print ("[", end=' ')
            for k in range (1,nc+1):
                print ("{0:1d}{1:1d}{2:1d}".format(i,j,k), end=' ')
            print ('] ', end=' ')
        print ('  ->  ', end=' ')
        for j in range (1,nc+1):
            print ("[", end=' ')
            for k in range (1,nc+1):
                print (model.symmetric_index_from_3d_array(elm=(i,j,k)), end=' ')
            print ('] ', end=' ')
        print ('')


.. parsed-literal::

    [ 111 112 113 ]  [ 121 122 123 ]  [ 131 132 133 ]    ->   [ 0 1 2 ]  [ 1 3 4 ]  [ 2 4 5 ]  
    [ 211 212 213 ]  [ 221 222 223 ]  [ 231 232 233 ]    ->   [ 1 3 4 ]  [ 3 6 7 ]  [ 4 7 8 ]  
    [ 311 312 313 ]  [ 321 322 323 ]  [ 331 332 333 ]    ->   [ 2 4 5 ]  [ 4 7 8 ]  [ 5 8 9 ]  


.. code:: ipython3

    def printResult(name, result, units):
        print ("{0:<10.10s}".format(name), end=' ')
        [print ("{0:10.3e}".format(x), end=' ') for x in result]
        print ("{0:<14.14s}".format(units))
    def printLabels(n):
        print ("{0:<15.15s}".format(''), end=' ')
        maxIdx = int(len(n)*(len(n)+1)*(len(n)+2)/6)
        [print ("[{0:3d}]{1:<5.5s}".format(idx, ''), end=' ') for idx in range(maxIdx)]
        print ()
    printLabels(n)
    try:
        printResult('d3Gdn3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d3gdn3(t,p,n), 'J/m^3')
        printResult('d4Gdn3dT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdn3dt(t,p,n), 'J/K-m^3')
        printResult('d4Gdn3dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d4gdn3dp(t,p,n), 'J/bar-m^3')
        printResult('d5Gdn3dT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn3dt2(t,p,n), 'J/K^2-m^3')
        printResult('d5Gdn3dTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn3dtdp(t,p,n), 'J/K-bar-m^3')
        printResult('d5Gdn3dP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d5gdn3dp2(t,p,n), 'J/bar^2-m^3')
        printResult('d6Gdn3dT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d6gdn3dt3(t,p,n), 'J/K^3-m^3')
        printResult('d6Gdn3dT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d6gdn3dt2dp(t,p,n), 'J/K^2-bar-m^3')
        printResult('d6Gdn3dTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d6gdn3dtdp2(t,p,n), 'J/K-bar^2-m^3')
        printResult('d6Gdn3dP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_d6gdn3dp3(t,p,n), 'J/bar^3-m^3')
    except AttributeError:
        pass
    try:
        printResult('d3Gdn3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d3gdn3(t,p,n), 'J/m^3')
        printResult('d4Gdn3dT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdn3dt(t,p,n), 'J/K-m^3')
        printResult('d4Gdn3dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d4gdn3dp(t,p,n), 'J/bar-m^3')
        printResult('d5Gdn3dT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn3dt2(t,p,n), 'J/K^2-m^3')
        printResult('d5Gdn3dTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn3dtdp(t,p,n), 'J/K-bar-m^3')
        printResult('d5Gdn3dP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d5gdn3dp2(t,p,n), 'J/bar^2-m^3')
        printResult('d6Gdn3dT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d6gdn3dt3(t,p,n), 'J/K^3-m^3')
        printResult('d6Gdn3dT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d6gdn3dt2dp(t,p,n), 'J/K^2-bar-m^3')
        printResult('d6Gdn3dTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d6gdn3dtdp2(t,p,n), 'J/K-bar^2-m^3')
        printResult('d6Gdn3dP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_d6gdn3dp3(t,p,n), 'J/bar^3-m^3')
    except AttributeError:
        pass


.. parsed-literal::

                    [  0]      [  1]      [  2]      [  3]      [  4]      [  5]      [  6]      [  7]      [  8]      [  9]      
    d3Gdn3     -8.506e+02  2.257e+02 -7.175e+02  5.764e+02  7.764e+02 -4.539e+02 -1.111e+04  2.712e+02  6.589e+03 -1.233e+04 J/m^3         
    d4Gdn3dT    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-m^3       
    d4Gdn3dP    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar-m^3     
    d5Gdn3dT2   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^2-m^3     
    d5Gdn3dTdP  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-bar-m^3   
    d5Gdn3dP2   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar^2-m^3   
    d6Gdn3dT3   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^3-m^3     
    d6Gdn3dT2d  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K^2-bar-m^3 
    d6Gdn3dTdP  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/K-bar^2-m^3 
    d6Gdn3dP3   0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00 J/bar^3-m^3   


Test and time the generated functions for Feldspar
--------------------------------------------------

Time the code

.. code:: ipython3

    try:
        %timeit Complex_Solution.cy_Orthopyroxene_Complex_Solution_g(t, p, n)
    except AttributeError:
        pass
    try:
        %timeit Complex_Solution.cy_Orthopyroxene_Complex_Solution_calib_g(t, p, n) 
    except AttributeError:
        pass


.. parsed-literal::

    8.22 µs ± 246 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)


Time the Rubicon wrapped Objective-C code

.. code:: ipython3

    from thermoengine import model as stdmodel
    modelDB = stdmodel.Database()
    CpxHC = modelDB.get_phase('Cpx')

.. code:: ipython3

    %timeit CpxHC.gibbs_energy(t,p,mol=np.array([1.1, 1.2, 1.3, 0.0, 0.0, 0.0, 0.0])) 


.. parsed-literal::

    69.2 µs ± 1.96 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)


Methods available only in the “Calib” versions of generated code
----------------------------------------------------------------

Execute the parameter value/metadata functions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions are only defined for the “calibration” model code
implementation:

.. code:: ipython3

    nparam = 0

.. code:: ipython3

    try:
        nparam = Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_number()
        names = Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_names()
        units = Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_units()
        values = Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_values()
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,nparam):
            print(fmt.format(names[i], values[i], Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_value(i), units[i]))
    except AttributeError:
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+00  1.000000e+00 bar       
    Fh         -1.380700e+04 -1.380700e+04 J/mol     
    Fs         -2.319000e+00 -2.319000e+00 J/K-mol   
    Fv         -5.878000e-02 -5.878000e-02 J/bar-mol 
    Hex        -7.824000e+03 -7.824000e+03 J/mol     
    Vex        -1.213000e-01 -1.213000e-01 J/bar-mol 
    Hx         -1.883000e+03 -1.883000e+03 J/mol     
    Vx          2.824000e-02  2.824000e-02 J/bar-mol 
    WhM2CaMg    3.163100e+04  3.163100e+04 J/mol     
    WvM2CaMg    3.347000e-02  3.347000e-02 J/bar-mol 
    WhM2CaFe    1.723800e+04  1.723800e+04 J/mol     
    WvM2CaFe    4.602000e-02  4.602000e-02 J/bar-mol 
    WhM1FeMg    8.368000e+03  8.368000e+03 J/mol     
    WvM1FeMg    1.412000e-02  1.412000e-02 J/bar-mol 
    WhM2FeMg    8.368000e+03  8.368000e+03 J/mol     
    WvM2FeMg    1.412000e-02  1.412000e-02 J/bar-mol 


Functions that allow modification of the array of parameter values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    try:
        values[1] = 100.0
        Complex_Solution.cy_Orthopyroxene_Complex_Solution_set_param_values(values)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,nparam):
            print(fmt.format(names[i], values[i], Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_value(i), units[i]))
    except (AttributeError, NameError):
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+02  1.000000e+02 bar       
    Fh         -1.380700e+04 -1.380700e+04 J/mol     
    Fs         -2.319000e+00 -2.319000e+00 J/K-mol   
    Fv         -5.878000e-02 -5.878000e-02 J/bar-mol 
    Hex        -7.824000e+03 -7.824000e+03 J/mol     
    Vex        -1.213000e-01 -1.213000e-01 J/bar-mol 
    Hx         -1.883000e+03 -1.883000e+03 J/mol     
    Vx          2.824000e-02  2.824000e-02 J/bar-mol 
    WhM2CaMg    3.163100e+04  3.163100e+04 J/mol     
    WvM2CaMg    3.347000e-02  3.347000e-02 J/bar-mol 
    WhM2CaFe    1.723800e+04  1.723800e+04 J/mol     
    WvM2CaFe    4.602000e-02  4.602000e-02 J/bar-mol 
    WhM1FeMg    8.368000e+03  8.368000e+03 J/mol     
    WvM1FeMg    1.412000e-02  1.412000e-02 J/bar-mol 
    WhM2FeMg    8.368000e+03  8.368000e+03 J/mol     
    WvM2FeMg    1.412000e-02  1.412000e-02 J/bar-mol 


Functions that allow modification of a particular parameter value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    try:
        Complex_Solution.cy_Orthopyroxene_Complex_Solution_set_param_value(1, 1.0)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,nparam):
            print(fmt.format(names[i], values[i], Complex_Solution.cy_Orthopyroxene_Complex_Solution_get_param_value(i), units[i]))
    except AttributeError:
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+02  1.000000e+00 bar       
    Fh         -1.380700e+04 -1.380700e+04 J/mol     
    Fs         -2.319000e+00 -2.319000e+00 J/K-mol   
    Fv         -5.878000e-02 -5.878000e-02 J/bar-mol 
    Hex        -7.824000e+03 -7.824000e+03 J/mol     
    Vex        -1.213000e-01 -1.213000e-01 J/bar-mol 
    Hx         -1.883000e+03 -1.883000e+03 J/mol     
    Vx          2.824000e-02  2.824000e-02 J/bar-mol 
    WhM2CaMg    3.163100e+04  3.163100e+04 J/mol     
    WvM2CaMg    3.347000e-02  3.347000e-02 J/bar-mol 
    WhM2CaFe    1.723800e+04  1.723800e+04 J/mol     
    WvM2CaFe    4.602000e-02  4.602000e-02 J/bar-mol 
    WhM1FeMg    8.368000e+03  8.368000e+03 J/mol     
    WvM1FeMg    1.412000e-02  1.412000e-02 J/bar-mol 
    WhM2FeMg    8.368000e+03  8.368000e+03 J/mol     
    WvM2FeMg    1.412000e-02  1.412000e-02 J/bar-mol 


Functions that evaluate parameter derivatives …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    try:
        fmt = "    {0:<10.10s} {1:13.6e}"
        for i in range(0, nparam):
            print ('Derivative with respect to parameter: ', names[i], ' of')
            print (fmt.format('G', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_g(t, p, n, i)))
            print (fmt.format('dGdT', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_dgdt(t, p, n, i)))
            print (fmt.format('dGdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_dgdp(t, p, n, i)))
            print (fmt.format('d2GdT2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d2gdt2(t, p, n, i)))
            print (fmt.format('d2GdTdP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d2gdtdp(t, p, n, i)))
            print (fmt.format('d2GdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d2gdp2(t, p, n, i)))
            print (fmt.format('d3GdT3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d3gdt3(t, p, n, i)))
            print (fmt.format('d3GdT2dP', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d3gdt2dp(t, p, n, i)))
            print (fmt.format('d3GdTdP2', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d3gdtdp2(t, p, n, i)))
            print (fmt.format('d3GdP3', Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_d3gdp3(t, p, n, i)))
    except (AttributeError, TypeError):
        pass


.. parsed-literal::

    Derivative with respect to parameter:  T_r  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3               nan
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3               nan
    Derivative with respect to parameter:  P_r  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3               nan
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3               nan
    Derivative with respect to parameter:  Fh  of
        G           4.048797e-01
        dGdT       -5.341382e-05
        dGdP        4.024796e-07
        d2GdT2      5.753133e-08
        d2GdTdP    -2.041586e-10
        d2GdP2     -1.897976e-13
        d3GdT3               nan
        d3GdT2dP    2.583639e-13
        d3GdTdP2    3.043538e-15
        d3GdP3               nan
    Derivative with respect to parameter:  Fs  of
        G          -8.097595e+02
        dGdT       -2.980521e-01
        dGdP       -8.049592e-04
        d2GdT2     -8.235011e-06
        d2GdTdP     5.837508e-09
        d2GdP2      3.795952e-10
        d3GdT3               nan
        d3GdT2dP    2.495297e-10
        d3GdTdP2   -5.517684e-12
        d3GdP3               nan
    Derivative with respect to parameter:  Fv  of
        G           4.048797e-01
        dGdT       -5.341382e-05
        dGdP        4.048801e-01
        d2GdT2      5.753133e-08
        d2GdTdP    -5.341403e-05
        d2GdP2      8.049590e-07
        d3GdT3               nan
        d3GdT2dP    1.117201e-07
        d3GdTdP2   -7.662544e-10
        d3GdP3               nan
    Derivative with respect to parameter:  Hex  of
        G           1.843090e-01
        dGdT       -1.137947e-04
        dGdP        8.574566e-07
        d2GdT2      1.225667e-07
        d2GdTdP    -4.349465e-10
        d2GdP2     -4.043514e-13
        d3GdT3               nan
        d3GdT2dP    5.504274e-13
        d3GdTdP2    6.484060e-15
        d3GdP3               nan
    Derivative with respect to parameter:  Vex  of
        G           1.843090e-01
        dGdT       -1.137947e-04
        dGdP        1.843099e-01
        d2GdT2      1.225667e-07
        d2GdTdP    -1.137951e-04
        d2GdP2      1.714913e-06
        d3GdT3               nan
        d3GdT2dP    2.380124e-07
        d3GdTdP2   -1.632455e-09
        d3GdP3               nan
    Derivative with respect to parameter:  Hx  of
        G           3.049122e-01
        dGdT       -5.240953e-05
        dGdP        3.949121e-07
        d2GdT2      7.198221e-08
        d2GdTdP    -3.173600e-10
        d2GdP2      6.956816e-13
        d3GdT3               nan
        d3GdT2dP    9.659313e-13
        d3GdTdP2    6.333638e-16
        d3GdP3               nan
    Derivative with respect to parameter:  Vx  of
        G           3.049122e-01
        dGdT       -5.240953e-05
        dGdP        3.049126e-01
        d2GdT2      7.198221e-08
        d2GdTdP    -5.240985e-05
        d2GdP2      7.898250e-07
        d3GdT3               nan
        d3GdT2dP    1.562180e-07
        d3GdTdP2   -1.454090e-09
        d3GdP3               nan
    Derivative with respect to parameter:  WhM2CaMg  of
        G           4.541294e-01
        dGdT        1.068276e-04
        dGdP       -8.049592e-07
        d2GdT2     -1.150627e-07
        d2GdTdP     4.083171e-10
        d2GdP2      3.795952e-13
        d3GdT3               nan
        d3GdT2dP   -5.167278e-13
        d3GdTdP2   -6.087076e-15
        d3GdP3               nan
    Derivative with respect to parameter:  WvM2CaMg  of
        G           4.541294e-01
        dGdT        1.068276e-04
        dGdP        4.541286e-01
        d2GdT2     -1.150627e-07
        d2GdTdP     1.068281e-04
        d2GdP2     -1.609918e-06
        d3GdT3               nan
        d3GdT2dP   -2.234402e-07
        d3GdTdP2    1.532509e-09
        d3GdP3               nan
    Derivative with respect to parameter:  WhM2CaFe  of
        G           3.764262e-01
        dGdT       -1.068276e-04
        dGdP        8.049592e-07
        d2GdT2      1.150627e-07
        d2GdTdP    -4.083171e-10
        d2GdP2     -3.795952e-13
        d3GdT3               nan
        d3GdT2dP    5.167278e-13
        d3GdTdP2    6.087076e-15
        d3GdP3               nan
    Derivative with respect to parameter:  WvM2CaFe  of
        G           3.764262e-01
        dGdT       -1.068276e-04
        dGdP        3.764270e-01
        d2GdT2      1.150627e-07
        d2GdTdP    -1.068281e-04
        d2GdP2      1.609918e-06
        d3GdT3               nan
        d3GdT2dP    2.234402e-07
        d3GdTdP2   -1.532509e-09
        d3GdP3               nan
    Derivative with respect to parameter:  WhM1FeMg  of
        G           5.071750e-01
        dGdT        1.104680e-04
        dGdP       -8.323900e-07
        d2GdT2     -1.345163e-07
        d2GdTdP     5.392714e-10
        d2GdP2     -4.893799e-13
        d3GdT3               nan
        d3GdT2dP   -1.246762e-12
        d3GdTdP2   -3.941558e-15
        d3GdP3               nan
    Derivative with respect to parameter:  WvM1FeMg  of
        G           5.071750e-01
        dGdT        1.104680e-04
        dGdP        5.071742e-01
        d2GdT2     -1.345163e-07
        d2GdTdP     1.104686e-04
        d2GdP2     -1.664780e-06
        d3GdT3               nan
        d3GdT2dP   -2.776530e-07
        d3GdTdP2    2.286975e-09
        d3GdP3               nan
    Derivative with respect to parameter:  WhM2FeMg  of
        G           1.163339e-01
        dGdT       -5.648974e-06
        dGdP        4.256570e-08
        d2GdT2     -9.448163e-09
        d2GdTdP     9.544848e-11
        d2GdP2     -9.019834e-13
        d3GdT3               nan
        d3GdT2dP   -6.851010e-13
        d3GdTdP2    2.674830e-15
        d3GdP3               nan
    Derivative with respect to parameter:  WvM2FeMg  of
        G           1.163339e-01
        dGdT       -5.648974e-06
        dGdP        1.163339e-01
        d2GdT2     -9.448163e-09
        d2GdTdP    -5.648879e-06
        d2GdP2      8.513050e-08
        d3GdT3               nan
        d3GdT2dP   -3.478314e-08
        d3GdTdP2    6.212044e-10
        d3GdP3               nan


Parameter derivatives of the chemical potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def printResult(name, result, units):
        print ("dmu[*]/d {0:<10.10s}".format(name), end=' ')
        [print ("{0:13.6e}".format(x), end=' ') for x in result]
        print ("{0:<12.12s}".format(units))
    def printLabels(n):
        print ("         {0:<18.18s}".format(''), end=' ')
        [print ("[{0:3d}]{1:<8.8s}".format(idx, ''), end=' ') for idx in range(len(n))]
        print ()
    try:
        printLabels(n)
        for i in range(0, nparam):
            result = Complex_Solution.cy_Orthopyroxene_Complex_Solution_dparam_dgdn(t,p,n, i)
            printResult(names[i], result, 'J/m^2/p-unit')
    except AttributeError:
        pass    


.. parsed-literal::

                                [  0]         [  1]         [  2]         
    dmu[*]/d T_r         0.000000e+00  0.000000e+00  0.000000e+00 J/m^2/p-unit
    dmu[*]/d P_r         0.000000e+00  0.000000e+00  0.000000e+00 J/m^2/p-unit
    dmu[*]/d Fh         -4.394940e-02  2.557182e-01  1.125863e-01 J/m^2/p-unit
    dmu[*]/d Fs          8.789879e+01 -5.114364e+02 -2.251726e+02 J/m^2/p-unit
    dmu[*]/d Fv         -4.394940e-02  2.557182e-01  1.125863e-01 J/m^2/p-unit
    dmu[*]/d Hex         2.268990e-03  7.547392e-02  7.018803e-02 J/m^2/p-unit
    dmu[*]/d Vex         2.268990e-03  7.547392e-02  7.018803e-02 J/m^2/p-unit
    dmu[*]/d Hx         -1.593030e-02  1.178342e-01  1.392573e-01 J/m^2/p-unit
    dmu[*]/d Vx         -1.593030e-02  1.178342e-01  1.392573e-01 J/m^2/p-unit
    dmu[*]/d WhM2CaMg    9.792966e-02 -1.402945e-01  3.959693e-01 J/m^2/p-unit
    dmu[*]/d WvM2CaMg    9.792966e-02 -1.402945e-01  3.959693e-01 J/m^2/p-unit
    dmu[*]/d WhM2CaFe    3.247158e-02  2.706957e-01  1.220968e-02 J/m^2/p-unit
    dmu[*]/d WvM2CaFe    3.247158e-02  2.706957e-01  1.220968e-02 J/m^2/p-unit
    dmu[*]/d WhM1FeMg    5.632443e-02  4.706427e-01 -9.196394e-02 J/m^2/p-unit
    dmu[*]/d WvM1FeMg    5.632443e-02  4.706427e-01 -9.196394e-02 J/m^2/p-unit
    dmu[*]/d WhM2FeMg   -3.372310e-02 -2.112597e-02  1.375234e-01 J/m^2/p-unit
    dmu[*]/d WvM2FeMg   -3.372310e-02 -2.112597e-02  1.375234e-01 J/m^2/p-unit


