
Berman SymPy Testing Theory and Code Generation
===============================================

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()

Standard State Properties - Structure of the Equations
------------------------------------------------------

There are three classes of terms: - Terms that apply over the whole of
:math:`T`-, :math:`P`-space, :math:`T_r \le T`, :math:`P_r \le P` -
Terms that apply over a specified range of :math:`T`-, :math:`P`-space,
:math:`(T_{r_\lambda},P_{r_\lambda}) \le (T,P) \le (T_\lambda,P_\lambda)`
- Terms that apply to a specific :math:`T_t` and :math:`P_t` and higher
:math:`T`, :math:`P`, :math:`T_t \le T`, :math:`P_t \le P`

Second-order phase transitions (:math:`lambda`-transitions) are an
example of the second type, as are order disorder transformations.
First-order phase transitions are an example of the third type.

Berman terms valid over the whole of T,P space
----------------------------------------------

| Parameters of the Berman (1988) model that pertain to the whole of
  :math:`T`-:math:`P` space: - :math:`\Delta{H_{T_r,P_r}}`, the enthalpy
  of formation from the elements at :math:`T_r` and :math:`P_r` -
  :math:`S_{T_r,P_r}`, the third law entropy at :math:`T_r` and
  :math:`P_r` - :math:`V_{T_r,P_r}`, the volume at :math:`T_r` and
  :math:`P_r` - :math:`k_0`, :math:`k_1`, :math:`k_2`, and :math:`k_3`,
  coefficients in the reference pressure heat capacity expression
| - :math:`v_1`, :math:`v_2`, :math:`v_3`, and :math:`v_4`, coefficients
  in Berman Equation of State

.. code:: ipython3

    Tr,Pr = sym.symbols('T_r P_r')
    k0,k1,k2,k3 = sym.symbols('k0 k1 k2 k3')
    VTrPr,v1,v2,v3,v4 = sym.symbols('V_TrPr v1 v2 v3 v4')
    STrPr,HTrPr = sym.symbols('S_TrPr H_TrPr')

Input variables:

.. code:: ipython3

    T,P = sym.symbols('T P')

Heat capacity at constant pressure. Berman and Brown (1985) form.

.. code:: ipython3

    CpPr = k0+k1/sym.sqrt(T)+k2/T**2+k3/T**3
    CpPr




.. math::

    k_{0} + \frac{k_{2}}{T^{2}} + \frac{k_{3}}{T^{3}} + \frac{k_{1}}{\sqrt{T}}



Apparent enthalpy of formation at :math:`T` and :math:`P_r`. The
quantity is apparent because the heat capacity of the substance is
integrated up in :math:`T` from :math:`T_r`, but the enthalpy of the
elements are left at :math:`T_r` and :math:`P_r`.

.. code:: ipython3

    HPr = HTrPr + sym.integrate(CpPr,(T,Tr,T))
    HPr




.. math::

    H_{TrPr} + 2 \sqrt{T} k_{1} + T k_{0} - 2 \sqrt{T_{r}} k_{1} - T_{r} k_{0} + \frac{k_{2}}{T_{r}} + \frac{k_{3}}{2 T_{r}^{2}} - \frac{k_{2}}{T} - \frac{k_{3}}{2 T^{2}}



Entropy at :math:`T` and :math:`P_r`.

.. code:: ipython3

    SPr = STrPr + sym.integrate(CpPr/T,(T,Tr,T))
    SPr




.. math::

    S_{TrPr} + k_{0} \log{\left (T \right )} - k_{0} \log{\left (T_{r} \right )} + \frac{k_{2}}{2 T_{r}^{2}} + \frac{k_{3}}{3 T_{r}^{3}} + \frac{2 k_{1}}{\sqrt{T_{r}}} - \frac{k_{2}}{2 T^{2}} - \frac{k_{3}}{3 T^{3}} - \frac{2 k_{1}}{\sqrt{T}}



Volume at :math:`T` and :math:`P` evaluated using the Berman (1988)
Equation of State.

.. code:: ipython3

    V = VTrPr*(1+v1*(P-Pr)+v2*(P-Pr)**2+v3*(T-Tr)+v4*(T-Tr)**2)
    V




.. math::

    V_{TrPr} \left(v_{1} \left(P - P_{r}\right) + v_{2} \left(P - P_{r}\right)^{2} + v_{3} \left(T - T_{r}\right) + v_{4} \left(T - T_{r}\right)^{2} + 1\right)



Apparent Gibbs Free Energy of formation at :math:`T` and :math:`P_r`.

.. code:: ipython3

    GPr = HPr - T*SPr
    sym.simplify(GPr)




.. math::

    H_{TrPr} - S_{TrPr} T + 4 \sqrt{T} k_{1} - T k_{0} \log{\left (T \right )} + T k_{0} \log{\left (T_{r} \right )} + T k_{0} - \frac{T k_{2}}{2 T_{r}^{2}} - \frac{T k_{3}}{3 T_{r}^{3}} - \frac{2 T k_{1}}{\sqrt{T_{r}}} - 2 \sqrt{T_{r}} k_{1} - T_{r} k_{0} + \frac{k_{2}}{T_{r}} + \frac{k_{3}}{2 T_{r}^{2}} - \frac{k_{2}}{2 T} - \frac{k_{3}}{6 T^{2}}



Apparent Gibbs Free Energy of formation at :math:`T` and :math:`P`.

.. code:: ipython3

    G = GPr + sym.integrate(V,(P,Pr,P))
    sym.simplify(G)




.. math::

    H_{TrPr} + \frac{P^{3} V_{TrPr} v_{2}}{3} - P^{2} P_{r} V_{TrPr} v_{2} + \frac{P^{2} V_{TrPr} v_{1}}{2} + P P_{r}^{2} V_{TrPr} v_{2} - P P_{r} V_{TrPr} v_{1} + P T^{2} V_{TrPr} v_{4} - 2 P T T_{r} V_{TrPr} v_{4} + P T V_{TrPr} v_{3} + P T_{r}^{2} V_{TrPr} v_{4} - P T_{r} V_{TrPr} v_{3} + P V_{TrPr} - \frac{P_{r}^{3} V_{TrPr} v_{2}}{3} + \frac{P_{r}^{2} V_{TrPr} v_{1}}{2} - P_{r} T^{2} V_{TrPr} v_{4} + 2 P_{r} T T_{r} V_{TrPr} v_{4} - P_{r} T V_{TrPr} v_{3} - P_{r} T_{r}^{2} V_{TrPr} v_{4} + P_{r} T_{r} V_{TrPr} v_{3} - P_{r} V_{TrPr} - S_{TrPr} T + 4 \sqrt{T} k_{1} - T k_{0} \log{\left (T \right )} + T k_{0} \log{\left (T_{r} \right )} + T k_{0} - \frac{T k_{2}}{2 T_{r}^{2}} - \frac{T k_{3}}{3 T_{r}^{3}} - \frac{2 T k_{1}}{\sqrt{T_{r}}} - 2 \sqrt{T_{r}} k_{1} - T_{r} k_{0} + \frac{k_{2}}{T_{r}} + \frac{k_{3}}{2 T_{r}^{2}} - \frac{k_{2}}{2 T} - \frac{k_{3}}{6 T^{2}}



Berman terms valid over an interval of T,P space
------------------------------------------------

Parameters of the Berman (1988) lambda transition model: - :math:`l_1`
and :math:`l_2`, coefficients in Berman (1988)’s :math:`lambda`-heat
capacity model - :math:`k_{\lambda}`,
:math:`\frac{{d{T_\lambda }}}{{dP}}` in Berman (1988)’s
:math:`lambda`-transition model - :math:`T_{\lambda,{P_r}}`, Temperature
of the :math:`\lambda`-transition at reference pressure -
:math:`T_{\lambda,{ref}}`, Temperature of the lower bound of the heat
capacity integral for the :math:`\lambda`-transition at reference
pressure

.. code:: ipython3

    l1,l2 = sym.symbols('l1 l2')
    kl = sym.symbols('k_lambda')
    TlPr, Tlref = sym.symbols('T_lambda_Pr T_lambda_ref')

:math:`T_{\lambda}` at :math:`P`.

.. code:: ipython3

    Tl = TlPr + kl*(P-Pr)
    Tl




.. math::

    T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)



Temperature difference between :math:`T_{lambda}` at :math:`P` and
:math:`P_r`

.. code:: ipython3

    td = TlPr - Tl
    td




.. math::

    - k_{\lambda} \left(P - P_{r}\right)



Reference temperature for lower limit of heat capacity integral.

.. code:: ipython3

    tr = Tlref - td
    tr




.. math::

    T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)



Heat capacity due to the :math:`\lambda`-transition at :math:`T` and
:math:`P`. Valid: :math:`T_r \le T \le T_{\lambda}`.

.. code:: ipython3

    Cpl = (T+td)*(l1+l2*(T+td))**2
    Cpl




.. math::

    \left(T - k_{\lambda} \left(P - P_{r}\right)\right) \left(l_{1} + l_{2} \left(T - k_{\lambda} \left(P - P_{r}\right)\right)\right)^{2}



Enthalpy due to the :math:`\lambda`-transition at :math:`T` and
:math:`P`. Valid: :math:`T_r \le T \le T_{\lambda}`.

.. code:: ipython3

    Hl = sym.integrate(Cpl,(T,tr,T))
    Hl




.. math::

    \frac{T^{4} l_{2}^{2}}{4} + T^{3} \left(- P k_{\lambda} l_{2}^{2} + P_{r} k_{\lambda} l_{2}^{2} + \frac{2 l_{1} l_{2}}{3}\right) + T^{2} \left(\frac{3 P^{2} k_{\lambda}^{2} l_{2}^{2}}{2} - 3 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 2 P k_{\lambda} l_{1} l_{2} + \frac{3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2}}{2} + 2 P_{r} k_{\lambda} l_{1} l_{2} + \frac{l_{1}^{2}}{2}\right) + T \left(- P^{3} k_{\lambda}^{3} l_{2}^{2} + 3 P^{2} P_{r} k_{\lambda}^{3} l_{2}^{2} + 2 P^{2} k_{\lambda}^{2} l_{1} l_{2} - 3 P P_{r}^{2} k_{\lambda}^{3} l_{2}^{2} - 4 P P_{r} k_{\lambda}^{2} l_{1} l_{2} - P k_{\lambda} l_{1}^{2} + P_{r}^{3} k_{\lambda}^{3} l_{2}^{2} + 2 P_{r}^{2} k_{\lambda}^{2} l_{1} l_{2} + P_{r} k_{\lambda} l_{1}^{2}\right) - \frac{l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{4}}{4} - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3} \left(- P k_{\lambda} l_{2}^{2} + P_{r} k_{\lambda} l_{2}^{2} + \frac{2 l_{1} l_{2}}{3}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(\frac{3 P^{2} k_{\lambda}^{2} l_{2}^{2}}{2} - 3 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 2 P k_{\lambda} l_{1} l_{2} + \frac{3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2}}{2} + 2 P_{r} k_{\lambda} l_{1} l_{2} + \frac{l_{1}^{2}}{2}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(- P^{3} k_{\lambda}^{3} l_{2}^{2} + 3 P^{2} P_{r} k_{\lambda}^{3} l_{2}^{2} + 2 P^{2} k_{\lambda}^{2} l_{1} l_{2} - 3 P P_{r}^{2} k_{\lambda}^{3} l_{2}^{2} - 4 P P_{r} k_{\lambda}^{2} l_{1} l_{2} - P k_{\lambda} l_{1}^{2} + P_{r}^{3} k_{\lambda}^{3} l_{2}^{2} + 2 P_{r}^{2} k_{\lambda}^{2} l_{1} l_{2} + P_{r} k_{\lambda} l_{1}^{2}\right)



Entropy due to the :math:`\lambda`-transition at :math:`T` and
:math:`P`. Valid: :math:`T_r \le T \le T_{\lambda}`.

.. code:: ipython3

    Sl = sym.integrate(Cpl/T,(T,tr,T))
    Sl




.. math::

    \frac{T^{3} l_{2}^{2}}{3} + T^{2} \left(- \frac{3 P k_{\lambda} l_{2}^{2}}{2} + \frac{3 P_{r} k_{\lambda} l_{2}^{2}}{2} + l_{1} l_{2}\right) + T \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right) - k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T \right )} + k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )} - \frac{l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3}}{3} - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(- \frac{3 P k_{\lambda} l_{2}^{2}}{2} + \frac{3 P_{r} k_{\lambda} l_{2}^{2}}{2} + l_{1} l_{2}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right)



Gibbs Free Energy due to the :math:`\lambda`-transition at :math:`T` and
:math:`P`. Valid: :math:`T_r \le T \le T_{\lambda}`.

.. code:: ipython3

    Gl = Hl - T*Sl
    Gl




.. math::

    \frac{T^{4} l_{2}^{2}}{4} + T^{3} \left(- P k_{\lambda} l_{2}^{2} + P_{r} k_{\lambda} l_{2}^{2} + \frac{2 l_{1} l_{2}}{3}\right) + T^{2} \left(\frac{3 P^{2} k_{\lambda}^{2} l_{2}^{2}}{2} - 3 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 2 P k_{\lambda} l_{1} l_{2} + \frac{3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2}}{2} + 2 P_{r} k_{\lambda} l_{1} l_{2} + \frac{l_{1}^{2}}{2}\right) - T \left(\frac{T^{3} l_{2}^{2}}{3} + T^{2} \left(- \frac{3 P k_{\lambda} l_{2}^{2}}{2} + \frac{3 P_{r} k_{\lambda} l_{2}^{2}}{2} + l_{1} l_{2}\right) + T \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right) - k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T \right )} + k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )} - \frac{l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3}}{3} - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(- \frac{3 P k_{\lambda} l_{2}^{2}}{2} + \frac{3 P_{r} k_{\lambda} l_{2}^{2}}{2} + l_{1} l_{2}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right)\right) + T \left(- P^{3} k_{\lambda}^{3} l_{2}^{2} + 3 P^{2} P_{r} k_{\lambda}^{3} l_{2}^{2} + 2 P^{2} k_{\lambda}^{2} l_{1} l_{2} - 3 P P_{r}^{2} k_{\lambda}^{3} l_{2}^{2} - 4 P P_{r} k_{\lambda}^{2} l_{1} l_{2} - P k_{\lambda} l_{1}^{2} + P_{r}^{3} k_{\lambda}^{3} l_{2}^{2} + 2 P_{r}^{2} k_{\lambda}^{2} l_{1} l_{2} + P_{r} k_{\lambda} l_{1}^{2}\right) - \frac{l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{4}}{4} - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3} \left(- P k_{\lambda} l_{2}^{2} + P_{r} k_{\lambda} l_{2}^{2} + \frac{2 l_{1} l_{2}}{3}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(\frac{3 P^{2} k_{\lambda}^{2} l_{2}^{2}}{2} - 3 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 2 P k_{\lambda} l_{1} l_{2} + \frac{3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2}}{2} + 2 P_{r} k_{\lambda} l_{1} l_{2} + \frac{l_{1}^{2}}{2}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(- P^{3} k_{\lambda}^{3} l_{2}^{2} + 3 P^{2} P_{r} k_{\lambda}^{3} l_{2}^{2} + 2 P^{2} k_{\lambda}^{2} l_{1} l_{2} - 3 P P_{r}^{2} k_{\lambda}^{3} l_{2}^{2} - 4 P P_{r} k_{\lambda}^{2} l_{1} l_{2} - P k_{\lambda} l_{1}^{2} + P_{r}^{3} k_{\lambda}^{3} l_{2}^{2} + 2 P_{r}^{2} k_{\lambda}^{2} l_{1} l_{2} + P_{r} k_{\lambda} l_{1}^{2}\right)



Gibbs Free Energy due to the :math:`\lambda`-transition at :math:`T` and
:math:`P` when :math:`k_{\lambda}` is zero. Valid:
:math:`T_r \le T \le T_{\lambda}`.

.. code:: ipython3

    sym.simplify(Gl.subs(kl, 0))




.. math::

    - \frac{T^{4} l_{2}^{2}}{12} - \frac{T^{3} l_{1} l_{2}}{3} - \frac{T^{2} l_{1}^{2}}{2} + \frac{T T_{\lambda ref}^{3} l_{2}^{2}}{3} + T T_{\lambda ref}^{2} l_{1} l_{2} + T T_{\lambda ref} l_{1}^{2} - \frac{T_{\lambda ref}^{4} l_{2}^{2}}{4} - \frac{2 T_{\lambda ref}^{3} l_{1} l_{2}}{3} - \frac{T_{\lambda ref}^{2} l_{1}^{2}}{2}



The volume of the :math:`\lambda`-transition is given by differentiation
of the Gibbs Free Energy:

.. code:: ipython3

    Vl = sym.simplify(Gl.diff(P))
    Vl




.. math::

    \frac{k_{\lambda} \left(T \left(2 k_{\lambda} \left(- P + P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} + \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(6 P^{2} k_{\lambda}^{2} l_{2}^{2} - 12 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 8 P k_{\lambda} l_{1} l_{2} + 6 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 8 P_{r} k_{\lambda} l_{1} l_{2} + 3 T^{2} l_{2}^{2} + 4 T l_{2} \left(- 3 P k_{\lambda} l_{2} + 3 P_{r} k_{\lambda} l_{2} + 2 l_{1}\right) - 4 k_{\lambda} l_{2} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right) \log{\left (T \right )} + 4 k_{\lambda} l_{2} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right) \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )} + 2 l_{1}^{2} - l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} - 2 l_{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(- 3 P k_{\lambda} l_{2} + 3 P_{r} k_{\lambda} l_{2} + 2 l_{1}\right) + 2 \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T \right )} - 2 \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )}\right)\right) + 2 \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(- T^{3} l_{2}^{2} + T^{2} l_{2} \left(3 P k_{\lambda} l_{2} - 3 P_{r} k_{\lambda} l_{2} - 2 l_{1}\right) - T \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right) - k_{\lambda} \left(- P^{3} k_{\lambda}^{2} l_{2}^{2} + 3 P^{2} P_{r} k_{\lambda}^{2} l_{2}^{2} + 2 P^{2} k_{\lambda} l_{1} l_{2} - 3 P P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} - 4 P P_{r} k_{\lambda} l_{1} l_{2} - P l_{1}^{2} + P_{r}^{3} k_{\lambda}^{2} l_{2}^{2} + 2 P_{r}^{2} k_{\lambda} l_{1} l_{2} + P_{r} l_{1}^{2}\right)\right)\right)}{2 \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)}



For :math:`T > T_{\lambda}`, we have the following contributions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    CpatTl = Cpl.subs(T,Tl)
    CpatTl




.. math::

    T_{\lambda Pr} \left(T_{\lambda Pr} l_{2} + l_{1}\right)^{2}



.. code:: ipython3

    HatTl = sym.simplify(Hl.subs(T,Tl))
    HatTl




.. math::

    \frac{T_{\lambda Pr}^{4} l_{2}^{2}}{4} + \frac{2 T_{\lambda Pr}^{3} l_{1} l_{2}}{3} + \frac{T_{\lambda Pr}^{2} l_{1}^{2}}{2} - \frac{T_{\lambda ref}^{4} l_{2}^{2}}{4} - \frac{2 T_{\lambda ref}^{3} l_{1} l_{2}}{3} - \frac{T_{\lambda ref}^{2} l_{1}^{2}}{2}



.. code:: ipython3

    SatTl = sym.simplify(Sl.subs(T,Tl))
    SatTl




.. math::

    - k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right) \right )} + k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )} + \frac{l_{2}^{2} \left(T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right)^{3}}{3} - \frac{l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3}}{3} + \frac{l_{2} \left(T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(- 3 P k_{\lambda} l_{2} + 3 P_{r} k_{\lambda} l_{2} + 2 l_{1}\right)}{2} + \frac{l_{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(3 P k_{\lambda} l_{2} - 3 P_{r} k_{\lambda} l_{2} - 2 l_{1}\right)}{2} + \left(T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right)



.. code:: ipython3

    _.subs(kl,0)




.. math::

    \frac{T_{\lambda Pr}^{3} l_{2}^{2}}{3} + T_{\lambda Pr}^{2} l_{1} l_{2} + T_{\lambda Pr} l_{1}^{2} - \frac{T_{\lambda ref}^{3} l_{2}^{2}}{3} - T_{\lambda ref}^{2} l_{1} l_{2} - T_{\lambda ref} l_{1}^{2}



.. code:: ipython3

    GatTl = sym.simplify(HatTl-T*SatTl)
    GatTl




.. math::

    \frac{T \left(6 k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right) \right )} - 6 k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )} - 2 l_{2}^{2} \left(T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right)^{3} + 2 l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3} - 3 l_{2} \left(T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(- 3 P k_{\lambda} l_{2} + 3 P_{r} k_{\lambda} l_{2} + 2 l_{1}\right) + 3 l_{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(- 3 P k_{\lambda} l_{2} + 3 P_{r} k_{\lambda} l_{2} + 2 l_{1}\right) - 6 \left(T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right) + 6 \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right)\right)}{6} + \frac{T_{\lambda Pr}^{4} l_{2}^{2}}{4} + \frac{2 T_{\lambda Pr}^{3} l_{1} l_{2}}{3} + \frac{T_{\lambda Pr}^{2} l_{1}^{2}}{2} - \frac{T_{\lambda ref}^{4} l_{2}^{2}}{4} - \frac{2 T_{\lambda ref}^{3} l_{1} l_{2}}{3} - \frac{T_{\lambda ref}^{2} l_{1}^{2}}{2}



.. code:: ipython3

    sym.simplify(_.subs(kl,0))




.. math::

    - \frac{T \left(T_{\lambda Pr}^{3} l_{2}^{2} + 3 T_{\lambda Pr}^{2} l_{1} l_{2} + 3 T_{\lambda Pr} l_{1}^{2} - T_{\lambda ref}^{3} l_{2}^{2} - 3 T_{\lambda ref}^{2} l_{1} l_{2} - 3 T_{\lambda ref} l_{1}^{2}\right)}{3} + \frac{T_{\lambda Pr}^{4} l_{2}^{2}}{4} + \frac{2 T_{\lambda Pr}^{3} l_{1} l_{2}}{3} + \frac{T_{\lambda Pr}^{2} l_{1}^{2}}{2} - \frac{T_{\lambda ref}^{4} l_{2}^{2}}{4} - \frac{2 T_{\lambda ref}^{3} l_{1} l_{2}}{3} - \frac{T_{\lambda ref}^{2} l_{1}^{2}}{2}



Berman terms valid at T :math:`\ge` :math:`T_t`
-----------------------------------------------

Berman parameters (in this case :math:`T_t` is equivalent to
:math:`T_{\lambda}`: - :math:`{{\Delta}_t}H`, First order enthalpy
contribution at :math:`T_{\lambda}`

.. code:: ipython3

    deltaHt = sym.symbols('H_t')

:math:`{{\Delta}_t}S`, First order enropy contribution at
:math:`T_{\lambda}`

.. code:: ipython3

    deltaSt = deltaHt/Tl
    deltaSt




.. math::

    \frac{H_{t}}{T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)}



Above :math:`T_{\lambda}` the following term is convenient to define
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    GaboveTl = -(T-Tl)*(deltaSt+Sl)
    GaboveTl




.. math::

    \left(- T + T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)\right) \left(\frac{H_{t}}{T_{\lambda Pr} + k_{\lambda} \left(P - P_{r}\right)} + \frac{T^{3} l_{2}^{2}}{3} + T^{2} \left(- \frac{3 P k_{\lambda} l_{2}^{2}}{2} + \frac{3 P_{r} k_{\lambda} l_{2}^{2}}{2} + l_{1} l_{2}\right) + T \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right) - k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T \right )} + k_{\lambda} \left(P - P_{r}\right) \left(- P k_{\lambda} l_{2} + P_{r} k_{\lambda} l_{2} + l_{1}\right)^{2} \log{\left (T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right) \right )} - \frac{l_{2}^{2} \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{3}}{3} - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right)^{2} \left(- \frac{3 P k_{\lambda} l_{2}^{2}}{2} + \frac{3 P_{r} k_{\lambda} l_{2}^{2}}{2} + l_{1} l_{2}\right) - \left(T_{\lambda ref} + k_{\lambda} \left(P - P_{r}\right)\right) \left(3 P^{2} k_{\lambda}^{2} l_{2}^{2} - 6 P P_{r} k_{\lambda}^{2} l_{2}^{2} - 4 P k_{\lambda} l_{1} l_{2} + 3 P_{r}^{2} k_{\lambda}^{2} l_{2}^{2} + 4 P_{r} k_{\lambda} l_{1} l_{2} + l_{1}^{2}\right)\right)



Berman terms valid over an interval of T,P space
------------------------------------------------

Parameters of the Berman (1988) order-disorder model: - :math:`d_0`,
:math:`d_1`, :math:`d_2`, :math:`d_3`, :math:`d_4`, :math:`d_5`,
order-disorder coefficients from the Berman (1988) model -
:math:`T_{D_{ref}}`, :math:`T_D`, minimum, maximum temperature of
ordering interval, :math:`T_{D_{ref}} \le T \le T_D`

.. code:: ipython3

    d0,d1,d2,d3,d4,d5 = sym.symbols('d0 d1 d2 d3 d4 d5')
    TD,TDref = sym.symbols('T_D T_D_ref')

Heat capacity of disorder at :math:`T` and :math:`P_r`

.. code:: ipython3

    CpDs = d0 + d1/sym.sqrt(T) + d2/T**2 + d3*T + d4*T**2
    CpDs




.. math::

    T^{2} d_{4} + T d_{3} + d_{0} + \frac{d_{2}}{T^{2}} + \frac{d_{1}}{\sqrt{T}}



Enthalpy of disorder at :math:`T` and :math:`P_r`

.. code:: ipython3

    HDs = sym.integrate(CpDs,(T,TDref,T))
    sym.collect(HDs,(d0,d1,d2,d3,d4))




.. math::

    d_{0} \left(T - T_{D ref}\right) + d_{1} \left(2 \sqrt{T} - 2 \sqrt{T_{D ref}}\right) + d_{2} \left(\frac{1}{T_{D ref}} - \frac{1}{T}\right) + d_{3} \left(\frac{T^{2}}{2} - \frac{T_{D ref}^{2}}{2}\right) + d_{4} \left(\frac{T^{3}}{3} - \frac{T_{D ref}^{3}}{3}\right)



Entropy of disorder at :math:`T` and :math:`P_r`

.. code:: ipython3

    SDs = sym.integrate(CpDs/T,(T,TDref,T))
    sym.collect(SDs,(d0,d1,d2,d3,d4))




.. math::

    d_{0} \left(2 \log{\left (\sqrt{T} \right )} - 2 \log{\left (\sqrt{T_{D ref}} \right )}\right) + d_{1} \left(\frac{2}{\sqrt{T_{D ref}}} - \frac{2}{\sqrt{T}}\right) + d_{2} \left(\frac{1}{2 T_{D ref}^{2}} - \frac{1}{2 T^{2}}\right) + d_{3} \left(T - T_{D ref}\right) + d_{4} \left(\frac{T^{2}}{2} - \frac{T_{D ref}^{2}}{2}\right)



Volume of disorder at :math:`T` and :math:`P`

.. code:: ipython3

    VDs = HDs/d5
    VDs




.. math::

    \frac{2 \sqrt{T} d_{1} + \frac{T^{3} d_{4}}{3} + \frac{T^{2} d_{3}}{2} + T d_{0} - 2 \sqrt{T_{D ref}} d_{1} - \frac{T_{D ref}^{3} d_{4}}{3} - \frac{T_{D ref}^{2} d_{3}}{2} - T_{D ref} d_{0} + \frac{d_{2}}{T_{D ref}} - \frac{d_{2}}{T}}{d_{5}}



Gibbs Free Energy of disorder at :math:`T` and :math:`P`

.. code:: ipython3

    GDs = HDs - T*SDs + VDs*(P-Pr)
    sym.collect(GDs,(d0,d1,d2,d3,d4))




.. math::

    - T \left(\frac{T^{2} d_{4}}{2} + T d_{3} - \frac{T_{D ref}^{2} d_{4}}{2} - T_{D ref} d_{3} + 2 d_{0} \log{\left (\sqrt{T} \right )} - 2 d_{0} \log{\left (\sqrt{T_{D ref}} \right )} + \frac{d_{2}}{2 T_{D ref}^{2}} + \frac{2 d_{1}}{\sqrt{T_{D ref}}} - \frac{d_{2}}{2 T^{2}} - \frac{2 d_{1}}{\sqrt{T}}\right) + d_{0} \left(T - T_{D ref}\right) + d_{1} \left(2 \sqrt{T} - 2 \sqrt{T_{D ref}}\right) + d_{2} \left(\frac{1}{T_{D ref}} - \frac{1}{T}\right) + d_{3} \left(\frac{T^{2}}{2} - \frac{T_{D ref}^{2}}{2}\right) + d_{4} \left(\frac{T^{3}}{3} - \frac{T_{D ref}^{3}}{3}\right) + \frac{\left(P - P_{r}\right) \left(2 \sqrt{T} d_{1} + \frac{T^{3} d_{4}}{3} + \frac{T^{2} d_{3}}{2} + T d_{0} - 2 \sqrt{T_{D ref}} d_{1} - \frac{T_{D ref}^{3} d_{4}}{3} - \frac{T_{D ref}^{2} d_{3}}{2} - T_{D ref} d_{0} + \frac{d_{2}}{T_{D ref}} - \frac{d_{2}}{T}\right)}{d_{5}}



Disordering Gibbs energy at :math:`T_D`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    GDsAtTD = GDs.subs(T,TD)
    GDsAtTD




.. math::

    2 \sqrt{T_{D}} d_{1} + \frac{T_{D}^{3} d_{4}}{3} + \frac{T_{D}^{2} d_{3}}{2} + T_{D} d_{0} - T_{D} \left(\frac{T_{D}^{2} d_{4}}{2} + T_{D} d_{3} - \frac{T_{D ref}^{2} d_{4}}{2} - T_{D ref} d_{3} + 2 d_{0} \log{\left (\sqrt{T_{D}} \right )} - 2 d_{0} \log{\left (\sqrt{T_{D ref}} \right )} + \frac{d_{2}}{2 T_{D ref}^{2}} + \frac{2 d_{1}}{\sqrt{T_{D ref}}} - \frac{d_{2}}{2 T_{D}^{2}} - \frac{2 d_{1}}{\sqrt{T_{D}}}\right) - 2 \sqrt{T_{D ref}} d_{1} - \frac{T_{D ref}^{3} d_{4}}{3} - \frac{T_{D ref}^{2} d_{3}}{2} - T_{D ref} d_{0} + \frac{\left(P - P_{r}\right) \left(2 \sqrt{T_{D}} d_{1} + \frac{T_{D}^{3} d_{4}}{3} + \frac{T_{D}^{2} d_{3}}{2} + T_{D} d_{0} - 2 \sqrt{T_{D ref}} d_{1} - \frac{T_{D ref}^{3} d_{4}}{3} - \frac{T_{D ref}^{2} d_{3}}{2} - T_{D ref} d_{0} + \frac{d_{2}}{T_{D ref}} - \frac{d_{2}}{T_{D}}\right)}{d_{5}} + \frac{d_{2}}{T_{D ref}} - \frac{d_{2}}{T_{D}}



Above :math:`T_D` the following term is convenient to define
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    GDsAboveTD = -(T-TD)*(SDs.subs(T,TD))
    GDsAboveTD




.. math::

    \left(- T + T_{D}\right) \left(\frac{T_{D}^{2} d_{4}}{2} + T_{D} d_{3} - \frac{T_{D ref}^{2} d_{4}}{2} - T_{D ref} d_{3} + 2 d_{0} \log{\left (\sqrt{T_{D}} \right )} - 2 d_{0} \log{\left (\sqrt{T_{D ref}} \right )} + \frac{d_{2}}{2 T_{D ref}^{2}} + \frac{2 d_{1}}{\sqrt{T_{D ref}}} - \frac{d_{2}}{2 T_{D}^{2}} - \frac{2 d_{1}}{\sqrt{T_{D}}}\right)



Use code printers to construct “C” package code
===============================================

Make the berman_calc.h and the berman_calib.h include files

.. code:: ipython3

    module = 'berman'
    params = ['H_TrPr', 'S_TrPr', 'k0', 'k1', 'k2', 'k3', 'V_TrPr', 'v1', 'v2', 'v3', 'v4', \
              'l1', 'l2', 'k_lambda', 'T_lambda_Pr', 'T_lambda_ref', 'H_t', \
              'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'T_D', 'T_D_ref', 'T_r', 'P_r']
    units = ['J/m', 'J/K-m', 'J/K-m', 'J/K^(1/2)-m', 'J-K/m', 'J-K^2/m', 'J/bar', '1/bar', '1/bar^2', '1/K', '1/K^2', \
            'J^(1/2)/K-m^(1/2)', 'J^(1/2)/K^2-m^(1/2)', 'K/bar', 'K', 'K', 'J/m', \
            'J/K-m', 'J/K^(1/2)-m', 'J-K/m', 'J/K^2-m', 'J/K^3-m', 'bar', 'K', 'K', 'K', 'bar']

Subclass the C99 code printer to expand :math:`x^2`, :math:`x^3`, and
:math:`x^4` as multiplication rather than pow(x, (double)n)

.. code:: ipython3

    from sympy.printing.ccode import C99CodePrinter
    class SubCodePrinter(C99CodePrinter):
        def _print_Pow(self, expr):
            if expr.exp.is_integer and expr.exp > 0 and expr.exp <= 4:
                result = ')*('.join([self._print(expr.base) for i in range(expr.exp)])
                return '((' + result + '))'
            else:
                return super()._print_Pow(expr)
    printer = SubCodePrinter()

.. code:: ipython3

    from sympy.printing.codeprinter import Assignment

.. code:: ipython3

    printer.doprint(Tr+2)




.. parsed-literal::

    'T_r + 2'



(1)
~~~

Function template for primary routines in berman_calc.h

.. code:: ipython3

    !mkdir -p working
    %cd working


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/working


.. code:: ipython3

    berman_calc_template = """\
    
    static double {module}_{func}(double T, double P) {{
        double result = {g_code};
        if ((l1 != 0.0) || (l2 != 0.0)) {{
          if ((({tlmin}) < T) && (({tlmax}) >= T)) {{
            result += {gl_code};  
          }} else if (T > ({tlmax})) {{
            result += {glmax_code};
          }}
        }}
        if ((H_t != 0.0) && (T >= ({tlmax}))) {{
          result += {glabove_code};
        }}
        if ((d0 != 0.0) || (d1 != 0.0) || (d2 != 0.0) || (d3 != 0.0) || (d4 != 0.0) || (d5 != 0.0)) {{
          if ((T_D_ref < T) && (T_D > T)) {{
            result += {gd_code};
          }} else if (T >= T_D) {{
            result += {gdmax_code};
            result += {gdabove_code};
          }}
        }}
        return result;
    }}
    \
    """

Primary routines in berman_calc.h

.. code:: ipython3

    f_list = ['g', 'dgdt', 'dgdp', 'd2gdt2', 'd2gdtdp', 'd2gdp2', 'd3gdt3', 'd3gdt2dp', 'd3gdtdp2', 'd3gdp3']
    G_matrix = sym.Matrix([G, G.diff(T), G.diff(P), G.diff(T,T), G.diff(T,P), G.diff(P,P), \
                           G.diff(T,T,T), G.diff(T,T,P), G.diff(T,P,P), G.diff(P,P,P)])
    Gl_matrix = sym.Matrix([Gl, Gl.diff(T), Gl.diff(P), Gl.diff(T,T), Gl.diff(T,P), Gl.diff(P,P), \
                            Gl.diff(T,T,T), Gl.diff(T,T,P), Gl.diff(T,P,P), Gl.diff(P,P,P)])
    GatTl_matrix = sym.Matrix([GatTl, GatTl.diff(T), GatTl.diff(P), GatTl.diff(T,T), GatTl.diff(T,P), \
                               GatTl.diff(P,P), GatTl.diff(T,T,T), GatTl.diff(T,T,P), GatTl.diff(T,P,P), \
                               GatTl.diff(P,P,P)])
    GaboveTl_matrix = sym.Matrix([GaboveTl, GaboveTl.diff(T), GaboveTl.diff(P), GaboveTl.diff(T,T), \
                                  GaboveTl.diff(T,P), GaboveTl.diff(P,P), GaboveTl.diff(T,T,T), \
                                  GaboveTl.diff(T,T,P), GaboveTl.diff(T,P,P), GaboveTl.diff(P,P,P)])
    GDs_matrix = sym.Matrix([GDs, GDs.diff(T), GDs.diff(P), GDs.diff(T,T), GDs.diff(T,P), GDs.diff(P,P), \
                             GDs.diff(T,T,T), GDs.diff(T,T,P), GDs.diff(T,P,P), GDs.diff(P,P,P)])
    GDsAtTD_matrix = sym.Matrix([GDsAtTD, GDsAtTD.diff(T), GDsAtTD.diff(P), GDsAtTD.diff(T,T), GDsAtTD.diff(T,P), \
                                 GDsAtTD.diff(P,P), GDsAtTD.diff(T,T,T), GDsAtTD.diff(T,T,P), GDsAtTD.diff(T,P,P), \
                                 GDsAtTD.diff(P,P,P)])
    GDsAboveTD_matrix = sym.Matrix([GDsAboveTD, GDsAboveTD.diff(T), GDsAboveTD.diff(P), GDsAboveTD.diff(T,T), \
                                    GDsAboveTD.diff(T,P), GDsAboveTD.diff(P,P), GDsAboveTD.diff(T,T,T), \
                                    GDsAboveTD.diff(T,T,P), GDsAboveTD.diff(T,P,P), GDsAboveTD.diff(P,P,P)])

.. code:: ipython3

    berman_calc = '#include <math.h>\n\n'
    
    for i in range(0,len(f_list)):
        berman_calc += berman_calc_template.format(\
                                  module=module,
                                  func=f_list[i],
                                  g_code=printer.doprint(G_matrix[i]), \
                                  tlmin=printer.doprint(tr), \
                                  tlmax=printer.doprint(Tl), \
                                  gl_code=printer.doprint(Gl_matrix[i]), \
                                  glmax_code=printer.doprint(GatTl_matrix[i]), \
                                  glabove_code=printer.doprint(GaboveTl_matrix[i]), \
                                  gd_code=printer.doprint(GDs_matrix[i]), \
                                  gdmax_code=printer.doprint(GDsAtTD_matrix[i]), \
                                  gdabove_code=printer.doprint(GDsAboveTD_matrix[i])
                                 )

Add convenience routiines to berman_calc.h

.. code:: ipython3

    convenience_template = """\
    
    static double {module}_s(double T, double P) {{
        double result = -{module}_dgdt(T, P);
        return result;
    }}
    
    static double {module}_v(double T, double P) {{
        double result = {module}_dgdp(T, P);
        return result;
    }}
    
    static double {module}_cv(double T, double P) {{
        double result = -T*{module}_d2gdt2(T, P);
        double dvdt = {module}_d2gdtdp(T, P);
        double dvdp = {module}_d2gdp2(T, P);
        result += T*dvdt*dvdt/dvdp;
        return result;
    }}
    
    static double {module}_cp(double T, double P) {{
        double result = -T*{module}_d2gdt2(T, P);
        return result;
    }}
    
    static double {module}_dcpdt(double T, double P) {{
        double result = -T*{module}_d3gdt3(T, P) - {module}_d2gdt2(T, P);
        return result;
    }}
    
    static double {module}_alpha(double T, double P) {{
        double result = {module}_d2gdtdp(T, P)/{module}_dgdp(T, P);
        return result;
    }}
    
    static double {module}_beta(double T, double P) {{
        double result = -{module}_d2gdp2(T, P)/{module}_dgdp(T, P);
        return result;
    }}
    
    static double {module}_K(double T, double P) {{
        double result = -{module}_dgdp(T, P)/{module}_d2gdp2(T, P);
        return result;
    }}
    
    static double {module}_Kp(double T, double P) {{
        double result = {module}_dgdp(T, P);
        result *= {module}_d3gdp3(T, P);
        result /= pow({module}_d2gdp2(T, P), 2.0);
        return result - 1.0;
    }}
    
    \
    """

.. code:: ipython3

    berman_calc += convenience_template.format(module=module)

Write contents to berman_calc.h

.. code:: ipython3

    #print(berman_calc)
    with open('berman_calc.h', 'w') as f:
        f.write(berman_calc)

(2)
~~~

Function template for primary routines in berman_calc.h

.. code:: ipython3

    berman_calib_template = """\
    
    #include <math.h>
    
    static double {module}_dparam_{func}(double T, double P, int index) {{
        double result = 0.0;
        switch (index) {{
        case 0:
            result += {g_code[0]};
            break;
        case 1:
            result += {g_code[1]};
            break;
        case 2:
            result += {g_code[2]};
            break;
        case 3:
            result += {g_code[3]};
            break;
        case 4:
            result += {g_code[4]};
            break;
        case 5:
            result += {g_code[5]};
            break;
        case 6:
            result += {g_code[6]};
            break;
        case 7:
            result += {g_code[7]};
            break;
        case 8:
            result += {g_code[8]};
            break;
        case 9:
            result += {g_code[9]};
            break;
        case 10:
            result += {g_code[10]};
            break;
        case 25:
            result += {g_code[25]};
            break;
        case 26:
            result += {g_code[26]};
            break;
        default:
            break;
        }}
        if ((l1 != 0.0) || (l2 != 0.0)) {{
          if ((({tlmin}) < T) && (({tlmax}) > T)) {{
            switch (index) {{
            case 11:
                result += {gl_code[11]};
                break;
            case 12:
                result += {gl_code[12]};
                break;
            case 13:
                result += {gl_code[13]};
                break;
            case 14:
                result += {gl_code[14]};
                break;
            case 15:
                result += {gl_code[15]};
                break;
            default:
                break;
            }}
          }} else if (T >= ({tlmax})) {{
            switch (index) {{
            case 11:
                result += {glmax_code[11]};
                break;
            case 12:
                result += {glmax_code[12]};
                break;
            case 13:
                result += {glmax_code[13]};
                break;
            case 14:
                result += {glmax_code[14]};
                break;
            case 15:
                result += {glmax_code[15]};
                break;
            default:
                break;
            }}
          }}
        }}
        if ((H_t != 0.0) && (T >= ({tlmax}))) {{
          switch (index) {{
            case 16:
                result += {glabove_code[16]};
                break;
            default:
                break;
            }}
        }}
        if ((d0 != 0.0) || (d1 != 0.0) || (d2 != 0.0) || (d3 != 0.0) || (d4 != 0.0) || (d5 != 0.0)) {{
          if ((T_D_ref < T) && (T_D > T)) {{
            switch (index) {{
            case 17:
                result += {gd_code[17]};
                break;
            case 18:
                result += {gd_code[18]};
                break;
            case 19:
                result += {gd_code[19]};
                break;
            case 20:
                result += {gd_code[20]};
                break;
            case 21:
                result += {gd_code[21]};
                break;
            case 22:
                result += {gd_code[22]};
                break;
            case 23:
                result += {gd_code[23]};
                break;
            case 24:
                result += {gd_code[24]};
                break;
            default:
                break;
            }}
          }} else if (T >= T_D) {{
            switch (index) {{
            case 17:
                result += {gdmax_code[17]};
                result += {gdabove_code[17]};
                break;
            case 18:
                result += {gdmax_code[18]};
                result += {gdabove_code[18]};
                break;
            case 19:
                result += {gdmax_code[19]};
                result += {gdabove_code[19]};
                break;
            case 20:
                result += {gdmax_code[20]};
                result += {gdabove_code[20]};
                break;
            case 21:
                result += {gdmax_code[21]};
                result += {gdabove_code[21]};
                break;
            case 22:
                result += {gdmax_code[22]};
                result += {gdabove_code[22]};
                break;
            case 23:
                result += {gdmax_code[23]};
                result += {gdabove_code[23]};
                break;
            case 24:
                result += {gdmax_code[24]};
                result += {gdabove_code[24]};
                break;
            default:
                break;
            }}
          }}
        }}
        return result;
    }}
    \
    """

.. code:: ipython3

    symparam = HTrPr, STrPr, k0, k1, k2, k3, VTrPr, v1, v2, v3, v4, l1, l2, kl, TlPr, Tlref, deltaHt, \
               d0, d1, d2, d3, d4, d5, TD, TDref, Tr, Pr

.. code:: ipython3

    G_param_jac = sym.Matrix(G_matrix).jacobian(symparam)
    Gl_param_jac = sym.Matrix(Gl_matrix).jacobian(symparam)
    GatTl_param_jac = sym.Matrix(GatTl_matrix).jacobian(symparam)
    GaboveTl_param_jac = sym.Matrix(GaboveTl_matrix).jacobian(symparam)
    GDs_param_jac = sym.Matrix(GDs_matrix).jacobian(symparam)
    GDsAtTD_param_jac = sym.Matrix(GDsAtTD_matrix).jacobian(symparam)
    GDsAboveTD_param_jac = sym.Matrix(GDsAboveTD_matrix).jacobian(symparam)

.. code:: ipython3

    berman_calib = ''
    
    for j in range(0,len(f_list)):
        G_jac_list = [ printer.doprint(G_param_jac[j,i]) for i in range(0, len(params)) ]
        Gl_jac_list = [ printer.doprint(Gl_param_jac[j,i]) for i in range(0, len(params)) ]
        GatTl_jac_list = [ printer.doprint(GatTl_param_jac[j,i]) for i in range(0, len(params)) ]
        GaboveTl_jac_list = [ printer.doprint(GaboveTl_param_jac[j,i]) for i in range(0, len(params)) ]
        GDs_jac_list = [ printer.doprint(GDs_param_jac[j,i]) for i in range(0, len(params)) ]
        GDsAtTD_jac_list = [ printer.doprint(GDsAtTD_param_jac[j,i]) for i in range(0, len(params)) ]
        GDsAboveTD_jac_list = [ printer.doprint(GDsAboveTD_param_jac[j,i]) for i in range(0, len(params)) ]
    
        berman_calib += berman_calib_template.format(\
                                  module=module,
                                  func=f_list[j],
                                  g_code=G_jac_list, \
                                  tlmin=printer.doprint(tr), \
                                  tlmax=printer.doprint(Tl), \
                                  gl_code=Gl_jac_list, \
                                  glmax_code=GatTl_jac_list, \
                                  glabove_code=GaboveTl_jac_list, \
                                  gd_code=GDs_jac_list, \
                                  gdmax_code=GDsAtTD_jac_list, \
                                  gdabove_code=GDsAboveTD_jac_list
                                 )

Add convenience routines to berman_calib.h

.. code:: ipython3

    extra_template = """\
    
    static int {module}_get_param_number(void) {{
        return {number_params};
    }}
    
    static const char *paramNames[{number_params}] = {names_params};
    static const char *paramUnits[{number_params}] = {units_params};
    
    static const char **{module}_get_param_names(void) {{
        return paramNames;
    }}
    
    static const char **{module}_get_param_units(void) {{
        return paramUnits;
    }}
    
    static void {module}_get_param_values(double *values) {{
        values[ 0] = {value_params[0]};
        values[ 1] = {value_params[1]};
        values[ 2] = {value_params[2]};
        values[ 3] = {value_params[3]};
        values[ 4] = {value_params[4]};
        values[ 5] = {value_params[5]};
        values[ 6] = {value_params[6]};
        values[ 7] = {value_params[7]};
        values[ 8] = {value_params[8]};
        values[ 9] = {value_params[9]};
        values[10] = {value_params[10]};
        values[11] = {value_params[11]};
        values[12] = {value_params[12]};
        values[13] = {value_params[13]};
        values[14] = {value_params[14]};
        values[15] = {value_params[15]};
        values[16] = {value_params[16]};
        values[17] = {value_params[17]};
        values[18] = {value_params[18]};
        values[19] = {value_params[19]};
        values[20] = {value_params[20]};
        values[21] = {value_params[21]};
        values[22] = {value_params[22]};
        values[23] = {value_params[23]};
        values[24] = {value_params[24]};
        values[25] = {value_params[25]};
        values[26] = {value_params[26]};
    }}
    
    static int {module}_set_param_values(double *values) {{
        {value_params[0]}  = values[ 0];
        {value_params[1]}  = values[ 1];
        {value_params[2]}  = values[ 2];
        {value_params[3]}  = values[ 3];
        {value_params[4]}  = values[ 4];
        {value_params[5]}  = values[ 5];
        {value_params[6]}  = values[ 6];
        {value_params[7]}  = values[ 7];
        {value_params[8]}  = values[ 8];
        {value_params[9]}  = values[ 9];
        {value_params[10]} = values[10];
        {value_params[11]} = values[11];
        {value_params[12]} = values[12];
        {value_params[13]} = values[13];
        {value_params[14]} = values[14];
        {value_params[15]} = values[15];
        {value_params[16]} = values[16];
        {value_params[17]} = values[17];
        {value_params[18]} = values[18];
        {value_params[19]} = values[19];
        {value_params[20]} = values[20];
        {value_params[21]} = values[21];
        {value_params[22]} = values[22];
        {value_params[23]} = values[23];
        {value_params[24]} = values[24];
        {value_params[25]} = values[25];
        {value_params[26]} = values[26];
        return 1;
    }}
    
    static double {module}_get_param_value(int index) {{
        double result = 0.0;
        switch (index) {{
        case 0:
            result = {value_params[0]};
            break;
        case 1:
            result = {value_params[1]};
            break;
        case 2:
            result = {value_params[2]};
            break;
        case 3:
            result = {value_params[3]};
            break;
        case 4:
            result = {value_params[4]};
            break;
        case 5:
            result = {value_params[5]};
            break;
        case 6:
            result = {value_params[6]};
            break;
        case 7:
            result = {value_params[7]};
            break;
        case 8:
            result = {value_params[8]};
            break;
        case 9:
            result = {value_params[9]};
            break;
        case 10:
            result = {value_params[10]};
            break;
        case 11:
            result = {value_params[11]};
            break;
        case 12:
            result = {value_params[12]};
            break;
        case 13:
            result = {value_params[13]};
            break;
        case 14:
            result = {value_params[14]};
            break;
        case 15:
            result = {value_params[15]};
            break;
        case 16:
            result = {value_params[16]};
            break;
        case 17:
            result = {value_params[17]};
            break;
        case 18:
            result = {value_params[18]};
            break;
        case 19:
            result = {value_params[19]};
            break;
        case 20:
            result = {value_params[20]};
            break;
        case 21:
            result = {value_params[21]};
            break;
        case 22:
            result = {value_params[22]};
            break;
        case 23:
            result = {value_params[23]};
            break;
        case 24:
            result = {value_params[24]};
            break;
        case 25:
            result = {value_params[25]};
            break;
        case 26:
            result = {value_params[26]};
            break;
        default:
            break;
        }}
        return result;
    }}
    
    static int {module}_set_param_value(int index, double value) {{
        int result = 1;
        switch (index) {{
        case 0:
            {value_params[0]}  = value;
            break;
        case 1:
            {value_params[1]}  = value;
            break;
        case 2:
            {value_params[2]}  = value;
            break;
        case 3:
            {value_params[3]}  = value;
            break;
        case 4:
            {value_params[4]}  = value;
            break;
        case 5:
            {value_params[5]}  = value;
            break;
        case 6:
            {value_params[6]}  = value;
            break;
        case 7:
            {value_params[7]}  = value;
            break;
        case 8:
            {value_params[8]}  = value;
            break;
        case 9:
            {value_params[9]}  = value;
            break;
        case 10:
            {value_params[10]} = value;
            break;
        case 11:
            {value_params[11]} = value;
            break;
        case 12:
            {value_params[12]} = value;
            break;
        case 13:
            {value_params[13]} = value;
            break;
        case 14:
            {value_params[14]} = value;
            break;
        case 15:
            {value_params[15]} = value;
            break;
        case 16:
            {value_params[16]} = value;
            break;
        case 17:
            {value_params[17]} = value;
            break;
        case 18:
            {value_params[18]} = value;
            break;
        case 19:
            {value_params[19]} = value;
            break;
        case 20:
            {value_params[20]} = value;
            break;
        case 21:
            {value_params[21]} = value;
            break;
        case 22:
            {value_params[22]} = value;
            break;
        case 23:
            {value_params[23]} = value;
            break;
        case 24:
            {value_params[24]} = value;
            break;
        case 25:
            {value_params[25]} = value;
            break;
        case 26:
            {value_params[26]} = value;
            break;
        default:
            result = 0;
            break;
        }}
        return result;
    }}
    
    \
    """

.. code:: ipython3

    import json
    berman_calib += extra_template.format(module=module, \
                                          number_params=len(params), \
                                          names_params=json.dumps(params).replace('[', '{').replace(']', '}'), \
                                          units_params=json.dumps(units).replace('[', '{').replace(']', '}'), \
                                          value_params=[ printer.doprint(symparam[i]) for i in range(0, len(params)) ])

.. code:: ipython3

    #print(berman_calib)
    with open('berman_calib.h', 'w') as f:
        f.write(berman_calib)

Get the Berman (1988) database as a Panda’s data frame
======================================================

.. code:: ipython3

    berman_df = pd.read_json('../berman_1988.json')
    berman_df['T_r'] = 298.15 
    berman_df['P_r'] = 1.0
    berman_df = berman_df[['Phase', 'Formula'] + params]
    berman_df.fillna(0, inplace=True)
    #print(berman_df.head())
    #print(list(berman_df))
    row = berman_df.iloc[42,:]
    print(row.Phase.title())
    print(row.Formula.title().replace("(1)","").replace("(", "").replace(")",""))
    print(params[0], row[params[0]])
    #print(printer.doprint(Assignment(symparam[0], row[params[0]])))
    print(row[2:])


.. parsed-literal::

    Potassium_Feldspar
    KAlSi3O8
    H_TrPr -3970790.781
    H_TrPr         -3.97079e+06
    S_TrPr              214.145
    k0                  381.372
    k1                 -1941.05
    k2                -12037252
    k3               1836425472
    V_TrPr               10.869
    v1              -1.8045e-06
    v2                        0
    v3              1.51451e-05
    v4                  5.5e-09
    l1                        0
    l2                        0
    k_lambda                  0
    T_lambda_Pr               0
    T_lambda_ref              0
    H_t                       0
    d0                  282.983
    d1                 -4831.38
    d2              3.62071e+06
    d3                 -0.15733
    d4                3.477e-05
    d5                   410630
    T_D                 1436.15
    T_D_ref              298.15
    T_r                  298.15
    P_r                       1
    Name: 47, dtype: object


.. code:: ipython3

    elements = pd.read_csv('../elements.csv')
    entropies = pd.read_csv('../entropies.csv')
    print(elements.head())


.. parsed-literal::

            Name Abbrv       MW
    0      dummy     !  0.00000
    1   hydrogen     H  1.00790
    2     helium    He  4.00260
    3    lithium    Li  6.94000
    4  beryllium    Be  9.01218


.. code:: ipython3

    def parse_formula(formula_string="Si(1)O(2)"):
        formula = row.Formula.title().replace('(',',').replace(')',',').split(',')
        mw = 0.0
        elmvector = np.zeros(120)
        for i in range(0,len(formula)-1,2):
            element = elements.loc[elements['Abbrv'] == formula[i]]
            ind = element.index.values[0]
            mw += element.MW.values[0]*float(formula[i+1])
            elmvector[ind] = float(formula[i+1])
        return mw, elmvector
    mw, elmvector = parse_formula(row.Formula)
    print (mw)
    print (elmvector)


.. parsed-literal::

    278.33524
    [0. 0. 0. 0. 0. 0. 0. 0. 8. 0. 0. 0. 0. 1. 3. 0. 0. 0. 0. 1. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]


.. code:: ipython3

    fast_c_template = """\
    
    static const int identifier = {git_identifier};
    static const double {param_names[0]} = {param_values[0]};
    static const double {param_names[1]} = {param_values[1]};
    static const double {param_names[2]} = {param_values[2]};
    static const double {param_names[3]} = {param_values[3]};
    static const double {param_names[4]} = {param_values[4]};
    static const double {param_names[5]} = {param_values[5]};
    static const double {param_names[6]} = {param_values[6]};
    static const double {param_names[7]} = {param_values[7]};
    static const double {param_names[8]} = {param_values[8]};
    static const double {param_names[9]} = {param_values[9]};
    static const double {param_names[10]} = {param_values[10]};
    static const double {param_names[11]} = {param_values[11]};
    static const double {param_names[12]} = {param_values[12]};
    static const double {param_names[13]} = {param_values[13]};
    static const double {param_names[14]} = {param_values[14]};
    static const double {param_names[15]} = {param_values[15]};
    static const double {param_names[16]} = {param_values[16]};
    static const double {param_names[17]} = {param_values[17]};
    static const double {param_names[18]} = {param_values[18]};
    static const double {param_names[19]} = {param_values[19]};
    static const double {param_names[20]} = {param_values[20]};
    static const double {param_names[21]} = {param_values[21]};
    static const double {param_names[22]} = {param_values[22]};
    static const double {param_names[23]} = {param_values[23]};
    static const double {param_names[24]} = {param_values[24]};
    static const double {param_names[25]} = {param_values[25]};
    static const double {param_names[26]} = {param_values[26]};
    
    #include "berman_calc.h"
    
    const int {phase}_{module}_identifier(void) {{
        return identifier;
    }}
    
    const char *{phase}_{module}_name(void) {{
        return "{phase}";
    }}
    
    const char *{phase}_{module}_formula(void) {{
        return "{formula}";
    }}
    
    const double {phase}_{module}_mw(void) {{
        return {mw};
    }}
    
    static const double elmformula[106] = {{
            {elmvector[0]},{elmvector[1]},{elmvector[2]},{elmvector[3]},{elmvector[4]},{elmvector[5]},
            {elmvector[6]},{elmvector[7]},{elmvector[8]},{elmvector[9]},{elmvector[10]},{elmvector[11]},
            {elmvector[12]},{elmvector[13]},{elmvector[14]},{elmvector[15]},{elmvector[16]},{elmvector[17]},
            {elmvector[18]},{elmvector[19]},{elmvector[20]},{elmvector[21]},{elmvector[22]},{elmvector[23]},
            {elmvector[24]},{elmvector[25]},{elmvector[26]},{elmvector[27]},{elmvector[28]},{elmvector[29]},
            {elmvector[30]},{elmvector[31]},{elmvector[32]},{elmvector[33]},{elmvector[34]},{elmvector[35]},
            {elmvector[36]},{elmvector[37]},{elmvector[38]},{elmvector[39]},{elmvector[40]},{elmvector[41]},
            {elmvector[42]},{elmvector[43]},{elmvector[44]},{elmvector[45]},{elmvector[46]},{elmvector[47]},
            {elmvector[48]},{elmvector[49]},{elmvector[50]},{elmvector[51]},{elmvector[52]},{elmvector[53]},
            {elmvector[54]},{elmvector[55]},{elmvector[56]},{elmvector[57]},{elmvector[58]},{elmvector[59]},
            {elmvector[60]},{elmvector[61]},{elmvector[62]},{elmvector[63]},{elmvector[64]},{elmvector[65]},
            {elmvector[66]},{elmvector[67]},{elmvector[68]},{elmvector[69]},{elmvector[70]},{elmvector[71]},
            {elmvector[72]},{elmvector[73]},{elmvector[74]},{elmvector[75]},{elmvector[76]},{elmvector[77]},
            {elmvector[78]},{elmvector[79]},{elmvector[80]},{elmvector[81]},{elmvector[82]},{elmvector[83]},
            {elmvector[84]},{elmvector[85]},{elmvector[86]},{elmvector[87]},{elmvector[88]},{elmvector[89]},
            {elmvector[90]},{elmvector[91]},{elmvector[92]},{elmvector[93]},{elmvector[94]},{elmvector[95]},
            {elmvector[96]},{elmvector[97]},{elmvector[98]},{elmvector[99]},{elmvector[100]},{elmvector[101]},
            {elmvector[102]},{elmvector[103]},{elmvector[104]},{elmvector[105]}
        }};
    
    const double *{phase}_{module}_elements(void) {{
        return elmformula;
    }}
    
    double {phase}_{module}_g(double T, double P) {{
        return {module}_g(T, P);
    }}
    
    double {phase}_{module}_dgdt(double T, double P) {{
        return {module}_dgdt(T, P);
    }}
    
    double {phase}_{module}_dgdp(double T, double P) {{
        return {module}_dgdp(T, P);
    }}
    
    double {phase}_{module}_d2gdt2(double T, double P) {{
        return {module}_d2gdt2(T, P);
    }}
    
    double {phase}_{module}_d2gdtdp(double T, double P) {{
        return {module}_d2gdtdp(T, P);
    }}
    
    double {phase}_{module}_d2gdp2(double T, double P) {{
        return {module}_d2gdp2(T, P);
    }}
    
    double {phase}_{module}_d3gdt3(double T, double P) {{
        return {module}_d3gdt3(T, P);
    }}
    
    double {phase}_{module}_d3gdt2dp(double T, double P) {{
        return {module}_d3gdt2dp(T, P);
    }}
    
    double {phase}_{module}_d3gdtdp2(double T, double P) {{
        return {module}_d3gdtdp2(T, P);
    }}
    
    double {phase}_{module}_d3gdp3(double T, double P) {{
        return {module}_d3gdp3(T, P);
    }}
    
    double {phase}_{module}_s(double T, double P) {{
        return {module}_s(T, P);
    }}
    
    double {phase}_{module}_v(double T, double P) {{
        return {module}_v(T, P);
    }}
    
    double {phase}_{module}_cv(double T, double P) {{
        return {module}_cv(T, P);
    }}
    
    double {phase}_{module}_cp(double T, double P) {{
        return {module}_cp(T, P);
    }}
    
    double {phase}_{module}_dcpdt(double T, double P) {{
        return {module}_dcpdt(T, P);
    }}
    
    double {phase}_{module}_alpha(double T, double P) {{
        return {module}_alpha(T, P);
    }}
    
    double {phase}_{module}_beta(double T, double P) {{
        return {module}_beta(T, P);
    }}
    
    double {phase}_{module}_K(double T, double P) {{
        return {module}_K(T, P);
    }}
    
    double {phase}_{module}_Kp(double T, double P) {{
        return {module}_Kp(T, P);
    }}
    
    \
    """

.. code:: ipython3

    calib_c_template = """\
    
    static int identifier = {git_identifier};
    static double {param_names[0]} = {param_values[0]};
    static double {param_names[1]} = {param_values[1]};
    static double {param_names[2]} = {param_values[2]};
    static double {param_names[3]} = {param_values[3]};
    static double {param_names[4]} = {param_values[4]};
    static double {param_names[5]} = {param_values[5]};
    static double {param_names[6]} = {param_values[6]};
    static double {param_names[7]} = {param_values[7]};
    static double {param_names[8]} = {param_values[8]};
    static double {param_names[9]} = {param_values[9]};
    static double {param_names[10]} = {param_values[10]};
    static double {param_names[11]} = {param_values[11]};
    static double {param_names[12]} = {param_values[12]};
    static double {param_names[13]} = {param_values[13]};
    static double {param_names[14]} = {param_values[14]};
    static double {param_names[15]} = {param_values[15]};
    static double {param_names[16]} = {param_values[16]};
    static double {param_names[17]} = {param_values[17]};
    static double {param_names[18]} = {param_values[18]};
    static double {param_names[19]} = {param_values[19]};
    static double {param_names[20]} = {param_values[20]};
    static double {param_names[21]} = {param_values[21]};
    static double {param_names[22]} = {param_values[22]};
    static double {param_names[23]} = {param_values[23]};
    static double {param_names[24]} = {param_values[24]};
    static double {param_names[25]} = {param_values[25]};
    static double {param_names[26]} = {param_values[26]};
    
    #include "berman_calc.h"
    #include "berman_calib.h"
    
    const int {phase}_{module}_calib_identifier(void) {{
        return identifier;
    }}
    
    const char *{phase}_{module}_calib_name(void) {{
        return "{phase}";
    }}
    
    const char *{phase}_{module}_calib_formula(void) {{
        return "{formula}";
    }}
    
    const double {phase}_{module}_calib_mw(void) {{
        return {mw};
    }}
    
    static const double elmformula[106] = {{
            {elmvector[0]},{elmvector[1]},{elmvector[2]},{elmvector[3]},{elmvector[4]},{elmvector[5]},
            {elmvector[6]},{elmvector[7]},{elmvector[8]},{elmvector[9]},{elmvector[10]},{elmvector[11]},
            {elmvector[12]},{elmvector[13]},{elmvector[14]},{elmvector[15]},{elmvector[16]},{elmvector[17]},
            {elmvector[18]},{elmvector[19]},{elmvector[20]},{elmvector[21]},{elmvector[22]},{elmvector[23]},
            {elmvector[24]},{elmvector[25]},{elmvector[26]},{elmvector[27]},{elmvector[28]},{elmvector[29]},
            {elmvector[30]},{elmvector[31]},{elmvector[32]},{elmvector[33]},{elmvector[34]},{elmvector[35]},
            {elmvector[36]},{elmvector[37]},{elmvector[38]},{elmvector[39]},{elmvector[40]},{elmvector[41]},
            {elmvector[42]},{elmvector[43]},{elmvector[44]},{elmvector[45]},{elmvector[46]},{elmvector[47]},
            {elmvector[48]},{elmvector[49]},{elmvector[50]},{elmvector[51]},{elmvector[52]},{elmvector[53]},
            {elmvector[54]},{elmvector[55]},{elmvector[56]},{elmvector[57]},{elmvector[58]},{elmvector[59]},
            {elmvector[60]},{elmvector[61]},{elmvector[62]},{elmvector[63]},{elmvector[64]},{elmvector[65]},
            {elmvector[66]},{elmvector[67]},{elmvector[68]},{elmvector[69]},{elmvector[70]},{elmvector[71]},
            {elmvector[72]},{elmvector[73]},{elmvector[74]},{elmvector[75]},{elmvector[76]},{elmvector[77]},
            {elmvector[78]},{elmvector[79]},{elmvector[80]},{elmvector[81]},{elmvector[82]},{elmvector[83]},
            {elmvector[84]},{elmvector[85]},{elmvector[86]},{elmvector[87]},{elmvector[88]},{elmvector[89]},
            {elmvector[90]},{elmvector[91]},{elmvector[92]},{elmvector[93]},{elmvector[94]},{elmvector[95]},
            {elmvector[96]},{elmvector[97]},{elmvector[98]},{elmvector[99]},{elmvector[100]},{elmvector[101]},
            {elmvector[102]},{elmvector[103]},{elmvector[104]},{elmvector[105]}
        }};
    
    const double *{phase}_{module}_calib_elements(void) {{
        return elmformula;
    }}
    
    double {phase}_{module}_calib_g(double T, double P) {{
        return {module}_g(T, P);
    }}
    
    double {phase}_{module}_calib_dgdt(double T, double P) {{
        return {module}_dgdt(T, P);
    }}
    
    double {phase}_{module}_calib_dgdp(double T, double P) {{
        return {module}_dgdp(T, P);
    }}
    
    double {phase}_{module}_calib_d2gdt2(double T, double P) {{
        return {module}_d2gdt2(T, P);
    }}
    
    double {phase}_{module}_calib_d2gdtdp(double T, double P) {{
        return {module}_d2gdtdp(T, P);
    }}
    
    double {phase}_{module}_calib_d2gdp2(double T, double P) {{
        return {module}_d2gdp2(T, P);
    }}
    
    double {phase}_{module}_calib_d3gdt3(double T, double P) {{
        return {module}_d3gdt3(T, P);
    }}
    
    double {phase}_{module}_calib_d3gdt2dp(double T, double P) {{
        return {module}_d3gdt2dp(T, P);
    }}
    
    double {phase}_{module}_calib_d3gdtdp2(double T, double P) {{
        return {module}_d3gdtdp2(T, P);
    }}
    
    double {phase}_{module}_calib_d3gdp3(double T, double P) {{
        return {module}_d3gdp3(T, P);
    }}
    
    double {phase}_{module}_calib_s(double T, double P) {{
        return {module}_s(T, P);
    }}
    
    double {phase}_{module}_calib_v(double T, double P) {{
        return {module}_v(T, P);
    }}
    
    double {phase}_{module}_calib_cv(double T, double P) {{
        return {module}_cv(T, P);
    }}
    
    double {phase}_{module}_calib_cp(double T, double P) {{
        return {module}_cp(T, P);
    }}
    
    double {phase}_{module}_calib_dcpdt(double T, double P) {{
        return {module}_dcpdt(T, P);
    }}
    
    double {phase}_{module}_calib_alpha(double T, double P) {{
        return {module}_alpha(T, P);
    }}
    
    double {phase}_{module}_calib_beta(double T, double P) {{
        return {module}_beta(T, P);
    }}
    
    double {phase}_{module}_calib_K(double T, double P) {{
        return {module}_K(T, P);
    }}
    
    double {phase}_{module}_calib_Kp(double T, double P) {{
        return {module}_Kp(T, P);
    }}
    
    int {phase}_{module}_get_param_number(void) {{
        return {module}_get_param_number();
    }}
    
    const char **{phase}_{module}_get_param_names(void) {{
        return {module}_get_param_names();
    }}
    
    const char **{phase}_{module}_get_param_units(void) {{
        return {module}_get_param_units();
    }}
    
    void {phase}_{module}_get_param_values(double *values) {{
        {module}_get_param_values(values);
    }}
    
    int {phase}_{module}_set_param_values(double *values) {{
        return {module}_set_param_values(values);
    }}
    
    double {phase}_{module}_get_param_value(int index) {{
        return {module}_get_param_value(index);
    }}
    
    int {phase}_{module}_set_param_value(int index, double value) {{
        return {module}_set_param_value(index, value);
    }}
    
    double {phase}_{module}_dparam_g(double T, double P, int index) {{
        return {module}_dparam_g(T, P, index);
    }}
    
    double {phase}_{module}_dparam_dgdt(double T, double P, int index) {{
        return {module}_dparam_dgdt(T, P, index);
    }}
    
    double {phase}_{module}_dparam_dgdp(double T, double P, int index) {{
        return {module}_dparam_dgdp(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d2gdt2(double T, double P, int index) {{
        return {module}_dparam_d2gdt2(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d2gdtdp(double T, double P, int index) {{
        return {module}_dparam_d2gdtdp(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d2gdp2(double T, double P, int index) {{
        return {module}_dparam_d2gdp2(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d3gdt3(double T, double P, int index) {{
        return {module}_dparam_d3gdt3(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d3gdt2dp(double T, double P, int index) {{
        return {module}_dparam_d3gdt2dp(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d3gdtdp2(double T, double P, int index) {{
        return {module}_dparam_d3gdtdp2(T, P, index);
    }}
    
    double {phase}_{module}_dparam_d3gdp3(double T, double P, int index) {{
        return {module}_dparam_d3gdp3(T, P, index);
    }}
    
    \
    """

.. code:: ipython3

    fast_h_template = """\
    
    const char *{phase}_{module}_name(void);
    const char *{phase}_{module}_formula(void);
    const double {phase}_{module}_mw(void);
    const double *{phase}_{module}_elements(void);
    
    double {phase}_{module}_g(double T, double P);
    double {phase}_{module}_dgdt(double T, double P);
    double {phase}_{module}_dgdp(double T, double P);
    double {phase}_{module}_d2gdt2(double T, double P);
    double {phase}_{module}_d2gdtdp(double T, double P);
    double {phase}_{module}_d2gdp2(double T, double P);
    double {phase}_{module}_d3gdt3(double T, double P);
    double {phase}_{module}_d3gdt2dp(double T, double P);
    double {phase}_{module}_d3gdtdp2(double T, double P);
    double {phase}_{module}_d3gdp3(double T, double P);
    
    double {phase}_{module}_s(double T, double P);
    double {phase}_{module}_v(double T, double P);
    double {phase}_{module}_cv(double T, double P);
    double {phase}_{module}_cp(double T, double P);
    double {phase}_{module}_dcpdt(double T, double P);
    double {phase}_{module}_alpha(double T, double P);
    double {phase}_{module}_beta(double T, double P);
    double {phase}_{module}_K(double T, double P);
    double {phase}_{module}_Kp(double T, double P);
    
    \
    """

.. code:: ipython3

    calib_h_template = """\
    
    const char *{phase}_{module}_calib_name(void);
    const char *{phase}_{module}_calib_formula(void);
    const double {phase}_{module}_calib_mw(void);
    const double *{phase}_{module}_calib_elements(void);
    
    double {phase}_{module}_calib_g(double T, double P);
    double {phase}_{module}_calib_dgdt(double T, double P);
    double {phase}_{module}_calib_dgdp(double T, double P);
    double {phase}_{module}_calib_d2gdt2(double T, double P);
    double {phase}_{module}_calib_d2gdtdp(double T, double P);
    double {phase}_{module}_calib_d2gdp2(double T, double P);
    double {phase}_{module}_calib_d3gdt3(double T, double P);
    double {phase}_{module}_calib_d3gdt2dp(double T, double P);
    double {phase}_{module}_calib_d3gdtdp2(double T, double P);
    double {phase}_{module}_calib_d3gdp3(double T, double P);
    
    double {phase}_{module}_calib_s(double T, double P);
    double {phase}_{module}_calib_v(double T, double P);
    double {phase}_{module}_calib_cv(double T, double P);
    double {phase}_{module}_calib_cp(double T, double P);
    double {phase}_{module}_calib_dcpdt(double T, double P);
    double {phase}_{module}_calib_alpha(double T, double P);
    double {phase}_{module}_calib_beta(double T, double P);
    double {phase}_{module}_calib_K(double T, double P);
    double {phase}_{module}_calib_Kp(double T, double P);
    
    int {phase}_{module}_get_param_number(void);
    const char *{phase}_{module}_get_param_names(void);
    const char *{phase}_{module}_get_param_units(void);
    void {phase}_{module}_get_param_values(double *values);
    int {phase}_{module}_set_param_values(double *values);
    double {phase}_{module}_get_param_value(int index);
    int {phase}_{module}_set_param_value(int index, double value);
    
    double {phase}_{module}_dparam_g(double T, double P, int index);
    double {phase}_{module}_dparam_dgdt(double T, double P, int index);
    double {phase}_{module}_dparam_dgdp(double T, double P, int index);
    double {phase}_{module}_dparam_d2gdt2(double T, double P, int index);
    double {phase}_{module}_dparam_d2gdtdp(double T, double P, int index);
    double {phase}_{module}_dparam_d2gdp2(double T, double P, int index);
    double {phase}_{module}_dparam_d3gdt3(double T, double P, int index);
    double {phase}_{module}_dparam_d3gdt2dp(double T, double P, int index);
    double {phase}_{module}_dparam_d3gdtdp2(double T, double P, int index);
    double {phase}_{module}_dparam_d3gdp3(double T, double P, int index);
    
    \
    """

Generate both fast computation and calibibration code for each phase in
the Berman (1988) database

.. code:: ipython3

    for i in range(0,berman_df.shape[0]): #range(0,1) [0,1,6,48,49,50] 
        row = berman_df.iloc[i,:]
        phase = row.Phase.title().replace("-","_")
        print (phase)
        formula = row.Formula.title().replace("(1)","").replace("(", "").replace(")","")
        mw, elmvector = parse_formula(row.Formula)
        berman_phase_c = ''
        berman_phase_h = ''
        berman_phase_c += fast_c_template.format( \
            module=module, \
            phase=phase, \
            formula=formula, \
            mw=mw, \
            elmvector=elmvector, \
            param_names=[ printer.doprint(symparam[i]) for i in range(0, len(params)) ], \
            param_values=row[2:], \
            git_identifier=1 \
            )
        berman_phase_h += fast_h_template.format(module=module, phase=phase)
        with open(phase + '_berman.h', 'w') as f:
            f.write(berman_phase_h)
        with open(phase + '_berman.c', 'w') as f:
            f.write(berman_phase_c)
        berman_phase_calib_c = ''
        berman_phase_calib_h = ''
        berman_phase_calib_c += calib_c_template.format( \
            module=module, \
            phase=phase, \
            formula=formula, \
            mw=mw, \
            elmvector=elmvector, \
            param_names=[ printer.doprint(symparam[i]) for i in range(0, len(params)) ], \
            param_values=row[2:], \
            git_identifier=1 \
            )
        berman_phase_calib_h += calib_h_template.format(module=module, phase=phase)
        with open(phase + '_berman_calib.h', 'w') as f:
            f.write(berman_phase_calib_h)
        with open(phase + '_berman_calib.c', 'w') as f:
            f.write(berman_phase_calib_c)


.. parsed-literal::

    Akermanite
    Albite
    Ca_Al_Pyroxene
    Calcite
    Chrysotile
    Clinochlore
    Coesite
    Cordierite
    Corundum
    Alpha_Cristobalite
    Beta_Cristobalite
    Diaspore
    High_Albite
    Diopside
    Dolomite
    Clinoenstatite
    Orthoenstatite
    Protoenstatite
    Fayalite
    Ferrosilite
    Forsterite
    Gehlenite
    Grossular
    Low_Albite
    Hematite
    Ilmenite
    Jadeite
    Kaolinite
    Kyanite
    Lawsonite
    Lime
    Magnesite
    Magnetite
    Margarite
    Almandine
    Meionite
    Merwinite
    Monticellite
    Muscovite
    Paragonite
    Periclase
    Phlogopite
    Potassium_Feldspar
    Sanidine
    Microcline
    Andalusite
    Prehnite
    Pyrope
    Pyrophyllite
    A_Quartz
    B_Quartz
    Rutile
    Sillimanite
    Sphene
    Spinel
    Talc
    Anorthite
    Tremolite
    Low_Tridymite
    High_Tridymite
    Wollastonite
    Pseudowollastonite
    Zoisite
    Clinozoisite
    Carbon_Dioxide
    Water
    Oxygen_Gas
    Anthophyllite
    Sulfur_Gas
    Hydrogen_Gas
    Antigorite
    Brucite


Load the Cython Jupyter magic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    %load_ext cython

Generate the C-code and python extension

::

    setuptools.extension.Extension
    self, name, sources, include_dirs=None, define_macros=None, undef_macros=None, library_dirs=None, libraries=None, runtime_library_dirs=None, extra_objects=None, extra_compile_args=None, extra_link_args=None, export_symbols=None, swig_opts=None, depends=None, language=None, optional=None, **kw

::

    -O0, -O1, -O2, -O3, -Ofast, -Os, -Oz, -Og, -O, -O4
    Specify which optimization level to use:

    -O0 Means “no optimization”: this level compiles the fastest and generates the most debuggable code.

    -O1 Somewhere between -O0 and -O2.

    -O2 Moderate level of optimization which enables most optimizations.

    -O3 Like -O2, except that it enables optimizations that take longer to perform or that may generate larger code (in an attempt to make the program run faster).

    -Ofast Enables all the optimizations from -O3 along with other aggressive optimizations that may violate strict compliance with language standards.

    -Os Like -O2 with extra optimizations to reduce code size.

    -Oz Like -Os (and thus -O2), but reduces code size further.

    -Og Like -O1. In future versions, this option might disable different optimizations in order to improve debuggability.

    -O Equivalent to -O2.

    -O4 and higher

.. code:: ipython3

    %%writefile berman.pyxbld
    import numpy
    
    #            module name specified by `%%cython_pyximport` magic
    #            |        just `modname + ".pyx"`
    #            |        |
    def make_ext(modname, pyxfilename):
        from setuptools.extension import Extension
        return Extension(modname,
                         sources=[pyxfilename, 'Akermanite_berman.c', 'Akermanite_berman_calib.c', \
                                  'A_Quartz_berman.c', 'A_Quartz_berman_calib.c', \
                                  'B_Quartz_berman.c', 'B_Quartz_berman_calib.c', \
                                 ],
                         include_dirs=['.', numpy.get_include()], extra_compile_args=['-O3'])


.. parsed-literal::

    Overwriting berman.pyxbld


::

    preallocate our output array
    cdef cnp.ndarray[cnp.double_t, ndim=1] dY = np.empty(y.size, dtype=np.double)
    now call the C function
    c_odes(<double *> y.data, <double *> dY.data)

.. code:: ipython3

    %%cython_pyximport berman
    import numpy as np
    cimport numpy as cnp # cimport gives us access to NumPy's C API
    
    # here we just replicate the function signature from the header
    cdef extern from "Akermanite_berman.h":
        double Akermanite_berman_g(double t, double p)
        double Akermanite_berman_s(double t, double p)
        double Akermanite_berman_cp(double t, double p)
        double Akermanite_berman_v(double t, double p)
    cdef extern from "Akermanite_berman_calib.h":
        double Akermanite_berman_calib_g(double t, double p)
        double Akermanite_berman_calib_s(double t, double p)
        double Akermanite_berman_calib_cp(double t, double p)
        double Akermanite_berman_calib_v(double t, double p)
    cdef extern from "A_Quartz_berman.h":
        double A_Quartz_berman_g(double t, double p)
        double A_Quartz_berman_s(double t, double p)
        double A_Quartz_berman_cp(double t, double p)
        double A_Quartz_berman_v(double t, double p)
    cdef extern from "A_Quartz_berman_calib.h":
        double A_Quartz_berman_calib_g(double t, double p)
        double A_Quartz_berman_calib_s(double t, double p)
        double A_Quartz_berman_calib_cp(double t, double p)
        double A_Quartz_berman_calib_v(double t, double p)
    cdef extern from "B_Quartz_berman.h":
        double B_Quartz_berman_g(double t, double p)
        double B_Quartz_berman_s(double t, double p)
        double B_Quartz_berman_cp(double t, double p)
        double B_Quartz_berman_v(double t, double p)
    cdef extern from "B_Quartz_berman_calib.h":
        double B_Quartz_berman_calib_g(double t, double p)
        double B_Quartz_berman_calib_s(double t, double p)
        double B_Quartz_berman_calib_cp(double t, double p)
        double B_Quartz_berman_calib_v(double t, double p)
    
    # here is the "wrapper" signature
    def cy_Akermanite_berman_g(double t, double p):
        result = Akermanite_berman_g(<double> t, <double> p)
        return result
    def cy_Akermanite_berman_s(double t, double p):
        result = Akermanite_berman_s(<double> t, <double> p)
        return result
    def cy_Akermanite_berman_cp(double t, double p):
        result = Akermanite_berman_cp(<double> t, <double> p)
        return result
    def cy_Akermanite_berman_v(double t, double p):
        result = Akermanite_berman_v(<double> t, <double> p)
        return result
    
    def cy_Akermanite_berman_calib_g(double t, double p):
        result = Akermanite_berman_calib_g(<double> t, <double> p)
        return result
    def cy_Akermanite_berman_calib_s(double t, double p):
        result = Akermanite_berman_calib_s(<double> t, <double> p)
        return result
    def cy_Akermanite_berman_calib_cp(double t, double p):
        result = Akermanite_berman_calib_cp(<double> t, <double> p)
        return result
    def cy_Akermanite_berman_calib_v(double t, double p):
        result = Akermanite_berman_calib_v(<double> t, <double> p)
        return result
    
    def cy_A_Quartz_berman_g(double t, double p):
        result = A_Quartz_berman_g(<double> t, <double> p)
        return result
    def cy_A_Quartz_berman_s(double t, double p):
        result = A_Quartz_berman_s(<double> t, <double> p)
        return result
    def cy_A_Quartz_berman_cp(double t, double p):
        result = A_Quartz_berman_cp(<double> t, <double> p)
        return result
    def cy_A_Quartz_berman_v(double t, double p):
        result = A_Quartz_berman_v(<double> t, <double> p)
        return result
    
    def cy_A_Quartz_berman_calib_g(double t, double p):
        result = A_Quartz_berman_calib_g(<double> t, <double> p)
        return result
    def cy_A_Quartz_berman_calib_s(double t, double p):
        result = A_Quartz_berman_calib_s(<double> t, <double> p)
        return result
    def cy_A_Quartz_berman_calib_cp(double t, double p):
        result = A_Quartz_berman_calib_cp(<double> t, <double> p)
        return result
    def cy_A_Quartz_berman_calib_v(double t, double p):
        result = A_Quartz_berman_calib_v(<double> t, <double> p)
        return result
    
    def cy_B_Quartz_berman_g(double t, double p):
        result = B_Quartz_berman_g(<double> t, <double> p)
        return result
    def cy_B_Quartz_berman_s(double t, double p):
        result = B_Quartz_berman_s(<double> t, <double> p)
        return result
    def cy_B_Quartz_berman_cp(double t, double p):
        result = B_Quartz_berman_cp(<double> t, <double> p)
        return result
    def cy_B_Quartz_berman_v(double t, double p):
        result = B_Quartz_berman_v(<double> t, <double> p)
        return result
    
    def cy_B_Quartz_berman_calib_g(double t, double p):
        result = B_Quartz_berman_calib_g(<double> t, <double> p)
        return result
    def cy_B_Quartz_berman_calib_s(double t, double p):
        result = B_Quartz_berman_calib_s(<double> t, <double> p)
        return result
    def cy_B_Quartz_berman_calib_cp(double t, double p):
        result = B_Quartz_berman_calib_cp(<double> t, <double> p)
        return result
    def cy_B_Quartz_berman_calib_v(double t, double p):
        result = B_Quartz_berman_calib_v(<double> t, <double> p)
        return result

.. code:: ipython3

    %cd ..


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen


Test and time the generated functions for Akermanite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    t = 1000.0
    p = 10000.0

.. code:: ipython3

    newCodeFast = cy_Akermanite_berman_g(t, p)
    newCodeCalib = cy_Akermanite_berman_calib_g(t, p)
    print (newCodeFast, newCodeCalib, newCodeFast-newCodeCalib)
    print (cy_Akermanite_berman_s(t, p))
    print (cy_Akermanite_berman_calib_s(t, p))


.. parsed-literal::

    -4105481.9950771444 -4105481.9950771444 0.0
    523.4611502049585
    523.4611502049585


Compare results with the Phases Python module (wrapped Objective-C code)

.. code:: ipython3

    from thermoengine import phases
    from thermoengine import model
    modelDB = model.Database()
    Akermanite =modelDB.get_phase('Ak')
    print(Akermanite.props['phase_name'])


.. parsed-literal::

    Akermanite


.. code:: ipython3

    oldCode = Akermanite.gibbs_energy(t, p)
    print (oldCode)
    print (oldCode-newCodeFast)
    print (Akermanite.entropy(t, p))


.. parsed-literal::

    -4105478.0080570094
    3.9870201349258423
    523.4518527550034


Instantiate the Objective-C code via Rubicon

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))
    AkermaniteRaw = ObjCClass('AkermaniteBerman')
    obj = AkermaniteRaw.alloc().init()
    print (obj.phaseName)
    print (obj.phaseFormula)


.. parsed-literal::

    Akermanite
    Ca2MgSi2O7


Time the fast code

.. code:: ipython3

    %timeit(cy_Akermanite_berman_g(t, p))


.. parsed-literal::

    125 ns ± 0.345 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)


Time the calibration code

.. code:: ipython3

    %timeit(cy_Akermanite_berman_calib_g(t, p))


.. parsed-literal::

    153 ns ± 0.679 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)


Time the Python+Rubicon wrapped Objective-C code

.. code:: ipython3

    %timeit(Akermanite.gibbs_energy(t, p))


.. parsed-literal::

    84 µs ± 1.16 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)


Time the Rubicon wrapped Objective-C code

.. code:: ipython3

    %timeit(obj.getGibbsFreeEnergyFromT_andP_(t, p))


.. parsed-literal::

    56.3 µs ± 972 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)


Test and time the generated functions for Quartz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    t = 1000.0
    p = 10000.0;

.. code:: ipython3

    QuartzRaw = ObjCClass('QuartzBerman')
    obj = QuartzRaw.alloc().init()
    print (obj.phaseName)
    print (obj.phaseFormula)
    QuartzRaw.disableQuartzCorrectionUsed()
    print (obj.isAlphaPhaseAtT_andP_(t, p))


.. parsed-literal::

    Quartz
    SiO2
    1


.. code:: ipython3

    new_g_quartz_Alpha = cy_A_Quartz_berman_g(t, p)
    new_s_quartz_Alpha = cy_A_Quartz_berman_s(t, p)
    new_cp_quartz_Alpha = cy_A_Quartz_berman_cp(t, p)
    new_v_quartz_Alpha = cy_A_Quartz_berman_v(t, p)
    new_g_quartz_Calib_Alpha = cy_A_Quartz_berman_calib_g(t, p)
    new_g_quartz_Beta = cy_B_Quartz_berman_g(t, p)
    new_s_quartz_Beta = cy_B_Quartz_berman_s(t, p)
    new_cp_quartz_Beta = cy_B_Quartz_berman_cp(t, p)
    new_v_quartz_Beta = cy_B_Quartz_berman_v(t, p)
    new_g_quartz_Calib_Beta = cy_B_Quartz_berman_calib_g(t, p)
    oldCode_g = obj.getGibbsFreeEnergyFromT_andP_(t, p)
    oldCode_s = obj.getEntropyFromT_andP_(t, p)
    oldCode_cp = obj.getHeatCapacityFromT_andP_(t, p)
    oldCode_v = obj.getVolumeFromT_andP_(t, p)

.. code:: ipython3

    print (new_g_quartz_Alpha, new_g_quartz_Calib_Alpha, new_g_quartz_Alpha-new_g_quartz_Calib_Alpha)
    print (new_g_quartz_Beta, new_g_quartz_Calib_Beta, new_g_quartz_Beta-new_g_quartz_Calib_Beta)
    print (new_g_quartz_Alpha-new_g_quartz_Beta)


.. parsed-literal::

    -958058.6858876313 -958058.6858876313 0.0
    -957890.7058826386 -957890.7058826386 0.0
    -167.98000499268528


.. code:: ipython3

    print (oldCode_g-new_g_quartz_Alpha, oldCode_g-new_g_quartz_Beta)


.. parsed-literal::

    8.858848060597666 -159.12115693208762


.. code:: ipython3

    print (oldCode_s-new_s_quartz_Alpha, oldCode_s-new_s_quartz_Beta)


.. parsed-literal::

    -0.0003208644001091443 -2.4145530879386428


.. code:: ipython3

    print (oldCode_cp-new_cp_quartz_Alpha, oldCode_cp-new_cp_quartz_Beta)


.. parsed-literal::

    0.01659417791225337 7.012753660230288


.. code:: ipython3

    oldCode_h = oldCode_g+t*oldCode_s
    new_h_quartz_Alpha = new_g_quartz_Alpha+t*new_s_quartz_Alpha
    new_h_quartz_Beta = new_g_quartz_Beta+t*new_s_quartz_Beta
    print (oldCode_h-new_h_quartz_Alpha, oldCode_h-new_h_quartz_Beta)


.. parsed-literal::

    8.537983660469763 -2573.674244870781


.. code:: ipython3

    print (oldCode_v-new_v_quartz_Alpha, oldCode_v-new_v_quartz_Beta)


.. parsed-literal::

    0.00234952007086342 -0.06335551468245093


.. code:: ipython3

    %timeit(cy_A_Quartz_berman_g(t, p))


.. parsed-literal::

    144 ns ± 0.184 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)


.. code:: ipython3

    %timeit(cy_A_Quartz_berman_calib_g(t, p))


.. parsed-literal::

    177 ns ± 1.4 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)


.. code:: ipython3

    %timeit(obj.getGibbsFreeEnergyFromT_andP_(t, p))


.. parsed-literal::

    56 µs ± 201 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)

