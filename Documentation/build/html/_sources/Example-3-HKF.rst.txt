Revised HKF - Standard State
============================

**This notebook develops code for the revised HKF
(Helgeson-Kirkham-Flowers) formalism for the thermodynamic properties of
aqueous species.**

The derivation follows that in Appendix B of

Shock EL, Oelkers EH, Johnson JW, Sverjensky DA, Helgeson HC (1992)
Calculation of the thermodynamic properties of aqueous species at high
pressures and temperatures. Journal of the Chemical Society Faraday
Transactions, 88(6), 803-826

and in

Tanger JC, Helgeson HC (1988) Calculation of the thermodynamic and
transport properties of aqueous species at high pressures and
temperatures: Revised equations of state for the standard partial molal
properties if ions and electrolytes. American Journal of Science, 288,
19-98

Required system packages and initialization

.. code:: ipython3

    import pandas as pd
    import numpy as np
    import sympy as sym
    sym.init_printing()

Import the coder module and retrieve sympy extensions for the Born
functions and for the Shock et al., (1992) “*g*” function:

.. code:: ipython3

    from thermoengine import coder
    from thermoengine.coder import B, Q, Y, U, N, X, dUdT, dUdP, dNdT, dNdP, dXdT, gSolvent

Create a standard state model instance …

.. code:: ipython3

    model = coder.StdStateModel()

… and retrieve sympy symbols for standard variables.

.. code:: ipython3

    T = model.get_symbol_for_t()
    P = model.get_symbol_for_p()
    Tr = model.get_symbol_for_tr()
    Pr = model.get_symbol_for_pr()

Equation of State (EOS)
-----------------------

Shock et al., 1992, eqn. B10:

:math:`{\bar v^o} = {a_1} + \frac{{{a_2}}}{{\Psi + P}} + \left( {{a_3} + \frac{{{a_4}}}{{\Psi + P}}} \right)\frac{1}{{T - \theta }} - \omega Q - \left( {B + 1} \right){\left( {\frac{{\partial \omega }}{{\partial P}}} \right)_T}`

| :math:`\Psi` has a value of 2600 bars as determined by Tanger and
  Helgeson (1988)
| :math:`\theta` has a value of 228 K as determined by Helgeson and
  Kirkham (1976)

:math:`\omega` is defined in terms of the charge on the ion (*z*) and a
mysterious function, *g*, as derived by Shock et al. (1992),building on
the work of Tanger and Helgeson (1988):

:math:`\omega = \eta z\left( {\frac{z}{{{r_{e,ref}} + \left| z \right|g}} - \frac{1}{{{r_{e,{H^ + }}} + g}}} \right)`

where :math:`\eta` is a conversion constant. The reference ionic radius
is:

:math:`{r_{e,ref}} = \frac{{{z^2}}}{{\frac{{{\omega _{ref}}}}{\eta } + \frac{z}{{{r_{e,{H^ + }}}}}}}`

which when inserted into the definition becomes

:math:`\omega = \left( {\frac{{{\omega _{ref}} + \frac{{z\eta }}{{{r_{e,{H^ + }}}}}}}{{1 + \frac{{\sqrt {{z^2}} }}{{{z^2}}}g\left( {\frac{{{\omega _{ref}}}}{\eta } + \frac{z}{{{r_{e,{H^ + }}}}}} \right)}} - \frac{{\eta z}}{{{r_{e,{H^ + }}} + g}}} \right)`

The above expression is cast into a slightly different form than the one
provided by Shock et al. (1992), but is otherwise identical. This form
is applicable to both charged and neutral species, as demonstrated
below.

.. code:: ipython3

    eta, rH, omega0, z = sym.symbols('eta rH omega0 z', real=True)

.. code:: ipython3

    omega_def = sym.Piecewise((omega0, sym.Eq(z,0)), \
                ((omega0 + z*eta/rH)/(1 + sym.Abs(z)*gSolvent(T,P)*(omega0/eta + z/rH)/z**2) - eta*z/(rH + gSolvent(T,P)), True))
    omega_def




.. math::

    \displaystyle \begin{cases} \omega_{0} & \text{for}\: z = 0 \\- \frac{\eta z}{rH + \operatorname{gSolvent}{\left(T,P \right)}} + \frac{\frac{\eta z}{rH} + \omega_{0}}{1 + \frac{\left(\frac{z}{rH} + \frac{\omega_{0}}{\eta}\right) \left|{z}\right| \operatorname{gSolvent}{\left(T,P \right)}}{z^{2}}} & \text{otherwise} \end{cases}



Hence,
:math:`{\bar v^o} = {a_1} + \frac{{{a_2}}}{{\Psi + P}} + \left( {{a_3} + \frac{{{a_4}}}{{\Psi + P}}} \right)\frac{1}{{T - \theta }} - \omega Q - \left( {B + 1} \right){\left( {\frac{{\partial \omega }}{{\partial P}}} \right)_T}`
is

.. code:: ipython3

    a1,a2,a3,a4 = sym.symbols('a1 a2 a3 a4')
    Psi,theta = sym.symbols('Psi theta')
    omega = sym.Function('omega')(T,P)

.. code:: ipython3

    vtp = a1 + a2/(Psi+P) + (a3 + a4/(Psi+P))/(T-theta) - omega*Q(T,P) - (B(T,P)+1)*omega.diff(P)
    vtp




.. math::

    \displaystyle a_{1} + \frac{a_{2}}{P + \Psi} - \left(B{\left(T,P \right)} + 1\right) \frac{\partial}{\partial P} \omega{\left(T,P \right)} - \omega{\left(T,P \right)} Q{\left(T,P \right)} + \frac{a_{3} + \frac{a_{4}}{P + \Psi}}{T - \theta}



Note that the derivative of $ -
:raw-latex:`\left[ {\left( {{B_{T,P}} + 1} \right){\omega _{T,P}} - \left( {{B_{T,{P_r}}} + 1} \right){\omega _{T,{P_r}}}} \right]`$

.. code:: ipython3

    -(1+B(T,P))*omega + (1+B(T,Pr))*(omega.subs(P,Pr))




.. math::

    \displaystyle \left(- B{\left(T,P \right)} - 1\right) \omega{\left(T,P \right)} + \left(B{\left(T,P_{r} \right)} + 1\right) \omega{\left(T,P_{r} \right)}



is

.. code:: ipython3

    (-(1+B(T,P))*omega + (1+B(T,Pr))*(omega.subs(P,Pr))).diff(P)




.. math::

    \displaystyle \left(- B{\left(T,P \right)} - 1\right) \frac{\partial}{\partial P} \omega{\left(T,P \right)} - \omega{\left(T,P \right)} Q{\left(T,P \right)}



So that :math:`\int_{{P_r}}^P {{{\bar v}^o}} dP` may be written:

.. code:: ipython3

    deltaG = sym.integrate(a1 + a2/(Psi+P) + (a3 + a4/(Psi+P))/(T-theta), (P,Pr,P)) - (1+B(T,P))*omega + (1+B(T,Pr))*omega.subs(P,Pr)
    deltaG




.. math::

    \displaystyle \frac{P \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - \frac{P_{r} \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - \left(B{\left(T,P \right)} + 1\right) \omega{\left(T,P \right)} + \left(B{\left(T,P_{r} \right)} + 1\right) \omega{\left(T,P_{r} \right)} + \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P + \Psi \right)}}{T - \theta} - \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P_{r} + \Psi \right)}}{T - \theta}



The second derivative of this volume integral gives the contribution to
the heat capacity,
i.e. :math:`\frac{{{\partial ^2}G}}{{\partial {T^2}}} = - \frac{{\partial S}}{{\partial T}} = - \frac{{{C_P}}}{T}`
so that

:math:`{{\bar c}_P}\left( {T,P} \right) - {{\bar c}_P}\left( {T,{P_r}} \right) = - T\frac{{{\partial ^2}\int_{{P_r}}^P {{{\bar v}^o}} dP}}{{\partial {T^2}}}`

.. code:: ipython3

    -T*deltaG.diff(T,2)




.. math::

    \displaystyle - T \left(- \frac{2 P a_{1}}{\left(T - \theta\right)^{2}} + \frac{2 P \left(T a_{1} - a_{1} \theta + a_{3}\right)}{\left(T - \theta\right)^{3}} + \frac{2 P_{r} a_{1}}{\left(T - \theta\right)^{2}} - \frac{2 P_{r} \left(T a_{1} - a_{1} \theta + a_{3}\right)}{\left(T - \theta\right)^{3}} - \frac{2 a_{2} \log{\left(P + \Psi \right)}}{\left(T - \theta\right)^{2}} + \frac{2 a_{2} \log{\left(P_{r} + \Psi \right)}}{\left(T - \theta\right)^{2}} - \left(B{\left(T,P \right)} + 1\right) \frac{\partial^{2}}{\partial T^{2}} \omega{\left(T,P \right)} + \left(B{\left(T,P_{r} \right)} + 1\right) \frac{\partial^{2}}{\partial T^{2}} \omega{\left(T,P_{r} \right)} - \omega{\left(T,P \right)} X{\left(T,P \right)} + \omega{\left(T,P_{r} \right)} X{\left(T,P_{r} \right)} - 2 Y{\left(T,P \right)} \frac{\partial}{\partial T} \omega{\left(T,P \right)} + 2 Y{\left(T,P_{r} \right)} \frac{\partial}{\partial T} \omega{\left(T,P_{r} \right)} + \frac{2 \left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P + \Psi \right)}}{\left(T - \theta\right)^{3}} - \frac{2 \left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P_{r} + \Psi \right)}}{\left(T - \theta\right)^{3}}\right)



Heat Capacity functions
-----------------------

Heat capacity is parameterized as (Tanger and Helgeson, 1988, eq A-4, p
78):

:math:`{{\bar c}_P} = {c_1} + \frac{{{c_2}}}{{{{\left( {T - \theta } \right)}^2}}} - \left[ {\frac{{2T}}{{{{\left( {T - \theta } \right)}^3}}}} \right]\left[ {{a_3}\left( {P - {P_r}} \right) + {a_4}\ln \left( {\frac{{\Psi + P}}{{\Psi + {P_r}}}} \right)} \right] + \omega TX + 2TY\frac{{\partial \omega }}{{\partial T}} - T\left( {\frac{1}{\varepsilon } - 1} \right)\frac{{{\partial ^2}\omega }}{{\partial {T^2}}}`

or equivalently, as the Born function is defined as
:math:`B = - \frac{1}{\varepsilon }`

:math:`{{\bar c}_P} = {c_1} + \frac{{{c_2}}}{{{{\left( {T - \theta } \right)}^2}}} - \left[ {\frac{{2T}}{{{{\left( {T - \theta } \right)}^3}}}} \right]\left[ {{a_3}\left( {P - {P_r}} \right) + {a_4}\ln \left( {\frac{{\Psi + P}}{{\Psi + {P_r}}}} \right)} \right] + \omega TX + 2TY\frac{{\partial \omega }}{{\partial T}} + T\left( {B + 1} \right)\frac{{{\partial ^2}\omega }}{{\partial {T^2}}}`

at the reference pressure this expression becomes

:math:`{\bar c_{{P_r}}} = {c_1} + \frac{{{c_2}}}{{{{\left( {T - \theta } \right)}^2}}} + {\omega _{{P_r}}}T{X_{T,{P_r}}} + 2T{Y_{T,{P_r}}}{\left. {\frac{{\partial \omega }}{{\partial T}}} \right|_{T,{P_r}}} + T\left( {{B_{T,{P_r}}} + 1} \right){\left. {\frac{{{\partial ^2}\omega }}{{\partial {T^2}}}} \right|_{T,{P_r}}}`

Note that the “Born” function terms are the equivalent of
:math:`T{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|_{T,{P_r}}}`,
so that the reference pressure heat capacity can be equivalently
written:

:math:`{{\bar c}_{{P_r}}} = {c_1} + \frac{{{c_2}}}{{{{\left( {T - \theta } \right)}^2}}} + T{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|_{T,{P_r}}}`

as demonstrated here:

.. code:: ipython3

    c1,c2 = sym.symbols('c1 c2')
    ctpr = c1 + c2/(T-theta)**2 + (T*((B(T,P)+1)*omega).diff(T,2)).subs(P,Pr)
    ctpr




.. math::

    \displaystyle T \left(\left(B{\left(T,P_{r} \right)} + 1\right) \frac{\partial^{2}}{\partial T^{2}} \omega{\left(T,P_{r} \right)} + \omega{\left(T,P_{r} \right)} X{\left(T,P_{r} \right)} + 2 Y{\left(T,P_{r} \right)} \frac{\partial}{\partial T} \omega{\left(T,P_{r} \right)}\right) + c_{1} + \frac{c_{2}}{\left(T - \theta\right)^{2}}



Gibbs free energy
-----------------

The above analysis gives a way to write the Gibbs free energy using the
identity :math:`dG = - SdT + VdP`, from which:

:math:`{G_{T,P}} = {G_{{T_r},{P_r}}} - \int_{{T_r}}^T {{S_{T,{P_r}}}dT} + \int_{{P_r}}^P {{V_{T,P}}dP}`

The entropy term is given by

:math:`{S_{T,{P_r}}} = {S_{{T_r},{P_r}}} + \int_{{T_r}}^T {\frac{{{{\bar c}_{T,{P_r}}}}}{T}dT}`

Combining expressions

:math:`{G_{T,P}} = {G_{{T_r},{P_r}}} - {S_{{T_r},{P_r}}}\left( {T - {T_r}} \right) - \int_{{T_r}}^T {\int_{{T_r}}^T {\frac{{{{\bar c}_{T,{P_r}}}}}{T}dT} dT} + \int_{{P_r}}^P {{V_{T,P}}dP}`

Now, using the above expression for the reference pressure heat capacity
we can write

:math:`{G_{T,P}} = {G_{{T_r},{P_r}}} - {S_{{T_r},{P_r}}}\left( {T - {T_r}} \right) - \int_{{T_r}}^T {\int_{{T_r}}^T {\frac{{\left( {{c_1} + \frac{{{c_2}}}{{{{\left( {T - \theta } \right)}^2}}} + T{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}} \right)}}{T}dT} dT} + \int_{{P_r}}^P {{V_{T,P}}dP}`

or

:math:`{G_{T,P}} = {G_{{T_r},{P_r}}} - {S_{{T_r},{P_r}}}\left( {T - {T_r}} \right) - \int_{{T_r}}^T {\int_{{T_r}}^T {\left[ {\frac{{{c_1}}}{T} + \frac{{{c_2}}}{{T{{\left( {T - \theta } \right)}^2}}} + {{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}} \right]dT} dT} + \int_{{P_r}}^P {{V_{T,P}}dP}`

which expands to

:math:`{G_{T,P}} = {G_{{T_r},{P_r}}} - {S_{{T_r},{P_r}}}\left( {T - {T_r}} \right) - \int_{{T_r}}^T {\int_{{T_r}}^T {\left[ {\frac{{{c_1}}}{T} + \frac{{{c_2}}}{{T{{\left( {T - \theta } \right)}^2}}}} \right]dT} dT} - \int_{{T_r}}^T {\int_{{T_r}}^T {{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}dT} dT} + \int_{{P_r}}^P {{V_{T,P}}dP}`

Note that
:math:`\int_{{T_r}}^T {\int_{{T_r}}^T {{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}dT} dT}`
evaluates to:

:math:`\int_{{T_r}}^T {\int_{{T_r}}^T {{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}dT} dT} = \int_{{T_r}}^T {\left[ {{{\left. {\frac{{\partial \left( {B + 1} \right)\omega }}{{\partial T}}} \right|}_{T,{P_r}}} - {{\left. {\frac{{\partial \left( {B + 1} \right)\omega }}{{\partial T}}} \right|}_{{T_r},{P_r}}}} \right]dT}`

and further to:

:math:`\int_{{T_r}}^T {\int_{{T_r}}^T {{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}dT} dT} = \left( {{B_{T,{P_r}}} + 1} \right){\omega _{T,P_r}} - \left( {{B_{{T_r},{P_r}}} + 1} \right){\omega _{{T_r},{P_r}}} - {\left. {\frac{{\partial \left( {B + 1} \right)\omega }}{{\partial T}}} \right|_{{T_r},{P_r}}}\left( {T - {T_r}} \right)`

and still further to:

:math:`\int_{{T_r}}^T {\int_{{T_r}}^T {{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}dT} dT} = \left( {{B_{T,{P_r}}} + 1} \right){\omega _{T,{P_r}}} - \left( {{B_{{T_r},{P_r}}} + 1} \right){\omega _{{T_r},{P_r}}} - {\left. {\frac{{\partial B}}{{\partial T}}} \right|_{{T_r},{P_r}}}{\omega _{{T_r},{P_r}}}\left( {T - {T_r}} \right) - {\left. {\left( {{B_{{T_r},{P_r}}} + 1} \right)\frac{{\partial \omega }}{{\partial T}}} \right|_{{T_r},{P_r}}}\left( {T - {T_r}} \right)`

recognizing that :math:`Y = \frac{{\partial B}}{{\partial T}}`, we have
finally

:math:`\int_{{T_r}}^T {\int_{{T_r}}^T {{{\left. {\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}} \right|}_{T,{P_r}}}dT} dT} = \left( {{B_{T,{P_r}}} + 1} \right){\omega _{T,{P_r}}} - \left( {{B_{{T_r},{P_r}}} + 1} \right){\omega _{{T_r},{P_r}}} - {Y_{{T_r},{P_r}}}\left( {T - {T_r}} \right){\omega _{{T_r},{P_r}}} - {\left. {\left( {{B_{{T_r},{P_r}}} + 1} \right)\left( {T - {T_r}} \right)\frac{{\partial \omega }}{{\partial T}}} \right|_{{T_r},{P_r}}}`

Note that:

$ - :raw-latex:`\int`\ *{{T_r}}^T {:raw-latex:`\int`*\ {{T_r}}^T
{{{:raw-latex:`\left`.
{:raw-latex:`\frac{{{\partial ^2}\left( {B + 1} \right)\omega }}{{\partial {T^2}}}`}
:raw-latex:`\right`\|}\ *{T,{P_r}}}dT} dT} = - :raw-latex:`\left`(
{{B*\ {T,{P_r}}} + 1}
:raw-latex:`\right`){:raw-latex:`\omega `\ *{T,{P_r}}} +
:raw-latex:`\left`( {{B*\ {{T_r},{P_r}}} + 1}
:raw-latex:`\right`){:raw-latex:`\omega `\ *{{T_r},{P_r}}} +
{Y*\ {{T_r},{P_r}}}:raw-latex:`\left`( {T - {T_r}}
:raw-latex:`\right`){:raw-latex:`\omega `\ *{{T_r},{P_r}}} +
{:raw-latex:`\left`. {:raw-latex:`\left`( {{B*\ {{T_r},{P_r}}} + 1}
:raw-latex:`\right`):raw-latex:`\left`( {T - {T_r}}
:raw-latex:`\right`):raw-latex:`\frac{{\partial \omega }}{{\partial T}}`}
:raw-latex:`\right`\|_{{T_r},{P_r}}}$

and

:math:`\int_{{P_r}}^P {{V_{T,P}}dP} = f\left( {T,\theta ,{a_1},{a_2},{a_3},{a_4},T,P,{P_r}} \right) - \left( {{B_{T,P}} + 1} \right){\omega _{T,P}} + \left( {{B_{T,{P_r}}} + 1} \right){\omega _{T,{P_r}}}`

Add to the above parameters the reference state enthalpy,
:math:`H_{ref}` and the reference state entropy, :math:`S_{ref}`

.. code:: ipython3

    Gref,Sref = sym.symbols('G_ref S_ref')

Derive an expression for the Gibbs free energy
----------------------------------------------

.. code:: ipython3

    gtp = Gref - Sref*(T-Tr) - sym.integrate(sym.integrate((c1 + c2/(T-theta)**2)/T, (T,Tr,T)), (T,Tr,T)) \
        - (B(T,Pr)+1)*omega.subs(P,Pr) + (B(Tr,Pr)+1)*omega.subs(P,Pr).subs(T,Tr) \
        + Y(Tr,Pr)*(T-Tr)*omega.subs(P,Pr).subs(T,Tr) + (B(Tr,Pr)+1)*(T-Tr)*omega.diff(T).subs(P,Pr) \
        + sym.integrate(a1 + a2/(Psi+P) + (a3 + a4/(Psi+P))/(T-theta), (P,Pr,P)) - (1+B(T,P))*omega + (1+B(T,Pr))*omega.subs(P,Pr)
    gtp




.. math::

    \displaystyle G_{ref} + \frac{P \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - \frac{P_{r} \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - S_{ref} \left(T - T_{r}\right) + \frac{T c_{2} \log{\left(T + \frac{- c_{1} \theta^{3} - 2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} - \frac{T \left(- T_{r} c_{1} \theta^{2} \log{\left(T_{r} \right)} - T_{r} c_{1} \theta^{2} - T_{r} c_{2} \log{\left(T_{r} \right)} + T_{r} c_{2} \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{1} \theta^{3} \log{\left(T_{r} \right)} + c_{1} \theta^{3} + c_{2} \theta \log{\left(T_{r} \right)} - c_{2} \theta \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{2} \theta\right)}{T_{r} \theta^{2} - \theta^{3}} - \frac{T_{r} c_{2} \log{\left(T_{r} + \frac{- c_{1} \theta^{3} - 2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} + \frac{T_{r} \left(- T_{r} c_{1} \theta^{2} \log{\left(T_{r} \right)} - T_{r} c_{1} \theta^{2} - T_{r} c_{2} \log{\left(T_{r} \right)} + T_{r} c_{2} \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{1} \theta^{3} \log{\left(T_{r} \right)} + c_{1} \theta^{3} + c_{2} \theta \log{\left(T_{r} \right)} - c_{2} \theta \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{2} \theta\right)}{T_{r} \theta^{2} - \theta^{3}} + \left(T - T_{r}\right) \left(B{\left(T_{r},P_{r} \right)} + 1\right) \frac{\partial}{\partial T} \omega{\left(T,P_{r} \right)} + \left(T - T_{r}\right) \omega{\left(T_{r},P_{r} \right)} Y{\left(T_{r},P_{r} \right)} - \left(B{\left(T,P \right)} + 1\right) \omega{\left(T,P \right)} + \left(B{\left(T_{r},P_{r} \right)} + 1\right) \omega{\left(T_{r},P_{r} \right)} + \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P + \Psi \right)}}{T - \theta} - \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P_{r} + \Psi \right)}}{T - \theta} - \frac{\left(T c_{1} \theta^{2} + T c_{2}\right) \log{\left(T + \frac{- c_{1} \theta^{3} - c_{2} \theta + \theta \left(c_{1} \theta^{2} + c_{2}\right)}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} + \frac{\left(T_{r} c_{1} \theta^{2} + T_{r} c_{2}\right) \log{\left(T_{r} + \frac{- c_{1} \theta^{3} - c_{2} \theta + \theta \left(c_{1} \theta^{2} + c_{2}\right)}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}}



At the reference pressure, :math:`\omega` is a constant,
:math:`\omega_0`, and its derivatives are zero too.

.. code:: ipython3

    gtp = gtp.subs(omega.subs(T,Tr).subs(P,Pr), omega0)
    gtp = gtp.subs(omega.diff(T).subs(P,Pr),0)
    gtp




.. math::

    \displaystyle G_{ref} + \frac{P \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - \frac{P_{r} \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - S_{ref} \left(T - T_{r}\right) + \frac{T c_{2} \log{\left(T + \frac{- c_{1} \theta^{3} - 2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} - \frac{T \left(- T_{r} c_{1} \theta^{2} \log{\left(T_{r} \right)} - T_{r} c_{1} \theta^{2} - T_{r} c_{2} \log{\left(T_{r} \right)} + T_{r} c_{2} \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{1} \theta^{3} \log{\left(T_{r} \right)} + c_{1} \theta^{3} + c_{2} \theta \log{\left(T_{r} \right)} - c_{2} \theta \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{2} \theta\right)}{T_{r} \theta^{2} - \theta^{3}} - \frac{T_{r} c_{2} \log{\left(T_{r} + \frac{- c_{1} \theta^{3} - 2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} + \frac{T_{r} \left(- T_{r} c_{1} \theta^{2} \log{\left(T_{r} \right)} - T_{r} c_{1} \theta^{2} - T_{r} c_{2} \log{\left(T_{r} \right)} + T_{r} c_{2} \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{1} \theta^{3} \log{\left(T_{r} \right)} + c_{1} \theta^{3} + c_{2} \theta \log{\left(T_{r} \right)} - c_{2} \theta \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{2} \theta\right)}{T_{r} \theta^{2} - \theta^{3}} + \omega_{0} \left(T - T_{r}\right) Y{\left(T_{r},P_{r} \right)} + \omega_{0} \left(B{\left(T_{r},P_{r} \right)} + 1\right) - \left(B{\left(T,P \right)} + 1\right) \omega{\left(T,P \right)} + \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P + \Psi \right)}}{T - \theta} - \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P_{r} + \Psi \right)}}{T - \theta} - \frac{\left(T c_{1} \theta^{2} + T c_{2}\right) \log{\left(T + \frac{- c_{1} \theta^{3} - c_{2} \theta + \theta \left(c_{1} \theta^{2} + c_{2}\right)}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} + \frac{\left(T_{r} c_{1} \theta^{2} + T_{r} c_{2}\right) \log{\left(T_{r} + \frac{- c_{1} \theta^{3} - c_{2} \theta + \theta \left(c_{1} \theta^{2} + c_{2}\right)}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}}



Substitute the definition of the omega function into the expression for
the Gibbs free energy

.. code:: ipython3

    gtp = gtp.subs(omega,omega_def)
    gtp




.. math::

    \displaystyle G_{ref} + \frac{P \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - \frac{P_{r} \left(T a_{1} - a_{1} \theta + a_{3}\right)}{T - \theta} - S_{ref} \left(T - T_{r}\right) + \frac{T c_{2} \log{\left(T + \frac{- c_{1} \theta^{3} - 2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} - \frac{T \left(- T_{r} c_{1} \theta^{2} \log{\left(T_{r} \right)} - T_{r} c_{1} \theta^{2} - T_{r} c_{2} \log{\left(T_{r} \right)} + T_{r} c_{2} \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{1} \theta^{3} \log{\left(T_{r} \right)} + c_{1} \theta^{3} + c_{2} \theta \log{\left(T_{r} \right)} - c_{2} \theta \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{2} \theta\right)}{T_{r} \theta^{2} - \theta^{3}} - \frac{T_{r} c_{2} \log{\left(T_{r} + \frac{- c_{1} \theta^{3} - 2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} + \frac{T_{r} \left(- T_{r} c_{1} \theta^{2} \log{\left(T_{r} \right)} - T_{r} c_{1} \theta^{2} - T_{r} c_{2} \log{\left(T_{r} \right)} + T_{r} c_{2} \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{1} \theta^{3} \log{\left(T_{r} \right)} + c_{1} \theta^{3} + c_{2} \theta \log{\left(T_{r} \right)} - c_{2} \theta \log{\left(T_{r} - \frac{c_{1} \theta^{3}}{c_{1} \theta^{2} + 2 c_{2}} - \frac{2 c_{2} \theta}{c_{1} \theta^{2} + 2 c_{2}} \right)} + c_{2} \theta\right)}{T_{r} \theta^{2} - \theta^{3}} + \omega_{0} \left(T - T_{r}\right) Y{\left(T_{r},P_{r} \right)} + \omega_{0} \left(B{\left(T_{r},P_{r} \right)} + 1\right) - \left(B{\left(T,P \right)} + 1\right) \left(\begin{cases} \omega_{0} & \text{for}\: z = 0 \\- \frac{\eta z}{rH + \operatorname{gSolvent}{\left(T,P \right)}} + \frac{\frac{\eta z}{rH} + \omega_{0}}{1 + \frac{\left(\frac{z}{rH} + \frac{\omega_{0}}{\eta}\right) \left|{z}\right| \operatorname{gSolvent}{\left(T,P \right)}}{z^{2}}} & \text{otherwise} \end{cases}\right) + \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P + \Psi \right)}}{T - \theta} - \frac{\left(T a_{2} - a_{2} \theta + a_{4}\right) \log{\left(P_{r} + \Psi \right)}}{T - \theta} - \frac{\left(T c_{1} \theta^{2} + T c_{2}\right) \log{\left(T + \frac{- c_{1} \theta^{3} - c_{2} \theta + \theta \left(c_{1} \theta^{2} + c_{2}\right)}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}} + \frac{\left(T_{r} c_{1} \theta^{2} + T_{r} c_{2}\right) \log{\left(T_{r} + \frac{- c_{1} \theta^{3} - c_{2} \theta + \theta \left(c_{1} \theta^{2} + c_{2}\right)}{c_{1} \theta^{2} + 2 c_{2}} \right)}}{\theta^{2}}



The “charge” contribution to the Gibbs free energy is

.. code:: ipython3

    gtp-gtp.subs(z,0)




.. math::

    \displaystyle \omega_{0} \left(B{\left(T,P \right)} + 1\right) - \left(B{\left(T,P \right)} + 1\right) \left(\begin{cases} \omega_{0} & \text{for}\: z = 0 \\- \frac{\eta z}{rH + \operatorname{gSolvent}{\left(T,P \right)}} + \frac{\frac{\eta z}{rH} + \omega_{0}}{1 + \frac{\left(\frac{z}{rH} + \frac{\omega_{0}}{\eta}\right) \left|{z}\right| \operatorname{gSolvent}{\left(T,P \right)}}{z^{2}}} & \text{otherwise} \end{cases}\right)



Add the derived expression for *G(T,P)* and its parameter list to the
model

.. code:: ipython3

    params = [('G_ref','J',Gref), ('S_ref','J/K',Sref), 
              ('a1','J/bar-m',a1), ('a2','J/bar^2-m',a2), ('a3','J/bar-m',a3),  ('a4','J/bar^2-m',a4), 
              ('c1','J/K-m',c1), ('c2','J-K/m',c2), 
              ('Psi', 'bar', Psi), ('eta', 'Å-J/mole', eta), ('rH', 'Å', rH), ('omega0','J',omega0), 
              ('theta','K',theta), ('z', '', z)]
    model.add_expression_to_model(gtp, params)

Code Print the Model, compile the code and link a Python module
---------------------------------------------------------------

Name the model class

.. code:: ipython3

    model.set_module_name('hkf')
    model.set_include_born_code(True)

Make a working sub-directory and move down into the directory. This is
done so that generated files will not clash between alternate model
configurations.

.. code:: ipython3

    model_working_dir = "aqueous"
    !mkdir -p {model_working_dir}
    %cd {model_working_dir}


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/aqueous


Define a parameter dictionary for a specific aqueous species. Parameters
from DEW.

.. code:: ipython3

    use_charged_species_example = True
    
    if use_charged_species_example:
        species     = 'Na+'
        formula     = 'Na(1)'
        DEW_g       = -62591.0     #       cal/mol
        DEW_h       = -57433.0     #       cal/mol
        DEW_s       = 13.96        #       cal/mol/K
        DEW_v       = -1.106808692 #       cm^3/mol
        DEW_cp      = 9.076552599  #       cal/mol/K
        DEW_omega   = 0.3306       # 10^-5 cal/mol
        DEW_z       = 1.0          #       charge
        DEW_a1      = 1.839        # 10    cal/mol/bar
        DEW_a2      = -2.285       # 10^-2 cal/mol
        DEW_a3      = 3.256        #       cal-K/mol/bar
        DEW_a4      = -2.726       # 10^-4 cal-K/mol
        DEW_c1      = 18.18        #       cal/mol/K
        DEW_c2      = -2.981       # 10^-4 cal-K/mol
        DEW_comment = "Shock & Helgeson (1988)"
    else:
        species     = 'NaCl'
        formula     = "Na(1)Cl(1)"
        DEW_g       = -92910.0     #       cal/mol
        DEW_h       = -96120.0     #       cal/mol
        DEW_s       = 28.00        #       cal/mol/K
        DEW_v       = 24.00006579  #       cm^3/mol
        DEW_cp      = 8.508360325  #       cal/mol/K
        DEW_omega   = -0.038       # 10^-5 cal/mol
        DEW_z       = 0.0          #       charge
        DEW_a1      = 5.0363       # 10    cal/mol/bar
        DEW_a2      = 4.7365       # 10^-2 cal/mol
        DEW_a3      = 3.4154       #       cal-K/mol/bar
        DEW_a4      = -2.9748      # 10^-4 cal-K/mol
        DEW_c1      = 10.8         #       cal/mol/K
        DEW_c2      = -1.3         # 10^-4 cal-K/mol
        DEW_comment = "Shock & Helgeson (1988)"

Scaling constants ..

.. code:: ipython3

    SCALEforOmega = 1.0e+5
    SCALEforA1    = 1.0e-1
    SCALEforA2    = 1.0e+2
    SCALEforA4    = 1.0e+4
    SCALEforC2    = 1.0e+4
    CALtoJOULES   = 4.184

.. code:: ipython3

    param_dict = {'G_ref':DEW_g*CALtoJOULES, 
                  'S_ref':DEW_s*CALtoJOULES, 
                  'a1':DEW_a1*CALtoJOULES*SCALEforA1, 
                  'a2':DEW_a2*CALtoJOULES*SCALEforA2, 
                  'a3':DEW_a3*CALtoJOULES,  
                  'a4':DEW_a4*CALtoJOULES*SCALEforA4, 
                  'c1':DEW_c1*CALtoJOULES, 
                  'c2':DEW_c2*CALtoJOULES*SCALEforC2, 
                  'omega0':DEW_omega*CALtoJOULES*SCALEforOmega, 
                  'theta':228.0,
                  'Psi':2600.0,
                  'eta':1.66027e5*CALtoJOULES,
                  'rH':3.082,
                  'z':DEW_z,
                  'T_r':298.15,
                  'P_r':1.0
                 }
    phase_name = "Na_1_cation"
    param_dict




.. parsed-literal::

    {'G_ref': -261880.744,
     'S_ref': 58.408640000000005,
     'a1': 0.7694376,
     'a2': -956.0440000000002,
     'a3': 13.623104,
     'a4': -114055.84000000001,
     'c1': 76.06512000000001,
     'c2': -124725.04000000001,
     'omega0': 138323.04,
     'theta': 228.0,
     'Psi': 2600.0,
     'eta': 694656.968,
     'rH': 3.082,
     'z': 1.0,
     'T_r': 298.15,
     'P_r': 1.0}



Note that the call to

::

   model.create_code_module(phase=phase_name, formula=formula, params=param_dict)

generates fast code with unmodifiable model parameters and
“calibration-” related functions. The call to:

::

   model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type='calib')

generates code suitable for model parameter calibration.

.. code:: ipython3

    #result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict)
    result = model.create_code_module(phase=phase_name, formula=formula, params=param_dict, module_type='calib')


.. parsed-literal::

    Creating (once only) generic fast model code file string
    Creating (once only) generic model calib code template include file string
    Creating (once only) generic model calib code template code file string
    Creating (once only) generic calib model code file string
    Creating include file ...
    ... done!
    Creating code file ...
    ... done
    Writing include file to working directory ...
    Writing code file to working directory ...
    Writing pyxbld file to working directory ...
    writing pyx file to working directory ...
    Compiling code and Python bindings ...
    Success! Import the module named  hkf


Import the new module and test the model
----------------------------------------

.. code:: ipython3

    import hkf
    %cd ..


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/Cython/Compiler/Main.py:369: FutureWarning: Cython directive 'language_level' not set, using 2 for now (Py2). This will change in a later release! File: /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen/aqueous/hkf.pyx
      tree = Parsing.p_module(s, pxd, full_module_name)


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Codegen


Evaluate functions at temperature (K) and pressure (bars)

.. code:: ipython3

    t = 1298.15
    p = 1001.0

Available in both “Fast” and “Calib” code versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execute the “fast” or “calibration” code metadata retrieval functions:

.. code:: ipython3

    try:
        print(hkf.cy_Na_1_cation_hkf_identifier())
        print(hkf.cy_Na_1_cation_hkf_name())
        print(hkf.cy_Na_1_cation_hkf_formula())
        print(hkf.cy_Na_1_cation_hkf_mw())
        print(hkf.cy_Na_1_cation_hkf_elements())
    except AttributeError:
        pass
    try:
        print(hkf.cy_Na_1_cation_hkf_calib_identifier())
        print(hkf.cy_Na_1_cation_hkf_calib_name())
        print(hkf.cy_Na_1_cation_hkf_calib_formula())
        print(hkf.cy_Na_1_cation_hkf_calib_mw())
        print(hkf.cy_Na_1_cation_hkf_calib_elements())
    except AttributeError:
        pass


.. parsed-literal::

    Wed Sep 23 09:56:24 2020
    Na_1_cation
    Na
    22.98977
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
     0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]


Execute the standard thermodynamic property retrieval functions:

.. code:: ipython3

    fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
    try:
        print(fmt.format('G', hkf.cy_Na_1_cation_hkf_g(t,p), 'J/m'))
        print(fmt.format('dGdT', hkf.cy_Na_1_cation_hkf_dgdt(t,p), 'J/K-m'))
        print(fmt.format('dGdP', hkf.cy_Na_1_cation_hkf_dgdp(t,p), 'J/bar-m'))
        print(fmt.format('d2GdP2', hkf.cy_Na_1_cation_hkf_d2gdt2(t,p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', hkf.cy_Na_1_cation_hkf_d2gdtdp(t,p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', hkf.cy_Na_1_cation_hkf_d2gdp2(t,p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', hkf.cy_Na_1_cation_hkf_d3gdt3(t,p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', hkf.cy_Na_1_cation_hkf_d3gdt2dp(t,p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', hkf.cy_Na_1_cation_hkf_d3gdtdp2(t,p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', hkf.cy_Na_1_cation_hkf_d3gdp3(t,p), 'J/bar^3-m'))
        print(fmt.format('S', hkf.cy_Na_1_cation_hkf_s(t,p), 'J/K-m'))
        print(fmt.format('V', hkf.cy_Na_1_cation_hkf_v(t,p), 'J/bar-m'))
        print(fmt.format('Cv', hkf.cy_Na_1_cation_hkf_cv(t,p), 'J/K-m'))
        print(fmt.format('Cp', hkf.cy_Na_1_cation_hkf_cp(t,p), 'J/K-m'))
        print(fmt.format('dCpdT', hkf.cy_Na_1_cation_hkf_dcpdt(t,p), 'J/K^2-m'))
        print(fmt.format('alpha', hkf.cy_Na_1_cation_hkf_alpha(t,p), '1/K'))
        print(fmt.format('beta', hkf.cy_Na_1_cation_hkf_beta(t,p), '1/bar'))
        print(fmt.format('K', hkf.cy_Na_1_cation_hkf_K(t,p), 'bar'))
        print(fmt.format('Kp', hkf.cy_Na_1_cation_hkf_Kp(t,p), ''))
    except AttributeError:
        pass
    try:
        print(fmt.format('G', hkf.cy_Na_1_cation_hkf_calib_g(t,p), 'J/m'))
        print(fmt.format('dGdT', hkf.cy_Na_1_cation_hkf_calib_dgdt(t,p), 'J/K-m'))
        print(fmt.format('dGdP', hkf.cy_Na_1_cation_hkf_calib_dgdp(t,p), 'J/bar-m'))
        print(fmt.format('d2GdP2', hkf.cy_Na_1_cation_hkf_calib_d2gdt2(t,p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', hkf.cy_Na_1_cation_hkf_calib_d2gdtdp(t,p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', hkf.cy_Na_1_cation_hkf_calib_d2gdp2(t,p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', hkf.cy_Na_1_cation_hkf_calib_d3gdt3(t,p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', hkf.cy_Na_1_cation_hkf_calib_d3gdt2dp(t,p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', hkf.cy_Na_1_cation_hkf_calib_d3gdtdp2(t,p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', hkf.cy_Na_1_cation_hkf_calib_d3gdp3(t,p), 'J/bar^3-m'))
        print(fmt.format('S', hkf.cy_Na_1_cation_hkf_calib_s(t,p), 'J/K-m'))
        print(fmt.format('V', hkf.cy_Na_1_cation_hkf_calib_v(t,p), 'J/bar-m'))
        print(fmt.format('Cv', hkf.cy_Na_1_cation_hkf_calib_cv(t,p), 'J/K-m'))
        print(fmt.format('Cp', hkf.cy_Na_1_cation_hkf_calib_cp(t,p), 'J/K-m'))
        print(fmt.format('dCpdT', hkf.cy_Na_1_cation_hkf_calib_dcpdt(t,p), 'J/K^2-m'))
        print(fmt.format('alpha', hkf.cy_Na_1_cation_hkf_calib_alpha(t,p), '1/K'))
        print(fmt.format('beta', hkf.cy_Na_1_cation_hkf_calib_beta(t,p), '1/bar'))
        print(fmt.format('K', hkf.cy_Na_1_cation_hkf_calib_K(t,p), 'bar'))
        print(fmt.format('Kp', hkf.cy_Na_1_cation_hkf_calib_Kp(t,p), ''))
    except AttributeError:
        pass


.. parsed-literal::

    G          -2.787254e+05 J/m       
    dGdT        1.769563e+01 J/K-m     
    dGdP       -1.553101e+02 J/bar-m   
    d2GdP2     -3.089866e-01 J/K^2-m   
    d2GdTdP    -1.692411e-01 J/K-bar-m 
    d2GdP2      3.750827e-01 J/bar^2-m 
    d3GdT3     -8.825633e-04 J/K^3-m   
    d3GdT2dP    5.486632e-04 J/K^2-bar-
    d3GdTdP2    2.847840e-04 J/K-bar^2-
    d3GdP3     -1.257095e-03 J/bar^3-m 
    S          -1.769563e+01 J/K-m     
    V          -1.553101e+02 J/bar-m   
    Cv          5.002420e+02 J/K-m     
    Cp          4.011110e+02 J/K-m     
    dCpdT       1.454686e+00 J/K^2-m   
    alpha       1.089698e-03 1/K       
    beta        2.415056e-03 1/bar     
    K           4.140690e+02 bar       
    Kp          3.877578e-01           


Available only in the “Calib” versions of generated code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Execute the parameter value/metadata functions.
| These functions are only defined for the “calibration” model code
  implementation:

.. code:: ipython3

    try:
        np = hkf.cy_Na_1_cation_hkf_get_param_number()
        names = hkf.cy_Na_1_cation_hkf_get_param_names()
        units = hkf.cy_Na_1_cation_hkf_get_param_units()
        values = hkf.cy_Na_1_cation_hkf_get_param_values()
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], hkf.cy_Na_1_cation_hkf_get_param_value(i), units[i]))
    except AttributeError:
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+00  1.000000e+00 bar       
    G_ref      -2.618807e+05 -2.618807e+05 J         
    S_ref       5.840864e+01  5.840864e+01 J/K       
    a1          7.694376e-01  7.694376e-01 J/bar-m   
    a2         -9.560440e+02 -9.560440e+02 J/bar^2-m 
    a3          1.362310e+01  1.362310e+01 J/bar-m   
    a4         -1.140558e+05 -1.140558e+05 J/bar^2-m 
    c1          7.606512e+01  7.606512e+01 J/K-m     
    c2         -1.247250e+05 -1.247250e+05 J-K/m     
    Psi         2.600000e+03  2.600000e+03 bar       
    eta         6.946570e+05  6.946570e+05 Å-J/mole  
    rH          3.082000e+00  3.082000e+00 Å         
    omega0      1.383230e+05  1.383230e+05 J         
    theta       2.280000e+02  2.280000e+02 K         
    z           1.000000e+00  1.000000e+00           


Test the functions that allow modification of the array of parameter
values

.. code:: ipython3

    try:
        values[1] = 100.0
        hkf.cy_Na_1_cation_hkf_set_param_values(values)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], hkf.cy_Na_1_cation_hkf_get_param_value(i), units[i]))
    except (AttributeError, NameError):
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+02  1.000000e+02 bar       
    G_ref      -2.618807e+05 -2.618807e+05 J         
    S_ref       5.840864e+01  5.840864e+01 J/K       
    a1          7.694376e-01  7.694376e-01 J/bar-m   
    a2         -9.560440e+02 -9.560440e+02 J/bar^2-m 
    a3          1.362310e+01  1.362310e+01 J/bar-m   
    a4         -1.140558e+05 -1.140558e+05 J/bar^2-m 
    c1          7.606512e+01  7.606512e+01 J/K-m     
    c2         -1.247250e+05 -1.247250e+05 J-K/m     
    Psi         2.600000e+03  2.600000e+03 bar       
    eta         6.946570e+05  6.946570e+05 Å-J/mole  
    rH          3.082000e+00  3.082000e+00 Å         
    omega0      1.383230e+05  1.383230e+05 J         
    theta       2.280000e+02  2.280000e+02 K         
    z           1.000000e+00  1.000000e+00           


Test the functions that allow modification of a particular parameter
value

.. code:: ipython3

    try:
        hkf.cy_Na_1_cation_hkf_set_param_value(1, 1.0)
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:<10.10s}"
        for i in range(0,np):
            print(fmt.format(names[i], values[i], hkf.cy_Na_1_cation_hkf_get_param_value(i), units[i]))
    except AttributeError:
        pass


.. parsed-literal::

    T_r         2.981500e+02  2.981500e+02 K         
    P_r         1.000000e+02  1.000000e+00 bar       
    G_ref      -2.618807e+05 -2.618807e+05 J         
    S_ref       5.840864e+01  5.840864e+01 J/K       
    a1          7.694376e-01  7.694376e-01 J/bar-m   
    a2         -9.560440e+02 -9.560440e+02 J/bar^2-m 
    a3          1.362310e+01  1.362310e+01 J/bar-m   
    a4         -1.140558e+05 -1.140558e+05 J/bar^2-m 
    c1          7.606512e+01  7.606512e+01 J/K-m     
    c2         -1.247250e+05 -1.247250e+05 J-K/m     
    Psi         2.600000e+03  2.600000e+03 bar       
    eta         6.946570e+05  6.946570e+05 Å-J/mole  
    rH          3.082000e+00  3.082000e+00 Å         
    omega0      1.383230e+05  1.383230e+05 J         
    theta       2.280000e+02  2.280000e+02 K         
    z           1.000000e+00  1.000000e+00           


Evaluate parameter derivatives …

.. code:: ipython3

    try:
        fmt = "    {0:<10.10s} {1:13.6e}"
        for i in range(0, np):
            print ('Derivative with respect to parameter: ', names[i], ' of')
            print (fmt.format('G', hkf.cy_Na_1_cation_hkf_dparam_g(t, p, i)))
            print (fmt.format('dGdT', hkf.cy_Na_1_cation_hkf_dparam_dgdt(t, p, i)))
            print (fmt.format('dGdP', hkf.cy_Na_1_cation_hkf_dparam_dgdp(t, p, i)))
            print (fmt.format('d2GdT2', hkf.cy_Na_1_cation_hkf_dparam_d2gdt2(t, p, i)))
            print (fmt.format('d2GdTdP', hkf.cy_Na_1_cation_hkf_dparam_d2gdtdp(t, p, i)))
            print (fmt.format('d2GdP2', hkf.cy_Na_1_cation_hkf_dparam_d2gdp2(t, p, i)))
            print (fmt.format('d3GdT3', hkf.cy_Na_1_cation_hkf_dparam_d3gdt3(t, p, i)))
            print (fmt.format('d3GdT2dP', hkf.cy_Na_1_cation_hkf_dparam_d3gdt2dp(t, p, i)))
            print (fmt.format('d3GdTdP2', hkf.cy_Na_1_cation_hkf_dparam_d3gdtdp2(t, p, i)))
            print (fmt.format('d3GdP3', hkf.cy_Na_1_cation_hkf_dparam_d3gdp3(t, p, i)))
    except (AttributeError, TypeError):
        pass


.. parsed-literal::

    Derivative with respect to parameter:  T_r  of
        G           3.519348e+02
        dGdT        2.935261e-01
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  P_r  of
        G           4.883484e-01
        dGdT        7.545207e-04
        dGdP        0.000000e+00
        d2GdT2      4.932880e-08
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3     -1.382857e-10
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  G_ref  of
        G           1.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  S_ref  of
        G          -1.000000e+03
        dGdT       -1.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  a1  of
        G           1.000000e+03
        dGdT        0.000000e+00
        dGdP        1.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  a2  of
        G           3.253156e-01
        dGdT        0.000000e+00
        dGdP        2.777006e-04
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2     -7.711764e-08
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      4.283124e-11
    Derivative with respect to parameter:  a3  of
        G           9.344484e-01
        dGdT       -8.731939e-04
        dGdP        9.344484e-04
        d2GdT2      1.631909e-06
        d2GdTdP    -8.731939e-07
        d2GdP2      0.000000e+00
        d3GdT3     -4.574805e-09
        d3GdT2dP    1.631909e-09
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  a4  of
        G           3.039907e-04
        dGdT       -2.840636e-07
        dGdP        2.594969e-07
        d2GdT2      5.308856e-10
        d2GdTdP    -2.424865e-10
        d2GdP2     -7.206246e-11
        d3GdT3     -1.488256e-12
        d3GdT2dP    4.531823e-13
        d3GdTdP2    6.733866e-14
        d3GdP3      4.002358e-14
    Derivative with respect to parameter:  c1  of
        G          -9.097068e+02
        dGdT       -1.471099e+00
        dGdP        0.000000e+00
        d2GdT2     -7.703270e-04
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      5.934037e-07
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  c2  of
        G          -3.121215e-02
        dGdT       -3.430487e-05
        dGdP        0.000000e+00
        d2GdT2     -6.726448e-10
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      1.775260e-12
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  Psi  of
        G           1.134530e-01
        dGdT       -1.063322e-05
        dGdP        8.194701e-05
        d2GdT2      1.987239e-08
        d2GdTdP    -7.680367e-09
        d2GdP2     -4.551347e-08
        d3GdT3     -5.570917e-11
        d3GdT2dP    1.435381e-11
        d3GdTdP2    4.265686e-12
        d3GdP3      3.791736e-11
    Derivative with respect to parameter:  eta  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  rH  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00
    Derivative with respect to parameter:  omega0  of
        G           7.694419e-01
        dGdT        1.328082e-03
        dGdP       -1.126329e-03
        d2GdT2     -1.810524e-06
        d2GdTdP    -1.223635e-06
        d2GdP2      2.711051e-06
        d3GdT3     -6.705944e-09
        d3GdT2dP    3.966748e-09
        d3GdTdP2    2.058888e-09
        d3GdP3     -9.087778e-09
    Derivative with respect to parameter:  theta  of
        G           6.970779e+01
        dGdT        7.319222e-02
        dGdP       -1.576139e-05
        d2GdT2      4.937114e-08
        d2GdTdP     2.945641e-08
        d2GdP2      7.680367e-09
        d3GdT3     -1.588061e-10
        d3GdT2dP   -8.257650e-11
        d3GdTdP2   -1.435381e-11
        d3GdP3     -4.265686e-12
    Derivative with respect to parameter:  z  of
        G           0.000000e+00
        dGdT        0.000000e+00
        dGdP        0.000000e+00
        d2GdT2      0.000000e+00
        d2GdTdP     0.000000e+00
        d2GdP2      0.000000e+00
        d3GdT3      0.000000e+00
        d3GdT2dP    0.000000e+00
        d3GdTdP2    0.000000e+00
        d3GdP3      0.000000e+00


Time execution of the code
--------------------------

.. code:: ipython3

    try:
        %timeit hkf.cy_Na_1_cation_hkf_calib_g(t,p) 
    except AttributeError:
        pass
    try:
        %timeit hkf.cy_Na_1_cation_hkf_g(t,p) 
    except AttributeError:
        pass


.. parsed-literal::

    166 µs ± 3.57 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)


Compare to Objective-C implementation of DEW model
--------------------------------------------------

.. code:: ipython3

    from ctypes import cdll
    from ctypes import util
    from rubicon.objc import ObjCClass, objc_method
    cdll.LoadLibrary(util.find_library('phaseobjc'))




.. parsed-literal::

    <CDLL '/usr/local/lib/libphaseobjc.dylib', handle 7fbfef904a70 at 0x7fc042695e50>



.. code:: ipython3

    DEWFluid = ObjCClass('DEWFluid')
    obj = DEWFluid.alloc().init()
    if use_charged_species_example:
        refModel = obj.endmembers.objectAtIndex_(25) # 25 is Na+
    else:
        refModel = obj.endmembers.objectAtIndex_(24) # 24 is NaCl
    refModel.phaseName




.. parsed-literal::

    'Na+'



.. code:: ipython3

    t = 1298.15
    p = 5001.0

.. code:: ipython3

    import math
    fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:13.6e} {4:6.2f}%"
    fmts = "{0:<10.10s} {1:13.6e}"
    try:
        x = hkf.cy_Na_1_cation_hkf_g(t,p)
        y = refModel.getGibbsFreeEnergyFromT_andP_(t,p)
        print(fmt.format('G', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_dgdt(t,p)
        y = -refModel.getEntropyFromT_andP_(t,p)
        print(fmt.format('dGdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_dgdp(t,p)
        y = refModel.getVolumeFromT_andP_(t,p)
        print(fmt.format('dGdP', x, y, x-y, 100.0*math.fabs((x-y)/y))) 
        x = hkf.cy_Na_1_cation_hkf_d2gdt2(t,p)
        print(fmts.format('d2GdT2', x))
        x = hkf.cy_Na_1_cation_hkf_d2gdtdp(t,p)
        y = refModel.getDvDtFromT_andP_(t,p)
        print(fmt.format('d2GdTdP', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_d2gdp2(t,p)
        y = refModel.getDvDpFromT_andP_(t,p)
        print(fmt.format('d2GdP2', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_d3gdt3(t,p)
        print(fmts.format('d3GdT3', x))
        x = hkf.cy_Na_1_cation_hkf_d3gdt2dp(t,p)
        y = refModel.getD2vDt2FromT_andP_(t,p)
        print(fmt.format('d3GdT2dP', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_d3gdtdp2(t,p)
        y = refModel.getD2vDtDpFromT_andP_(t,p)
        print(fmt.format('d3GdTdP2', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_d3gdp3(t,p)
        y = refModel.getD2vDp2FromT_andP_(t,p)
        print(fmt.format('d3GdP3', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_s(t,p)
        y = refModel.getEntropyFromT_andP_(t,p)
        print(fmt.format('S', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_v(t,p)
        y = refModel.getVolumeFromT_andP_(t,p)
        print(fmt.format('V', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_cv(t,p)
        print(fmts.format('Cv', x))
        x = hkf.cy_Na_1_cation_hkf_cp(t,p)
        y = refModel.getHeatCapacityFromT_andP_(t,p)
        print(fmt.format('Cp', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_dcpdt(t,p)
        y = refModel.getDcpDtFromT_andP_(t,p)
        print(fmt.format('dCpdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_alpha(t,p)
        print(fmts.format('alpha', x))
        x = hkf.cy_Na_1_cation_hkf_beta(t,p)
        print(fmts.format('beta', x))
        x = hkf.cy_Na_1_cation_hkf_K(t,p)
        print(fmts.format('K', x))
        x = hkf.cy_Na_1_cation_hkf_Kp(t,p)
        print(fmts.format('Kp', x))
    except AttributeError:
        pass
    try:
        x = hkf.cy_Na_1_cation_hkf_calib_g(t,p)
        y = refModel.getGibbsFreeEnergyFromT_andP_(t,p)
        print(fmt.format('G', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_dgdt(t,p)
        y = -refModel.getEntropyFromT_andP_(t,p)
        print(fmt.format('dGdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_dgdp(t,p)
        y = refModel.getVolumeFromT_andP_(t,p)
        print(fmt.format('dGdP', x, y, x-y, 100.0*math.fabs((x-y)/y))) 
        x = hkf.cy_Na_1_cation_hkf_calib_d2gdt2(t,p)
        print(fmts.format('d2GdT2', x))
        x = hkf.cy_Na_1_cation_hkf_calib_d2gdtdp(t,p)
        y = refModel.getDvDtFromT_andP_(t,p)
        print(fmt.format('d2GdTdP', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_d2gdp2(t,p)
        y = refModel.getDvDpFromT_andP_(t,p)
        print(fmt.format('d2GdP2', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_d3gdt3(t,p)
        print(fmts.format('d3GdT3', x))
        x = hkf.cy_Na_1_cation_hkf_calib_d3gdt2dp(t,p)
        y = refModel.getD2vDt2FromT_andP_(t,p)
        print(fmt.format('d3GdT2dP', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_d3gdtdp2(t,p)
        y = refModel.getD2vDtDpFromT_andP_(t,p)
        print(fmt.format('d3GdTdP2', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_d3gdp3(t,p)
        y = refModel.getD2vDp2FromT_andP_(t,p)
        print(fmt.format('d3GdP3', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_s(t,p)
        y = refModel.getEntropyFromT_andP_(t,p)
        print(fmt.format('S', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_v(t,p)
        y = refModel.getVolumeFromT_andP_(t,p)
        print(fmt.format('V', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_cv(t,p)
        print(fmts.format('Cv', x))
        x = hkf.cy_Na_1_cation_hkf_calib_cp(t,p)
        y = refModel.getHeatCapacityFromT_andP_(t,p)
        print(fmt.format('Cp', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_dcpdt(t,p)
        y = refModel.getDcpDtFromT_andP_(t,p)
        print(fmt.format('dCpdT', x, y, x-y, 100.0*math.fabs((x-y)/y)))
        x = hkf.cy_Na_1_cation_hkf_calib_alpha(t,p)
        print(fmts.format('alpha', x))
        x = hkf.cy_Na_1_cation_hkf_calib_beta(t,p)
        print(fmts.format('beta', x))
        x = hkf.cy_Na_1_cation_hkf_calib_K(t,p)
        print(fmts.format('K', x))
        x = hkf.cy_Na_1_cation_hkf_calib_Kp(t,p)
        print(fmts.format('Kp', x))
    except AttributeError:
        pass


.. parsed-literal::

    G          -3.732904e+05 -3.690727e+05 -4.217695e+03   1.14%
    dGdT       -1.354628e+02 -1.311818e+02 -4.280958e+00   3.26%
    dGdP       -2.288937e+00 -2.288937e+00  0.000000e+00   0.00%
    d2GdT2     -2.985789e-02
    d2GdTdP    -7.511966e-03 -7.511966e-03 -1.734723e-18   0.00%
    d2GdP2      1.143502e-03  1.143502e-03  2.168404e-19   0.00%
    d3GdT3     -8.917339e-05
    d3GdT2dP   -3.366689e-06 -3.476772e-06  1.100835e-07   3.17%
    d3GdTdP2    3.553811e-06  3.561965e-06 -8.154136e-09   0.23%
    d3GdP3     -6.291056e-07 -6.232179e-07 -5.887680e-09   0.94%
    S           1.354628e+02  1.311818e+02  4.280996e+00   3.26%
    V          -2.288937e+00 -2.288943e+00  5.626684e-06   0.00%
    Cv          1.028212e+02
    Cp          3.876002e+01  3.876000e+01  2.184210e-05   0.00%
    dCpdT       1.456183e-01  1.450226e-01  5.957707e-04   0.41%
    alpha       3.281857e-03
    beta        4.995779e-04
    K           2.001690e+03
    Kp          1.012431e-01


