
Calculating Equilibrium Phase Assemblages in Systems Subject to Generalized Constraints on System Chemical Potentials
=====================================================================================================================

Mark Ghiorso^1^

**1** OFM Research, (Seattle, WA, United States)

# 1. Introduction {#MPSection:261389CE-F82D-4E5F-FE73-6C36D5950E5B}

# 2. Generalized thermodynamic potentials {#MPSection:E7CF61A7-68DC-4778-883A-D29F41BD8D0E}

The differential form of the Gibbs free energy (*G*) of a system may be written (Prigogene and Defay, 1954):

<div id="MPEquationElement:6E22730B-B717-4B8C-B6B8-4D057A4646E0">

dG = - SdT + VdP + \\sum\\limits\_i\^c {{\\mu \_i}d{n\_i}}

</div>

where *T*, *P,* and *n~i~* are independent variables denoting the temperature, pressure and number of moles of the *i^th^* thermodynamic component in a system of *c*-components. The entropy of the system is given by *S* ( = - \\frac{{\\partial G}}{{\\partial T}} ), the volume by *V* ( = - \\frac{{\\partial G}}{{\\partial P}} ) and the chemical potential of the *i^th^* component by \\mu~i~ ( = - \\frac{{\\partial G}}{{\\partial {n\_i}}} ). The integral form of <span class="kind elementIndex">**Equation 1** </span> is:

<div id="MPEquationElement:DE19C888-9CE1-4F06-8E66-A814195E784E">

G = \\sum\\limits\_i\^c {{\\mu \_i}{n\_i}}

</div>

The equilibrium state of a system is characterized by a minimum in the Gibbs free energy at specified temperature, pressure and bulk composition. Under these conditions the system evolves from a disequilibrium to a equilibrium state under the restriction that the system boundary is closed to mass transfer (constant ;*n~i~*), but is open to heat exchange with its surroundings and to any volume change of the system associated with phase formation. Khorzhinskii (ref) generalized the concept of thermodynamic equilibrium to “open” systems subject to mass exchange with their surroundings for the specific case that the constancy of component mole number may be exchanged with constancy of the equivalent chemical potential, e.g. equilibrium in a system open to exchange of a fugitive component like H~2~O may be uniquely determined if the chemical potential of H~2~\[Evaluating reaction stoichiometry in magmatic systems evolving under generalized thermodynamic constraints: Examples comparing isothermal and …\](\#Ghiorso1987)O is somehow externally imposed upon the system. The thermodynamic potential that is minimal under these “open” system conditions is the “Korzhinskii potential” ;, which may be defined by a suitable Legendre transformation of <span class="kind elementIndex">**Equation 2** </span> , e.g., for the example involving H~2~O:

<div id="MPEquationElement:FDBB0D6C-309E-44BB-949D-00FB4253B16B">

L = G - \\frac{{\\partial G}}{{\\partial {n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}}}{n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}} = G - {\\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}{n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}

</div>

While the usual application of the Khorzhinskii potential is to the characterization of equilibrium in systems constrained by externally fixed chemical potentials, the theory equally applies to the case where a chemical potential, or some linear combination of chemical potentials, is constrained by the demand that a particular phase is present in the system. For example, if the system is forced to be in equilibrium with a pure H~2~O fluid phase at some temperature and pressure, then <span class="kind elementIndex">**Equation 3** </span> also defines the thermodynamic state function that is minimal at equilibrium for these conditions. Note, that the amount of fluid is not known *a priori* before minimizing *L* subject to fixed *T*, *P*, {n\_{i \\ne {{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}, and {\\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}; {n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}is a *derived* quantity.

In general, for a phase of fixed stoichiometry, indexed as \\phi*~j~*, whose composition can be expressed using the system components *vis-a-vis* a reaction of the form:

<div id="MPEquationElement:3F3248E1-CD38-4AE7-99AA-243895E7C03D">

\\sum\\limits\_k\^c {{c\_{{\\phi \_j}k}}} componen{t\_k} \\Leftrightarrow phas{e\_{{\\phi \_j}}}

</div>

the chemical potential of the phase is:

<div id="MPEquationElement:46636B9E-E958-419D-9138-D31FDDFD6F5D">

{\\mu \_{{\\phi \_j}}} = \\sum\\limits\_k\^c {{c\_{{\\phi \_j}k}}{\\mu \_k}}

</div>

The Khorzhinskii potential that guarantees the presence of phase \\phi*~j~* in the assemblage is:

<div id="MPEquationElement:CCCD231B-39AA-4B0B-AFB6-24C093A64B6C">

L = G - {n\_{{\\phi \_j}}}{\\mu \_{{\\phi \_j}}}

</div>

where {n\_{{\\phi \_j}}} is the mole number of phase {\\phi \_j} in the system. This quantity will not be known *a priori* and must be determined in the course of finding the minimum of *L*. <span class="kind elementIndex">**Equation 6** </span> may be expanded to:

<div id="MPEquationElement:48B7CE38-8B04-4077-C680-D7F4BEE71930">

L = G - {n\_{{\\phi \_j}}}\\left( {\\sum\\limits\_k\^c {{c\_{{\\phi \_j}k}}{\\mu \_k}} } \\right)

</div>

For a collection of *f*-phases, which are assumed to be present in the system, <span class="kind elementIndex">**Equation 6** </span> is generalized to:

<div id="MPEquationElement:DFC5BE25-953A-4A29-8765-3A8FA4F642AC">

L = G - \\sum\\limits\_j\^f {{n\_{{\\phi \_j}}}{\\mu \_{{\\phi \_j}}}}

</div>

and <span class="kind elementIndex">**Equation 7** </span> becomes

<div id="MPEquationElement:F2DFC594-CBD9-4D42-9840-0BD135FF3C53">

L = G - \\sum\\limits\_j\^f {{n\_{{\\phi \_j}}}\\sum\\limits\_k\^c {{c\_{{\\phi \_j}k}}{\\mu \_k}} }

</div>

Equation <span class="kind elementIndex">**Equation 9** </span> is the generalized Khorzhinskii potential for a thermodynamic system containing a specified collection of *f*-phases.

<div id="MPFootnotesElement:6B2F0337-9F2E-4B3F-B49F-EAF21924A6E4">

<div id="MPFootnote:BA043F1B-A7C9-4758-9088-1D06B627C1B4-contents">

A solution phase consisting of two or more stoichiometric end member components may be treated by writing a statement of <span class="kind elementIndex">**Equation 4** </span> for each end member; the amount of each end member (i.e. <span class="kind elementIndex">**Equation 6** </span> ) thereafter determines the composition of the solution phase. See further discussion later in the text.

</div>

</div>

# 3. Constraints on the minimum {#MPSection:736862DD-BD17-4212-AA84-E13BBFAA3A77}

In order to determine the equilibrium phase assemblage governed by the potential of ;<span class="kind elementIndex">**Equation 9** </span>, temperature, pressure and the mass of all system components *that are not allowed to vary by phase present constraints* must be specified. In addition, constraints on chemical potentials, i.e., <span class="kind elementIndex">**Equation 5** </span>\[Chemical mass transfer in magmatic processes. i. thermodynamic relations and numerical algorithms\](\#Ghiorso1985), must be imposed in order to guarantee that the stipulated phase assemblage is present at equilibrium. Without the phase present constraints, the specification of a constraint matrix that guarantees constancy of system bulk composition is straightforward, and results in a linear system of equality constraints on the mole numbers of phases in the equilibrium assemblage ;. For this phase-present scenario, these constraints become more involved.

Consider a system of *p*-phases at some specified temperature and pressure, with *f* of those *p* phases specified or forced to be present in the system. The number of “extra” or non-specified phases is e \\equiv p - f \\ge 0. For the moment, we assume that the identity of the *e* phases in the system is known1. At fixed temperature and pressure, the Gibbs phase rule guarantees that c \\ge p = f + e.

There always exists a linear mapping between the mole numbers of phases present in the system and the mole of the *k^th^* system component:

<div id="MPEquationElement:33A3894F-E6C2-4FFA-9880-B419E84AF395">

{n\_k} = \\sum\\limits\_j\^p {{n\_{{\\phi \_j}}}c\_{k{\\phi \_j}}\^\*}

</div>

This mapping can be extended to all components and succinctly written in matrix form:

<div id="MPEquationElement:7F5060CF-6092-41F3-DDF5-B1AB004616F9">

{{\\bf{n}}\_c} = \\left\[ {\\begin{array}{\*{20}{c}} {{n\_1}}\\\\ \\vdots \\\\ {{n\_k}}\\\\ \\vdots \\\\ {{n\_c}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} {c\_{1{\\phi \_1}}\^\*}& \\cdots &{c\_{1{\\phi \_j}}\^\*}& \\cdots &{c\_{1{\\phi \_p}}\^\*}\\\\ \\vdots & \\ddots & \\vdots & {\\mathinner{\\mkern2mu\\raise1pt\\hbox{.}\\mkern2mu \\raise4pt\\hbox{.}\\mkern2mu\\raise7pt\\hbox{.}\\mkern1mu}} & \\vdots \\\\ {c\_{k{\\phi \_1}}\^\*}& \\cdots &{c\_{k{\\phi \_j}}\^\*}& \\cdots &{c\_{k{\\phi \_p}}\^\*}\\\\ \\vdots & {\\mathinner{\\mkern2mu\\raise1pt\\hbox{.}\\mkern2mu \\raise4pt\\hbox{.}\\mkern2mu\\raise7pt\\hbox{.}\\mkern1mu}} & \\vdots & \\ddots & \\vdots \\\\ {c\_{c{\\phi \_1}}\^\*}& \\cdots &{c\_{c{\\phi \_j}}\^\*}& \\cdots &{c\_{c{\\phi \_p}}\^\*} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{n\_{{\\phi \_1}}}}\\\\ \\vdots \\\\ {{n\_{{\\phi \_k}}}}\\\\ \\vdots \\\\ {{n\_{{\\phi \_p}}}} \\end{array}} \\right\] = {\\bf{C}}{{\\bf{n}}\_{\\phi}}

</div>

where the column dimension of **C** is greater than or equal to its row dimension. The matrix **C** may be partitioned as \\left\[ {{{\\bf{C}}\_f}|{{\\bf{C}}\_e}} \\right\] , and <span class="kind elementIndex">**Equation 11** </span> can be re-written

<div id="MPEquationElement:5310787F-DB88-4416-DB57-A6919AB1D781">

{\\bf{n}\_c} = {{\\bf{C}}\_f}{{\\bf{n}}\_f} + {{\\bf{C}}\_e}{{\\bf{n}}\_e}

</div>

Note that {{\\bf{n}}\_\\phi }, the vector of mole numbers of all phases in the system, is partitioned as {\\left\[ {\\begin{array}{\*{20}{c}} {{{\\bf{n}}\_f}}&{{{\\bf{n}}\_e}} \\end{array}} \\right\]\^T} . The matrix **C**~*f*~ may have rank less than *f*. Any number of situations can generate this rank deficiency, including imposing a collection of fixed phases whose compositions are redundant, such as multiple structural forms of the same stoichiometry. The *transpose* of the matrix **C**~*f*~ may be rendered into the product of three matrices using Singular Value Decomposition (Lawson and Hanson, 1974), as:

<div id="MPEquationElement:3504A45D-6B6E-4D72-ACE1-B77CB660F8F7">

{\\bf{C}}\_f\^T = {{\\bf{U}}\_f}{{\\bf{S}}\_f}{\\bf{V}}\_f\^T

</div>

where **U**~*f*~ and **V**~*f*~ are square orthogonal matrices and **S**~*f*~ is a diagonal-dominant rectangular matrix; the number of non-zero diagonal elements of **S***~f~* is equal to the rank, *r*, of **C**~*f*~. The columns of **V**~*f*~, which are ordered to pair with the diagonal elements of **S**~*f*~, partition into an *r* column range space and f - r column null space:

<div id="MPEquationElement:321FDC7D-3B55-49DD-F61A-84AD01924C77">

\\left\[ {\\begin{array}{\*{20}{c}} {{v\_{11}}}& \\cdots &{{v\_{1r}}}&{{v\_{1r + 1}}}& \\cdots &{{v\_{1f}}}\\\\ \\vdots & \\cdots & \\vdots & \\vdots & \\cdots & \\vdots \\\\ {{v\_{k1}}}& \\cdots &{{v\_{kr}}}&{{v\_{kr + 1}}}& \\cdots &{{v\_{kf}}}\\\\ \\vdots & \\cdots & \\vdots & \\vdots & \\cdots & \\vdots \\\\ {{v\_{c1}}}& \\cdots &{{v\_{cr}}}&{{v\_{cr + 1}}}& \\cdots &{{v\_{cf}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} {{{\\bf{v}}\_1}}& \\cdots &{{{\\bf{v}}\_r}}&{{{\\bf{v}}\_{r + 1}}}& \\cdots &{{{\\bf{v}}\_f}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} {{{\\bf{V}}\_{\\left. f \\right|r}}}&{{{\\bf{V}}\_{\\left. f \\right|f}}} \\end{array}} \\right\]

</div>

The range space defines the set of linear combinations of system component mole numbers that are free to assume any value in order to accommodate the forced equilibration of the system with the *f* fixed phases. The null space defines those linear combinations of system component mole numbers are fixed — i.e., that are independent of any mass transfer associated with bringing the system into equilibrium with the imposed assemblage. If **b**~*c*~ denotes a vector of mole numbers of system components, then if **b***~c~* is substituted for the left hand side of <span class="kind elementIndex">**Equation 12** </span>,

<div id="MPEquationElement:EA1DAEF2-CBEE-48B7-82C7-A2F2F10E86BC">

{\\bf{b}\_c} = {{\\bf{C}}\_f}{{\\bf{n}}\_f} + {{\\bf{C}}\_e}{{\\bf{n}}\_e}

</div>

and both sides of the expression are multiplied by {\\bf{V}}\_{\\left. f \\right|f}\^T, the result is a set of constraints that permit partial mass exchange according to the stoichiometry of the fixed phase constraints, yet maintain initial system composition in the null space of those phase constraints:

<div id="MPEquationElement:9492D404-9C9A-4350-C66C-487C1B1942BF">

{\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{b}}\_{\\bf{c}}} = {\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{C}}\_f}{{\\bf{n}}\_f} + {\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{C}}\_e}{{\\bf{n}}\_e}

</div>

At this point a brief concrete example is welcome and illustrative.

<div id="MPFootnotesElement:D2D01F1B-C50C-4A62-AE69-B36CDFEE84D6">

<div id="MPFootnote:2490AB11-C356-40DF-DB8A-CFEA82B678B2-contents">

<span class="citation">Ghiorso 2013)</span> describes algorithms for estimating the phase identity and approximate phase compositions for systems where the final phase assemblage requires estimation.

</div>

<div id="MPFootnote:1131E832-416D-4A64-FC96-74CA761AE615-contents">

If a phase is a solution, there is an entry in the vector for each end member component of that solution.

</div>

</div>

## 3.1. Example of phase constraints {#MPSection:B3C87549-4E0E-449C-A324-EC09BD655C23}

Consider a chemical system with thermodynamic components SiO~2~, Al~2~O~3~, CaO, Na~2~O and K~2~O. For this system we are interested in phase equilibria at some specified temperature and pressure assuming that both quartz (SiO~2~) and corundum (Al~2~O~3~) are present in the system. In addition, the system is known to contain a silicate liquid. We do not know the composition of the silicate liquid, but we have a thermodynamic model that describes its Gibbs free energy of this solution in terms of components SiO~2~, Al~2~O~3~, CaSiO~3~, Na~2~SiO~3~, KAlSiO~4~. The system is free to exchange both SiO~2~ and Al~2~O~3~ with its surroundings to alter the initially specified bulk composition, given by

<div id="MPEquationElement:59783F1B-FF14-4CA7-CFC5-7073E0E88D81">

{{\\bf{b}}\_c} = {\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_2}}}}&{{b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}}&{{b\_{{\\rm{CaO}}}}}&{{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}&{{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}} \\end{array}} \\right\]\^T}

</div>

in order to bring the system to equilibrium with quartz and corundum at the specified temperature and pressure.

<span class="kind elementIndex">**Equation 15** </span> may be written for this case as:

<div id="MPEquationElement:CEC405BD-17B1-4FBF-87C8-BE118E683854">

\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_2}}}}\\\\ {{b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}}\\\\ {{b\_{{\\rm{CaO}}}}}\\\\ {{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 1&0\\\\ 0&1\\\\ 0&0\\\\ 0&0\\\\ 0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{n\_{Qz}}}\\\\ {{n\_{Cr}}} \\end{array}} \\right\] + \\left\[ {\\begin{array}{\*{20}{c}} 1&0&1&1&1\\\\ 0&1&0&0&{\\frac{1}{2}}\\\\ 0&0&1&0&0\\\\ 0&0&0&1&0\\\\ 0&0&0&0&{\\frac{1}{2}} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}} \\end{array}} \\right\]

</div>

and <span class="kind elementIndex">**Equation 13** </span> becomes:

<div id="MPEquationElement:E2D33DA1-0982-4A5F-DC28-4C631C515E0C">

{\\bf{C}}\_f\^T = \\left\[ {\\begin{array}{\*{20}{c}} 1&0&0&0&0\\\\ 0&1&0&0&0 \\end{array}} \\right\] = {{\\bf{U}}\_f}{{\\bf{S}}\_f}{\\bf{V}}\_f\^T = \\left\[ {\\begin{array}{\*{20}{c}} 0&1\\\\ 1&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 1&0&0&0&0\\\\ 0&1&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 0&1&0&0&0\\\\ 1&0&0&0&0\\\\ 0&0&0&0&1\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]

</div>

from which

<div id="MPEquationElement:48C2AF36-F98D-46E0-8F6B-35A5D203B200">

V\_{\\left. f \\right|f}\^T = \\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&1\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]

</div>

Multiplying <span class="kind elementIndex">**Equation 17** </span> by <span class="kind elementIndex">**Equation 19** </span> gives

<div id="MPEquationElement:79A70C20-7037-4509-8A1A-85074ED3EE43">

\\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&1\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_2}}}}\\\\ {{b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}}\\\\ {{b\_{{\\rm{CaO}}}}}\\\\ {{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&1\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 1&0\\\\ 0&1\\\\ 0&0\\\\ 0&0\\\\ 0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{n\_{Qz}}}\\\\ {{n\_{Cr}}} \\end{array}} \\right\] + \\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&1\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 1&0&1&1&1\\\\ 0&1&0&0&{\\frac{1}{2}}\\\\ 0&0&1&0&0\\\\ 0&0&0&1&0\\\\ 0&0&0&0&{\\frac{1}{2}} \\end{array}} \\right\]

</div>

which simplified to

<div id="MPEquationElement:BDABB612-8384-45B2-EC16-0B998BAC47F5">

\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{\\rm{CaO}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 0&0\\\\ 0&0\\\\ 0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{n\_{Qz}}}\\\\ {{n\_{Cr}}} \\end{array}} \\right\] + \\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&{\\frac{1}{2}}\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{}} \\end{array}} \\right\]

</div>

<span class="kind elementIndex">**Equation 22** </span> is a concrete example of <span class="kind elementIndex">**Equation 16** </span>. Note that the final set of constraints do not involve the initial concentrations of SiO~2~ or Al~2~O~3~ in the system, which is consistent with the open system nature of the phase present constraints.

To complete this example, note the the Khorzhinskii potential, <span class="kind elementIndex">**Equation 9** </span>, would be written:

<div id="MPEquationElement:4E0DEE58-BE21-4138-F6BE-DBF5E68FC512">

L\\left( {T,P,n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq},n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq},{n\_{Qz}},{n\_{Cr}}} \\right) = G - {n\_{Qz}}{\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}} - {n\_{Cr}}{\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}

</div>

and that to compute thermodynamic equilibrium, this function would be minimized subject to constant temperature, pressure, and

<div id="MPEquationElement:0449FDED-5E78-465D-DD8C-A31F89E5D742">

\\begin{array}{c} \\left\[ {\\begin{array}{\*{20}{c}} 0&0\\\\ 0&0\\\\ 0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{n\_{Qz}}}\\\\ {{n\_{Cr}}} \\end{array}} \\right\] + \\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&{\\frac{1}{2}}\\\\ 0&0&0&1&0\\\\ 0&0&1&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{}} \\end{array}} \\right\] - \\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{\\rm{CaO}}}}} \\end{array}} \\right\] = 0\\\\ \\mu \_{Qz}\^o - {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}} = 0\\\\ \\mu \_{\_{Cr}}\^o - {\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}} = 0\\\\ n\_{Qz} = n\_{Fix}\\\\ n\_{Cr} = n\_{Fix} \\end{array}

</div>

The last two constraints must be specified because the system is now open to mass transfer of SiO~2~ and Al~2~O~3~, and although the chemical potentials of quartz and corundum fix the amounts of both components in the liquid phase, the amount of quartz and corundum present in the system is arbitrary. If a constraint is not imposed on these quantities, the Khorzhinskii potential would attempt to seek a minimal value by adding an infinite amount of quartz and corundum to the system. The additional constraints prevent this from happening. Finally, note that because the liquid is an omnicomponent phase, the system chemical potentials can be written in terms of liquid components:

<div id="MPEquationElement:45C36947-1A6B-47E2-D97B-ADE11D29DE38">

\\begin{array}{c} {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}} = \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\\\ {\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}} = \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}\\\\ {\\mu \_{{\\rm{CaO}}}} = \\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\\\ {\\mu \_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} = \\mu \_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\\\ {\\mu \_{{{\\rm{K}}\_2}{\\rm{O}}}} = 2\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq} - 2\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} \\end{array}

</div>

In the next section, we consider numerical methods of performing this minimization.

# 4. Minimizing the generalized potential {#MPSection:6045E840-FBB7-42D4-CF25-DB09CCE21CCD}

<span class="kind elementIndex">**Equation 9** </span>is a non-linear function of *T*, *P*, and {{\\bf{n}}\_\\phi } . At constant temperature and pressure, an equilibrium system may be found by minimizing this equation with respect to all {n\_{{\\phi \_j}}} subject to all relevant statements of <span class="kind elementIndex">**Equation 5** </span> and to <span class="kind elementIndex">**Equation 16** </span>. This problem is one of non-linear optimization subject to non-linear constraints. There are numerous methods for computing numerical solutions to this problem. One of the computationally fastest and proven methods (Ghiorso & Kelemen, 1987) is based on a second-order Newton iterative technique that employs explicit second-order derivatives.

We transform <span class="kind elementIndex">**Equation 9** </span> to its Lagrangian form

<div id="MPEquationElement:CC4A876E-C613-4F91-D9DA-BC1B1C7EFC94">

\\Lambda = G - \\sum\\limits\_j\^f {{n\_{{\\phi \_j}}}\\sum\\limits\_k\^c {{c\_{{\\phi \_j}k}}{\\mu \_k}} } - {\\left\[ {\\begin{array}{\*{20}{c}} {{\\lambda \_{{n\_1}}}}\\\\ \\vdots \\\\ {{\\lambda \_{{n\_{c - f}}}}} \\end{array}} \\right\]\^T}\\left( {{\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{C}}\_f}{{\\bf{n}}\_f} + {\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{C}}\_e}{{\\bf{n}}\_e} - {\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{b}}\_c}} \\right) - \\sum\\limits\_j\^f {{\\lambda \_{{\\phi \_j}}}\\left( {\\sum\\limits\_k\^c {{c\_{{\\phi \_j}k}}{\\mu \_k}} - {\\mu \_{{\\phi \_j}}}} \\right)} - {\\left\[ {\\begin{array}{\*{20}{c}} {{\\lambda \_{{n\_{{\\phi \_1}}}}}}\\\\ \\vdots \\\\ {{\\lambda \_{{n\_{{\\phi \_f}}}}}} \\end{array}} \\right\]\^T}\\left( {{{\\bf{n}}\_f} - {{\\bf{n}}\_{fixed}}} \\right)

</div>

by introducing Lagrange multipliers, {{\\lambda \_{{n\_1}}}} to {{\\lambda \_{{n\_{c - f}}}}} and {{\\lambda \_{{\\phi \_1}}}} to {{\\lambda \_{{\\phi \_f}}}} whose values are never negative and should tend towards zero as the minimum of \\Lambda is achieved. The third, fourth and fifth terms of <span class="kind elementIndex">**Equation 26** </span> derive from the equality constraints. These equality constraints can be rearranged into a zero vector of length *c-r+f*:

<div id="MPEquationElement:BB5977AF-7869-456F-CF04-5074710F8899">

\\left\[ {\\begin{array}{\*{20}{c}} {{\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{C}}\_f}{{\\bf{n}}\_f} + {\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{C}}\_e}{{\\bf{n}}\_e} - {\\bf{V}}\_{\\left. f \\right|f}\^T{{\\bf{b}}\_c}}\\\\ {\\sum\\limits\_k\^c {{c\_{{\\phi \_1}k}}{\\mu \_k}} - {\\mu \_{{\\phi \_1}}}}\\\\ \\vdots \\\\ {\\sum\\limits\_k\^c {{c\_{{\\phi \_f}k}}{\\mu \_k}} - {\\mu \_{{\\phi \_f}}}}\]\]\\\\ {{{\\bf{n}}\_f} - {{\\bf{n}}\_{fixed}}} \\end{array}} \\right\]

</div>

This vector can be differentiated with respect to the mole numbers of phases in the system, {{\\bf{n}}\_\\phi }, to define a “linearized” constraint matrix, **A**,

<div id="MPEquationElement:A5233121-2AD6-44A3-EBCF-1CAD4F0464F8">

{\\bf{A}} = \\left\[ {\\begin{array}{\*{20}{c}} {{{\\bf{C}}\^T}{{\\bf{V}}\_{\\left. f \\right|f}}}\\\\ {\\sum\\limits\_k\^c {{c\_{{\\phi \_1}k}}\\frac{{\\partial {\\mu \_k}}}{{\\partial {{\\bf{n}}\_\\phi }}}} - \\frac{{\\partial {\\mu \_{{\\phi \_1}}}}}{{\\partial {{\\bf{n}}\_\\phi }}}}\\\\ \\vdots \\\\ {\\sum\\limits\_k\^c {{c\_{{\\phi \_f}k}}\\frac{{\\partial {\\mu \_k}}}{{\\partial {{\\bf{n}}\_\\phi }}}} - \\frac{{\\partial {\\mu \_{{\\phi \_f}}}}}{{\\partial {{\\bf{n}}\_\\phi }}}}\\\\ \\begin{array}{\*{20}{c}} {{{\\bf{0}}\_{c \\times c}}}&{{{\\bf{I}}\_{f \\times f}}} \\end{array} \\end{array}} \\right\]

</div>

The dimensions of **A** are *c-r+f* rows and *p* columns. A matrix **Z** can be computed from **A** using either orthogonal or singular value decomposition that has the property of projecting the optimal parameters, {{\\bf{n}}\_\\phi }, into the “linearized” null space of the constraints. The elements of **Z** are themselves functions of {{\\bf{n}}\_\\phi } because the last *f* rows of **A** depend explicitly on {{\\bf{n}}\_\\phi }. This is what makes the constraints non-linear in the optimal parameters. Nevertheless, **A** and **Z** can be applied to embody the non-linear constraints in an approximate fashion in order to develop an algorithm to compute the minimum of \\Lambda.

Choose some initial guess for the minimum of \\Lambda that satisfies the equality constraints and call it {\\bf{n}}\_{\_\\phi }\^i. Define a vector of differences, \\Delta {\\bf{n}}\_{\_\\phi }\^i, so that a better estimate for the minimum is given by {\\bf{n}}\_\\phi \^{i + 1} = \\Delta {\\bf{n}}\_{\_\\phi }\^i + {\\bf{n}}\_\\phi \^i.\\Delta {\\bf{n}}\_{\_\\phi }\^ Applying the null space projection operator obtained from the linearized equality constraints,

<div id="MPEquationElement:6AA33256-6D0B-43EC-8801-1FCB29A7442B">

{\\bf{Zn}}\_\\phi \^{i + 1} = {\\bf{Z}}\\Delta {\\bf{n}}\_{\_\\phi }\^i + {\\bf{Zn}}\_\\phi \^i

</div>

\\Lambdamay now be expanded in a Taylor series, as:

<div id="MPEquationElement:790EE42F-8130-4A01-FB4E-58800DC0E855">

{\\Lambda \_{i + 1}} = {\\Lambda \_i} + {{\\bf{g}}\^T}{\\bf{Z}}\\Delta {\\bf{n}}\_\\phi \^i + {\\left( {\\Delta {\\bf{n}}\_\\phi \^i} \\right)\^T}{{\\bf{Z}}\^T}{\\bf{WZ}}\\Delta {\\bf{n}}\_\\phi \^i + \\cdots

</div>

In <span class="kind elementIndex">**Equation 30** </span> **g** is the gradient of *L* with respect to {\\bf{n}}\_\\phi evaluated at {\\bf{n}}\_\\phi \^i and **W** is a matrix of second derivatives of \\Lambda with respect to {\\bf{n}}\_\\phi again evaluated at {\\bf{n}}\_\\phi \^i. This matrix is known as the Wronskian. Truncating the Taylor expansion after the quadratic term, differentiating the result with respect to {\\bf{n}}\_\\phi \^{i + 1} and setting the resultant expresser to zero, results in:

<div id="MPEquationElement:35C8C7EA-4886-4EE6-EA59-92C2590799C2">

{{\\bf{Z}}\^T}{\\bf{WZ}}\\Delta {\\bf{n}}\_\\phi \^i = - {{\\bf{Z}}\^T}{\\bf{g}}

</div>

Evaluation of **W** require estimates of the Lagrange multipliers. These estimates are provided by solving the system of equations:

<div id="MPEquationElement:076A4500-9B28-41FE-90F0-C30284C08C27">

{\\bf{g}} = {{\\bf{A}}\^T}\\left\[ {\\begin{array}{\*{20}{c}} {{\\lambda \_{{n\_1}}}}\\\\ \\vdots \\\\ {{\\lambda \_{{n\_{c - f}}}}}\\\\ {{\\lambda \_{{\\phi \_1}}}}\\\\ \\vdots \\\\ {{\\lambda \_{{\\phi \_f}}}}\\\\ {{\\lambda \_{{n\_{{\\phi \_1}}}}}}\\\\ \\vdots \\\\ {{\\lambda \_{{n\_{{\\phi \_f}}}}}} \\end{array}} \\right\]

</div>

The numerical solution to <span class="kind elementIndex">**Equation 31** </span>, \\Delta {\\bf{n}}\_{\_\\phi }\^i is the quadratic approximation to the minimum of \\Lambda. \\Delta {\\bf{n}}\_{\_\\phi }\^i + {\\bf{n}}\_\\phi \^i is not expected to yield the optimal solution, and in general \\Delta {\\bf{n}}\_{\_\\phi }\^i not result in an incremental approximation to the optimum that is feasible in light of the constraints. This is because the constraints have been linearized in order to compute \\Delta {\\bf{n}}\_{\_\\phi }\^i, and that linearization is exact only for infinitesimal values of \\Delta {\\bf{n}}\_{\_\\phi }\^i. Furthermore, the quadratic approximation to the minimum may require scaling to better locate a more optimal approximation to the mimimm, and usually this is accomplished by taking some step length, *s*, along the search direction, as:

<div id="MPEquationElement:0A878B2F-6027-4853-EFB2-E7AD7F232648">

{\\bf{\\tilde n}}\_\\phi \^{i + 1} = s\\Delta {\\bf{n}}\_\\phi \^i + {\\bf{n}}\_\\phi \^i + \\Delta {\\bf{n}}\_\\phi \^{corr}

</div>

where \\Delta {\\bf{n}}\_\\phi \^{corr} is a correction term that adjusts the result ( {\\bf{\\tilde n}}\_\\phi \^{i + 1} ) to maintain feasibility of the constraints, i.e. <span class="kind elementIndex">**Equation 27** </span>. The correction term is generally small and must be computed iteratively.

<div id="MPFootnotesElement:6FAE36B6-4E9E-42A4-F21F-6D2355083D5F">

</div>

# 5. Case Studies {#MPSection:128B24A0-DD5B-4E0D-9B46-92BDD2CF4B3C}

In this section we present a few illustrations of solving various phase-constrained equilibrium calculations using the methods developed above.

## 5.1. Pure phase constraints {#MPSection:27B956EC-6DAC-4346-CEA5-7BF8226F9953}

We first consider the simple case explored previously involving two pure phases, quartz and corundum, in equilibrium with a silicate liquid. The Gibbs free energy of this system is given by

<div id="MPEquationElement:427C80A9-4BC1-434E-BF39-ABD703588828">

G = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq} + {n\_{Qz}}\\mu \_{Qz}\^o + {n\_{Cr}}\\mu \_{Cr}\^o

</div>

and from <span class="kind elementIndex">**Equation 23** </span> the Khorzhinskii potential is:

<div id="MPEquationElement:0CA130F1-480A-4CD3-8586-7D0CE3A4539C">

\\begin{array}{c} L = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\\\ + {n\_{Qz}}\\left( {\\mu \_{Qz}\^o - {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}}} \\right) + {n\_{Cr}}\\left( {\\mu \_{Cr}\^o - {\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right) \\end{array}

</div>

where a distinction is made between the imposed chemical potential of SiO~2~ and Al~2~O~3~, i..e {\\mu \_{Qz}\^o} and {\\mu \_{Cr}\^o}, and the system quantities. In other words, the Khorzhinskii potential is defined in such a way that values can be computed for the general disequilibrium case. Of course, at equilibrium, the chemical potential of silica in the system is the same as that of quartz, by design. Similarly for alumina and corundum. Because the liquid is an omnicomponent phase, without loss of generality and for both the equilibrium and disequilibrium case, the identities in <span class="kind elementIndex">**Equation 25** </span> hold, and <span class="kind elementIndex">**Equation 35** </span> becomes:

<div id="MPEquationElement:F297C776-129E-4242-E714-8A2ED0D93368">

\\begin{array}{c} L = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\\\ + {n\_{Qz}}\\left( {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) + {n\_{Cr}}\\left( {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) \\end{array}

</div>

We now have a working for the Khorzhinskii potential.

In order to construct the Lagrangian associated with this potential, we first simply the equality constraint vector derived previously ( <span class="kind elementIndex">**Equation 24** </span> ) by multiplying out the matrices and substituting liquid potentials for those of the system:

<div id="MPEquationElement:4D2AC1F9-11F5-48AF-9779-A103EE710260">

\\left\[ {\\begin{array}{\*{20}{c}} {\\frac{1}{2}n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq} - {b\_{{{\\rm{K}}\_2}{\\rm{O}}}}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} - {b\_{{\\rm{CaO}}}}}\\\\ {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ n\_{Qz} - 1\\\\ n\_{Cr} - 1 \\end{array}} \\right\]

</div>

The matrix **A** is the derivative of this vector with respect to n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq},n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq},{n\_{Qz}},{n\_{Cr}}:

<div id="MPEquationElement:1C348E70-93C1-4344-C528-97A9DB82973C">

{\\bf{A}} = \\left\[ {\\begin{array}{\*{20}{c}} 0&0&0&0&{\\frac{1}{2}}&0&0\\\\ 0&0&0&1&0&0&0\\\\ 0&0&1&0&0&0&0\\\\ { - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&0\\\\ { - \\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&0\\\\ 0&0&0&0&0&1&0\\\\ 0&0&0&0&0&0&1 \\end{array}} \\right\]

</div>

and the gradient, **g**, of <span class="kind elementIndex">**Equation 36** </span> is:

<div id="MPEquationElement:7BDAED01-2278-4F11-D0AA-1CF3EDF6183B">

{\\bf{g}} = \\left\[ {\\begin{array}{\*{20}{c}} {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}} - {n\_{Cr}}\\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}\\\\ {\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}} - {n\_{Cr}}\\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}\\\\ {\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {\\mu \_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\end{array}} \\right\]

</div>

The Lagrange multiplier may now be obtained by solving <span class="kind elementIndex">**Equation 32** </span> using both <span class="kind elementIndex">**Equation 38** </span> and <span class="kind elementIndex">**Equation 39** </span> evaluated at an initial guess to the solution. How is this initial guess obtained? By specifying an initial system bulk composition, **b***~c~*, and assigning that composition entirely to the liquid phase, at which point concentrations of n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} and n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} are adjusted to satisfy the equalities \\mu \_{Qz}\^o = \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} and \\mu \_{Cr}\^o = \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}. Once the multipliers are computed the Lagrangian function,

<div id="MPEquationElement:4E2B7ECA-2540-4E21-EACA-BCC921E50288">

\\begin{array}{c} \\Lambda = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}\\\\ + {n\_{Qz}}\\left( {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) + {n\_{Cr}}\\left( {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - {\\lambda \_{{{\\rm{K}}\_2}{\\rm{O}}}}\\left( {\\frac{1}{2}n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq} - {b\_{{{\\rm{K}}\_2}{\\rm{O}}}}} \\right) - {\\lambda \_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}\\left( {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}} \\right)\\\\ { - {\\lambda \_{{\\rm{CaO}}}}\\left( {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} - {b\_{{\\rm{CaO}}}}} \\right) - {\\lambda \_{{\\rm{Qz}}}}\\left( {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) - {\\lambda \_{{\\rm{Cr}}}}\\left( {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - {\\lambda \_{{n\_{{\\rm{Qz}}}}}}\\left( {{n\_{Qz}} - 1} \\right) - {\\lambda \_{{n\_{{\\rm{Cr}}}}}}\\left( {{n\_{Cr}} - 1} \\right)} \\end{array}

</div>

may be evaluated and it’s Wronskian computed by taking the second partial derivatives with respect to n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}, n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}, n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}, n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}, n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}, {n\_{Qz}}, and {n\_{Cr}}. For this purpose, <span class="kind elementIndex">**Equation 40** </span> is written more compactly as

<div id="MPEquationElement:33EFC217-8764-4FD0-F47D-5874D3561334">

\\begin{array}{c} \\Lambda = {G\^{liq}} - \\lambda \_{1:f}\^T\\left( {{\\bf{V}}\_{\\left. f \\right|f}\^T{\\bf{Cn}} - {\\bf{V}}\_{\\left. f \\right|f}\^T{\\bf{b}}} \\right) + \\left( {{n\_{Qz}} - {\\lambda \_{{\\rm{Qz}}}}} \\right)\\left( {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) + \\left( {{n\_{Cr}} - {\\lambda \_{{\\rm{Cr}}}}} \\right)\\left( {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right)\\\\ - {\\lambda \_{{n\_{{\\rm{Qz}}}}}}\\left( {{n\_{Qz}} - 1} \\right) - {\\lambda \_{{n\_{{\\rm{Cr}}}}}}\\left( {{n\_{Cr}} - 1} \\right) \\end{array}

</div>

and its first derivative with respect to *n~i~* is

<div id="MPEquationElement:5C54AEEA-2EDE-4E23-D3A9-5B0F9EA8321C">

\\begin{array}{c} \\frac{{\\partial \\Lambda }}{{\\partial {n\_i}}} = \\frac{{\\partial {G\^{liq}}}}{{\\partial {n\_i}}} - {\\left( {{{\\bf{C}}\^T}{{\\bf{V}}\_{\\left. f \\right|f}}{\\lambda \_{1:f}}} \\right)\_i} + \\frac{{\\partial {n\_{Qz}}}}{{\\partial {n\_i}}}\\left( {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) - \\left( {{n\_{Qz}} - {\\lambda \_{{\\rm{Qz}}}}} \\right)\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}}} + \\frac{{\\partial {n\_{Cr}}}}{{\\partial {n\_i}}}\\left( {\\mu \_{Cr}\^o - \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - \\left( {{n\_{Cr}} - {\\lambda \_{{\\rm{Cr}}}}} \\right)\\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial {n\_i}}}\\\\ - {\\delta \_{{n\_i},{n\_{Qz}}}}{\\lambda \_{{n\_{Qz}}}} - {\\delta \_{{n\_i},{n\_{Cr}}}}{\\lambda \_{{n\_{Cr}}}} \\end{array}

</div>

which yields a second derivative

<div id="MPEquationElement:9AB7B596-467B-4C00-C417-95BFC9EA2A29">

\\frac{{{\\partial \^2}\\Lambda }}{{\\partial {n\_i}\\partial {n\_j}}} = \\frac{{{\\partial \^2}{G\^{liq}}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{\\partial {n\_{Qz}}}}{{\\partial {n\_i}}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_j}}} - \\frac{{\\partial {n\_{Qz}}}}{{\\partial {n\_j}}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}}} - \\left( {{n\_{Qz}} - {\\lambda \_{{\\rm{Qz}}}}} \\right)\\frac{{{\\partial \^2}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{\\partial {n\_{Cr}}}}{{\\partial {n\_i}}}\\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial {n\_j}}} - \\frac{{\\partial {n\_{Cr}}}}{{\\partial {n\_j}}}\\frac{{\\partial \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial {n\_i}}} - \\left( {{n\_{Cr}} - {\\lambda \_{{\\rm{Cr}}}}} \\right)\\frac{{{\\partial \^2}\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}

</div>

The elements of the Wronskian are given by <span class="kind elementIndex">**Equation 43** </span>.

In order to compute a search direction according to <span class="kind elementIndex">**Equation 31** </span> the orthogonal null space projection matrix, **Z**, must be computed from **A**. Previously we obtained a similar matrix to project the equality constraints by employing Singular Value Decomposition. That numerical technique is computationally costly, and its use is acceptable because the projection only needed to be applied once for the problem at hand. However, **A** must be evaluated and decomposed numerous times in the course of minimizing the Lagrangian, so the fastest possible numerical algorithm should be employed. The preferred method is RQ decomposition, which factors the matrix into the product:

<div id="MPEquationElement:9EB740CE-2A66-4AD1-D5DD-581B818B8B39">

{\\bf{A}} = {\\bf{RQ}} = \\left\[ {\\begin{array}{\*{20}{c}} {{{\\bf{R}}\_1}\_1}&{\\bf{0}} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{{\\bf{Q}}\_1}}\\\\ {{{\\bf{Q}}\_2}} \\end{array}} \\right\]

</div>

where the matrix **Q** is orthonormal ( {{\\bf{Q}}\^T}{\\bf{Q}} = {\\bf{I}} ) and **R** is upper triangular. From <span class="kind elementIndex">**Equation 44** </span> it follows that {{\\bf{Q}}\_2}\^T{\\bf{A}} = {\\bf{0}} , and consequently that {{\\bf{Q}}\_2}\^T is the **Z** matrix that we seek.

## 5.2. Fixed chemical potential of oxygen and water {#MPSection:99975B8E-4502-4E2F-AD38-2A5CC95F68F4}

Next we consider a silicate liquid in the system SiO~2~-Fe~2~O~3~-Fe~2~SiO~4~-H~2~O. The problem is similar to the one in the previous section, except that it is assumed that an O~2~-H~2~O fluid phase is not saturated. The Gibbs free energy of this system is written as:

<div id="MPEquationElement:539341CA-CB1E-4B4E-9493-E4CFBD8FB3E8">

G = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} + n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}

</div>

from which the Khorzhinskii potential is formally given by

<div id="MPEquationElement:E4CED649-DDCB-4F19-B7EF-28E291C67ADC">

L = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} + n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq} - n\_{{{\\rm{O}}\_2}}\^{sys}\\mu \_{{{\\rm{O}}\_2}}\^{sys} - n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}

</div>

Because the system is not saturated with either pure O~2~ gas, pure H~2~O fluid, or a mixed H~2~O-O~2~ fluid, mass balance considerations demand that n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq} = n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} and \\frac{1}{2}n\_{{{\\rm{O}}\_2}}\^{sys} = n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}. The equality constraints corresponding to this scenario can be expressed as:

<div id="MPEquationElement:846F733B-A2DF-458F-AE04-3E49F065CE4B">

\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_{\\rm{2}}}}}}\\\\ {{b\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_{\\rm{3}}}}}}\\\\ {{b\_{{\\rm{FeO}}}}}\\\\ {{b\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 1&0&1&0&0&0\\\\ 0&1&0&0&2&0\\\\ 0&0&2&0&{ - 4}&0\\\\ 0&0&0&1&0&1\\\\ \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}}\\\\ {n\_{{{\\rm{O}}\_2}}\^{sys}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}} \\end{array}} \\right\]

</div>

where the last two rows account for the additional mass balance constraints. The fixed constraints may be decomposed using SVD as in <span class="kind elementIndex">**Equation 19** </span>

<div id="MPEquationElement:63ABD1AD-20F8-42D1-C03C-BE80D9737E0A">

{\\bf{C}}\_f\^T = \\left\[ {\\begin{array}{\*{20}{c}} 0&2&{ - 4}&0\\\\ 0&0&0&1 \\end{array}} \\right\] = {{\\bf{U}}\_f}{{\\bf{S}}\_f}{\\bf{V}}\_f\^T = \\left\[ {\\begin{array}{\*{20}{c}} 1&0\\\\ 0&1 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {2\\sqrt 5 }&0&0&0\\\\ 0&1&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{1}{{\\sqrt 5 }}}&{ - \\frac{2}{{\\sqrt 5 }}}&0\\\\ 0&0&0&1\\\\ 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{1}{{\\sqrt 5 }}}&0\\\\ 1&0&0&0 \\end{array}} \\right\]

</div>

yielding

<div id="MPEquationElement:35BE27B9-E23B-4EC5-EDF2-4201DDA18872">

V\_{\\left. f \\right|f}\^T = \\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{1}{{\\sqrt 5 }}}&0\\\\ 1&0&0&0 \\end{array}} \\right\]

</div>

Multiplication of the constraint equation by this null-space projection operator gives

<div id="MPEquationElement:52A10275-C691-47CD-E265-E001C2263390">

\\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{1}{{\\sqrt 5 }}}&0\\\\ 1&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_{\\rm{2}}}}}}\\\\ {{b\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_{\\rm{3}}}}}}\\\\ {{b\_{{\\rm{FeO}}}}}\\\\ {{b\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{1}{{\\sqrt 5 }}}&0\\\\ 1&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 1&0&1&0&0&0\\\\ 0&1&0&0&2&0\\\\ 0&0&2&0&{ - 4}&0\\\\ 0&0&0&1&0&1 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}}\\\\ {n\_{{{\\rm{O}}\_2}}\^{sys}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}} \\end{array}} \\right\]

</div>

and finally the projected bulk composition constraint matrix that we seek:

<div id="MPEquationElement:AA5D2915-5818-491A-F436-DCE30B61AAF7">

\\left\[ {\\begin{array}{\*{20}{c}} {\\frac{2}{{\\sqrt 5 }}{b\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_{\\rm{3}}}}} + \\frac{1}{{\\sqrt 5 }}{b\_{{\\rm{FeO}}}}}\\\\ {{b\_{{\\rm{Si}}{{\\rm{O}}\_{\\rm{2}}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{2}{{\\sqrt 5 }}}&0&0&0\\\\ 1&0&1&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}}\\\\ {n\_{{{\\rm{O}}\_2}}\^{sys}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}} \\end{array}} \\right\]

</div>

The complete constraint vector is constructed by stacking this projection along with the fixed external constraints, e.g. <span class="kind elementIndex">**Equation 24** </span>

<div id="MPEquationElement:E459E444-1CC9-4B02-AD4B-36C49DEF345C">

\\begin{array}{\*{20}{c}} {\\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{2}{{\\sqrt 5 }}}&0&0&0\\\\ 1&0&1&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}}\\\\ {n\_{{{\\rm{O}}\_2}}\^{sys}}\\\\ {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}} \\end{array}} \\right\] - \\left\[ {\\begin{array}{\*{20}{c}} {\\frac{2}{{\\sqrt 5 }}{b\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_{\\rm{3}}}}} + \\frac{1}{{\\sqrt 5 }}{b\_{{\\rm{FeO}}}}}\\\\ {{b\_{{\\rm{Si}}{{\\rm{O}}\_{\\rm{2}}}}}} \\end{array}} \\right\] = {\\bf{0}}}\\\\ {\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} - \\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq} = 0}\\\\ {\\mu \_{{{\\rm{O}}\_2}}\^{sys} - \\mu \_{{{\\rm{O}}\_2}}\^{liq} = \\mu \_{{{\\rm{O}}\_2}}\^{sys} - 2\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - 2\\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_3}}\^{liq} + 2\\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} = 0}\\\\ n\_{{{\\rm{O}}\_2}}\^{sys} + 2n\_{{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} - 2n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - 2n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} = 0\\\\ n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} - n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq} = 0 \\end{array}

</div>

The last two rows of constraints reflect the fact that there is only one phase in the system, liquid. From <span class="kind elementIndex">**Equation 52** </span> , the **A** matrix may be derived

<div id="MPEquationElement:F96AF3BC-E6EF-41BC-D2B7-7CA443AB2721">

{\\bf{A}} = \\left\[ {\\begin{array}{\*{20}{c}} 0&{\\frac{2}{{\\sqrt 5 }}}&{\\frac{2}{{\\sqrt 5 }}}&0&0&0\\\\ 1&0&1&0&0&0\\\\ { - \\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}{{\\partial n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}{{\\partial n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}{{\\partial n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}}&0&{\\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{sys}}}{{\\partial n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{sys}}}}\\\\ { - \\frac{{\\partial \\mu \_{{{\\rm{O}}\_{\\rm{2}}}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{{\\rm{O}}\_{\\rm{2}}}}\^{liq}}}{{\\partial n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{{\\rm{O}}\_{\\rm{2}}}}\^{liq}}}{{\\partial n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{{\\rm{O}}\_{\\rm{2}}}}\^{liq}}}{{\\partial n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}}&{\\frac{{\\partial \\mu \_{{{\\rm{O}}\_{\\rm{2}}}}\^{sys}}}{{\\partial n\_{{{\\rm{O}}\_2}}\^{sys}}}}&0\\\\ { - 2}&{ - 2}&2&0&1&0\\\\ 0&0&0&{ - 1}&0&1 \\end{array}} \\right\]

</div>

Note that the entries for the **A** matrix in the fourth row may be computed using the definition of the liquid chemical potential of oxygen:

<div id="MPEquationElement:0B2E8025-ED86-43D8-90C4-AFD43412BF16">

\\frac{{\\partial \\mu \_{{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_i\^{liq}}} = 2\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_i\^{liq}}} + 2\\frac{{\\partial \\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{{\\rm{O}}\_3}}\^{liq}}}{{\\partial n\_i\^{liq}}} - 2\\frac{{\\partial \\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}}{{\\partial n\_i\^{liq}}}

</div>

The two derivatives, {\\frac{{\\partial \\mu \_{{{\\rm{O}}\_{\\rm{2}}}}\^{sys}}}{{\\partial n\_{{{\\rm{O}}\_2}}\^{sys}}}} and {\\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{sys}}}{{\\partial n\_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{sys}}}} are zero if the imposed chemical potentials are fixed.

The gradient of the Khorzhinskii potential is computed from <span class="kind elementIndex">**Equation 46** </span> as

<div id="MPEquationElement:3C91031F-9803-430E-A52C-FE86632AE386">

{\\bf{g}} = \\left\[ {\\begin{array}{\*{20}{c}} {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {\\mu \_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {\\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}}\\\\ { - \\mu \_{{{\\rm{O}}\_2}}\^{sys} - n\_{{{\\rm{O}}\_2}}\^{sys}\\frac{{\\partial \\mu \_{{{\\rm{O}}\_2}}\^{sys}}}{{\\partial n\_{{{\\rm{O}}\_2}}\^{sys}}}}\\\\ { - \\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} - n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}\\frac{{\\partial \\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}}}{{\\partial n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}}}} \\end{array}} \\right\]

</div>

and the Lagrangian is given by:

<div id="MPEquationElement:5BB3D1EE-BCF0-4933-D73D-A4AD452EFF8D">

\\begin{array}{c} \\Lambda = n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq}\\mu \_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq}\\mu \_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} + n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq} - n\_{{{\\rm{O}}\_2}}\^{sys}\\mu \_{{{\\rm{O}}\_2}}\^{sys} - n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}\\\\ - {\\lambda \_{{\\rm{F}}{{\\rm{e}}\_{Tot}}}}\\left( {\\frac{2}{{\\sqrt 5 }}n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} + \\frac{2}{{\\sqrt 5 }}n\_{{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} - \\frac{2}{{\\sqrt 5 }}{b\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}} - \\frac{1}{{\\sqrt 5 }}n\_{{\\rm{FeO}}}\^{liq}} \\right) - {\\lambda \_{{\\rm{Si}}{{\\rm{O}}\_2}}}\\left( {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - {b\_{{\\rm{Si}}{{\\rm{O}}\_2}}}} \\right)\\\\ - {\\lambda \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O,}}\\mu }}\\left( {\\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{sys} - \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}} \\right) - {\\lambda \_{{{\\rm{O}}\_2},\\mu }}\\left( {\\mu \_{{{\\rm{O}}\_2}}\^{sys} - \\mu \_{{{\\rm{O}}\_2}}\^{liq}} \\right)\\\\ - {\\lambda \_{{{\\rm{O}}\_2},n}}\\left( {n\_{{{\\rm{O}}\_2}}\^{sys} + 2n\_{{\\rm{F}}{{\\rm{e}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_4}}\^{liq} - 2n\_{{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}\^{liq} - 2n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) - {\\lambda \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O,}}n}}\\left( {n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} - n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{liq}} \\right) \\end{array}

</div>

from which

<div id="MPEquationElement:26A9526D-B5CA-426A-A2C8-4DD456D29296">

\\frac{{\\partial \\Lambda }}{{\\partial {n\_i}}} = \\mu \_i\^{liq} - \\frac{{\\partial n\_{{{\\rm{O}}\_2}}\^{sys}}}{{\\partial {n\_i}}}\\mu \_{{{\\rm{O}}\_2}}\^{sys} - {\\delta \_{i,{{\\rm{H}}\_2}{\\rm{O}}}}\\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} - {\\lambda \_{{\\rm{F}}{{\\rm{e}}\_{Tot}}}}\\left( {\\frac{2}{{\\sqrt 5 }}{\\delta \_{i,{\\rm{F}}{{\\rm{e}}\_2}{{\\rm{O}}\_3}}} + \\frac{2}{{\\sqrt 5 }}{\\delta \_{i,{\\rm{F}}{{\\rm{e}}\_2}{\\rm{Si}}{{\\rm{O}}\_4}}}} \\right) - {\\lambda \_{{\\rm{Si}}{{\\rm{O}}\_2}}}{\\delta \_{i,{\\rm{Si}}{{\\rm{O}}\_2}}} + {\\lambda \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O,}}\\mu }}\\frac{{\\partial \\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}{{\\partial {n\_i}}} + {\\lambda \_{{{\\rm{O}}\_2},\\mu }}\\frac{{\\partial \\mu \_{{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}}}

</div>

and

<div id="MPEquationElement:DA281CF5-CE28-4D23-9EA4-14195FC70830">

\\frac{{\\partial \\Lambda }}{{\\partial n\_{{{\\rm{O}}\_2}}\^{sys}}} = - \\mu \_{{{\\rm{O}}\_2}}\^{sys} - {\\lambda \_{{{\\rm{O}}\_2},n}}

</div>

and

<div id="MPEquationElement:617332EE-4C4D-49C7-A82A-8566EBC225E7">

\\frac{{\\partial \\Lambda }}{{\\partial n\_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys}}} = - \\mu \_{{{\\rm{H}}\_2}{\\rm{O}}}\^{sys} - {\\lambda \_{{{\\rm{H}}\_2}{\\rm{O}},n}}

</div>

<div id="MPEquationElement:572F6FA4-6EDD-4673-C681-9580EFCE816A">

\\frac{{{\\partial \^2}\\Lambda }}{{\\partial {n\_i}\\partial {n\_j}}} = \\frac{{\\partial \\mu \_i\^{liq}}}{{\\partial {n\_j}}} + {\\lambda \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O,}}\\mu }}\\frac{{{\\partial \^2}\\mu \_{{{\\rm{H}}\_{\\rm{2}}}{\\rm{O}}}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}} + {\\lambda \_{{{\\rm{O}}\_2},\\mu }}\\frac{{{\\partial \^2}\\mu \_{{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}

</div>

where all oxygen and water derivatives are zero.

## 5.3. Fixed saturation with a solid solution {#MPSection:6C8BC2EE-DE8C-4B1C-A362-7B9B4A59AE8F}

We next consider the case of one pure phase, quartz, and a solid solution, feldspar, in equilibrium with a silicate liquid. The Gibbs free energy of this system is given by

<div id="MPEquationElement:C1065AAB-ED83-4084-B272-521AA034A0A9">

G = {G\^{liq}}\\left( {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq},n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}} \\right) + {n\_{Qz}}\\mu \_{Qz}\^o + {G\^{plag}}\\left( {n\_{Ab}\^{plag},n\_{An}\^{plag},n\_{Sn}\^{plag}} \\right) + {G\^{alk}}\\left( {n\_{Ab}\^{alk},n\_{An}\^{alk},n\_{Sn}\^{alk}} \\right)

</div>

and the Khorzhinskii potential by

<div id="MPEquationElement:437E2000-839B-4FFB-8EA5-451F96148926">

\\begin{array}{c} L = {G\^{liq}}\\left( {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq},n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}} \\right) + {n\_{Qz}}\\mu \_{Qz}\^o + {G\^{plag}}\\left( {n\_{Ab}\^{plag},n\_{An}\^{plag},n\_{Sn}\^{plag}} \\right) + {G\^{alk}}\\left( {n\_{Ab}\^{alk},n\_{An}\^{alk},n\_{Sn}\^{alk}} \\right)\\\\ { - {n\_{Qz}}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\Omega \_{Ab}\^{liq} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\Omega \_{An}\^{liq} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\Omega \_{Sn}\^{liq}} \\end{array}

</div>

where

<div id="MPEquationElement:097C97FD-5B1F-4BAA-C271-2B18E255EFDC">

\\begin{array}{l} \\Omega \_{Ab}\^{liq} = \\frac{5}{2}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + \\frac{1}{2}\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} + \\frac{1}{2}\\mu \_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}\\\\ \\Omega \_{An}\^{liq} = \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + \\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} + \\frac{1}{2}\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}\\\\ \\Omega \_{Sn}\^{liq} = 2\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + \\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq} \\end{array}

</div>

The equality constraint matrix for this problem is

<div id="MPEquationElement:B41655A5-BBB9-4A1B-A06A-3506C9B300D8">

\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_2}}}}\\\\ {{b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}}\\\\ {{b\_{{\\rm{CaO}}}}}\\\\ {{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 1&0&1&1&1&1&3&2&3&3&2&3\\\\ 0&1&0&0&{\\frac{1}{2}}&0&{\\frac{1}{2}}&1&{\\frac{1}{2}}&{\\frac{1}{2}}&1&{\\frac{1}{2}}\\\\ 0&0&1&0&0&0&0&1&0&0&1&0\\\\ 0&0&0&1&0&0&{\\frac{1}{2}}&0&0&{\\frac{1}{2}}&0&0\\\\ 0&0&0&0&{\\frac{1}{2}}&0&0&0&{\\frac{1}{2}}&0&0&{\\frac{1}{2}} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {{n\_{Qz}}}\\\\ {n\_{Ab}\^{plag}}\\\\ {n\_{An}\^{plag}}\\\\ {n\_{Sn}\^{plag}}\\\\ {n\_{Ab}\^{alk}}\\\\ {n\_{An}\^{alk}}\\\\ {n\_{Sn}\^{alk}} \\end{array}} \\right\]

</div>

and the decomposition of the “fixed” constraint partition of the constraint matrix is given by:

<div id="MPEquationElement:124E1C3D-DC37-4FB1-AAE5-8CA8DA6D903D">

{\\bf{C}}\_f\^T = \\left\[ {\\begin{array}{\*{20}{c}} 1&0&0&0&0\\\\ 3&{\\frac{1}{2}}&0&{\\frac{1}{2}}&0\\\\ 2&1&1&0&0\\\\ 3&{\\frac{1}{2}}&0&0&{\\frac{1}{2}}\\\\ 3&{\\frac{1}{2}}&0&{\\frac{1}{2}}&0\\\\ 2&1&1&0&0\\\\ 3&{\\frac{1}{2}}&0&0&{\\frac{1}{2}} \\end{array}} \\right\] = {{\\bf{U}}\_f}{{\\bf{S}}\_f}{\\bf{V}}\_f\^T = {{\\bf{U}}\_f}\\left\[ {\\begin{array}{\*{20}{c}} x&0&0&0&0\\\\ 0&x&0&0&0\\\\ 0&0&x&0&0\\\\ 0&0&0&x&0\\\\ 0&0&0&0&0\\\\ 0&0&0&0&0\\\\ 0&0&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} x&x&x&x&x\\\\ x&x&x&x&x\\\\ 0&0&0&x&x\\\\ x&x&x&x&x\\\\ 0&{ - \\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}} \\end{array}} \\right\]

</div>

where *x* denotes some non-zero entry. The null space projection operator applied to the contain matrix, <span class="kind elementIndex">**Equation 63** </span>, gives:

<div id="MPEquationElement:E668E6FE-998E-43D3-8FE5-3CF403D519E7">

\\left\[ {\\begin{array}{\*{20}{c}} 0&{ - \\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {{b\_{{\\rm{Si}}{{\\rm{O}}\_2}}}}\\\\ {{b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}}\\\\ {{b\_{{\\rm{CaO}}}}}\\\\ {{b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}}}\\\\ {{b\_{{{\\rm{K}}\_2}{\\rm{O}}}}} \\end{array}} \\right\] = \\left\[ {\\begin{array}{\*{20}{c}} 0&{ - \\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} 1&0&1&1&1&1&3&2&3&3&2&3\\\\ 0&1&0&0&{\\frac{1}{2}}&0&{\\frac{1}{2}}&1&{\\frac{1}{2}}&{\\frac{1}{2}}&1&{\\frac{1}{2}}\\\\ 0&0&1&0&0&0&0&1&0&0&1&0\\\\ 0&0&0&1&0&0&{\\frac{1}{2}}&0&0&{\\frac{1}{2}}&0&0\\\\ 0&0&0&0&{\\frac{1}{2}}&0&0&0&{\\frac{1}{2}}&0&0&{\\frac{1}{2}} \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {{n\_{Qz}}}\\\\ {n\_{Ab}\^{plag}}\\\\ {n\_{An}\^{plag}}\\\\ {n\_{Sn}\^{plag}}\\\\ {n\_{Ab}\^{alk}}\\\\ {n\_{An}\^{alk}}\\\\ {n\_{Sn}\^{alk}} \\end{array}} \\right\]

</div>

which simplifies to

<div id="MPEquationElement:98CC42D5-E1CD-485B-C7BE-92EE0E219AF5">

\\frac{1}{2}\\left( {{b\_{{\\rm{CaO}}}} + {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} + {b\_{{{\\rm{K}}\_2}{\\rm{O}}}} - {b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right) = \\left\[ {\\begin{array}{\*{20}{c}} 0&{ - \\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}}&0&0&0&0&0&0&0&0 \\end{array}} \\right\]\\left\[ {\\begin{array}{\*{20}{c}} {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}\\\\ {n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}\\\\ {{n\_{Qz}}}\\\\ {n\_{Ab}\^{plag}}\\\\ {n\_{An}\^{plag}}\\\\ {n\_{Sn}\^{plag}}\\\\ {n\_{Ab}\^{alk}}\\\\ {n\_{An}\^{alk}}\\\\ {n\_{Sn}\^{alk}} \\end{array}} \\right\]

</div>

or simply the one constraint:

<div id="MPEquationElement:019566F8-F9DA-41B0-9BE2-07AD25AA1972">

\\frac{1}{2}\\left( {{b\_{{\\rm{CaO}}}} + {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} + {b\_{{{\\rm{K}}\_2}{\\rm{O}}}} - {b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right) = \\frac{1}{2}\\left( {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right)

</div>

The complete set of constraints is now written:

<div id="MPEquationElement:F3DD925F-F86F-4E98-FE6E-7188D80EAACE">

\\begin{array}{c} \\frac{1}{2}\\left( {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - \\frac{1}{2}\\left( {{b\_{{\\rm{CaO}}}} + {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} + {b\_{{{\\rm{K}}\_2}{\\rm{O}}}} - {b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right) = 0\\\\ \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{qtz} - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} = 0\\\\ \\mu \_{Ab}\^{plag} - {\\Omega \_{Ab}\^{liq}} = 0\\\\ \\mu \_{An}\^{plag} - {\\Omega \_{An}\^{liq}} = 0\\\\ \\mu \_{Sn}\^{plag} - {\\Omega \_{Sn}\^{liq}} = 0\\\\ \\mu \_{Ab}\^{alk} - {\\Omega \_{Ab}\^{liq}} = 0\\\\ \\mu \_{An}\^{alk} - {\\Omega \_{An}\^{liq}} = 0\\\\ \\mu \_{Sn}\^{alk} - {\\Omega \_{Sn}\^{liq}} = 0\\\\ {n\_{Qz}} = {n\_{Fix}}\\\\ n\_{Ab}\^{plag} + n\_{An}\^{plag} + n\_{Sn}\^{plag} = {n\_{Fix}}\\\\ n\_{Ab}\^{alk} + n\_{An}\^{alk} + n\_{Sn}\^{alk} = {n\_{Fix}}\\\\ n\_{Ab}\^{plag} = {n\_{Fix}}X\_{Ab}\^{plag}\\\\ n\_{An}\^{plag} = {n\_{Fix}}X\_{An}\^{plag}\\\\ n\_{Ab}\^{alk} = {n\_{Fix}}X\_{Ab}\^{alk}\\\\ n\_{An}\^{alk} = {n\_{Fix}}X\_{An}\^{alk} \\end{array}

</div>

<div id="MPEquationElement:333B9059-0577-45D0-A2C9-B751545A9002">

(equation)

</div>

In <span class="kind elementIndex">**Equation 68** </span> there are 15 constraints, but the problem only has 12 unknowns. This apparent discrepancy is resolved by recognizing that three of these constraints are redundant. The second constraint on the chemical potential of sanidine in alkali feldspar may be removed. The Darken equation provides a definition of the chemical potentials of feldspar components:

<div id="MPEquationElement:99D765E1-DB48-4490-D3E1-3E412DCEFE63">

\\begin{array}{l} \\mu \_{Ab}\^{feld} = {{\\hat G}\^{feld}} + \\left( {1 - X\_{Ab}\^{feld}} \\right)\\frac{{\\partial {{\\hat G}\^{feld}}}}{{\\partial X\_{Ab}\^{feld}}} - X\_{An}\^{feld}\\frac{{\\partial {{\\hat G}\^{feld}}}}{{\\partial X\_{An}\^{feld}}}\\\\ \\mu \_{An}\^{feld} = {{\\hat G}\^{feld}} - X\_{Ab}\^{feld}\\frac{{\\partial {{\\hat G}\^{feld}}}}{{\\partial X\_{Ab}\^{feld}}} + \\left( {1 - X\_{An}\^{feld}} \\right)\\frac{{\\partial {{\\hat G}\^{feld}}}}{{\\partial X\_{An}\^{feld}}}\\\\ \\mu \_{Sn}\^{feld} = {{\\hat G}\^{feld}} - X\_{Ab}\^{feld}\\frac{{\\partial {{\\hat G}\^{feld}}}}{{\\partial X\_{Ab}\^{feld}}} - X\_{An}\^{feld}\\frac{{\\partial {{\\hat G}\^{feld}}}}{{\\partial X\_{An}\^{feld}}} \\end{array}

</div>

from which if follows that

<div id="MPEquationElement:5DD68DB2-D8A1-44CF-CD94-A1A046C0A7AC">

\\begin{array}{l} \\mu \_{Ab}\^{plag} = \\mu \_{Sn}\^{plag} + \\frac{{\\partial {{\\hat G}\^{plag}}}}{{\\partial X\_{Ab}\^{plag}}} = \\mu \_{Ab}\^{alk} = \\mu \_{Sn}\^{alk} + \\frac{{\\partial {{\\hat G}\^{alk}}}}{{\\partial X\_{Ab}\^{alk}}}\\\\ \\mu \_{An}\^{plag} = \\mu \_{Sn}\^{plag} + \\frac{{\\partial {{\\hat G}\^{plag}}}}{{\\partial X\_{An}\^{plag}}} = \\mu \_{An}\^{alk} = \\mu \_{Sn}\^{alk} + \\frac{{\\partial {{\\hat G}\^{alk}}}}{{\\partial X\_{An}\^{alk}}} \\end{array}

</div>

Simplification gives \\frac{{\\partial {{\\hat G}\^{plag}}}}{{\\partial X\_{Ab}\^{plag}}} = \\frac{{\\partial {{\\hat G}\^{alk}}}}{{\\partial X\_{Ab}\^{alk}}}or \\frac{{\\partial {{\\hat G}\^{plag}}}}{{\\partial X\_{An}\^{plag}}} = \\frac{{\\partial {{\\hat G}\^{alk}}}}{{\\partial X\_{An}\^{alk}}} , and consequently the equality of chemical potentials of Ab and An in two coexisting feldspars implies \\mu \_{Sn}\^{plag} = \\mu \_{Sn}\^{alk} . The constraint \\mu \_{Sn}\^{alk} - \\left( {2\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + \\frac{1}{2}\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}} \\right) = 0 may be removed because it is redundant to \\mu \_{Sn}\^{plag} - \\left( {2\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} + \\frac{1}{2}\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}} \\right) = 0 . The constraint n\_{Ab}\^{plag} + n\_{An}\^{plag} + n\_{Sn}\^{plag} = {n\_{Fix}} is also redundant, since when it is combined with n\_{Ab}\^{plag} = {n\_{Fix}}X\_{Ab}\^{plag} and n\_{An}\^{plag} = {n\_{Fix}}X\_{An}\^{plag} as {n\_{Fix}}X\_{Ab}\^{plag} + {n\_{Fix}}X\_{An}\^{plag} + n\_{Sn}\^{plag} = {n\_{Fix}} , we obtain a known identity: X\_{Ab}\^{plag} + X\_{An}\^{plag} + X\_{Sn}\^{plag} = 1 . A similar argument also eliminates the constraint: n\_{Ab}\^{alk} + n\_{An}\^{alk} + n\_{Sn}\^{alk} = {n\_{Fix}} . <span class="kind elementIndex">**Equation 68** </span>may now be simplified without loss to the twelve constraints:

<div id="MPEquationElement:CC452263-DC28-43B0-D814-E610EC195889">

\\begin{array}{c} \\frac{1}{2}\\left( {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - \\frac{1}{2}\\left( {{b\_{{\\rm{CaO}}}} + {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} + {b\_{{{\\rm{K}}\_2}{\\rm{O}}}} - {b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right) = 0\\\\ \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{qtz} - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} = 0\\\\ \\mu \_{Ab}\^{plag} - {\\Omega \_{Ab}\^{liq}} = 0\\\\ \\mu \_{An}\^{plag} - {\\Omega \_{An}\^{liq}} = 0\\\\ \\mu \_{Sn}\^{plag} - {\\Omega \_{Sn}\^{liq}} = 0\\\\ \\mu \_{Ab}\^{alk} - {\\Omega \_{Ab}\^{liq}} = 0\\\\ \\mu \_{An}\^{alk} - {\\Omega \_{An}\^{liq}} = 0\\\\ {n\_{Qz}} = {n\_{Fix}}\\\\ n\_{Ab}\^{plag} = {n\_{Fix}}X\_{Ab}\^{plag}\\\\ n\_{An}\^{plag} = {n\_{Fix}}X\_{An}\^{plag}\\\\ n\_{Ab}\^{alk} = {n\_{Fix}}X\_{Ab}\^{alk}\\\\ n\_{An}\^{alk} = {n\_{Fix}}X\_{An}\^{alk} \\end{array}

</div>

Even without formally writing down definitions of the **A** matrix, gradient vector or the Lagrangian, if follows that if the composition of the feldspars is specified as a constraint as in <span class="kind elementIndex">**Equation 72** </span> then the solution is uniquely determined without the need for potential minimization.

A generalized solution is more revealing. Suppose we seek to find the composition of a silicate liquid in equilibrium with quartz and two feldspars, but make no stipulation as to the composition of the feldspar. Under these conditions, <span class="kind elementIndex">**Equation 72** </span> transforms to:

<div id="MPEquationElement:8F3D1EAD-2C88-40A6-CFE9-1FE6CE92CD1D">

\\begin{array}{c} \\frac{1}{2}\\left( {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - \\frac{1}{2}\\left( {{b\_{{\\rm{CaO}}}} + {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} + {b\_{{{\\rm{K}}\_2}{\\rm{O}}}} - {b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right) = 0\\\\ \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{qtz} - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} = 0\\\\ \\mu \_{Ab}\^{plag} - {\\Omega \_{Ab}\^{liq}} = 0\\\\ \\mu \_{An}\^{plag} - {\\Omega \_{An}\^{liq}} = 0\\\\ \\mu \_{Sn}\^{plag} - {\\Omega \_{Sn}\^{liq}} = 0\\\\ \\mu \_{Ab}\^{alk} - {\\Omega \_{Ab}\^{liq}} = 0\\\\ \\mu \_{An}\^{alk} - {\\Omega \_{An}\^{liq}} = 0\\\\ {n\_{Qz}} = {n\_{Fix}} \\end{array}

</div>

and liberates four unknowns from the constraint set. The **A** matrix of this constraint set may now be written:

<div id="MPEquationElement:C478D83A-8A82-4560-C036-445B177C0FA5">

{\\bf{A}} = \\left\[ {\\begin{array}{\*{20}{c}} 0&{ - \\frac{1}{2}}&{\\frac{1}{2}}&{\\frac{1}{2}}&0&0&0&0&0&0&0&0\\\\ { - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&0&0&0&0&0&0\\\\ { - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&{\\frac{{\\partial \\mu \_{Ab}\^{plag}}}{{\\partial n\_{Ab}\^{plag}}}}&{\\frac{{\\partial \\mu \_{Ab}\^{plag}}}{{\\partial n\_{An}\^{plag}}}}&{\\frac{{\\partial \\mu \_{Ab}\^{plag}}}{{\\partial n\_{Sn}\^{plag}}}}&0&0&0\\\\ { - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&{\\frac{{\\partial \\mu \_{An}\^{plag}}}{{\\partial n\_{Ab}\^{plag}}}}&{\\frac{{\\partial \\mu \_{An}\^{plag}}}{{\\partial n\_{An}\^{plag}}}}&{\\frac{{\\partial \\mu \_{An}\^{plag}}}{{\\partial n\_{Sn}\^{plag}}}}&0&0&0\\\\ { - \\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&{\\frac{{\\partial \\mu \_{Sn}\^{plag}}}{{\\partial n\_{Ab}\^{plag}}}}&{\\frac{{\\partial \\mu \_{Sn}\^{plag}}}{{\\partial n\_{An}\^{plag}}}}&{\\frac{{\\partial \\mu \_{Sn}\^{plag}}}{{\\partial n\_{Sn}\^{plag}}}}&0&0&0\\\\ { - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&0&0&0&{\\frac{{\\partial \\mu \_{Ab}\^{alk}}}{{\\partial n\_{Ab}\^{alk}}}}&{\\frac{{\\partial \\mu \_{Ab}\^{alk}}}{{\\partial n\_{An}\^{alk}}}}&{\\frac{{\\partial \\mu \_{Ab}\^{alk}}}{{\\partial n\_{Sn}\^{alk}}}}\\\\ { - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}&{ - \\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}&0&0&0&0&{\\frac{{\\partial \\mu \_{An}\^{alk}}}{{\\partial n\_{Ab}\^{alk}}}}&{\\frac{{\\partial \\mu \_{An}\^{alk}}}{{\\partial n\_{An}\^{alk}}}}&{\\frac{{\\partial \\mu \_{An}\^{alk}}}{{\\partial n\_{Sn}\^{alk}}}}\\\\ 0&0&0&0&0&1&0&0&0&0&0&0 \\end{array}} \\right\]

</div>

The gradient of the Khorzhinskii potential is computed from <span class="kind elementIndex">**Equation 62** </span> as

<div id="MPEquationElement:11E59D4D-FDAA-4A57-FBF2-7AE2A1756919">

{\\bf{g}} = \\left\[ {\\begin{array}{\*{20}{c}} {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}}\\\\ {\\mu \_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}}}}\\\\ {\\mu \_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq}}}}\\\\ {\\mu \_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq}}}}\\\\ {\\mu \_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq} - {n\_{Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Ab}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{An}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\frac{{\\partial \\Omega \_{Sn}\^{liq}}}{{\\partial n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}}}}\\\\ {\\mu \_{Qz}\^o - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}\\\\ {\\mu \_{Ab}\^{plag} - \\Omega \_{Ab}\^{liq}}\\\\ {\\mu \_{An}\^{plag} - \\Omega \_{An}\^{liq}}\\\\ {\\mu \_{Sn}\^{plag} - \\Omega \_{Sn}\^{liq}}\\\\ {\\mu \_{Ab}\^{alk} - \\Omega \_{Ab}\^{liq}}\\\\ {\\mu \_{An}\^{alk} - \\Omega \_{An}\^{liq}}\\\\ {\\mu \_{Sn}\^{alk} - \\Omega \_{Sn}\^{liq}} \\end{array}} \\right\]

</div>

The Lagrangian function is:

<div id="MPEquationElement:A42C8045-0353-4ABD-BAC5-B5AE3B79EC52">

\\begin{array}{c} \\Lambda = {G\^{liq}}\\left( {n\_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq},n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{N}}{{\\rm{a}}\_{\\rm{2}}}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq},n\_{{\\rm{KAlSi}}{{\\rm{O}}\_4}}\^{liq}} \\right) + {n\_{Qz}}\\mu \_{Qz}\^o + {G\^{plag}}\\left( {n\_{Ab}\^{plag},n\_{An}\^{plag},n\_{Sn}\^{plag}} \\right) + {G\^{alk}}\\left( {n\_{Ab}\^{alk},n\_{An}\^{alk},n\_{Sn}\^{alk}} \\right)\\\\ - {n\_{Qz}}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\Omega \_{Ab}\^{liq} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\Omega \_{An}\^{liq} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\Omega \_{Sn}\^{liq}\\\\ - {\\lambda \_b}\\left\[ {\\frac{1}{2}\\left( {n\_{{\\rm{CaSi}}{{\\rm{O}}\_3}}\^{liq} + n\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{Si}}{{\\rm{O}}\_3}}\^{liq} - n\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}\^{liq}} \\right) - \\frac{1}{2}\\left( {{b\_{{\\rm{CaO}}}} + {b\_{{\\rm{N}}{{\\rm{a}}\_2}{\\rm{O}}}} + {b\_{{{\\rm{K}}\_2}{\\rm{O}}}} - {b\_{{\\rm{A}}{{\\rm{l}}\_2}{{\\rm{O}}\_3}}}} \\right)} \\right\] - {\\lambda \_{Qz}}\\left( {\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{qtz} - \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}} \\right) - \\lambda \_{Ab}\^{plag}\\left( {\\mu \_{Ab}\^{plag} - \\Omega \_{Ab}\^{liq}} \\right)\\\\ - \\lambda \_{An}\^{plag}\\left( {\\mu \_{An}\^{plag} - \\Omega \_{An}\^{liq}} \\right) - \\lambda \_{Sn}\^{plag}\\left( {\\mu \_{Sn}\^{plag} - \\Omega \_{Sn}\^{liq}} \\right) - \\lambda \_{Ab}\^{alk}\\left( {\\mu \_{Ab}\^{alk} - \\Omega \_{Ab}\^{liq}} \\right) - \\lambda \_{An}\^{alk}\\left( {\\mu \_{An}\^{alk} - \\Omega \_{An}\^{liq}} \\right) - {\\lambda \_{{n\_{Qz}}}}\\left( {{n\_{Qz}} - {n\_{Fix}}} \\right) \\end{array}

</div>

from which the elements of the Wronskian can be computed:

<div id="MPEquationElement:7058114C-DC63-4498-B92D-5B9E49C39A69">

\\begin{array}{c} \\frac{{{\\partial \^2}\\Lambda }}{{\\partial {n\_i}\\partial {n\_j}}} = \\frac{{{\\partial \^2}{G\^{liq}}}}{{\\partial {n\_i}\\partial {n\_j}}} + \\frac{{{\\partial \^2}{G\^{plag}}}}{{\\partial {n\_i}\\partial {n\_j}}} + \\frac{{{\\partial \^2}{G\^{alk}}}}{{\\partial {n\_i}\\partial {n\_j}}} - {\\delta \_{i,Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_j}}} - {\\delta \_{j,Qz}}\\frac{{\\partial \\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}}} - \\left( {{n\_{Qz}} - {\\lambda \_{Qz}}} \\right)\\frac{{{\\partial \^2}\\mu \_{{\\rm{Si}}{{\\rm{O}}\_2}}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}\\\\ - \\left( {{\\delta \_{i,plag - Ab}} + {\\delta \_{i,alk - Ab}}} \\right)\\frac{{\\Omega \_{Ab}\^{liq}}}{{\\partial {n\_j}}} - \\left( {{\\delta \_{j,plag - Ab}} + {\\delta \_{j,alk - Ab}}} \\right)\\frac{{\\Omega \_{Ab}\^{liq}}}{{\\partial {n\_i}}} - \\left( {n\_{Ab}\^{plag} + n\_{Ab}\^{alk}} \\right)\\frac{{{\\partial \^2}\\Omega \_{Ab}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}\\\\ - \\left( {{\\delta \_{i,plag - An}} + {\\delta \_{i,alk - An}}} \\right)\\frac{{\\Omega \_{An}\^{liq}}}{{\\partial {n\_j}}} - \\left( {{\\delta \_{j,plag - An}} + {\\delta \_{j,alk - An}}} \\right)\\frac{{\\Omega \_{An}\^{liq}}}{{\\partial {n\_i}}} - \\left( {n\_{An}\^{plag} + n\_{An}\^{alk}} \\right)\\frac{{{\\partial \^2}\\Omega \_{An}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}\\\\ - \\left( {{\\delta \_{i,plag - Sn}} + {\\delta \_{i,alk - Sn}}} \\right)\\frac{{\\Omega \_{Sn}\^{liq}}}{{\\partial {n\_j}}} - \\left( {{\\delta \_{j,plag - Sn}} + {\\delta \_{j,alk - Sn}}} \\right)\\frac{{\\Omega \_{Sn}\^{liq}}}{{\\partial {n\_i}}} - \\left( {n\_{Sn}\^{plag} + n\_{Sn}\^{alk}} \\right)\\frac{{{\\partial \^2}\\Omega \_{Sn}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}\\\\ - \\lambda \_{Ab}\^{plag}\\left( {\\frac{{{\\partial \^2}\\mu \_{Ab}\^{plag}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{{\\partial \^2}\\Omega \_{Ab}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}} \\right) - \\lambda \_{An}\^{plag}\\left( {\\frac{{{\\partial \^2}\\mu \_{An}\^{plag}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{{\\partial \^2}\\Omega \_{An}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}} \\right) - \\lambda \_{Sn}\^{plag}\\left( {\\frac{{{\\partial \^2}\\mu \_{Sn}\^{plag}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{{\\partial \^2}\\Omega \_{Sn}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}} \\right)\\\\ - \\lambda \_{Ab}\^{alk}\\left( {\\frac{{{\\partial \^2}\\mu \_{Ab}\^{alk}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{{\\partial \^2}\\Omega \_{Ab}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}} \\right) - \\lambda \_{An}\^{alk}\\left( {\\frac{{{\\partial \^2}\\mu \_{An}\^{alk}}}{{\\partial {n\_i}\\partial {n\_j}}} - \\frac{{{\\partial \^2}\\Omega \_{An}\^{liq}}}{{\\partial {n\_i}\\partial {n\_j}}}} \\right) \\end{array}

</div>

# 6. Bibliography {#MPSection:DD6D405C-05BC-482E-DC3F-62C9120E8FE3}

<div id="MPBibliographyElement:5C8B992A-1764-4F74-C684-0F123A8B7662">

</div>
