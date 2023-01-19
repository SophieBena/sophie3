Examples (Jupyter notebooks)
****************************

The **Notebooks** folder contains Jupyter notebooks that illustrate various features of the ENKI software infrastructure. 
   
The notebooks are organized in subfolders: 

* `Codegen`_
* `DEW (Deep Earth Water)`_
* `Deprecated`_
* `Development`_
* `MELTS-pMELTS`_
* `Pure-Phases`_
* `Solutions`_
* `Utilities`_

**The links below access static versions only.** You can access executable versions from the Notebooks folder.

Codegen  
=======  
These notebooks, which generate code for the example notebooks, are in development and subject to change.

- Berman-SymPy-Testing

  Demonstrates how to create C-code for calculating standard state properties using the model framework of Berman (1988).

- Simple-Solution-SymPy-Testing

  Demonstrates how to create C-code for calculating n-component solution properties using a simple configurational entropy model and excess enthalpies defined by an asymmetric sub-regular solution with ternary terms. Generates code for a feldspar solid solution as an example implementation.

- Non-convergent-Ordering-Solution

- Convergent-Ordering-Solution

  Not yet implemented.


DEW (Deep Earth Water)  
======================  
These notebooks implement various calculations involving DEW, the Deep Earth Water model (Sverjensky et al., 2014, Geochimica et Cosmochimica Acta 129, 125-145).  

- :doc:`DEW-QuartzSolubility<DEW-QuartzSolubility>`

  Shows how to calculate quartz solubility in aqueous solutions using the Deep Earth Water model.  

- :doc:`DEW-Standard-State<DEW-Standard-State>`

  Shows how to calculate standard state properties of ions and neutral species using parameters from the Deep Earth Water model.  

- :doc:`DEWFluid<DEWFluid>`

  Illustrates how to calculate a variety of thermodynamic properties given a specified fluid composition using the Deep Earth Water model.

Deprecated  
==========  
These Jupyter notebooks illustrate features of the ENKI library or calculational methods that are deprecated. 

- :doc:`Quartz-Berman-oldstyle<Quartz-Berman-oldstyle>`

  Illustrates calculation of the thermodynamic properties of quartz using the Berman (1988) database and direct access methods to the underlying code framework. Use of the old-style direct access methods is discouraged in favor of the Python wrapper methods illustrated in the other notebooks.  

- :doc:`Quartz-Holland-Powell-oldstyle<Quartz-Holland-Powell-oldstyle>`

  Illustrates calculation of the thermodynamic properties of quartz using the Holland and Powell (2000) database and direct access methods to the underlying code framework. Use of the old-style direct access methods is discouraged in favor of the Python wrapper methods illustrated in the other notebooks.

- :doc:`DEW-Muscovite and DEW-muscovite-fast<DEW-Muscovite and DEW-muscovite-fast>`

  Demonstrate solubility calculations involving muscovite and Cl-rich fluids at elevated pressure and temperature. The file **SverjenskyEtAl1991-Fig1a.png** is required for their execution.


Development  
===========  
The notebooks in development illustrate new features and capabilities of ENKI software infrastructure.

Equilibrium calculators
+++++++++++++++++++++++++

- Equilibrate-extension-pure-phases

- Equilibrate-extension-solutions

- Equilibrate-extension-oxygen-water

Support a paper in progress that develops theory for calculating phase equilibria in systems subject to various constraints, e.g., given that the following minerals coexist with a silicate liquid at some T, P condition, at equilibrium, what is the composition of the liquid?
 
The manuscript is located in the folder **ms-Generalized-Equilibrium**.

MELTS-DEW integration
+++++++++++++++++++++

- MELTS-DEW 

Is a prototype for calculations involving the coupled MELTS and DEW solution models.


Thermodynamic Database Calibration
++++++++++++++++++++++++++++++++++

The remaining notebooks in the Development folder deal with Bayesian methods for calibration of pure-component thermodynamic databases. These methods and code are under active development. Data files in the folder **phase-reversal-data** are required for execution of these notebooks.

MELTS-pMELTS  
============  
These notebooks implement MELTS and pMELTS model calculations. See the `MELTS website <http://melts.ofm-research.org>`_ for further information about which version of MELTS/pMELTS to use.  

- :doc:`MELTS-v1.0.2-equilibrium<MELTS-v1.0.2-equilibrium>`

  Illustrates equilibrium crystallization calculations over a sequence of temperatures using the original calibration of rhyolite-MELTS.

- :doc:`MELTS-v1.0.2-fractionation<MELTS-v1.0.2-fractionation>`

  Illustrates fractional crystallization calculations over a sequence of temperatures using the original calibration of rhyolite-MELTS.

- :doc:`MELTS-v1.1.0-equilibrium<MELTS-v1.1.0-equilibrium>`

  Illustrates equilibrium crystallization calculations over a sequence of temperatures using the H\ :sub:`2`\ O-CO\ :sub:`2`\  fluid model of Ghiorso and Gualda, adjusted to recover the ternary minimum in the two-feldspar-quartz-fluid saturated "ternary" system.

- :doc:`MELTS-v1.2.0-equilibrium<MELTS-v1.2.0-equilibrium>`

  Illustrates equilibrium crystallization calculations over a sequence of temperatures using the H\ :sub:`2`\ O-CO\ :sub:`2`\  fluid model of Ghiorso and Gualda. This model should not be used for the two-feldspar-quartz-fluid saturated "ternary" system.

- :doc:`pMELTS-v5.6.1-melting<pMELTS-v5.6.1-melting>`

  Illustrates equilibrium partial melting calculations over a sequence of temperatures and pressures using the pMELTS model on a peridotite bulk composition.

- :doc:`pMELTS-v5.6.1-adiabatic<pMELTS-v5.6.1-adiabatic>`

  Illustrates equilibrium partial melting calculations over a sequence of pressures under adiabatic constraints using the pMELTS model on a peridotite bulk composition.

Pure-Phases  
===========  
These notebooks illustrate calculation of thermodynamic properties of pure component phases using established thermodynamic databases.

- :doc:`Quartz-Berman<Quartz-Berman>`

  Illustrates calculation of the thermodynamic properties of quartz using the Berman (1988) database and the Python ThermoEngine package.

- :doc:`Forsterite-Stixrude<Forsterite-Stixrude>`

  Illustrates calculation of the thermodynamic properties of forsterite using the Stixrude-Lithgow-Bertelloni (2012) database and the Python ThermoEngine package.

- :doc:`Phase-Diagram<Phase-Diagram>`

  Illustrates construction of the alumino-silicate phase diagram using the Berman (1988) database and the Python ThermoEngine package.

- :doc:`Compare-Phases<Compare-Phases>`

  Illustrates how to compare thermodynamic properties of quartz calculated from three different databases: Berman (1988), Holland and Powell (2000) and Stixrude-Lithgow-Bertelloni (2012), using the Python ThermoEngine package. This notebook requires the file **Ghiorso-cp-quartz.txt**. 

- :doc:`Plot-Reaction<Plot-Reaction>`

  Illustrates how to calculate and plot a univariant curve for the reaction forsterite + quartz -> enstatite using the Berman (1988) database and the Python ThermoEngine package.

- :doc:`Plot-Reaction-Stixrude<Plot-Reaction-Stixrude>`

  Illustrates how to calculate and plot a univariant curve for the reaction forsterite -> magnesio-wadsleyite using the Stixrude-Lithgow-Bertelloni (2012) database and the Python ThermoEngine package.

- :doc:`Water<Water>`

  Illustrates calculation of the thermodynamic properties of water using the Simulated Water Interpolated Model (SWIM) and direct access methods to the underlying code framework.

Solutions  
=========  
These notebooks illustrate calculation of the thermodynamic properties of non-aqueous solutions.

- :doc:`Clinopyroxene-MELTS<Clinopyroxene-MELTS>`

  Illustrates calculation of the thermodynamic properties of clinopyroxene solid solutions (as utilized in MELTS) utilizing ThermoEngine Python wrappers.  

- :doc:`Feldspar-MELTS<Feldspar-MELTS>`

  Illustrates calculation of the thermodynamic properties of feldspar solid solutions (as utilized in MELTS) utilizing ThermoEngine Python wrappers.  

- :doc:`Pyroxene-Geothermometer<Pyroxene-Geothermometer>` 

  Illustrates construction of a two-pyroxene geothermometer using thermodynamic properties from teh MELTS package as exposed by ThermoEngine Python wrappers. Calculations rely on Python optimization tools to compute the equilibrium temperature.


Utilities  
=========  
These Jupyter notebooks are designed to perform once only or seldom required operations.
