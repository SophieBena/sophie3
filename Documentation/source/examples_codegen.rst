Codegen Examples (Jupyter notebooks)
************************************

The **Notebooks/Codegen** folder contains Jupyter notebooks that generate code for the example notebooks. They are in development and subject to change.

    Modules used:  
        * Coder
        * Coder Templates       

**The links below access static versions only.** You can access executable versions from the Notebooks folder.

❇️ :doc:`Example 1: Berman Standard State Code Generator<Example-1-Berman-std-state>`
  
Demonstrates how to construct the Berman standard state property model for the Gibbs  free energy of a stoichiomtric phase, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. This notebook generates files in a working directory called *working*.   
  
❇️ :doc:`Example 2: Berman Standard State Code Generator + Birch-Murnaghan<Example-2-Berman-plus-BM>`  
  
Demonstrates how to construct a "Berman-like" model for standard state properties of stoichiometric phases using a Berman-Brown heat capacity model and a third order Birch-Murnaghan Equation of State model. "C-" code for fast calculation and for calibration of model parameters is generated. This notebook generates files in a working directory called *working.*  
  
❇️ :doc:`Example 3: Revised HKF - Standard State<Example-3-HKF>`  
    
Demonstrates how to construct a model for the standard state properties of aqeuous species using the revised Helgeson, Kirkham and Flowers standard state formulation, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. This notebook generates files in a working directory called *aqueous*. This notebook depends on the file *slop16_v3_1.dict*, which is generated from the standard SLOP file, *slop16_v3_1.dat*, using the utility notebook **Utility-Read-SLOP-file**. 
 
❇️ :doc:`Example 4: Helmholtz energy (Stixrude - Debye integral)<Example-4-Stixrude-Debye>`  
  
Demonstrates how to construct a Stixrude standard state property model for the Helmholtz free energy of a stoichiometric phase, using a Debye model for the entropy.  "C-" code for fast calculation of the Gibbs free energy and for calibration of model parameters is generated. This notebook generates files in a working directory called *working*.  

❇️ Example 5: Helgeson, Delany, Nesbitt and Bird Standard State Code Generator

Demonstrates how to construct a model for the standard state properties of minerals and gases using the Helgeson, Delany, Nesbitt and Bird standard state formulation and database.  A "C-" coded implementation for fast calculation and for calibration of model parameters is produced.  This notebook generates files in a working directory called *hdnb*. This notebook depends on the file *slop16_v3_1.dict*, which is generated from the standard SLOP file, *slop16_v3_1.dat*, using the utility notebook **Utility-Read-SLOP-file**.
  
For more information, see the executable notebook (Example-5-HDNB).

❇️ Example 6: Helgeson, Delany, Nesbitt and Bird Standard State Code Generator

Demonstrates how to construct a model for the standard state properties of quartz using the Helgeson, Delany, Nesbitt and Bird standard state formulation and database, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters.  This notebook generates files in a working directory called *hdnb*. This notebook depends on the file *slop16_v3_1.dict*, which is generated from the standard SLOP file, *slop16_v3_1.dat*, using the utility notebook **Utility-Read-SLOP-file**.
  
For more information, see the executable notebook (Example-6-HDNB-Quartz).

❇️ :doc:`Example 7: Simple Solution SymPy Code Generation<Example-7-Simple-Solution>`  
  
Demonstrates use of the Coder package to construct an asymmetric regular solution with ternary terms, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. The example illustrated is the feldspar solid solution model used in MELTS. This notebook generates files in a working directory called *working*.  
  
❇️ :doc:`Example 8: How to import a coder generated module into ThermoEngine<Example-8-Importing-modules>`  
  
Demonstrates how to import Coder modules into the ThermoEngine package, illustrating method calls for both pure (stoichiometric) and solution phases.  This notebook requires files previously written to a working directory called *working* by example notebooks 1 and 7.
  
❇️ :doc:`Example 9: Complex Solution SymPy Code Generation<Example-9-Complex-Solution>`  
  
Demonstrates use of the Coder package to construct a thermodynamic model for a reciprocal solution with ordering parameters, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. The example illustrated is the orthopyroxene solid solution model calibrated by Sack and Ghiorso (1994). This notebook generates files in a working directory called *working*. 
  
❇️ Example 10: Speciation Solution SymPy Code Generation

Demonstrates use of the Coder package to construct a thermodynamic model for a speciated ideal gas solution. The formulation uses ordering parameters and produces a "C-" coded implementation for calculation and calibration of model parameters. The example illustrated is a simplified ideal gas model based on code for a solar condensation model being developed by Grayson Boyer. The speciated solution model is derived from the complex solution model class. It utilizes a more robust numerical method to solve for species mole fractions (ordering parameters) that allows for extreme differences in the abundance of species in the equilibrium state. This notebook generates files in a working directory called *working*.
  
For more information, see the executable notebook (Example-10-gas-Speciation-Solution).

❇️ Example 11: Speciation Solution SymPy Code Generation

For details, see the executable notebook (Example-11-non-ideal-Speciation-Solution).

❇️ :doc:`Example 12: Standard Water Intergrated Model<Example-12-SWIM>`  
  
Demonstrates SWIM accessed using Cython wrappers to C/C++ code. 

❇️ :doc:`Example 13: Sulfide liquid model of Kress<Example-13-Sulfide_liquid>`  
  
Demonstrates use of the ThermoEngine sulfide_liquid module, which implements the thermodynamic model of Kress for liquids in the system O-S-Fe-Ni-Cu. Two methods of accessing the methods of this module are illustrated.
