# Codegen Notebooks

This directory contains Jupyter notebooks and supporting files that illustrate use of the coder module of the ThermoEngine package.  Coder contains methods that permit thermodynamic models to be generated that describe the properties of both stoichiometric compounds and solutions.  Coder uses the SymPy Python package for symbolic mathematics and for computer code generation.  The notebooks illustrate how to specify the structure of thermodynamic models and how to automatically generate computer code (e.g., C-code), wrap that in a Python module, and utilize the generated module for property calculation and for model calibration. 

The notebooks are divided into *Examples*, *Demos*, *Utilities* and *Documentation*, the latter are found in the Documentation sub-folder. Deprecated notebooks are located in a correspondingly labeled sub-folder. 

Start with the example notebooks if you are new to the coder module

### Examples of use cases

- [**Example-1-Berman-std-state**](./Example-1-Berman-std-state.ipynb)  
Demonstrates how to construct the Berman standard state property model for the Gibbs free energy of a stoichiomtric phase, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. This notebook generates files in a working directory called *working*.
- [**Example-2-Berman-plus-BM**](./Example-2-Berman-plus-BM.ipynb)   
Demonstrates how to construct a "Berman-like" model for standard state properties of stoichiometric phases using a Berman-Brown heat capacity model and a third order Birch-Murnaghan Equation of State model. "C-" code for fast calculation and for calibration of model parameters is generated. This notebook generates files in a working directory called *working*.
- [**Example-3-HKF**](./Example-3-HKF.ipynb)   
Demonstrates how to construct a model for the standard state properties of aqeuous species using the revised Helgeson, Kirkham and Flowers standard state formulation, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters.  This notebook generates files in a working directory called *aqueous*. This notebook depends on the file *slop16_v3_1.dict*, which is generated from the standard SLOP file, *slop16_v3_1.dat* using the utility notebook [**Utility-Read-SLOP-file**](./Utility-Read-SLOP-file.ipynb).
- [**Example-4-Stixrude-Debye**](./Example-4-Stixrude-Debye.ipynb)   
Demonstrates how to construct a Stixrude standard state property model for the Helmholtz free energy of a stoichiometric phase, using a Debye model for the entropy.  "C-" code for fast calculation of the Gibbs free energy and for calibration of model parameters is generated. This notebook generates files in a working directory called *working*.
- [**Example-5-HDNB**](./Example-5-HDNB.ipynb)   
Demonstrates how to construct a model for the standard state properties of minerals and gases using the Helgeson, Delany, Nesbitt and Bird standard state formulation and database.  A "C-" coded implementation for fast calculation and for calibration of model parameters is produced.  This notebook generates files in a working directory called *hdnb*. This notebook depends on the file *slop16_v3_1.dict*, which is generated from the standard SLOP file, *slop16_v3_1.dat* using the utility notebook [**Utility-Read-SLOP-file**](./Utility-Read-SLOP-file.ipynb).
- [**Example-6-HDNB-Quartz**](./Example-6-HDNB-Quartz.ipynb)   
Demonstrates how to construct a model for the standard state properties of quartz using the Helgeson, Delany, Nesbitt and Bird standard state formulation and database, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters.  This notebook generates files in a working directory called *hdnb*. This notebook depends on the file *slop16_v3_1.dict*, which is generated from the standard SLOP file, *slop16_v3_1.dat* using the utility notebook [Utility-Read-SLOP-file](./Utility-Read-SLOP-file.ipynb).
- [**Example-7-Simple-Solution**](./Example-7-Simple-Solution.ipynb)  
Demonstrates use of the coder package to construct an asymmetric regular solution with ternary terms, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. The example illustrated is the feldspar solid solution model used in MELTS. This notebook generates files in a working directory called *working*.
- [**Example-8-Importing-modules**](./Example-8-Importing-modules.ipynb)  
Demonstrates how to import coder modules into the ThermoEngine package, illustrating method calls for both pure (stoichiometric) and solution phases.  This notebook requires files previously written to a working directory called *working* by example notebooks one and seven. 
- [**Example-9-Complex-Solution**](./Example-9-Complex-Solution.ipynb)  
Demonstrates use of the coder package to construct a thermodynamic model for a reciprocal solution with ordering parameters, and from that model produce a "C-" coded implementation for fast calculation and for calibration of model parameters. The example illustrated is the orthopyroxene solid solution model calibrated by Sack and Ghiorso (1994). This notebook generates files in a working directory called *working*. 
- [**Example-10-gas-Speciation-Solution**](./Example-10-gas-Speciation-Solution.ipynb)  
Demonstrates use of the coder package to construct a thermodynamic model for a speciated ideal gas solution.  The formulation uses ordering parameters, and produces a "C-" coded implementation for calculation and calibration of model parameters. The example illustrated is a simplified ideal gas model based on code for a solar condensation model being developed by Grayson Boyer.  The speciated solution model is derived from the complex solution model class. It utilizes a more robust numerical method to solve for species mole fractions (ordering parameters) that allows for extreme differences in the abundance of species in the equilibrium state.  This notebook generates files in a working directory called *working*.  
- [**Example-13-Sulfide_liquid**](./Example-13-Sulfide_liquid.ipynb)  
Demonstrates use of the ThermoEngine sulfide_liquid module, which implememnts the thermodynamic model of Kress for liquids in the system O-S-Fe-Ni-Cu. Two methods of accessing the methods of this module are illustrated. 
### Demo and Utility Notebooks
- [**Demo-Access-Module**](./Demo-Access-Module.ipynb)  
Demonstrates how to access a previously generated and compiled module from an independent notebook.
- [**Utility-Read-SLOP-file**](./Utility-Read-SLOP-file.ipynb)  
Reads a SLOP file (*slop16_v3_1.dat*) and generates a pickled dictionary whose keys are mineral, gas and aqueous species datatypes contained in the original file, and whose values are Pandas databases of parameter values.
- [**azero**](./azero.ipynb)   
Reference notebook that contains thermodynamic identities for the Debye-Hückel azero functions and its temperature, pressure and compositional derivatives.
- [**ibar**](./ibar.ipynb)   
Reference notebook that contains thermodynamic identities for the true ionic strtength function and its compositional derivative.
- [**Mu-excess**](./Mu-excess.ipynb)   
Reference notebook that contains thermodynamic identities for the Debye-Hückel excess chemical potential function and its temperature, pressure and compositional derivatives.

### Documentation Notebooks
- [**Born_functions**](./Documentation/Born_functions.ipynb)   
Reference notebook that contains thermodynamic identities for the Born functions and their T,P derivatives.
- [**DH_functions**](./Documentation/DH_functions.ipynb)   
Reference notebook that contains thermodynamic identities for the Debye-Hückel functions and their T,P derivatives.
- [**Helmholtz-to-Gibbs**](./Documentation/Helmholtz-to-Gibbs.ipynb)   
Reference notebook that contains thermodynamic identities for converting Helholtz functions and their T,V derivatives into Gibbs functions and their T,P derivatives.
- [**StdState-piecewise**](./Documentation/StdState-piecewise.ipynb)   
Reference notebook that contains an illustration of how to use the SymPy *Piecewise* function to construct thermodynamic identities for models that contain mathematical terms that are applicable over restricted domains, e.g. "lambda-type" heat capacity functions.
- [**Theory-Complex-Solutions-docs**](./Documentation/Theory-Complex-Solutions-docs.ipynb)   
Reference notebook that contains thermodynamic relations for computing equilibrium properties of solutions that are formulated with one or more ordering parameters.  The ordering parameters are implicit functions of composition, temperature and pressure.

### Deprecated Notebooks and Files

These notebooks and dependent files are no longer supported and were developed mainly for testing code in the coder module. Caveat emptor.  

- [**berman_1988.json**](./Deprecated/berman_1988.json)  
Parameter file for the Berman (1988) database in JSON format.
- [**Berman-SymPy-Testing**](./Deprecated/Berman-SymPy-Testing.ipynb)  
Demonstrates how to create C-code for calculating standard state properties using the model framework of Berman (1988).
- [**Born_SymPy_and_HKF**](./Deprecated/Born_SymPy_and_HKF.ipynb)  
Notebook for testing Born function implementation in SymPy.
- [**elements.csv**](./Deprecated/elements.csv)  
Comma-separated-value file of element names, abbreviations and molecular weights.
- [**entropies.csv**](./Deprecated/entropies.csv)  
Comma-separated-value file of third law entropies of the elements used in converting Berman database reference state conditions from the conventional mode to that adopted by Helgeson et al. (1978). 
- [**Simple-Solution-SymPy-Testing**](#simsol)  
Demonstrates how to create C-code for calculating n-component solution properties using a simple configurational entropy model and excess enthalpies defined by an asymmetric sub-regular solution with ternary terms. Generates code for a feldspar solid solution as an example implementation.
- [**SWIM-BornAndTests**](./Deprecated/SWIM-BornAndTests.ipynb)  
Notebook for testing SWIM water code and Born functions.
- [**SymPy-Deriv-print**](./Deprecated/SymPy-Deriv-print.ipynb)  
Notebook for testing implementation of a SymPy derivative printer.
- [**Development-Non-convergent-Ordering-Solution**](./Development-Non-convergent-Ordering-Solution.ipynb)  
Incomplete implementation.  Illustration, using the Sack and Ghiorso model for multicomponent pyroxenes, of the generation of a thermodynamic solution model that requires the computation of cation site occupancy between symmetrically non-equivalent sites, necessitating the computation of internal (homogeneous) equilibrium.  Code generation is not implemented.

