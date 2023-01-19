.. PhaseObjC documentation master file, created by
   sphinx-quickstart on Mon Jun 19 20:37:00 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ThermoEngine documentation 
**************************
Welcome to the ThermoEngine documentation! 

This reference material describes the Python modules, and their methods, available in ThermoEngine. It also provides details about the libPhaseObjC compiled library from which the methods are called.

Python modules
==============
This section describes the Python modules and how to call Python functions that (1) compute thermodynamic properties of phases, and (2) calculate equilibrium phase assemblages. 

The classes and functions described call methods from the :ref:`library-label`. To permit Python to interact with this library, the modules use Rubicon.  

Phases
------

.. toctree::
   :maxdepth: 2

   phases

The **Phases module** of ThermoEngine documents Python functions and classes for computing standard thermodynamic calculations utilizing the Berman, Holland and Powell, or Stixrude-Lithgow-Bertelloni endmember databases, and calculations based on solution properties utilized by MELTS. There are many helper functions available in this module that assist in the calculation of pseudosections, univariant equilibria and the construction of phase diagrams.

    Examples:  
        * :doc:`Codegen examples<examples_codegen>`
        * :doc:`DEW (Deep Earth Water) examples<examples_DEW>`
        * :doc:`Equilibrate examples<examples_equilibrate>`               
        * :doc:`Pure Phases examples<examples_pure-phases>` 
        * :doc:`Solutions examples<examples_solutions>` 


Equilibrate
-----------

.. toctree::
   :maxdepth: 2

   equilibrate

The **Equilibrate module** of ThermoEngine documents Python functions and classes for computing equilibrium phase assemblages with focus on MELTS calculations.

    Examples:  
        * :doc:`Equilibrate examples<examples_equilibrate>`  
        * :doc:`MELTS-pMELTS examples<examples_MELTS-pMELTS>`  


Coder, Coder Templates
-----

.. toctree::
   :maxdepth: 2

   coder
   coder_templates

The **Coder module** of ThermoEngine documents Python functions that generate thermodynamic models for pure component stoichiometric compounds and for solutions. The module provides methods that generate computer code for both fast computation and for calibration of the specified models. The Coder module relies upon another module, **Coder Templates**, whose methods are generally not directly called by users.

    Examples:  
        * :doc:`Codegen examples<examples_codegen>` 
        * :doc:`Equilibrate examples<examples_equilibrate>`  
        
Calibrate
---------

.. toctree::
   :maxdepth: 2

   calibrate

The **Calibrate module** of ThermoEngine documents Python functions that support calibration of thermodynamic models.  The models may be imported from the Model module or may be generated using the Coder module.


Core
----

.. toctree::
   :maxdepth: 2

   core

The **Core module** contains Python support functions shared by the Phases and Equilibrate modules. Most users do not need to call these functions directly.


Model
-----

.. toctree::
   :maxdepth: 2

   model 

The **Model module** of ThermoEngine implements a Python interface with the Phase objective-C classes as well as the infrastructure for pure phase thermodynamic calibration. The module contains methods that allow for loading and selection of built-in thermodynamic databases.

    Examples:  
        * :doc:`Codegen examples<examples_codegen>`
        * :doc:`Equilibrate examples<examples_equilibrate>`  
        * :doc:`Pure Phases examples<examples_pure-phases>`
        * :doc:`Solutions examples<examples_solutions>`


.. _library-label:
   
libPhaseObjC compiled library
=============================
This section covers phase property and equilibrium calculators that are available in the libPhaseObjC compiled library. 

Except for protocol declarations, each class name is represented in the repository by a header file (.h) and usually by an associated Objective-C source file (.m).

.. toctree::
   :maxdepth: 1
   
   EquilibrateCalcs.md
   PhaseClassesDoc.md

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
