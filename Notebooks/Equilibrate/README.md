# Equilibrate notebooks
Jupyter notebooks that illustrate features and capabilities of the Equilibrate class, which is defined in the equilibrate module.  The Equilibrate class computes equilibrium phase assemblages given a collection of phases, a bulk composition, and specified values of two thermodynamic variables, e.g., (T,P), (S,P), (T,V) or (S,V).  The class can also compute equilibrium in systems open to mass transfer under conditions of specified chemical potential constraints, e.g., fixed oxygen fugacity, fixed saturation state or fluid or mineral phases. 

### Equilibrate Python Module

- **Equilibrate class**
- **MELTS class** See MELTS-pMELTS notebook folder

### Equilibrium example notebooks

Notebooks that demonostrate use of the Equilibrate class for 

- **1-G-Equilibrate** Gibbs free energy minimization; closed system; rhyolite-MELTS
- **2-K-O2-Equilibrate** Gibbs free energy minimization; open system; fixed chemical potential of oxygen using the empirical constraint of Kress and Carmichael (1991); rhyolite-MELTS
- **3-K-O2-Equilibrate** Khorzhinskii potential minimization in an open system; fixed chemical potential of water
- **4-H-Equilibrate** Enthalpy minimization; closed system; fixed entropy and pressure 
- **5-A-Equilibrate** Helmholtz energy minimization; closed system; fixed temperature and volume
- **6-E-Equilibrate** Internal energy energy minimization; closed system; fixed eentropy and volume
- **7a-K-Quartz-Equilibrate** Khorzhinskii potential minimization in an open system; fixed chemical potential of silica (quartz)
- **7b-K-Cor-Equilibrate** Khorzhinskii potential minimization in an open system; fixed chemical potential of alumina (corundum)
- **8-K-Qtz-Cor-Equilibrate** Khorzhinskii potential minimization in an open system; fixed chemical potential of silica (quartz) and alumina (corundum)
- **9-K-Qtz-Fld-Equilibrate** Khorzhinskii potential minimization in an open system; fixed chemical potential of silica (quartz) and feldspar saturation
- **9-Olivine-loop** Calculation of the olivine solid-liquid phase loop

### Development notebooks

Notebooks that illustrate features of the Equilibrate class that are under development.

- **9-Testing-Equilibrate** Detailed tests to assist in understanding the algoithm underlying notebook 9
- **9a-Pseudo-phase-MELTS** Algorithmic details for omnicomponent pseudo-phase generation using MELTS as a basis
- **9b-Pseudo-phase-MELTS** Algorithmic details for omnicomponent pseudo-phase generation using MELTS as a basis
- **9c-Pseudo-phase-Stixrude** Algorithmic details for omnicomponent pseudo-phase generation using Stixrude as a basis
- **9d-Pseudo-phase-Affinity** Notebook to test returns from Objective-C versions of affinity_and_comp

### Equilibrium class theory notebooks
In the *Generalized Equilibrium Paper* folder, the three notebooks

- **Equilibrate-extension-pure-phases**
- **Equilibrate-extension-solutions**
- **Equilibrate-extension-oxygen-water**

support a paper in progress that develops theory for calculating phase equilibria in systems subject to various constraints, e.g. given that the following minerals coexist with a silicate liquid at some T, P condition, at equilibrium, what is the composition of the liquid?

The manuscript in progress accompanying these notebooks may be found in the notebook **Generalized-Equilibrium**. 

### Equilibrium class testing and deprecated notebooks
These notebooks are located in the folder *Deprecated*.

- **1-Testing-Equilibrate** Gibbs free energy minimization with garnet
- **9-Testing-Equilibrate** Testing details of minimization algorithm with liquid+feldspar
- **MELTS-v1.0.2-equilibrium** Objective-C implementation of MELTS (see MELTS notebook folder)
- Other miscelaneous notebooks 


