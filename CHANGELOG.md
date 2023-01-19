# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0-beta-10] - 2021-10-14 (wolf)

### Added
- New method for equilibrating subsolidus phase assemblages (absent an omnicomponent phase); method is thoroughly tested for Stixrude database, but does not currently apply to others (which posess complicating reciprocal solutions).
- Sizable rapid developer test suite for working locally and for CI pipeline; Ensures that all new code returns expected values for all public interfaces in newly expanded codebase and tests some of the legacy codebase.
- Provides simple (extensible) interfaces for defining and equilibrating systems; provides swappable techniques for equilibration method.
- Defines phase samples as phase models evaluated at particular T, P, X conditions; key for allowing immiscible coexisting samples from a single phase.

### Changed
- CI/CD yaml file is modified to run both developer AND notebook tests, and report combined coverage to coveralls.
- Equilibration of System object only offers non-omnicomponent algorithm (Gridded solution in place of MELTS).

## [1.0.0-beta-9] - 2021-04-02 (ghiorso)

### Added 
- Cluster infrastructure support for Knative web services
- Notebook version of prototype algorithm for rhyolite-MELTS geobarometer and liquidus surface phase diagram calculator
- Environment and runner tags added to CI/CD yaml files
- ENKI production server cluster on Google Cloud added to Gitlab managed clusters; installed Prometheus for monitoring

### Changed
- CI/CD yaml file for GitLab to use Kaniko instead of Docker for container builds
- ENKI production server by installing the Knative server component (required upgrading default node pool to 4 compute nodes and reconfiguring it to autoscale)
- Kubernetes software upgrade on all cluster nodes to latest stable version (required for Knative)

### Fixed
- Docker build files to correctly construct images when building in a branch other than master

## [1.0.0-beta-8] - 2021-02-12 (ghiorso, spiegelman, johnson)
### Added 
- Methods and enhancements to ComplexSoln class in the coder module
- Added shell script for running the CI generated Docker container locally

### Changed
- Updated README and Cluster documentation
- Updated package documentation
- Updated CI scripts and other build/cluster related scripts; forced production server to always use latest container image
- Updated Stixrude Debye model in coder package
- Updated code to address an issue with Ubuntu Rubicon
- Other miscellaneous minor updates

## [1.0.0-beta-7] - 2020-08-20 (wolf)
### Added
- Bug-fix for sulfide_liquid module by adding cmake dependency
- Simplify integration of MELTS Fe-Ni alloy model
  - Liquid & Solid Alloy phases available as 'MtlL' and 'MtlS'
  - demonstrate fO2-dependent metal segregation in Equilibrate/Metal-Segregation NB
- Add work in development for sub-solidus ('sunken-hull') equilibration method
  - code under development in Notebooks/Equilibrate/develop
  - to be incorporated when fully functional

## [1.0.0-beta-6] - 2020-07-15 (ghiorso)
### Added
- A new module to ThermoEngine (sulfide_liquid) that implements the thermodynamic model of Kress for liquids in the system O-S-Fe-Ni-Cu
- Jupyter notebook in Codegen folder to test this module

### Changed
- Testing scripts to include the additional notebook

## [1.0.0-beta-5] - 2020-06-25 (ghiorso)
### Added
- Jupyter notebook to provide edge case tests for Phase class affinity_and_comp() method
- Objective-C derived class for reciprocal solutions of the Stixrude thermodynamic database. The derived class provides a method for the calculation of affinity and composition from specified chemical potentials that utilizes both component species and dependent species in the calculation. This addition extends the range of permitted compositions for the affected solutions.
- Jupyter notebook to test new affinity and composition method for Stixrude reciprocal solutions

### Changed
- Algorithm for the calculation of affinity and composition from specified chemical potentials for Stixrude reciprocal solutions (opx, cpx, garnet)
- Revised configurational entropy model for cpx and garnet in the Stixrude database (see notes in testing notebook)

### Fixed
- Objective-C implementation of getAffinityAndComposition... for Berman database (also MELTS) solution models OpxBerman, CpxBerman, OlivineBerman. Added bound constraints on derived composition to catch machine underflow.

## [1.0.0-beta-4] - 2020-06-23 (ghiorso)
### Added
- Jupyter notebook to demonstrate calculation of the olivine liquid-solid phase loop

### Changed
- Default behavior of the Database class in the model module. When a coder module generated phase is imported, "H2O" and "Liquid" phases are no longer included as extra phases in the returned database.

### Fixed
- Bugs are corrected related to the calculation of simple system phase relations using the Equilibrate class. These fixes apply to systems where the number of components in the omnicomponent phase is fewer than the number of distinct elements used to describe the composition of the phase. For example, calculation of the olivine phase loop in the system Mg2SiO4-Fe2SiO4.

## [1.0.0-beta-3] - 2020-05-22 (ghiorso)
### Added
- Constrained oxygen chemical potential open system calculations to the equilibrate Python module
- Jupyter notebooks demonstrating and testing the above
- Access to the phases module specialized affinity calculators that are available in objective-C coded solution phases

### Changed
- Notebooks/README in the Equilibrate notebook folder
- Shell and CI scripts pertaining to code testing and coverage, adding references to additional open-system equilibrium calculations

### Fixed
- Core module calls to Objective-C vector, matrix and tensor storage classes to speed up access to the underlying data structures
- Tuned default numerical threshold and convergence constants in the Equilibrate Python class

## [1.0.0-beta-2] - 2020-04-24 (ghiorso)
### Added
- CI/CD enhancements for Docker image generation, code testing, and code coverage when pushing to the master branch
- Code coverage is exposed at https://coveralls.io

### Changed
- README, Cluster/README, MAINTENANCE, CONTRIBUTING and CHANGELOG documents

### Fixed
- Coding errors in ThermoEngine Python package within the coder and equilibrate modules

## [1.0.0-beta-1] - 2020-03-19 (ghiorso)
### Added
- CHANGELOG document
- Experimental Kubernetes cluster on Google Cloud

### Changed
- CONTRIBUTING document  

*standard changelog tags: Added, Changed, Deprecated, Removed, Fixed, Security*
