# Maintenance Notes
This document provides maintenance notes for the repository. The content is useful for owners and maintainers of the repository and supported services. It is not directly relevant to code contributors or general users.

## Contents

- [Git Management](#git-management)
- [Associated repositories](#associated-repositories)
- [Associated accounts](#associated-accounts)
- [Repository structure](#respository-structure)
- [Continuous Integration and Deployment (CI/CD)](#continuous-integration-and-deployment-cicd)
- [Settings](#settings)

## Git management
### General
- The use of git (and GitLab-specific) features are integral to the management of this project
- The forking workflow employed by the ENKI-project relies upon merge requests from forked repos to push the project forward. As detailed in contributing.md document, each developer forks the main repo, does work, and then posts a merge request to merge changes upstream into the main ENKI repository.
### Merge Revert info for maintainers
- For maintainers, the forking workflow requires review and acceptance of merge requests. If a merge request is accepted prematurely or is found to cause a major problem, the merge can be reverted (though this has consequences, detailed below, which are tolerable but best avoided whenever possible).
- In order for any material from this merge to be later accepted and incorporated into the repo (after fixing issues that prompted the revert), maintainers must first revert the merge-revert. This revert-revert request is itself issued as a separate merge request. 
- **Reverting the revert is critical, otherwise git would dutifully ignore all material covered by the first revert.**. 
- Code covered by a merge-revert is forever ear-marked to be ignored by the repo, until the revert-revert command is given at GitLab. This is a non-destructive method general to git that yields a somewhat messy but accurate repo history, guaranteeing that no history of the event is lost.
- see example resolving this issue over period from Oct 7-11, 2021, commits:
  - Original Merge: eaa1c8fa7c57932725a0849fef5d3819745f4ed6
  - Revert Merge: 2fe1bc839e62c63e81cb55920744ca95d5765cec
  - Revert-revert Merge: 55ec9d289cf7246c8b9ffe4f3e63bc62e02dd870

## Associated repositories
- ENKI-portal/**ThermoEngine**
  - Main code repository (this repository)
- ENKI-portal/**jupyterlab-gitlab**
  - JupyterLab extension that provides access to the ENKI-portal group public repositories
- ENKI-portal/**jupyterlab_enkiintro**
  - JupyterLab extension that provides a splash screen on login to the ENKI cluster server
- ENKI-portal/**jupyterlab_shared**
  - JupyterLab extension that provides a tab for access to ThermoEngine notebooks
- ENKI-portal/**jupyterhub_custom**
  - JupyterHub template extensions that provide a login screen for the ENKI cluster server
- ENKI-portal/**ENKI-cluster-config**
  - Google Cloud cluster configuration and deployment, automatically triggered by tagged uploads to ThermoEngine
  - Deprecated
- ENKI-portal/**Theoretical_Geochemistry**
  - Textbook on theoretical geochemistry authored by H.C. Helgeson
- ENKI-portal/**SulfLiq**
  - Python package that implements the thermodynamic model of Victor Kress for liquids in the system O-S-Fe-Ni-Cu

## Associated accounts
Accounts on the following platforms are required:
- Google Cloud
  - A Google Cloud account must be created at https://console.cloud.google.com/ and associated with a mechanism for billing and a project name for cluster creation. It is used by the shell scripts in **ThermoEngine**/Cluster and by the CI scripts in **ENKI-cluster-config** for interaction of GitLab with a Google Cloud Kubernetes cluster.
- npm
  - An npm account is required to deploy GitLab extensions that are uploaded to https://www.npmjs.com.  An access token is generated from this account on the website, and that token is defined as an environment variable in **Settings->CI/CD/Variables** on GitLab.com for the repository that performs an upload during the deploy phase of the CI/CD.  
    - The token name is `NPM_TOKEN`.
    - The token should be available to all environments and masked from display in the pipeline scripts.
  - Used by **jupyterlab-gitlab**, **jupyterlab_enkiintro**, and **jupyterlab_shared**.
- PyPi
  - A PyPi account is required to deploy Python packages to PyPi (https://pypi.org) using twine. It is used by the CI/CD scripts to deploy Python packages.  Two tokens should be defined as environment variables in the repository at **Settings->CI/CD/Variable**.  
    - The token names are `TWINE_PASSWORD` and `TWINE_USERNAME`.
    - These tokens are generated at the PyPi website by creating a token for login. The username is __token__ and the password is the generated token value.
  - Used by **jupyterlab-gitlab**.
- Coveralls
  - A Coveralls account is required to perform a code coverage analysis of a repository. The account is created on https://coveralls.io, and the repository is linked at that site.  Once those tasks are completed, links are generated that describe coverage.
  - Use your GitLab account to log in/register.
  - Once a repository is selected for coverage at the Coveralls site, a token will be generated. This token is used to upload coverage data from the GitLab repository using CI/CD. The token value should be set as a masked CI/CD environment variable in the relevant repository on GitLab (Settings -> CI/CD -> Variables) using the name `COVERALLS_REPO_TOKEN`.
  - Used by **ThermoEngine**.
- Zenodo
  - A Zenodo account (https://zenodo.org) is required to generate DOIs that can be used to tag software in a repository.
  - Used by **ThermoEngine** and **Theoretical_Geochemistry**.
- ReadTheDocs
  - A ReadTheDocs account (https://readthedocs.org) is required for automated builds of documentation.
  - Used by **ThermoEngine**.

## Respository structure
- Two required branches, master and Documentation, both protected
- optional dev branch (key for forked personal repos and empowered with CI/CD)
- README file
- GNU AGPLv3 public license; source must be included and exposed in any subsequent use of the code; otherwise, code is open source.
- CHANGELOG documents every major push to master.
- CONTRIBUTING describes the workflow for contributing to the code base.
- Issues
  - Service Desk is turned on so that issues may be emailed to the repository.
  - Issue templates are stored in .gitlab/issue_templates.
- Packages -> Container Registry is the GitLab container registry that holds Docker images for Kubernetes cluster or standalone deployment.
- Wiki is a wiki for the repository.
- Google Cloud Cluster linked through Operations -> Kubernetes
  - Production server deploying software from this repository.
  - Created and maintained externally on the Google Cloud platform and integrated into the repository.
  - Prometheus installed as a GitLab-managed application to enable metrics and logs to be gathered and visualized by GitLab.

## Continuous Integration and Deployment (CI/CD)
- CI/CD integration is scripted in the *.gitlab-ci.yml* file. There are four stages:
  - build (executed on push to master/develop)
    - Create Docker container image on upload to the master/develop branch; uses caching; pushes to GitLab container registry associated with the repository; tagged staging.
    - Create runtime (no Jupyter interface; only libraries installed) Docker container image on upload to the master/develop branch; uses caching; pushes to GitLab container registry associated with the repository; tagged runtime.
  - test (executed on push to master/develop)
    - Executes tests on the staging container image using pytest called from a shell script called `test-image.sh`.
    - Test results are used to assess code coverage.
    - Executes image builds for Knative web service apps based on the runtime image; tagged as latest or as denoting the function/app being implemented.
  - release (manually executed on a push to master)
    - Copies the staging image generated under the build stage to a version tagged master, and pushes that image to the container repository.
    - Installs Knative web services onto the attached Kubernetes cluster.
  - deploy
    - Pages (executed on push to Documentation branch)
    - Copies Sphinx-generated documentation into the Pages website, which is exposed at https://enki-portal.gitlab.io/.
  

## Settings
- General
    - Visibility, project features, permissions
      - Public
    - Badges
      - Set for status of ReadTheDocs upload (see below)
      - Pipeline status
      - Code coverage
      - Zenodo DOI
    - Service Desk
      - Set to On
      - Note email address for issues
- Members
  - Special members
    - Owners: Wolf and Ghiorso
    - Maintainers: Johnson (for Documentation branch access)
  - ENKI-portal group
    - These members: Ghiorso, Wolf, Johnson
- Webhooks
  - A ReadTheDocs hook is defined so that on pushing the Documentation branch, ReadTheDocs generates Sphinx-based documentation.
- Repository
  - Protected branches
    - Documentation is set up as a protected branch (master is by default). Only maintainers or owners can merge into or push a protected branch.
- CI/CD
  - Variables
    - Environment variables are set here for CI/CD builds.
- Operations
  - Kubernetes
    - Integration of external Google Cloud Kubernetes cluster
    - GitLab managed apps: Prometheus, GitLab runner
- Pages
  - Pages are served (using SSL certificates). These html pages serve a copy of locally generated documentation, which mirrors documentation automatically generated via a web hook to ReadThe Docs.
