# ThermoEngine README
The ThermoEngine repository contains Python packages, C code and header files, and C++ and Objective-C class implementations that allow thermodynamic properties of minerals, fluids, and melts to be estimated. The repository also includes phase equilibrium calculators: generic as well as those that implement the MELTS, pMELTS, MELTS+DEW, and Stixrude-Lithgow-Bertelloni thermodynamic model/data collections. Examples in Jupyter notebooks demonstrate how to call the class methods using Python wrappers.

This README explains how to access a pre-built version of the package, compile the C and Objective-C code, and build the Python wrappers. For API documentation, see the [ThermoEngine documentation](https://enki-portal.gitlab.io/ThermoEngine/).      



For details to build and deploy the ThermoEngine package on Google Cloud and to create a minikube Kubernetes cluster on your local computer, please see the README in the Cluster folder.

For descriptions of Jupyter notebooks that illustrate various uses and capabilities of the ThermoEngine package, please see the READMEs in the Notebooks folder.

If you are interested in contributing to this software project, please consult [How to Contribute to ThermoEngine](CONTRIBUTING.md).

Repository configuration and maintenance, and the relationship of this repository to others in the ENKI-portal group, are described in [Maintenance Notes](MAINTENANCE.md).

This software is released under the GNU Affero General Public License ([GNU AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html)).

## Contents  
- [Accessing a pre-built implementation of this software package](#Accessing-a-pre-built-implementation-of-this-software-package)
  * [Cloud server](#Cloud-server)  
  * [Running a container image locally](#Running-a-container-image-locally)  

- [Building and installing the code and testers](#building-and-installing-the-code-and-testers)
  * [Prerequisites](#prerequisites)  
  * [Compilation](#compilation)  
  * [Installation](#installation)   
  * [Xcode builds](#xcode-builds)   

- [Obtaining a Docker image](#obtaining-a-docker-image)  
- [Building the Python package](#building-the-python-package)  
- [Finding examples (Jupyter notebooks)](#finding-examples-jupyter-notebooks)  

## Accessing a pre-built implementation of this software package

For  a video presentation of the procedures described below, see [Accessing the ENKI software ecosystem]( https://www.youtube.com/watch?v=jUshguhFpjk) on our YouTube channel.
<center><div class="videowrapper">
  <iframe  src="https://www.youtube.com/embed/jUshguhFpjk" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe></div></center>

### Cloud server
The ThermoEngine package is accessible from a Google Cloud JupyterHub server. Point your browser to [ENKI-portal](https://server.enki-portal.org/hub/login), and log in with your GitLab credentials. If you lack credentials, [register for a free GitLab account](#https://gitlab.com/users/sign_up) to obtain a username and password. 

### Running a container image locally
The ThermoEngine package is available as a pre-built [Docker container image](#docker-image). You can run this image locally to use ThermoEngine software, to execute Jupyter notebooks, and to develop new applications.  All this can be accomplished with minimal system configuration. 

#### Setup
1. If Docker is not already installed on your system, download and install [Docker Desktop](#https://www.docker.com/get-started).

2. Make sure that [Git](#https://git-scm.com) is installed on your system. (Most likely, Git is already installed. Macs and Linux machines come preinstalled with Git.) 

2. Clone the ThermoEngine repository to a convenient location on your computer: Open a terminal window, navigate to a location to store the package, and execute the command:
    ```
    git clone https://gitlab.com/ENKI-portal/ThermoEngine.git
    ```

#### Run the image  
1. Change directory down into the root of the repository (e.g., on Linux and Mac):
    ```
    cd ThermoEngine
    ```
2. At this file location, run the shell script:
    ```
    ./run_docker_locally.sh
    ```
    The first time you execute this script, the complete Docker image is fetched from the server. This might take some time. Subsequent invocations of the script fetch updates only and should complete much faster. The script provides a URL that you paste into a browser. This URL gives access to the local Docker container running JupyterLab with an installed ThermoEngine package and its Python/dynamic library dependencies. The container will have read/write access to the files in the ThermoEngine root and directory tree.  

    The shell script referenced above accepts additional arguments.  Pass "help" as an argument to see options.  

#### Keep the repository updated
* Update the cloned version of this repository on a regular basis. The following command, executed from the root repository directory, will download changed files and keep your local version in sync with the container image: 
    ```
    git pull
    ```

## Building and installing the code and testers
These instructions pertain to a Macintosh system with Xcode installed. For an Ubuntu-based Linux system, follow the workflow in the Docker file included in the Cluster folder of this project.

### Prerequisites
The build requires that the `gsl` library is installed.  You can obtain this library from Homebrew (Macintosh) or download it directly from [GNU](https://www.gnu.org/software/gsl/). The build assumes that the `gsl` library is installed into **/usr/local**.

>NOTE: The gsl library downloaded using pip or conda does not function (perhaps because they place the library in the wrong location). Homebrew places the gsl library in the correct location.

### Compilation
You can compile the Objective-C library in the root project directory by executing the command:
```
make
```
### Installation
You can install the library into **/usr/local/lib** (where Python can locate it) with the command:
```
make install
```
You can install the **thermoengine** Python interface package with the command:
```
make pyinstall
```

### Testing and code coverage
You can run tests with the command:

```
make tests
```

### Xcode builds
You can also use Xcode to build the package. Launch Xcode, and open the ThermoEngine.workspace project. Schema are provided to build dynamic libraries and test routines.  

## Obtaining a Docker image
A fully configured and working executable of the project software is available as a Docker image. The image is named `registry.gitlab.com/enki-portal/thermoengine:master`.
Please consult the README in the Cluster folder for instructions on deploying this image.

## Building the Python package
Accessing the Objective-C library from Python requires installing the Rubicon-ObjC package. You can easily install this package if you have the Anaconda (conda) software package manager.  

If Anaconda is installed and if the root environment is based upon Python 2.7, you must create a new environment based on Python 3.7+ in order to install the ThermoEngine package. For example, execute these commands from a terminal window to create a new environment (python38) and to activate that Python 3.8 environment:
```
conda create -n python38 python=3.8 anaconda
source activate python38
```
Any Python 3.7+ environment works with Rubicon. Simply substitute the correct version numbers in the above commands.
  If you want to switch back to the original Python 2.7 environment, execute this command:
```
source activate root
```
Execute the following to install the latest version of Rubicon-ObjC 2.10 from [pypi](https://pypi.python.org/pypi):
```
pip install rubicon-objc==0.2.10
```
 Do not install a later version of Rubicon-Objc, as these are incompatible with the ThermoEngine Python package. Once installed, the package is available to Jupyter notebooks running the Python 3.7+ kernel on your local machine.  

 If you have a traditional Python 3.7+ installation, the above `pip` command will also install Rubicon, but additional Python dependencies may require manual installation.

## Finding examples (Jupyter notebooks)  
The Notebooks folder contains numerous Jupyter notebooks that demonstrate how to call the class methods using Python wrappers. For more details, see the [ThermoEngine documentation](https://enki-portal.gitlab.io/ThermoEngine/examples.html).
