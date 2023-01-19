# Workflow to duplicate the Python package installation of the Docker notebook image container

Use this workflow to mirror the cluster production environment on your Desktop. The commands below are appropriate for Macintosh Intel.

1. In a clean environment install a conda miniforge distribution from:
    ```
    https://github.com/conda-forge/miniforge
    ```
1. From this website, download 
    ```
    Miniforge3-MacOSX-x86_64.sh 
    ```
    a shell script, and execute that script
    ```
    bash Miniforge3-MacOSX-x86_64.sh
    ```
    answering all questions in the affirmative. This script will install Python version 3.9.7.
    The default installation puts Python packages in ~/miniforge3
    ```
    ls ./miniforge3/
    ```
    gives the output:
    ```
    LICENSE.txt	conda-meta	envs		include		pkgs		shell
    bin		condabin	etc		lib		share		ssl
    ```
1. To verify that all went as expected, execute:
    ```
    conda list
    ```
    which gives the output
    ```
    # packages in environment at /Users/ghiorso/miniforge3:
    #
    # Name                    Version                   Build  Channel
    brotlipy                  0.7.0           py39h89e85a6_1001    conda-forge
    ca-certificates           2021.10.8            h033912b_0    conda-forge
    certifi                   2021.10.8        py39h6e9494a_0    conda-forge
    cffi                      1.14.6           py39he338e87_1    conda-forge
    chardet                   4.0.0            py39h6e9494a_1    conda-forge
    charset-normalizer        2.0.0              pyhd8ed1ab_0    conda-forge
    colorama                  0.4.4              pyh9f0ad1d_0    conda-forge
    conda                     4.10.3           py39h6e9494a_2    conda-forge
    conda-package-handling    1.7.3            py39h89e85a6_0    conda-forge
    cryptography              3.4.8            py39ha2c9959_0    conda-forge
    idna                      3.1                pyhd3deb0d_0    conda-forge
    libcxx                    12.0.1               habf9029_0    conda-forge
    libffi                    3.4.2                he49afe7_4    conda-forge
    libzlib                   1.2.11            h9173be1_1013    conda-forge
    ncurses                   6.2                  h2e338ed_4    conda-forge
    openssl                   1.1.1l               h0d85af4_0    conda-forge
    pip                       21.3               pyhd8ed1ab_0    conda-forge
    pycosat                   0.6.3           py39h89e85a6_1006    conda-forge
    pycparser                 2.20               pyh9f0ad1d_2    conda-forge
    pyopenssl                 21.0.0             pyhd8ed1ab_0    conda-forge
    pysocks                   1.7.1            py39h6e9494a_3    conda-forge
    python                    3.9.7           h1248fe1_3_cpython    conda-forge
    python_abi                3.9                      2_cp39    conda-forge
    readline                  8.1                  h05e3726_0    conda-forge
    requests                  2.26.0             pyhd8ed1ab_0    conda-forge
    ruamel_yaml               0.15.80         py39h89e85a6_1004    conda-forge
    setuptools                58.2.0           py39h6e9494a_0    conda-forge
    six                       1.16.0             pyh6c4a22f_0    conda-forge
    sqlite                    3.36.0               h23a322b_2    conda-forge
    tk                        8.6.11               h5dbffcc_1    conda-forge
    tqdm                      4.62.3             pyhd8ed1ab_0    conda-forge
    tzdata                    2021c                he74cb21_0    conda-forge
    urllib3                   1.26.7             pyhd8ed1ab_0    conda-forge
    wheel                     0.37.0             pyhd8ed1ab_1    conda-forge
    xz                        5.2.5                haf1e3a3_1    conda-forge
    yaml                      0.2.5                haf1e3a3_0    conda-forge
    zlib                      1.2.11            h9173be1_1013    conda-forge
    
    ===============================
    ```
    Note that the package source is conda-forge.
1. Next, update miniconda
    ```
    conda update —all
    ```
1. Then, install the **mamba** package manager
    ```
    conda install mamba -n base -c conda-forge
    ```
1. Use **mamba** to install the Jupyter environment
    ```
    mamba install notebook jupyterhub jupyterlab
    ```
1. Add pacakges required for datascience
    ```
    mamba install altair beautifulsoup4 bokeh bottleneck cloudpickle conda-forge::blas=*=openblas cython \
    dask dill h5py ipympl ipywidgets matplotlib-base numba numexpr pandas patsy protobuf pytables scikit-image \
    scikit-learn scipy seaborn sqlalchemy statsmodels sympy widgetsnbextension xlrd
    ```
1. Add packages available from conda-forge that support ThermoEngine
    ```
    mamba install deprecation openpyxl pytest numdifftools cmake elasticsearch elasticsearch-dsl
    ```
1. Move to the root of the ThermoEngine repository and execute
    ```
    make
    make install
    make pyinstall
    ```
    The last make command will install rubicon-objc and thermoengine using pip. Pip will announce that
    ```
    ERROR: pip's dependency resolver does not currently take into account all the packages that are installed.
    This behaviour is the source of the following dependency conflicts.
    numba 0.54.1 requires numpy<1.21,>=1.17, but you have numpy 1.21.4 which is incompatible.
    ``` 
    but this inconsistency does not seem to affect execution of notebooks or tests.
