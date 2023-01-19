# Welcome to VESIcal on the ENKI server!

## What files are where?

### Manuscript

**VESIcal Part I: An open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility
in silicate melts** 
You can read and interact with the VESIcal Part I manuscript as a jupyter notebook. Just open the
manuscript.ipynb file, and start reading and coding. The manuscript is peer-reviewed and published
in Earth and Space Science. You can cite the manuscript as:

Iacovino K, Matthews, S, Wieser PE, Moore GM, and Bégué F (2021) VESIcal Part I: An open-source
thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts. Earth and
Space Science, doi: 10.1029/2020EA001584

### Manuscript supplementary files

All supplementary files and notebooks are included in the /Supplement folder.	

### Tutorials

The /Tutorials folder contains jupyter notebooks and datasets (xlsx files) related to
existing tutorials. Some tutorials have companion videos that live on the VESIcal YouTube.
Others refer to tutorials hosted on VESIcal's ReadTheDocs page and should be self explanatory.
Links to the YouTube videos can also be found on our ReadTheDocs page:
https://vesical.readthedocs.io/en/latest/index.html

## How to cite the VESIcal code

To cite computations done using VESIcal, please cite the Part I manuscript, the VESIcal version
number, as well as the model(s) used. Note that if a model was not specified during calculations,
the default model of MagmaSat was used and should be cited as “MagmaSat, Ghiorso and Gualda
(2015)”. For example: “Calculations were performed using VESIcal (v. 1.0.1; Iacovino et al., 2021)
with the models of Shishkina et al. (2014) and Dixon(1997, “VolatileCalc”).” The web-app always
runs on the most up-to-date version of the VESIcal code, but it is best practice to note if the
web-app was used (“Calculations were performed using the VESIcal web-app (v. 1.0.1; Iacovino et
al., 2021)...”). We also encourage users to be as explicit as possible as to the conditions used
for modelling. This includes stating the pressure, temperature, volatile concentration, and bulk
magma composition used in modelling. In the best case, VESIcal users will provide their code
(e.g., as a Jupyter Notebook or .py file) along with their publication such that it can be easily
replicated.

## More VESIcal resources

- Code repository: https://github.com/kaylai/VESIcal
- Documentation: https://vesical.readthedocs.io/en/latest/index.html
- Manuscript at Earth and Space Science: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001584
- Anvil web-app: https://vesical.anvil.app/
- Snapshot of code repository at time of manuscript publication: https://zenodo.org/record/5095382#.YVnSSUbMLUJ

## Get in touch

If you find issues with the VESIcal code or would like to request a certain feature, please do so
using GitHub's Issues feature. On VESIcal's GitHub page, select Issues in the menu bar toward the
top, and create a new issue. If you are reporting a bug, please provide a minimum working example
of the code you are trying to execute as well as how/where you are running the code (locally, on
ENKI, on anvil?), so that we can best assist in debugging.

If you would like help getting started or with a particular case, feel free to get in touch at
kayla.iacovino@nasa.gov.
