#!/bin/bash
echo $1
echo $2
docker run -v $2 -e COVERALLS_REPO_TOKEN=$3 --user root --rm $1 bash -c '
printenv &&
pip install nbval pytest-cov wurlitzer coveralls &&
git clone https://gitlab.com/ENKI-portal/ThermoEngine.git &&
cd ThermoEngine/Notebooks &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Compare-Phases.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Forsterite-Stixrude.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Phase-Diagram.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Plot-Reaction.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Quartz-Berman.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Water.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Water+For-EH.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Solutions/Feldspar-ss-Berman.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/1-G-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/2-KandC-O2-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/3-K-H2O-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/4-H-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/5-A-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/6-E-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/7a-K-Quartz-Equilibrate.ipynb || true &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/7b-K-Cor-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/8-K-Qtz-Cor-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/9-K-Qtz-Fld-Equilibrate.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/9-Olivine-loop.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval MELTS-pMELTS/MELTS-v1.0.2-equilibrium.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval MELTS-pMELTS/MELTS-v1.0.2-fractionation.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval MELTS-pMELTS/MELTS-v1.1.0-equilibrium.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval MELTS-pMELTS/MELTS-v1.2.0-equilibrium.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval MELTS-pMELTS/pMELTS-v5.6.1-melting.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval MELTS-pMELTS/pMELTS-v5.6.1-adiabatic.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-1-Berman-std-state.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-3-HKF.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-4-Stixrude-Debye.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-7-Simple-Solution.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-9-Complex-Solution.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-10-gas-Speciation-Solution.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-12-SWIM.ipynb &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-13-Sulfide_liquid.ipynb  || true &&
coverage html &&
mv htmlcov /mnt &&
ls -al /mnt &&
coveralls || echo "WARNING: Coveralls report failed, likely b/c repo not supported. To get coveralls report, need to setup repo support on coveralls (see https://docs.coveralls.io/supported-ci-services for details)."
'
