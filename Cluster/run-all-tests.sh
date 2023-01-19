#!/bin/bash

echo $1
echo $2
echo "CI_REPOSITORY_URL in run-all-tests.sh:"
echo $4
echo $5
echo "------------"
# docker run -v $2 -e COVERALLS_REPO_TOKEN=$3 CI_REPOSITORY_URL=$4 --user root --rm $1 bash -c '
docker run -v $2 -e COVERALLS_REPO_TOKEN=$3 -e GIT_REPO=$4 -e GIT_BRANCH=$5 --user root --rm $1 bash -c '
echo "CI_REPOSITORY_URL inside docker:" &&
echo $CI_REPOSITORY_URL &&
echo "-------------" &&
pip install nbval pytest-cov wurlitzer coveralls &&
# git clone $CI_REPOSITORY_URL &&
git clone --branch $GIT_BRANCH $GIT_REPO.git &&
echo "==== run devtests ====" &&
cd ThermoEngine &&
export PROJECT_DIR=$(pwd) &&
cd $PROJECT_DIR/thermoengine/thermoengine/test &&
mkdir -p htmlcov &&
pytest --cov=thermoengine --cov-append --cov-report= . &&
echo "===== find .coverage files =====" &&
ls -al &&
coverage html &&
mv htmlcov $PROJECT_DIR/htmlcov-devtests &&
mv .coverage $PROJECT_DIR/.coverage.devtests &&
echo "==== run nbtests ====" &&
cd $PROJECT_DIR/Notebooks &&
mkdir -p htmlcov &&
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
(pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/7a-K-Quartz-Equilibrate.ipynb || true) &&
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
(pytest --cov=thermoengine --cov-append --cov-report= --nbval Codegen/Example-13-Sulfide_liquid.ipynb  || true) &&
coverage html &&
mv htmlcov $PROJECT_DIR/htmlcov-nbtests &&
mv .coverage $PROJECT_DIR/.coverage.nbtests &&
echo "==== combine coverage reports ====" &&
cd $PROJECT_DIR &&
mkdir -p htmlcov &&
ls -al &&
ls -al htmlcov-nbtests &&
ls -al htmlcov-devtests &&
coverage combine &&
ls -al htmlcov &&
coverage html &&
ls -al htmlcov &&
(coveralls || echo "WARNING: Coveralls report failed, likely b/c repo not supported. To get coveralls report, need to setup repo support on coveralls (see https://docs.coveralls.io/supported-ci-services for details).") &&
mv htmlcov-nbtests /mnt &&
mv htmlcov-devtests /mnt &&
mv htmlcov /mnt &&
ls -al /mnt
'
