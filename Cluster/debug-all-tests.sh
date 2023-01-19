#!/bin/bash

echo $1
echo $2
docker run -v $2 -e COVERALLS_REPO_TOKEN=$3 --user root --rm $1 bash -c '
pip install nbval pytest-cov wurlitzer coveralls &&
git clone https://gitlab.com/ENKI-portal/ThermoEngine.git &&
echo "==== run devtests ====" &&
cd ThermoEngine &&
export PROJECT_DIR=$(pwd) &&
cd $PROJECT_DIR/thermoengine/thermoengine/test &&
mkdir -p htmlcov &&
pytest --cov=thermoengine --cov-append --cov-report= test_samples.py &&
echo "===== find .coverage files =====" &&
ls -al &&
coverage html &&
mv htmlcov $PROJECT_DIR/htmlcov-devtests &&
mv .coverage $PROJECT_DIR/.coverage.devtests &&
echo "==== run nbtests ====" &&
cd $PROJECT_DIR/Notebooks &&
mkdir -p htmlcov &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Compare-Phases.ipynb &&
(pytest --cov=thermoengine --cov-append --cov-report= --nbval Equilibrate/7a-K-Quartz-Equilibrate.ipynb || true) &&
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