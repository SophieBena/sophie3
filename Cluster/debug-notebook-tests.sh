#!/bin/bash

echo $1
echo $2
docker run -v $2 -e COVERALLS_REPO_TOKEN=$3 --user root --rm $1 bash -c '
pip install nbval pytest-cov wurlitzer coveralls &&
git clone https://gitlab.com/ENKI-portal/ThermoEngine.git &&
cd ThermoEngine/Notebooks &&
mkdir -p htmlcov &&
pytest --cov=thermoengine --cov-append --cov-report= --nbval Pure-Phases/Compare-Phases.ipynb &&
coverage html &&
mv htmlcov /mnt
'