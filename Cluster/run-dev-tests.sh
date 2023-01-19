#!/bin/bash
echo $1
echo $2
docker run -v $2 -e COVERALLS_REPO_TOKEN=$3 --user root --rm $1 bash -c '
printenv &&
pip install nbval pytest-cov wurlitzer coveralls &&
git clone https://gitlab.com/ENKI-portal/ThermoEngine.git &&
cd ThermoEngine/thermoengine/thermoengine/test &&
mkdir -p htmlcov &&
pytest --cov=thermoengine --cov-append --cov-report= . &&
coverage html &&
mv htmlcov /mnt &&
ls -al /mnt &&
coveralls || echo "WARNING: Coveralls report failed, likely b/c repo not supported. To get coveralls report, need to setup repo support on coveralls (see https://docs.coveralls.io/supported-ci-services for details)."
'
