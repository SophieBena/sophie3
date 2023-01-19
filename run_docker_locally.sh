#!/bin/bash
if [ -z "$1" ]
then 
	echo "Accessing a Jupyter Lab ThermoEngine container rooted to this directory ..."
	docker pull registry.gitlab.com/enki-portal/thermoengine:master
	docker run -p 8888:8888 \
               --env JUPYTER_ENABLE_LAB=yes \
               --user root \
               -e GRANT_SUDO=yes \
               -v $PWD:/home/jovyan/ThermoEngine \
               registry.gitlab.com/enki-portal/thermoengine:master start-notebook.sh
 elif [ "$1" = "term" ]
 then
 	echo "Running an interactive shell into a ThermoEngine container ..."
 	docker pull registry.gitlab.com/enki-portal/thermoengine:master
    docker run --user root \
               -it \
               -v $PWD:/home/jovyan/ThermoEngine \
               --rm registry.gitlab.com/enki-portal/thermoengine:master bash
 elif [ "$1" = "stop" ]
 then
 	echo "Removing running docker containers ..."
 	if [ ! -z "$(docker ps -q)" ]
 	then
 		docker kill $(docker ps -q)
 	fi
 else
 	echo "Usage:"
 	echo "  <no argument> - Run latest ThermoEngine docker container from GitLab."
 	echo "  term - Run a terminal shell in the latest ThermoEngine docker container from GitLab."
 	echo "  stop - Remove any running docker containers from your system."
 	echo "  help - This message."
 fi