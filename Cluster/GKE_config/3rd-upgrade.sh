#!/bin/bash
echo Upgrade jupyterhub for enki cluster using Helm and config-*.yaml
echo Upgrade to add nbgitpuller, https encryption and gitlab authentication ...
helm upgrade jhub jupyterhub/jupyterhub \
  --values config-upgrade.yaml
