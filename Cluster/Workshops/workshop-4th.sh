#!/bin/bash
echo Upgrade jupyterhub for enki workshop cluster using Helm and config-*.yaml
echo Upgrade to add nbgitpuller, https encryption and gitlab authentication ...
helm upgrade jwork jupyterhub/jupyterhub \
  --version=0.8.2 \
  --values config-workshop.yaml \
  --tiller-namespace workshop
