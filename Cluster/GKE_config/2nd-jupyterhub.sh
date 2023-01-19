#!/bin/bash
echo Setup jupyterhub for enki cluster using Helm and config.yaml
echo Get the Jupyterhub Helm chart and update ...
helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/
helm repo update
echo Now build the JupyterHub cluster according to the config.yaml file ...
RELEASE=jhub
NAMESPACE=jhub
helm upgrade --cleanup-on-fail \
  --install $RELEASE jupyterhub/jupyterhub \
  --namespace $NAMESPACE  \
  --create-namespace \
  --version=0.11.1 \
  --values config.yaml
echo Pods should be in running mode ...
kubectl get pod --namespace jhub
echo Find the IP address of the cluster ...
kubectl get service --namespace jhub
echo JupyterHub is running