#!/bin/bash
echo Setup jupyterhub for enki workshop cluster using Helm and config.yaml
echo Get the Jupyterhub Helm chart and update ...
helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/ --tiller-namespace workshop 
helm update --tiller-namespace workshop 
echo Now build the JupyterHub cluster according to the config.yaml file ...
RELEASE=jwork
NAMESPACE=jwork
helm upgrade --install $RELEASE jupyterhub/jupyterhub \
  --namespace $NAMESPACE  \
  --version=0.8.2 \
  --values config.yaml \
  --tiller-namespace workshop 
echo Pods should be in running mode ...
kubectl get pod --namespace jwork
echo Find the IP address of the cluster ...
kubectl get service --namespace jwork
echo JupyterHub is running