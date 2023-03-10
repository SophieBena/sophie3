#!/bin/bash
echo Setup helm for enki jupyterhub cluster using Google Cloud and Kubernetes
echo Service account and tiller ...
kubectl create namespace workshop
kubectl --namespace workshop create serviceaccount tiller
echo Give accout appropriate permissions ...
kubectl create clusterrolebinding tiller --clusterrole cluster-admin --serviceaccount=workshop:tiller --namespace workshop
echo Initialize Helm and Tiller ...
helm init --tiller-namespace workshop --service-account tiller --wait
echo Secure Tiller ...
kubectl patch deployment tiller-deploy --namespace=workshop --type=json --patch='[{"op": "add", "path": "/spec/template/spec/containers/0/command", "value": ["/tiller", "--listen=localhost:44134"]}]'
echo Check Helm is functioning ...
helm version --tiller-namespace workshop
echo Helm installation done