#!/bin/bash
echo Setup enki jupyterhub workshop cluster using Google Cloud and Kubernetes 
echo Creating standard workshop cluster with two nodes ...
gcloud container clusters create \
  --machine-type n1-standard-2 \
  --num-nodes 2 \
  --zone europe-west1-b \
  --cluster-version latest \
  enkiworkshop
echo Test the cluster ...
kubectl get node
echo Setup workshop cluster administration ...
kubectl create clusterrolebinding workshop-admin-binding \
  --clusterrole=workshop-admin \
  --user=ghiorso@gmail.com
echo Create user node pool ...
gcloud beta container node-pools create user-pool \
  --machine-type n1-standard-2 \
  --num-nodes 0 \
  --enable-autoscaling \
  --min-nodes 0 \
  --max-nodes 6 \
  --node-labels hub.jupyter.org/node-purpose=user \
  --node-taints hub.jupyter.org_dedicated=user:NoSchedule \
  --zone europe-west1-b \
  --cluster enkiworkshop
echo Cluster initialization done.