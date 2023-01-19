#!/bin/bash
echo Setup enki jupyterhub cluster using Google Cloud and Kubernetes 
echo Creating standard cluster with two nodes ...
gcloud container clusters create \
  --machine-type n1-standard-2 \
  --enable-autoscaling \
  --min-nodes=2 \
  --max-nodes=4 \
  --zone us-west1-a \
  --cluster-version latest \
  enkiserver
echo Test the cluster ...
kubectl get node
echo Setup cluster administration ...
kubectl create clusterrolebinding cluster-admin-binding \
  --clusterrole=cluster-admin \
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
  --zone us-west1-a \
  --cluster enkiserver
echo Cluster initialization done.