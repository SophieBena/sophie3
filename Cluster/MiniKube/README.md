# Additional notes on cluster configuration
These are notes that pertain to running a kubernetes cluster on your laptop or desktop machine.  This procedure is useful in testing kubernetes cluster configuration files prior to deployment in a production environment like AWS or GKS. The workflow described here uses MiniKube and deploys a JupyterHub server using a Helm chart. 

The instructions described below rely on Homebrew, which is an installation package for Macintosh OS. If you are using Linux or Windows, then an appropriate substitute for Homebrew must be found. 

## Minikube and preliminaries
Minikube is a "single-user" kubernetes "cluster" that runs on your local machine. Open a terminal window to run the commands shown below. 

First, kubectl, the kubernetes manager, must be installed:
```
brew install kubectl
```
MiniKube may be installed using brew:
```
brew install minikube
```
and requires a hypervisor to run:"
```
brew install hyperkit
```
Start minikube with:
```
minikube start --cpus 4 --memory 8192 --vm-driver hyperkit
```
## Helm
Helm is the package manager for kubernetes. Install as follows:
```
curl https://raw.githubusercontent.com/helm/helm/master/scripts/get-helm-3 | bash
```
## JupyterHub
JupyterHub is the notebook server that will run in your MiniKube. First, add its installation chart to Helm:
```
helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/
helm repo update
```
Release and namespace are both set to jhub when Helm installs JupyterHub on your MiniKube cluster:
```
helm upgrade --cleanup-on-fail \
  --install jhub jupyterhub/jupyterhub \
  --namespace jhub \
  --create-namespace \
  --version=0.10.2 \
  --values config.yaml
```
The *config.yaml* file should contain helm chart values to configure the cluster, including the Docker container image for the SingleUser pods. Helm chart 0.10.2 installs JupyterHub 1.2.1.  The Docker image specified in the config file is the latest image of the ThermoEngine container. 

Check on the status of the cluster:
```
kubectl --namespace=jhub get all
```
Access JupyterHub as:
```
kubectl --namespace=jhub port-forward deployment/proxy 8000
```
as http://127.0.0.1:8000 with username: *enki* and password: *enki*.

## Prometheus monitoring using Grafana
Your MiniKube cluster can be monitored using Prometheus and Grafana. Install as follows: 
```
helm repo add prometheus-community https://prometheus-community.github.io/helm-charts
```
```
helm install prometheus prometheus-community/kube-prometheus-stack
```
Access Grafana in the cluster by setting up port forwarding:
```
kubectl port-forward deployment/prometheus-grafana 3000
```
as http://127.0.0.1:3000 with username: *admin* and password: *prom-operator*.  

### Grafana dashboard
A Grafana dashboard for JupyterHub can be downloaded from:
```
https://github.com/yuvipanda/jupyterhub-grafana
```
dependencies (jsonnet) can be installed with brew; json file for dashboard can be generated in place and upoaded in the Grafana UI.  Follow instructions in the repository README:
```
brew install jsonnet
```
Move to repository top level directory and run the command (in this directory):
```
jsonnet -J vendor jupyterhub.jsonnet > jupyterhub.json
```
Upload the jupyterhub.json dashboard file to Grafana.  

Grafana doesn't populate the 'hub' variable properly by default. You'll need to:
1. Go to your Dashboard
1. Go to settings (gear icon in top right)
1. Select the 'hub' variable
1. Click the 'Update' button
## On GKE
Set up namespace for prometheus:
```
kubectl create ns monitoring
```
Install:
```
helm install prometheus prometheus-community/kube-prometheus-stack --namespace monitoring
```
Check its status by running:
```
kubectl --namespace monitoring get pods -l "release=prometheus"
```
To monitor using Grafana:
- Go to Cloud Platform -> Kubernetes Egine -> Services and Ingress
- Click on prometheus-grafana to show service details
- scroll yo bottom and click on PORT FORWARDING
- Click RUN IN CLOUD SHELL, and return to shell and type return.
- Use the Web review icon at the top right of the shell frame to start web preview; Choose preview on port 8080
- login to the grafana instance (admin, prom-operator)
- when done, contrl-C in Cloud shell to disable forwarding

## Issues
- JupyterHub Helm Chart 0.10.2 does not function with GitLab login, but works standalone
