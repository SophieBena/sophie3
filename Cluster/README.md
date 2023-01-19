# ENKI Cluster README

## Contents
- [Automatic cluster updates (CI/CD actions)](#automatic-cluster-updates-(ci/cd-actions))  
- [Generating and configuring the cluster on Google Cloud](#generating-and-configuring-the-cluster-on-google-cloud)  

  * [Preliminaries](#preliminaries)   
  * [1️⃣ Create the cluster](#1️⃣-create-the-cluster)  
  * [2️⃣ Install JupyterHub using Helm](#2️⃣-install-jupyterhub-using-helm)  
  * [3️⃣ Update DNS, allow GitLab access, and reconfigure JupyterHub for public access](#3️⃣-update-dns-allow-gitlab-access-and-reconfigure-jupyterhub-for-public-access)  
  * [4️⃣ Integrate the cluster into the GitLab repository](#4️⃣-integrate-the-cluster-into-the-gitlab-repository)  
  * [5️⃣ Install Knative](#5️⃣-install-knative)  
  * [Tearing down the cluster](#tearing-down-the-cluster)

- [Running the Docker image locally](#running-the-docker-image-locally)  
- [JupyterHub cluster implementation reference](#jupyterhub-cluster-implementation-reference)

## Automatic cluster updates (CI/CD actions)
The following CI/CD actions update the ENKI server on Google Cloud as a consequence of commits to the repository. For these actions to succeed, the cluster must exist, be properly configured, and be integrated with the repository. See [Generating and configuring the cluster on Google Cloud](#generating-and-configuring-the-cluster-on-google-cloud).

- The Docker container image used by the ENKI server is generated on a commit to the master or develop branch of this project using CI commands in .gitlab-ci.yml. It is built from the Cluster/Dockerfile.  The generated  image is publicly accessible and can be run locally, as described in [Running the Docker image locally](#running-the-docker-image-locally). The image *enki-portal/thermoengine:master* is automatically loaded on every compute pod instantiated by the cloud server.
- The GitLab CI/CD pipeline process automatically generates a staging image (tagged as *staging*), and retagging this image to the release image (tagged as *master*) must be performed manually. Visit the pipeline dashboard to execute this process.
- The GitLab CI/CD pipeline also automatically generates a minimalist image (tagged as *runtime*). This image lacks the Jupyter interface tools of the *staging* or *master* image, but is a suitable base for building Python scripts such as Knative web service functions and apps.
- The GitLab CI/CD pipeline builds application Docker container images derived from *runtime* that implement Knative server functions and apps. See the *Knative* subdirectory of this folder. These container images can be deployed to the production server by manually executing the release stage on the pipeline dashboard.

  See Settings > Packages > Container Registry for location and inventory of generated Docker images. 

## Generating and configuring the cluster on Google Cloud
You can configure the cluster from your desktop computer utilizing a cloned version of this repository, which the following steps assume is also named ThermoEngine. You can also accomplish the configuration by copying referenced yaml files from this repository to the Google Cloud shell.

The following steps assume that you have some familiarity with Google Cloud Engine.  

### Preliminaries
1. If you are installing from the desktop, open a terminal window and do the following:
    1. cd into the Cluster/GKE_config directory of the repository.
    1. **(Once only)** [Install Google Cloud resources](https://cloud.google.com/sdk/docs/install) ("gcloud") onto the desktop machine.
    1. **(Once only)** Install *kubectl*, the command line interface to Kubernetes:  
        ```
        gcloud components install kubectl
        ```
1. On either the desktop or Google Cloud shell, install Helm as follows:
    1.  **(Once only)** Install Helm:
        ```
        curl https://raw.githubusercontent.com/helm/helm/master/scripts/get-helm-3 | bash
        ```
    1.  Verify the Helm installation:
        ```
        helm list
        ```

### 1️⃣ Create the cluster
1. Set up the ENKI JupyterHub cluster using Google Cloud and Kubernetes:
    ```
    gcloud container clusters create \
      --machine-type n1-standard-2 \
      --enable-autoscaling \
      --min-nodes=2 \
      --max-nodes=4 \
      --zone us-west1-a \
      --cluster-version latest \
      enkiserver
    ```
2. Test that the cluster exists:
    ```
    kubectl get node
    ```
3. Set up cluster administration (Substitute administrator email address as appropriate):
    ```
    kubectl create clusterrolebinding cluster-admin-binding \
      --clusterrole=cluster-admin \
      --user=administrator@example.com
    ```
4. Create a user node pool that is expandable to six nodes:
    ```
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
    ```

### 2️⃣ Install JupyterHub using Helm
1. Install an initial configuration of JupyterHub on the cluster using Helm and the provided *config.yaml* file: 

    1. Note that the secret token defined in this file must appear in all yaml files used to build and upgrade the cluster. You can generate this token using:  
       ```openssl rand -hex 32```
    1. Get the Jupyterhub Helm chart and update:
       ```
       helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/
       helm repo update
       ```
    1. Build the JupyterHub cluster according to the config.yaml file:
       ```
       helm upgrade --cleanup-on-fail \
         --install jhub jupyterhub/jupyterhub \
         --namespace jhub  \
         --create-namespace \
         --version=0.11.1 \
         --values config.yaml
       ```
       The cluster *namespace* and *install* release are both named *jhub*, but you can specify another name. The JupyterHub version number is the latest stable version as of 04/2021. See the [Helm chart repository](https://jupyterhub.github.io/helm-chart/) for later versions.  

1. Check that the pods are in running mode:
    ```
    kubectl get pod --namespace jhub
    ```
1. Find the IP address of the proxy server (Execute this command repeatedly until one is assigned):
    ```
    kubectl get service --namespace jhub
    ```
    At this stage the cluster permits non-secure login for a single user (*enki*) who has admin privileges, with password *ENKI070719*

### 3️⃣ Update DNS, allow GitLab access, and reconfigure JupyterHub for public access
1. Update domain DNS (at your domain service provider) so that an *a-record* for your domain points to the generated IP address for the proxy server.
1. **(Once only)** At your GitLab user profile web portal, generate an oauth token for a new application (enkiserver) that allows GitLab credentials to be used for cluster login.
1. Edit *config-encrypt.yaml* to reflect oauth token credentials, your domain name, and the correct image TAG.
1. Make sure that the secret token used in *config-encrypt.yaml* is the same as the secret token used in *config.yaml* as discussed above.
1. Upgrade the cluster to add https access, to enable read-only pulls of the ThermoEngine repository, and to implement GitLab authentication for login:
    ```
    helm upgrade jhub jupyterhub/jupyterhub --values config-encrypt.yaml
    ```

    The Helm commands reference yaml files in the Cluster directory. Specify alternate JupyterHub configuration options by editing these files. See the [Zero to JupyterHub documentation](https://zero-to-jupyterhub.readthedocs.io/en/latest/).

### 4️⃣ Integrate the cluster into the GitLab repository

1. On the web portal to your repository at GitLab.com, navigate to the **Operations > Kubernetes** page, and click **Integrate with a cluster certificate.**  
1. Select the **Connect existing cluster** tab.  

1. Fill in your cluster name, follow the instructions to fill in your **CA Certificate** and **Service Token**, and make sure that **GitLab-managed cluster** is checked.  Click **Add Kubernetes cluster.**  
1. Once the cluster is integrated, do the following:  
    1. Select the **Applications** tab.  
    1. Install *Prometheus* and *GitLab Runner* as GitLab-managed cluster apps.  
    1. Optionally, install *elasticsearch* to facilitate analysis of cluster logs.  

### 5️⃣ Install Knative
Knative is used to distribute ENKI web services.  Do not use GitLab to install the managed app version of Knative. Instead, follow this procedure to install a Knative server from the desktop or from Google Cloud shell:
1. Ensure that your default node pool has autoscale turned on with a minimum of two nodes and a maximum of at least four nodes.
1. Ensure that your cluster nodes are running the latest stable version of GKE.
1. Install Knative Serving (The commands below install v0.21) with the Istio networking layer:
    ```
    kubectl apply --filename https://github.com/knative/serving/releases/download/v0.21.0/serving-crds.yaml
    kubectl apply --filename https://github.com/knative/serving/releases/download/v0.21.0/serving-core.yaml
    kubectl apply --filename https://github.com/knative/net-istio/releases/download/v0.21.0/istio.yaml
    kubectl apply --filename https://github.com/knative/net-istio/releases/download/v0.21.0/net-istio.yaml
    ```
1. Fetch the IP address of the external Istio gateway:
    ```
    kubectl --namespace istio-system get service istio-ingressgateway
    ```
1. Configure this address with your DNS provider using a wildcard CNAME entry:
    ```
    # Here knative.example.com is the domain suffix for your cluster
    *.knative.example.com == A 35.233.41.212
    ```
1. Once your DNS is configured, direct Knative to use that domain:
    ```
    # Replace knative.example.com with your domain suffix
    kubectl patch configmap/config-domain \
      --namespace knative-serving \
      --type merge \
      --patch '{"data":{"knative.example.com":""}}'
    ```
   Knative Serving is now installed.  
1. Monitor components using:
    ```
    kubectl get pods --namespace knative-serving
    ```
Optional Knative installations (e.g., TLS certificates, Knative Eventing) are documented at https://knative.dev/docs/install/any-kubernetes-cluster/.

### Tearing down the cluster
To tear down the cluster:
1. At the gitlab.com repository portal, navigate to the **Operations > Kubernetes** page, select the cluster, expand the **Advanced Settings** section, and execute *Remove integration and resources*.  
2. Once the integration and associated resources are removed from GitLab, go to a desktop terminal window or Google Cloud shell, and run the commands:
    ```
    helm delete jhub --purge
    kubectl delete namespace jhub
    gcloud container clusters delete enkiserver --zone=us-west1-a
    ```

## Running the Docker image locally   

The following steps assume that the latest automatically generated Docker image at ENKI-portal/ThermoEngine is being accessed. If you prefer not to use the automatically generated image, build the Docker image and upload it to your GitLab repository.
1. Make sure that Docker is installed on your machine.
1. At a terminal on your local machine, log in using GitLab credentials.
    ```
    docker login registry.gitlab.com
    ``` 
1. At a terminal on your local machine, run the command:  
    ```
    docker run -p 8888:8888 --env JUPYTER_ENABLE_LAB=yes --user root -e GRANT_SUDO=yes registry.gitlab.com/enki-portal/thermoengine:latest start-notebook.sh
    ```  
    This command launches JupyterLab in the container. The user has root permission with no password asked for sudo.  
1. Access the notebook in a browser pointed to the URL printed to your terminal, such as:  
    ``` http://127.0.0.1:8888/?token=5193b5714271b44006edb619fb26cfd203cd113a998eff8b```
1. If the container doesn't exit once you exit JupyerLab, kill the executing container with the terminal command:  
    ```
    docker kill $(docker ps -q)
    ```
See notes on other options for running the image at https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html.

You can obtain an interactive shell into the container with:  
```
docker run --user root -it --rm registry.gitlab.com/enki-portal/thermoengine:latest bash
```
    
# JupyterHub cluster implementation reference
The [JupyterHub cluster implementation reference](https://zero-to-jupyterhub.readthedocs.io/en/latest/) is available with the *Zero to JupyterHub with Kubernetes* documentation.