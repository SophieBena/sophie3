# Notes on building Knative functions 
- Builds are managed by gitlab-ci.yml file in repository root
- Dockerfile for minimal ThermoEngine container
    + Built from Minimal Jupyter notebook Dockerfile (without notebook)
    + And, non-interactive components of the principal Dockerfile
