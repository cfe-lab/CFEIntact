---
title: Installing CFEIntact
---

To streamline the installation process and avoid potential software conflicts, CFEIntact is distributed as a Docker container. Docker packages all required components inside a container, ensuring consistency across different systems.

## Install Docker

- **Windows / macOS**:

  1. Visit the [Docker Desktop download page](https://www.docker.com/products/docker-desktop/) and download the installer for your operating system.
  2. Follow the installation instructions provided on Docker's website.
  3. Once installed, launch Docker Desktop and ensure it is running.

- **Linux**:

  1. Follow the [Docker Linux installation guide](https://docs.docker.com/engine/install/) for your distribution (for example, Ubuntu or CentOS).

  2. Start the Docker service using:

     ```shell
     sudo systemctl start docker
     ```

  3. Optionally, enable Docker to start at boot:

     ```shell
     sudo systemctl enable docker
     ```

## Pull the Docker Image for CFEIntact

With Docker installed and running, you can now download the pre-configured Docker image for CFEIntact:

- Open your terminal or command prompt and execute the following command:

  ```shell
  docker pull cfelab/cfeintact
  ```

  This command downloads the latest version of the CFEIntact image from Docker Hub.

- To verify the installation, run:

  ```shell
  docker run --rm cfelab/cfeintact version
  ```

  You should see an output similar to:

  ```
  1.23.1
  ```

---

Next: [data preparation](data_prep.html).
