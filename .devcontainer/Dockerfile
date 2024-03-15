
# Use the official image as a parent image
FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noniteractive

# Set the working directory in the container
WORKDIR /workspaces/cfeintact

RUN apt-get update
RUN    apt-get install -y python3 python3-pip git-all \
    && python3 -m pip install --upgrade pip

# Copy the rest of your app's source code from your host to your image filesystem.
COPY . /workspaces/cfeintact
