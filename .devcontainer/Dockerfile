
# Use the official image as a parent image
FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noniteractive

RUN apt-get update
RUN apt-get install -y python3 python3-pip git-all mafft ncbi-blast+
RUN python3 -m pip install pytest
