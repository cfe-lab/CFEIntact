
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y mafft ncbi-blast+
RUN apt-get install -y python3 python3-pip

COPY . /tmp/CFEIntact

RUN python3 -m pip install /tmp/CFEIntact

WORKDIR /w

ENTRYPOINT ["cfeintact"]
