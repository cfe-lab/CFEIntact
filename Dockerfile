
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noniteractive

RUN apt-get update
RUN apt-get install -y mafft ncbi-blast+
RUN apt-get install -y python3 python3-pip
RUN apt-get install -y git

WORKDIR /workspaces/CFEIntact
COPY . .

RUN python3 -m pip install .

ENTRYPOINT ["cfeintact"]
