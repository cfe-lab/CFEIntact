
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y ncbi-blast+
RUN apt-get install -y python3 python3-pip

COPY . /tmp/CFEIntact

RUN python3 -m pip install --break-system-packages /tmp/CFEIntact

RUN rm -rf /tmp/CFEIntact

WORKDIR /w

ENTRYPOINT ["cfeintact"]
