
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y --no-install-recommends ncbi-blast+
RUN apt-get install -y --no-install-recommends python3 python3-pip

COPY . /tmp/CFEIntact

RUN python3 -m pip install --break-system-packages /tmp/CFEIntact

RUN rm -rf /tmp/CFEIntact
RUN rm -rf /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt

WORKDIR /w

ENTRYPOINT ["cfeintact"]
