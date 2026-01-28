
#
# ---- builder: install the project into a venv ----
#

FROM debian:bookworm-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Python + tools needed to install your package into a venv
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      python3 python3-pip python3-venv

# Create a self-contained venv for runtime
RUN python3 -m venv /opt/venv

# Make sure the venv is used by default
ENV PATH="/opt/venv/bin:${PATH}"

RUN pip install --upgrade pip setuptools wheel

# Copy your project and install it into the venv
COPY . /tmp/CFEIntact
RUN pip install /tmp/CFEIntact

RUN rm -rf /tmp/CFEIntact /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt /var/lib/dpkg /root/.cache

#
# ---- runtime: only what we need to run ----
#

FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime deps: BLAST + python runtime (so the venv has a compatible interpreter)
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ncbi-blast+ \
      python3 \
 && rm -rf /tmp/CFEIntact /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt /var/lib/dpkg /root/.cache

# Copy the prebuilt venv from the builder
COPY --from=builder /opt/venv /opt/venv

# Make sure the venv is used by default
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /w

ENTRYPOINT ["cfeintact"]
