
#
# ---- builder: install the project into a venv ----
#

FROM debian:bookworm-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Tools needed to install uv
RUN apt-get update
RUN apt-get install -y --no-install-recommends curl tar ca-certificates

# Install uv
ENV HOME=/opt/uv-home
RUN curl -LsSf https://astral.sh/uv/install.sh -o /tmp/uv-install.sh
RUN sh /tmp/uv-install.sh

ENV UV_PROJECT_ENVIRONMENT=/opt/venv
ENV PATH="/opt/uv-home/.local/bin:${PATH}"

# Copy your project and install it into the venv
COPY . /tmp/CFEIntact
RUN uv --project /tmp/CFEIntact sync --no-editable

# Cleanup
RUN rm -rf /tmp /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt /var/lib/dpkg /root/.cache /opt/uv-home/.cache

#
# ---- runtime: only what we need to run ----
#

FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime deps: BLAST + python runtime (so the venv has a compatible interpreter)
RUN apt-get update \
 && apt-get install -y --no-install-recommends ncbi-blast+ ncdu \
 && rm -rf /tmp /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt /var/lib/dpkg /root/.cache /opt/uv-home/.cache \
 && mkdir -p /tmp

# Copy the prebuilt venv from the builder
COPY --from=builder /opt/venv /opt/venv
COPY --from=builder /opt/uv-home /opt/uv-home

# Make sure the venv is used by default
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /w

# Use the venv-installed entrypoint (most reliable)
ENTRYPOINT ["cfeintact"]
