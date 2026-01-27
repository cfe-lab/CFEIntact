#! /bin/sh

set -xe

# Install necessary system packages.
apt-get install -y git curl tar

# Install python's uv package manager.
curl -LsSf https://astral.sh/uv/install.sh > /tmp/install_uv.sh
sh /tmp/install_uv.sh

uv sync --all-etras --refresh
uv run cfeintact version
