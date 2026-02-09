#!/usr/bin/env bash
set -euo pipefail

# Download the packaged quickstart demo bundle without cloning the full repo.
ASSET_URL="https://raw.githubusercontent.com/YunlongCao/MR_MOSS/main/inst/extdata/quickstart_demo_assets.zip"
OUT_DIR="${1:-quickstart_demo_assets}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"
wget -O quickstart_demo_assets.zip "${ASSET_URL}"
unzip -o quickstart_demo_assets.zip
echo "Downloaded and extracted to: $(pwd)/quickstart_demo"
