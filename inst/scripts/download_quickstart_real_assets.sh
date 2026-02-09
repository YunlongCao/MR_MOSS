#!/usr/bin/env bash
set -euo pipefail

BASE_URL="https://raw.githubusercontent.com/YunlongCao/MR_MOSS/main/inst/extdata"
OUT_DIR="${1:-quickstart_real_assets}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

wget -O quickstart_real_assets.zip.part-aa "${BASE_URL}/quickstart_real_assets.zip.part-aa"
wget -O quickstart_real_assets.zip.part-ab "${BASE_URL}/quickstart_real_assets.zip.part-ab"
cat quickstart_real_assets.zip.part-aa quickstart_real_assets.zip.part-ab > quickstart_real_assets.zip
unzip -o quickstart_real_assets.zip

echo "Extracted to: $(pwd)/quickstart_real"
