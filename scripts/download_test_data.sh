#!/usr/bin/env bash
set -euo pipefail

# Create the local example-data directory.
mkdir -p data/example

# Download the public HapMap3 archive used by the local toy workflow.
echo "Downloading public ADMIXTURE HapMap3 sample archive for reference/testing..."
curl -L \
  --fail \
  --output data/example/hapmap3-files.tar.gz \
  https://dalexander.github.io/admixture/hapmap3-files.tar.gz

# Extract PLINK bed/bim/fam files for the example run.
tar -xzf data/example/hapmap3-files.tar.gz -C data/example

# Print the expected files so users can verify the download quickly.
echo "Test data ready:"
echo "  data/example/hapmap3.bed"
echo "  data/example/hapmap3.bim"
echo "  data/example/hapmap3.fam"
echo "  data/example/hapmap3-files.tar.gz"
