#!/usr/bin/env bash
set -euo pipefail

# Create local tool directories that can be copied or added to PATH on HPC.
mkdir -p software/bin software/downloads

# Download PLINK2 Linux binary for HPC deployments.
echo "Downloading Linux x86_64 PLINK2..."
curl -L \
  --fail \
  --output software/downloads/plink2_linux_x86_64.zip \
  https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip

unzip -o software/downloads/plink2_linux_x86_64.zip -d software/bin

# Download ADMIXTURE only for deployments that prefer a manual/site-managed binary.
echo "Downloading optional Linux x86_64 ADMIXTURE..."
curl -L \
  --fail \
  --output software/downloads/admixture_linux-1.4.0.tar.gz \
  https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz

tar -xzf software/downloads/admixture_linux-1.4.0.tar.gz -C software/downloads

# Install ADMIXTURE into the same local bin directory when the archive layout matches.
if [ -f software/downloads/admixture_linux-1.4.0/admixture ]; then
  cp software/downloads/admixture_linux-1.4.0/admixture software/bin/admixture
  chmod +x software/bin/admixture
fi

# Print the PATH line users can paste into an HPC job or shell session.
echo "Downloaded tools to software/bin."
echo "Add this directory to PATH on HPC if you are not using modules or conda:"
echo "  export PATH=\"$PWD/software/bin:\$PATH\""
