#!/usr/bin/env bash
set -e

# === Step 0: check conda environment ===
if conda info --envs | grep -q "^scDock"; then
  echo "✅ Conda environment 'scDock' already exists."
else
  echo "=== Creating conda environment: scDock ==="
  conda create -y -n scDock python=3.10 r-base=4.1.3 -c conda-forge
fi

# activate environment
eval "$(conda shell.bash hook)"
conda activate scDock

echo "=== Step 1: install system packages ==="
sudo apt update
sudo apt install -y \
  build-essential \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  zlib1g-dev \
  libncurses-dev \
  libbz2-dev \
  liblzma-dev \
  libpcre2-dev \
  libhdf5-dev \
  libgsl-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libglpk-dev \
  pkg-config \
  gfortran

echo "=== Step 2: install conda packages ==="
conda install -y -c conda-forge \
  zlib \
  r-curl \
  r-httr \
  r-hdf5r \
  r-mass \
  r-ggplot2 \
  r-igraph \
  r-leidenbase \
  r-devtools \
  r-nloptr \
  rdkit \
  hdf5 \
  requests \
  r-remotes \
  r-matrix \
  r-dplyr \
  r-cowplot \
  r-data.table \
  r-reticulate \
  r-fastdummies \
  r-ggrepel \
  r-ggridges \
  r-miniui \
  r-patchwork \
  plotly \
  r-png \
  r-scattermore \
  r-sctransform \
  shiny \
  r-rspectra \
  r-leiden \
  r-irlba \
  r-yaml \
  r-seurat \
  r-biocmanager \
  r-nmf \
  r-circlize \
  r-lme4 \
  r-pbkrtest \
  r-car \
  r-rstatix \
  r-ggpubr

conda install -y -c bioconda \
  bioconductor-biobase \
  bioconductor-complexheatmap \
  bioconductor-biocneighbors \
  presto

echo "=== Step 3: install Python pip modules ==="
python -m pip install --upgrade pip
python -m pip install \
  git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3 \
  umap-learn

echo "=== Step 4: install R GitHub packages ==="
Rscript -e 'if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")'
Rscript -e 'devtools::install_github("jinworks/CellChat")'
Rscript -e 'devtools::install_github("chloelulu/scMayoMap")'
Rscript -e 'devtools::install_github("immunogenomics/presto")'

echo "=== Step 5: check R packages for scMayoMap ==="
Rscript -e 'pkgs <- c("ggplot2","dplyr","tidyr","tibble","reshape2"); new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]; if(length(new.pkgs)) install.packages(new.pkgs, repos="https://cloud.r-project.org")'

echo "✅ Install complete！"
