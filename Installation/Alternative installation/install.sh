#!/bin/bash

set -e

ENV_NAME="scDock"
YAML_FILE="environment.yml"

create_env() {

    echo "===== Step 1: Create conda environment ====="

    if [ -z "$CONDA_EXE" ]; then
        echo "Conda executable not found."
        exit 1
    fi

    eval "$("$CONDA_EXE" shell.bash hook)"

    CONDA_BIN="$(dirname "$CONDA_EXE")"

    if [ -x "$CONDA_BIN/mamba" ]; then
        echo "Using mamba..."
        "$CONDA_BIN/mamba" env create -f "$YAML_FILE" -y
    else
        echo "Using conda..."
        conda env create -f "$YAML_FILE" -y
    fi

    conda activate "$ENV_NAME"

    echo "Environment: $CONDA_DEFAULT_ENV"
}

install_python() {

    echo "===== Step 2: Install Python packages ====="

    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"

    python -m pip install --no-deps \
        pynndescent==0.6.0 \
        git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3.git@a62a5d98116f1590183b58a9ad732b997cf2579c
}

install_r() {

    echo "===== Step 3: Install R packages ====="

    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"

    ## remotes

    Rscript -e '
    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes", repos="https://cloud.r-project.org")
    '

    ## NMF

    echo "Installing NMF..."
    Rscript -e 'remotes::install_version("NMF", version = "0.28", repos = "https://cloud.r-project.org")'

    echo "Installing CellChat..."
    Rscript -e 'remotes::install_github("jinworks/CellChat@75253cd0c9e68410e6e721a6d3a0419a1d7e358f", upgrade="never", dependencies=FALSE)'

    echo "Installing scMayoMap..."
    Rscript -e 'remotes::install_github("chloelulu/scMayoMap@993e81a12edf45a80e58a846f7184ef32ee6e488", upgrade="never", dependencies=FALSE)'

    echo "Installing presto..."
    Rscript -e 'remotes::install_github("immunogenomics/presto@a24772a135c7895a8183b007376050556c60a05b", upgrade="never", dependencies=FALSE)'

    echo "Installing harmony..."
    Rscript -e 'remotes::install_github("immunogenomics/harmony@df19af23ae0639bd6ea2da63898f973f08c85862", upgrade="never", dependencies=FALSE)'
}

run_check() {

    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"

    bash env_check.sh
}

usage() {

    echo "Usage:"
    echo "  bash install.sh env"
    echo "  bash install.sh python"
    echo "  bash install.sh r"
    echo "  bash install.sh check"
    echo "  bash install.sh all"
}

case "${1:-all}" in

    env)
        create_env
        ;;

    python)
        install_python
        ;;

    r)
        install_r
        ;;

    check)
        run_check
        ;;

    all)
        create_env
        install_python
        install_r
        run_check
        ;;

    *)
        usage
        exit 1
        ;;

esac
