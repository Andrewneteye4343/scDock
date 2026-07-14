#!/bin/bash

FAILED=0

eval "$(conda shell.bash hook)"
conda activate scDock

echo "===== Python ====="

if python3 --version; then
    echo "Python OK"
else
    echo "Python FAILED"
    FAILED=1
fi

python3 <<EOF
import sys

try:
    import AutoDockTools
except ModuleNotFoundError:
    print("Missing Python package: AutoDockTools")
    sys.exit(1)

try:
    import umap
except ModuleNotFoundError:
    print("Missing Python package: umap")
    sys.exit(1)
EOF

if [ $? -eq 0 ]; then
    echo "Python packages OK"
else
    echo "Python packages FAILED"
    FAILED=1
fi

echo
echo "===== R ====="

if command -v Rscript >/dev/null 2>&1; then
    Rscript - <<EOF
pkgs <- c(
    "Seurat",
    "CellChat",
    "scMayoMap",
    "presto",
    "harmony",
    "NMF"
)

failed <- FALSE

for (p in pkgs){
    if (!requireNamespace(p, quietly = TRUE)){
        cat("Missing R package:", p, "\n")
        failed <- TRUE
    }
}

if (failed) quit(status = 1)
EOF

    if [ $? -eq 0 ]; then
        echo "R packages OK"
    else
        echo "R packages FAILED"
        FAILED=1
    fi
else
    echo "Rscript not found"
    echo "R packages FAILED"
    FAILED=1
fi

echo
echo "===== Versions ====="

python3 <<EOF
import importlib.metadata

packages = [
    "autodocktools-py3",
    "umap-learn"
]

for p in packages:
    try:
        print(f"{p:<20} {importlib.metadata.version(p)}")
    except importlib.metadata.PackageNotFoundError:
        print(f"{p:<20} NOT INSTALLED")
EOF

if command -v Rscript >/dev/null 2>&1; then
    Rscript - <<EOF
pkgs <- c(
    "CellChat",
    "scMayoMap",
    "presto",
    "harmony",
    "NMF"
)

for (p in pkgs){
    if (requireNamespace(p, quietly = TRUE)){
        cat(sprintf("%-12s %s\n",
            p,
            as.character(packageVersion(p))
        ))
    }
    else{
        cat(sprintf("%-12s NOT INSTALLED\n", p))
    }
}
EOF
fi

echo
echo "================================"

if [ $FAILED -eq 0 ]; then
    echo "Environment verification PASSED"
    echo "================================"
    exit 0
else
    echo "Environment verification FAILED"
    echo "================================"
    exit 1
fi
