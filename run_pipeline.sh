#!/bin/bash

# Check for 3 arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: ./run_pipeline.sh <h5ad_file> <reference_variant> <figures_dir>"
  exit 1
fi

# Run normalize.py
echo "Running normalize.py..."
python3 normalize.py "$1" "$3" output.h5ad

# Run scoring.py
echo "Running scoring.py..."
python3 scoring2.py output.h5ad "$2" "$3" output.csv

# Run Cluster.py
echo "Running Cluster.py..."
python3 Cluster.py output.h5ad output.csv "$3"

rm output.h5ad
rm output.csv

echo "Pipeline complete."
