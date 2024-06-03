#!/bin/bash

for script in src/generate_*.sh; do
    echo "----------------------------------"
    echo "Running $script..."
    echo "----------------------------------"
    bash "$script"
    echo
done
