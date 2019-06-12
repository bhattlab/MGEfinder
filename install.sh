#!/usr/bin/env bash

echo "Deactivating any active conda environment..."
source deactivate;

echo "Removing mgefinder environment if already installed..."
conda env remove -n mgefinder --yes

echo "Installing mgefinder environment..."
conda env create -f envs/environment.yml

echo "Installation Complete."
echo ""
echo "Before running mgefinder, activate the mgefinder environment with"
echo "> source activate mgefinder"
echo ""
echo "You can then run mgefinder by typing"
echo "> mgefinder [command]"
echo "Send any questions to mdurrant@stanford.edu"
