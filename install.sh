#!/usr/bin/env bash

echo "Deactivating any active conda environment..."
source deactivate;

echo "Removing mustache environment if already installed..."
conda env remove -n mustache --yes

echo "Installing mustache environment..."
conda env create -f envs/environment.yml

echo "Installation Complete."
echo ""
echo "Before running mustache, activate the mustache environment with"
echo "> source activate mustache"
echo ""
echo "You can then run mustache by typing"
echo "> mustache [command]"
echo "Send any questions to mdurrant@stanford.edu"
