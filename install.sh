#!/usr/bin/env bash

echo "Removing mgefinder environment if already installed..."
conda env remove -n mgefinder --yes

echo "Installing mgefinder environment..."
conda env create -f conda.yaml

echo "Installation Complete."
echo ""
echo "Before running mgefinder, activate the mgefinder environment with"
echo "> conda activate mgefinder"
echo ""
echo "You can then run mgefinder by typing"
echo "> mgefinder [command]"
echo "Send any questions to mdurrant@stanford.edu"
