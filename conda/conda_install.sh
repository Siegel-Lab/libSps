#!/bin/bash

# build conda package
conda build libsps -c conda-forge 

conda install --use-local libsps