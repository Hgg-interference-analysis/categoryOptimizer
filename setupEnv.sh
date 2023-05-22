#!/usr/bin/env bash

DIR=$1;

conda config --prepend envs_dirs ${DIR}/.conda/envs
conda config --prepend pkgs_dirs ${DIR}/.conda/pkgs
