# Diphoton pT Boundary Optimizer

This branch is for developements specific to H-->gamma gamma interference analysis. In particular, for optimizing the diphoton pT boundaries. This was achieved by modifying [Neil Raymond Schroeder / diphoton mva optimizer](https://gitlab.cern.ch/nschroed/diphoton-mva-optimizer) framework.

## Features

This code uses the scipy.optimize.minimize and makes use of the Nelder-Mead minimization algorithm. This algorithm is considered to be "greedy" and will try to find the nearest minimum using step sizes determined by the local shape of the loss function. To combat the unusually complex shape of the loss function, a monte carlo approach is used and for a solid set of boundaries it is suggested to run at least 1000 minimizations.

## Installation

the installation is fairly straight forward, and requires scipy.
```
git clone ssh://git@gitlab.cern.ch:7999/amkrishn/diphoton-mva-optimizer.git
python -m pip install scipy
```

PS: I used CERNBox and SWAN to work on this project. If you wish to do the same, follow the instructions given [here](https://swan.docs.cern.ch/swan/create_proj/) to create a SWAN project from this repository. Then you can use the SWAN terminal to run the framework.

One could also follow the instructions given [here](https://gitlab.cern.ch/bmarzocc/diphoton-mva-optimizer/-/tree/dev_for_hgg_interference/#installation) to setup a conda environment on lxplus to install this framework.

## Getting Started

To Run the framework you'll need a few things.
First, you'll need signal and background samples of the processes you're optimizing over.
These samples must have a diphoton invariant mass variable called `CMS_hgg_mass`, the diphoton invariant mass `dipho_mva`, and an event weight `weight`. There are additional tools available for use if you are using a data-driven dataset that has additional weights.

A config file must be provided in the format
```
type    treeName    path    legendEntry
```
and an example can be seen in `config/ul18.cfg`. Please note that the config file is a tab separated file (\t) and if you do not use a tab, the program will crash.

With samples in hand, and prepared, you can now run the framework


## Running

running the framework looks is done using the `optimize.py` tool and an example is given below:

```
./optimize.py -i <path/to/config/file> -o <output-name> --num-iter=<number of iterations> --num-bounds=<number of bounds> --lumi-scale <lumi16> <lumi17> <lumi18> --lumi-scale-bkg <lumibkg16> <lumibkg17> <lumibkg18> --log <path/to/log/file> -b <b1> <b2> ... <bN>
```

The options are as follows:
- `--inputFile`,`-i`: path to config file
- `--output`,`-o`: string to give to outputs
- `--num-bounds`,`-n`: number of categories to use
- `--lumi-scale`: list (separation is 1 space) of lumis (16 17 18) to use in minimization (if only minimizing one year, set the others to 0)
- `--lumi-scale-bkg`: list (separation is 1 space) of background lumis (16 17 18) to use in minimization (if only minimizing one year, set the others to 0)
- `--bkg-scale`: (advanced) scale for background weights
- `--xcheck`: plots variables for each file provided in the config file
- `--plot`: makes a stack plot (will add boundaries if provided)
- `--boundaries`,`-b`: list (separation is 1 space) of boundaries to use as the initial set (or to plot if `--plot` is provided)
- `--log`: path to the log file that will be the output of this framework. 

