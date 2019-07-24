## glm-haskell (v 0.1.0.0)

# Introduction
A version of R's lme4 functionality written using hmatrix and sparse-linear-algebra (with Polysemy for IO, error and logging effects).

# Installation
This package requires two external libraries:
- [nlopt](https://nlopt.readthedocs.io/en/latest/), for non-linear optimization
- [cholmod](https://developer.nvidia.com/cholmod) (for fast sparse Cholesky decomposition) 
which is part of Tim Davis's open-source [suite-sparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) set of libraries.

Both are open source and available via, apt-get or homebrew.  On a mac, for example, the following should suffice:

brew install nlopt
brew install suite-sparse

on a linux system: 

sudo apt-get update
sudo apt-get install libnlopt-dev
sudo apt-get install libsuitesparse-dev
