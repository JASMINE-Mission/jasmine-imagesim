[![pytest](https://github.com/JASMINE-Mission/jasmine-imagesim/actions/workflows/pytest.yml/badge.svg?branch=develop)](https://github.com/JASMINE-Mission/jasmine-imagesim/tree/develop)
[![build](https://github.com/JASMINE-Mission/jasmine-imagesim/actions/workflows/build.yml/badge.svg?branch=develop)](https://github.com/JASMINE-Mission/jasmine-imagesim/actions/workflows/build.yml)
[![Maintainability](https://api.codeclimate.com/v1/badges/6c68173a4da5ef19621c/maintainability)](https://codeclimate.com/github/JASMINE-Mission/jasmine-imagesim/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/6c68173a4da5ef19621c/test_coverage)](https://codeclimate.com/github/JASMINE-Mission/jasmine-imagesim/test_coverage)

# jasmine-imagesim (jis)
This is an image simulator package developed for the JASMINE project.

## pixsim
Detector simulator in a subpixel level. Pixsim uses a GPU parallel computing to boost the computation.

## photonsim
This module calculates PSF taking the pupil pattern, wavefront error, spectral response, 

and spectral property of the target into account. This module also calculates

the attitude control error.

# Structure
This package basically consists of 'bin', 'eval', 'src', test directories.
- bin: Executable examples for pixel simulations are stored.
- eval: examples for evaluation.
- src: Modules are stored. The modules can be used for general purposes.
- test: unit tests (pytest) and manual tests

# Requirements
- python3 >= 3.7
- numpy >= 1.17
- cuda >= 7.5 (at least)

# Tested

- CUDA 10.1, RTX 2080 MAX-Q, Python 3.7.3,  Ubuntu 18.04.3 LTS
- CUDA 10.0, TESLA V100, Python 3.7.4, Ubuntu 18.04.5 LTS
- CUDA 9.1, Titan X, Python 3.7.0, Ubuntu 18.04.1 LTS
- CUDA 9.1, TESLA V100, Python 3.7.0, Ubuntu 18.04.1 LTS
- CUDA 9.1, GTX 1080Ti, Python 3.7.3, Linux Mint 19
- CUDA 7.5, GTX Titan, Python 3.6.8, Ubuntu 14.04

# Installation

```
(for global) python setup.py install
(for local) python setup.py install --user
```

## setting for pixsim
Add a line below in your .bash_profile.

```
export CPLUS_INCLUDE_PATH="[jis-home]/src/jis/pixsim/include:$CPLUS_INCLUDE_PATH"
```

For a c shell-based environment, add the below in your .cshrc or .tcsrh etc.
```
setenv CPLUS_INCLUDE_PATH [jis-home]/src/jis/pixsim/include
```
or copy .h files to your "include" directory.


The error such as 
```
fatal error: pixlight_custom.h: No such file or directory\n     #include "pixlight_custom.h"\n                                 ^\ncompilation terminated.\n']
```
indicates that you do not set CPLUS_INCLUDE_PATH properly.

# Examples and Test

## Examples

- modify some parameter files.
  - [jis-home]/params/templates/det.json
    - rename intrapix.dirname ("/home/kamizuka/jasmine-imagesim/" -> [jis-home])
  - [jis-home]/params/templates/variability.json
    - rename dirname(s) ("/home/kawahara/jasmine-imagesim/" -> [jis-home])
- cd [jis-home]/bin
- execute test_light.sh

## Unit test

```
cd test
pytest
```

## other test

There are many manual tests in test directory
- test.sh (not maintained)


# Usage
- To use modules in the jis package, you can call them by typing 'from jis import xxxx' in your script.
- To use executables, you can execute files in [jis-home]/bin/mkimage.py


# Sphinx documents
- To generate sphinx documents, perform

```
apt install python3-sphinx
sphinx-build -b html ./docs ./docs/_build/
```

Then, check docs/_build/index.html with your browser.
