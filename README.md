# jasmine-imagesim (jis)
This is an image simulator package developed for the JASMINE project.

## pixsim
combining detector characteristics. To use pixsim you need to specify include directory like

setenv CPLUS_INCLUDE_PATH /home/kawahara/jasmine-imagesim/src/jis/pixsim/include

or copy .h files to your "include" directory.

## photonsim
This module calculates PSF taking the pupil pattern, wavefront error, spectral response, 

and spectral property of the target into account. This module also calculates

the attitude control error.

# Structure
This package basically consists of 'bin' and 'src' directories.
- bin: Executables are stored.
- src: Modules are stored. The modules can be used for general purposes.

# Requirements
- python 3.8 (3.7 might be ok)

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


# Test
- cd [jis-home]/bin
- bash test.sh

# Usage
- To use modules in the jis package, you can call them by typing 'from jis import xxxx' in your script.
- To use executables, you can execute files in [jis-home]/bin/xxxx.py.
