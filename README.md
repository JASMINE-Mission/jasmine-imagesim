# jasmine-imagesim (jis)
This is an image simulator package developed for the JASMINE project.

## pixsim
combining detector characteristics 

## (Kataza-san code)
trajectory, optics, ...

# Structure
This package basically consists of 'bin' and 'src' directories.
-bin: Executables are stored.
-src: Modules are stored. The modules can be used for general purposes.

# Requirements
- python 3.8 (3.7 might be ok)

# Installation

```
(for global) python setup.py install
(for local) python setup.py install --user
```

# Test
- cd [jis-home]/bin
- bash test.sh

# Usage
- To use modules in the jis package, you can call them by typing 'from jis import xxxx' in your script.
- To use executables, you can execute files in [jis-home]/bin/xxxx.py.
