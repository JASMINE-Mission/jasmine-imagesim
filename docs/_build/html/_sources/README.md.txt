# how to generate sphinx documents

'''
sphinx-apidoc -f -o docs/bin bin
sphinx-apidoc -f -o docs/src src
'''

then

'''
cd docs
make clean
make html
'''

