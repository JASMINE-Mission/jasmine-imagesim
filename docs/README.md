# how to generate sphinx documents

'''

cd <JIS_HOME>

sphinx-apidoc -f -o docs/bin bin

sphinx-apidoc -f -o docs/src src

'''

then

'''

cd docs

make clean

make html

'''

