echo "Light version of test.sh for debugging"
python mksimace.py  -n 100000 -t 7.0 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0001.fits 
python mksimace.py  -n 100000 -t 7.0 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0002.fits
python mkpixcube.py -v 1 -w 1 -p 0.025 -x ../photonsim/data/ace_1_0001.fits -y ../photonsim/data/ace_1_0002.fits -n 520 -s 7.0 -f 1

#python mksimace.py  -n 100000 -t 2400 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0001.fits 
#python mksimace.py  -n 100000 -t 2400 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0002.fits
#python mkpixcube.py -v 1 -w 1 -p 0.025 -x ../photonsim/data/ace_1_0001.fits -y ../photonsim/data/ace_1_0002.fits -n 520 -s 7.0 -f 342
