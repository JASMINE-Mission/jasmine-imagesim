echo "Light version of test.sh for debugging"

python ../../bin/mksimace.py  -n 100000 -t 2400 -e ../../photonsim/data/ace_001.json -m ../../photonsim/data/ace_1_0001.fits 
python ../../bin/mksimace.py  -n 100000 -t 2400 -e ../../photonsim/data/ace_001.json -m ../../photonsim/data/ace_1_0002.fits
python ../../bin/mkplanetlc.py -s 12.5 -f 170 -p ../../params/templates/planet.json -o ../../photonsim/data/lc.fits -m

python mkpixcube_uniform.py -l ../../photonsim/data/lc.fits -v 0.276 -w 0.276 -p 0.423 -x ../../photonsim/data/ace_1_0001.fits -y ../../photonsim/data/ace_1_0002.fits -n 520 -s 12.5 -f 170 -d 0.0 --det ../../params/templates/detk.json -o pixcube_uni_v1.h5  --psf ../../photonsim/data/psf400-2-125-917.fits.gz -m
