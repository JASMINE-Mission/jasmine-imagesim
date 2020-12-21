echo "Light version of test.sh for debugging"

#python ../../bin/mksimace.py  -n 1500000 -t 15.0 -e ../../photonsim/data/ace_001.json -m ../../photonsim/data/ace_1_0001.fits 
#python ../../bin/mksimace.py  -n 1500000 -t 15.0 -e ../../photonsim/data/ace_001.json -m ../../photonsim/data/ace_1_0002.fits

python mkpixcube_sampling.py -v 0.276 -w 0.276 -p 0.423 -x ../../photonsim/data/ace_1_0001.fits -y ../../photonsim/data/ace_1_0002.fits -n 520 -s 12.5 -f 1 -d 0.0 --det ../../params/templates/detk.json -o pixcube_samp1.h5  --psf ../../photonsim/data/psf400-2-125-917.fits.gz -m -u 1000

python mkpixcube_sampling.py -v 0.276 -w 0.276 -p 0.423 -x ../../photonsim/data/ace_1_0001.fits -y ../../photonsim/data/ace_1_0002.fits -n 520 -s 12.5 -f 1 -d 0.0 --det ../../params/templates/detk.json -o pixcube_samp2.h5  --psf ../../photonsim/data/psf400-2-125-917.fits.gz -m -u 100

python mkpixcube_sampling.py -v 0.276 -w 0.276 -p 0.423 -x ../../photonsim/data/ace_1_0001.fits -y ../../photonsim/data/ace_1_0002.fits -n 520 -s 12.5 -f 1 -d 0.0 --det ../../params/templates/detk.json -o pixcube_samp3.h5  --psf ../../photonsim/data/psf400-2-125-917.fits.gz -m -u 10

#python mkpixcube_sampling.py -v 0.276 -w 0.276 -p 0.423 -x ../../photonsim/data/ace_1_0001.fits -y ../../photonsim/data/ace_1_0002.fits -n 520 -s 12.5 -f 1 -d 0.0 --det ../../params/templates/detk.json -o pixcube_samp0.h5  --psf ../../photonsim/data/psf400-2-125-917.fits.gz -m -u 1

