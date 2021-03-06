echo "test mkaperture"
python mkaperture.py -t ../photonsim/data/Tel400a.json -a ../photonsim/data/apert400a.fits.gz
echo "test mkresponse"
python mkresponse.py -t ../photonsim/data/Tel400a.json -s ../photonsim/data/Obj2-125.json -d ../photonsim/data/Det01.json -r ../photonsim/data/spc400-2-125-00.json
echo "test random wfe gen"
python rnd_wfe.py -n 8 -e 50 -o 25 -z ../photonsim/data/wfe01.json
echo "test wfe"
python mkwfe.py  -t ../photonsim/data/Tel400a.json -e ../photonsim/data/wfe01.json -m ../photonsim/data/wfe01-400.fits.gz
echo "test psf"
python mkpsf.py -a ../photonsim/data/apert400a.fits.gz -w ../photonsim/data/wfe01-400.fits.gz -s ../photonsim/data/spc400-2-125-00.json -c ../photonsim/data/Ctl917.json -p ../photonsim/data/psf400-2-125-917.fits.gz
echo "test ace"

#python mksimace.py  -n 1048576 -t 15 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0001.fits -p ../photonsim/data/ace_1_1.png
#python mksimace.py  -n 1048576 -t 15 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0002.fits

python mksimace.py  -n 104850 -t 15 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0001.fits -p ../photonsim/data/ace_1_1.png
python mksimace.py  -n 104850 -t 15 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0002.fits

#python mkace2d.py  -v 1 -w 1 -p 0.025 -x ../photonsim/data/ace_1_0001.fits -y ../photonsim/data/ace_1_0002.fits -n 520 -o ../photonsim/data/ace2d_001_00010002.fits

echo "test mkpixcube"
python mkpixcube.py -v 0.276 -w 0.276 -p 0.423 -x ../photonsim/data/ace_1_0001.fits -y ../photonsim/data/ace_1_0002.fits -n 520 -s 12.5 -f 1 -d 0.0 --det ../params/templates/det.json --psf ../photonsim/data/psf400-2-125-917.fits.gz -o pixcube_v1.h5 -m

echo "test mkimage"
python mkimage.py --pd ../params/templates/ --starplate star_plate.csv --var variability.json --det det.json --tel tel.json --ace ace_001.json --ctl ctl.json --format platefits --overwrite

echo "finish."
