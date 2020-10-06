echo "test mkaperture"
python mkaperture.py -t ../photonsim/data/Tel400a.json -a ../photonsim/data/apert400a.fits.gz
echo "test mkresponse"
python mkresponse.py -t ../photonsim/data/Tel400a.json -s ../photonsim/data/Obj2-125.json -d ../photonsim/data/Det01.json -r ../photonsim/data/spc400-2-125-00.json
echo "test wfe"
python mkwfe.py  -t ../photonsim/data/Tel400a.json -e ../photonsim/data/wfe01.json -m ../photonsim/data/wfe01-400.fits.gz
echo "test psf"
python mkpsf.py -a ../photonsim/data/apert400a.fits.gz -w ../photonsim/data/wfe01-400.fits.gz -s ../photonsim/data/spc400-2-125-00.json -c ../photonsim/data/Ctl917.json -p ../photonsim/data/psf400-2-125-917.fits.gz
echo "finish."
