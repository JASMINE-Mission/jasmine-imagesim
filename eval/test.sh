echo "test mkaperture"
python mkaperture.py -t ../original_codes/photonsim/data/Tel400a.json -a ../original_codes/photonsim/data/apert400a.fits.gz
echo "test mkresponse"
python mkresponse.py -t ../original_codes/photonsim/data/Tel400a.json -s ../original_codes/photonsim/data/Obj2-125.json -d ../original_codes/photonsim/data/Det01.json -r ../original_codes/photonsim/data/spc400-2-125-00.json
#echo "test random wfe gen"
#python rnd_wfe.py -n 8 -e 50 -o 25 -z ../original_codes/photonsim/data/wfe01.json
echo "test wfe"
python mkwfe.py  -t ../original_codes/photonsim/data/Tel400a.json -e ../original_codes/photonsim/data/wfe01.json -m ../original_codes/photonsim/data/wfe01-400.fits.gz
echo "test psf"
python mkpsf.py -a ../original_codes/photonsim/data/apert400a.fits.gz -w ../original_codes/photonsim/data/wfe01-400.fits.gz -s ../original_codes/photonsim/data/spc400-2-125-00.json -c ../original_codes/photonsim/data/Ctl917.json -p ../original_codes/photonsim/data/psf400-2-125-917.fits.gz
echo "test ace"

