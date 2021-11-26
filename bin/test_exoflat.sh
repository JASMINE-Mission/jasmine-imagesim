echo "Light version of test.sh for debugging"

python mkimage.py --pd ../params/exoflat/ --starplate star_plate.csv --det det.json --tel tel.json --ace ace_001.json --ctl ctl.json --format platefits --overwrite --var variability.json --dft drift.json
