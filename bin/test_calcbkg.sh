echo "BKGD PSF IMAGE"

python mkbkgimg.py --pd ../params/bkgimg/ --starplate star_plate.csv --det det.json --tel tel.json --ace ace_001.json --ctl ctl.json --format platefits --overwrite --var variability.json
