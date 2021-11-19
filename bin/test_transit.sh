echo "Light version of test.sh for debugging"

python mkimage.py --pd ../params/templates/ --starplate star_plate.csv --det det.json --tel tel.json --ace ace_001.json --ctl ctl_planet.json --format platefits --overwrite --var variability.json --dft drift.json
