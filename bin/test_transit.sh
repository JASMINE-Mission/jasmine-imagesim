echo "Light version of test.sh for debugging"

python mkimage_transit.py --pd ../params/templates/ --starplate star_plate.csv --det detk.json --tel tel.json --ace ace_001.json --ctl ctlp.json --var variability.json --format platefits --overwrite --dft drift.json
