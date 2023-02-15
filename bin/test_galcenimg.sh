echo "Sample to make galcen images"
python galcenimg.py --pd ../params/templates/ --starplate star_plate.csv --det det.json --tel tel.json --ace ace_001.json --ctl ctl.json --format platefits --overwrite --var variability.json
