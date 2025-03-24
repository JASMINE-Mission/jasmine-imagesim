#!/bin/bash -x
TEMPLATE='../../params/templates/'

python ../../bin/mkimage.py \
  --starplate ${TEMPLATE}/star_plate.csv \
  --det ${TEMPLATE}/det.json \
  --tel ${TEMPLATE}/tel.json \
  --ace ${TEMPLATE}/ace_001.json \
  --ctl ./params/ctl.json \
  --format platefits \
  --od ./out --overwrite
