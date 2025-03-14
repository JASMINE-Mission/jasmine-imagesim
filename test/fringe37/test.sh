#!/bin/bash -x


python ../../bin/mkimage.py \
  --starplate ./params/star_plate.csv \
  --det ./params/det.json \
  --tel ./params/tel.json \
  --ace ./params/ace_001.json \
  --ctl ./params/ctl.json \
  --format platefits \
  --od ./out --overwrite
