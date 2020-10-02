# jasmine-imagesim/photonsim 

## はじめに

これはJASMINEで得られる画像のシミュレーションツールの説明です。

## 共通事項

* 入力パラメータとしてファイルから与える時は、json 形式を使う

* 画像は全て fits 形式。

## ツール

1. Aperture mask 生成 : 
望遠鏡の直径、中心・スパイダー遮蔽から 入射瞳面の波面強度マップを作成。
これは、単純な円形開口であれば開口内で1で開口外では0となる画像。
```sh
../bin/mkaperture.py -t ../data/Tel400a.json -a ../data/apert400a.fits.gz
```
入力 望遠鏡のパラメータ        data/Tel400a.json  
出力 Aperture mask のFITS画像  data/apert400a.fits.gz

1. 波長依存性 :
天体の波長依存性、光学系の波長特性、検出器の波長感度特性から
波長ごとの光子数を求める。
```sh
../bin/mkresponse.py -t ../data/Tel400a.json -s ../data/Obj2-125.json -d ../data/Det01.json -r ../data/spc400-2-125-00.json
```
入力 望遠鏡のパラメータ data/Tel400a.json  
入力 天体のパラメータ   data/Obj2-125.json  
入力 検出器パラメータ   data/Det01.json   
出力 光子の波長分布     data/spc400-2-125-00.json

1. 波面収差マップ作成 :
波面収差を Zernike 展開したときの、各項の係数を与え、そこから
波面収差マップを作成。
```sh
../bin/mkwfe.py  -t ../data/Tel400a.json -e ../data/wfe01.json -m ../data/wfe01-400.fits.gz
```
入力 望遠鏡のパラメータ data/Tel400a.json
入力 WFEのパラメータ    data/wfe01.json
出力 WFE map            data/wfe01-400.fits.gz

1. 実効PSF作成:
Aperture mask、波長依存性、波面収差から PSFを作成する。これは
焦点面における、ある時刻のPSFを、波長範囲内を光子数平均をとって
求めたもので、計算セルごとの光子数/秒を与える。時間的な揺らぎは
考慮しておらず、結果は時間平均値である。
```sh
../bin/mkpsf.py -a ../data/apert400a.fits.gz -w ../data/wfe01-400.fits.gz -s ../data/spc400-2-125-00.json -c ../data/Ctl917.json -p ../data/psf400-2-125-917.fits.gz
```
入力 Aperture mask        data/apert400a.fits.gz 
入力 WFE map              data/wfe01-400.fits.gz
入力 光子の波長分布       data/spc400-2-125-00.json
入力 コントールパラメータ data/Ctl917.json
出力 実効PSF              data/psf400-2-125-917.fits.gz


