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
mkaperture.py -t ../photonsim/data/Tel400a.json -a ../photonsim/data/apert400a.fits.gz
```
入力 望遠鏡のパラメータ        data/Tel400a.json  
出力 Aperture mask のFITS画像  data/apert400a.fits.gz

1. 波長依存性 :
天体の波長依存性、光学系の波長特性、検出器の波長感度特性から
波長ごとの光子数を求める。
```sh
mkresponse.py -t ../photonsim/data/Tel400a.json -s ../photonsim/data/Obj2-125.json -d ../photonsim/data/Det01.json -r ../photonsim/data/spc400-2-125-00.json
```
入力 望遠鏡のパラメータ data/Tel400a.json  
入力 天体のパラメータ   data/Obj2-125.json  
入力 検出器パラメータ   data/Det01.json   
出力 光子の波長分布     data/spc400-2-125-00.json

1. 波面収差の Zernike 係数生成
wave front error のために zernike で n=2..(nmax-1) の 強度データを作成する。
長さのスケールに用いる波長は 1.4 に固定です。
n が奇数のときは、強度を 1.4/zodd とし、偶数の時は 1.4/zeven とする。
そして、 m=0 ではこの強度をそのまま用い、 m がゼロでなければ、
この強度を正の m の項に に cos(θ) 倍、 m が負の項に sin(θ) 倍したものを
用いることにする。ここで θ は 0 から 2π の一様乱数である。
```sh
rnd_wfe.py -n 8 -e 50 -o 25 -z ../photonsim/data/wfe01.json
```
   -n nmax      Zernike max n  
   -e z_even    Strength of Zernike poly. at even odder is 1.4/z_even  
   -o z_odd     Strength of Zernike poly. at odd  odder is 1.4/z_odd   
   -z wfe.json  Generated Zernike polynomials strength  

1. 波面収差マップ作成 :
波面収差を Zernike 展開したときの、各項の係数を与え、そこから
波面収差マップを作成。
```sh
mkwfe.py  -t ../photonsim/data/Tel400a.json -e ../photonsim/data/wfe01.json -m ../photonsim/data/wfe01-400.fits.gz
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
mkpsf.py -a ../photonsim/data/apert400a.fits.gz -w ../photonsim/data/wfe01-400.fits.gz -s ../photonsim/data/spc400-2-125-00.json -c ../photonsim/data/Ctl917.json -p ../photonsim/data/psf400-2-125-917.fits.gz
```
入力  Aperture mask        data/apert400a.fits.gz  
入力  WFE map              data/wfe01-400.fits.gz  
入力  光子の波長分布       data/spc400-2-125-00.json  
入力  コントールパラメータ data/Ctl917.json  
出力  実効PSF              data/psf400-2-125-917.fits.gz  

1. Attitude Control Error のシミュレーション:
一次元の姿勢揺らぎデータを乱数によって作成する。姿勢ゆらぎのパワースペクトル密度のモデルを与えておく必要がある。このモデルについては別資料(position2020B.pdf)を参考にすること。結果は一次元の FITS画像で保存。また、グラフ出力も可能。
```sh
mksimace.py  -n 1048576 -t 15 -e ../photonsim/data/ace_001.json -m ../photonsim/data/ace_1_0001.fits -p ../photonsim/data/ace_1_1.png
```
  -n N         Length of output N steps  
  -t T         total time  
  -e ace.json  atitude control error parameter file  
  -m ace.fits  output atitude control error  
  -p plot.png  plot output  

1. ACE シミュレーション から 2次元マップ:
```sh
mkace2d.py -v 1 -w 1 -p 0.025 -x ../photonsim/data/ace_1_0001.fits -y ../photonsim/data/ace_1_0002.fits -n 520 -o ../photonsim/data/ace2d_001_00010002.fits
```
   -x aX.fits   X-axis simulated ACE file  
   -v xs        aX.fits is scaled by xs  
   -y aY.fits   Y-axis simulated ACE file  
   -w ys        aY.fits is scaled by ys  
   -n N         output image is N x N   
   -p ps        pixel scale of Output image  
   -o out.fits  output image file  

1. PSFの gaussian convolution 
この計算の実体は scipy の ndimage.gaussian_filter です。なので、モジュールとしては作成していません。
```sh
convolv.py -p ../photonsim/data/psf400-2-125-917.fits.gz  -c ../photonsim/data/ace_152.json -o ../photonsim/data/psf400-2-125-917-152.fits.gz
```
   -p psf.fits  
   -c param.json  
   -o output.fits  
