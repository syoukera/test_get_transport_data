# test_get_transport_data
chemkin純正コードのcklib.fとtplib.fを用いて，拡散係数と熱伝導率，定圧比熱を返すプログラムです．

## 使い方
このコード群は二つの機能に分かれています．1のステップについて，GRI-mechを用いる場合には不要です．（GRIの出力ファイルをlinkフォルダに既に配置しているので）
1. 反応機構と熱力学データおよび輸送データからcklinkとtplinkへの変換
2. cklinkとtplinkを用いてCFDの任意の条件から拡散係数と熱伝導率の取得

### 動作確認した環境
- Ubuntu 18.04 LTS
- gfortranが使用できること

### 1. cklinkとtplinkへの変換
必要なファイルは以下，現状ではGRI-mechを置いています．
- data/mech
- data/thermdat
- data/trandat

手順は以下の通りです．
1. 上に記述したような形で，反応機構をmechに，熱力学データをthermdatに，輸送係数データをtrndatとして保存する
2. 以下のように変換のスクリプトを実行
```bash
bash interp.sh
```
3. linkフォルダ以下にcklinkとtplink（バイナリファイル）が出力されることを確認
4. outputフォルダckoutとtpout（テキストファイル）が出力されているので，詳しく知りたい際にはそちらを参照

### 2. 拡散係数と熱伝導率の取得
サンプルプログラム（main.f90）を実行する場合．
上の手順を全て実行した後の状態を想定（linkフォルダ以下にcklinkとtplinkがある状態）

1. コンパイルを実行
```bash
bash compile.sh
```
2. コードの実行
```bash
./get_transport
```
3. 標準出力に化学種ごとの拡散係数，熱伝導率，定圧比熱が出力される．

### 備考
- 拡散係数の算出には，平均化された(Averaged)ものを使用している．より厳密な他成分系の拡散係数も算出可能だが，今回は使用していない．
- 定圧比熱はtplibではなく，cklibから計算できるが，ついでなので一緒に出力するようにした．
