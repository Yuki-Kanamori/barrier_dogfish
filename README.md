(1) メッシュの作成@make_spatial_info
fine_scale==FALSEの場合はINLA::inla.mesh.create()を
fine_scale==TRUEの場合はINLA::inla.nonconvex.hull()とINLA::inla.mesh.create()を使う
-> メッシュを作成する時にouterを入れるか入れないかが違う

(2) spde@make_mesh
anisotropic_spde = INLA::inla.spde2.matern(anisotropic_mesh, alpha=2)

(3) A matrix@make_spatial_info
fine_scale==FALSEの場合はA matrixを手作り
fine_scale==TRUEの場合は A_is = INLA::inla.spde.make.A()をしてdgCMatrix型をdgTMatrix型に変換



本来のTmbDataに入っているもの
TmbData$Xconfig_zcp [2×1×1] これがうまく作成されていない．make_data()でX1config_cpとX2config_cpを作っているところ, TmbData$X_gtp [100×25×1]knot×year×cate?, TmbData$X_itp [2548×25×1]DG×year×cate?

TmbDataの中身
ObsModel_ezが20と10で違うからチェックが必要 ->確認．モデルの違い
a_glの値
X_itp・X_gtp・Q_ikの3次元目が0になっている（本当は1になるはず） ->revised
X_itp・X_gtpは3次元目は0が入る（データを入れていないから） ->revised
Q_ikも（catchabilityを考慮しないなら）0が入るはず ->revised

コンパイルで失敗
make_model()で動くものでも，個別で動かすとコンパイルでエラーが出る
Parametersの中身の確認<-make_parameters(#7)がおかしいはず
gamma1_ctpの3次元目が0になっている．3次元目の値は全部0が入る
lambda1_kがNAになっている．0が入る
*Beta_rho1_fが0になっている．1が入る
*Epsilon_rho1_fが0になっている．1が入る
log_sigmaXi1_cpの2次元目が0になっている．1×1で値は0が入る
Xiinput1_scpの3次元目が0になっている．3次元目は1で0が入る
gamma2_ctpの3次元目が0になっている．3次元目の値は全部0が入る
lambda2_kがNAになっている．0が入る
*Beta_rho2_fが0になっている．1が入る
*Epsilon_rho2_fが0になっている．1が入る
log_sigmaXi2_cpの2次元目が0になっている．1×1で値は0が入る
*delta_iがNAになっている．何の数字が入る？（0.07583）
Xiinput2_scpの3次元目が0になっている．3次元目は1で0が入る

Mapの中身
L_beta1_z, L_beta_2z, Beta_mean_1c, Beta_mean_2cが余分に入っている -> 削除するコードを追加

Obj <- TMB::MakeADFun(data=TmbData, parameters=Parameters, hessian=FALSE, map=Map, random=Random, inner.method="newton", DLL=Version)でエラー: 型 "i" のオブジェクトからはスロット ("NULL") を得ることは出来ません 
引数を一つずつ調べる
①TmbDataの中身（make_data_yk()）
spdeがない
Ais_ijが0×2になっている．DG×2となるはず
Ais_xが0になっている．DGの行数分だけ1が入る
Ags_ijが0×2になっている．knot×2となるはず
Ais_xが0になっている．knotの数分だけ1が入る

loc_g = Extrapolation_List$Data_Extrap[ which(Extrapolation_List$Area_km2_x>0), c('E_km','N_km') ]
スパース行列がよく分からん　dgTMatrix
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    0    0
[2,]    0    0    0    1    1
[3,]    0    1    1    0    1
[4,]    0    0    1    0    0
[5,]    1    0    0    1    1
!!! 1行1列->2行1列->・・・->1行2列->2行2列といった順番で見ている !!!
 gtest@i
[1] 4 2 2 3 1 4 1 2 4
0じゃないところを行方向で調べている．行数は0始まり．1列目は5行目に1が，2列目は3行目に1が，3列目は3と4行目に1が，，，という感じ
 gtest@j
[1] 0 1 2 2 3 3 4 4 4（1,2,3,3,4,4,5,5,5）
0じゃないところを列方向で調べている．1列目は5行目に1が，2列目は3行目に1が，3列目は3と4行目に1が（2,2），4列目は2と5行目に1が（3,3），5列目は2と3と5行目に1が（4,4,4）



②Parameters（make_model()とmake_parameters()）
*Beta_rho1_fが0になっている．1が入る <-RhoConfigによって違う
*Epsilon_rho1_fが0になっている．1が入る  <-RhoConfigによって違う
*Beta_rho2_fが0になっている．1が入る  <-RhoConfigによって違う
*Epsilon_rho2_fが0になっている．1が入る  <-RhoConfigによって違う
*delta_iがNAになっている．何の数字が入る？（0.07583）<-make_parameters()は合っていそう．make_map()で何かしている．make_map()じゃなくてやっぱりmake_parametersっぽいけど，関数を使わずに動かすと何も値が入らない．なぜ？乱数が入るように修正
③Map（make_model()とmake_map()）
