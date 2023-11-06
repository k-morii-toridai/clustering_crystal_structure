from .clustermanager import ClusterManager
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist 
import polars as pl
import pandas as pd
import numpy as np


import itertools
import os
import pandas as pd
import numpy as np
import glob
import re
'''
def cal_distances(cluster_manager: ClusterManager,reference=1e-8,target_atoms=['Si1','O1'],chunksize=10000):
        if cluster_manager.target_combination_df is None:
            cluster_manager.calculate_self_distance_file()
        target_files = pl.from_pandas(cluster_manager.target_combination_files)
        chunksize = chunksize if len(target_files) > chunksize else len(target_files)
        result = list()
        for fram in target_files.iter_slices(n_rows=chunksize):
            results_dis = np.zeros(len(fram))
            target_cluster_coordinates = fram.select(
                        pl.all().apply(lambda file_i: pl.scan_csv(file_i))
                    )
            for pickup_atom in target_atoms:
                pickup_coordinates = target_cluster_coordinates.select(
                                pl.all().apply(lambda x: x.filter(pl.col('atom') == pickup_atom).select(['x', 'y', 'z']).collect().to_numpy())
                            )
                pickup_coordinates = pickup_coordinates.to_numpy()
                A = np.stack(pickup_coordinates[:,0],0)
                B = np.stack(pickup_coordinates[:,1],0)
                if A.shape==B.shape:
                    n = A.shape[1]
                    # ユークリッド距離
                    d = [cdist(ai, bi) for ai, bi in zip(A, B)]
                    # 線形割当問題の解
                    assignment = [linear_sum_assignment(di) for di in d]
                    # コスト
                    distance = np.array([di[assignmenti].sum() / n for di,assignmenti in zip(d, assignment)])
                else:
                    distance = np.array([float('nan')]*A.shape[0])
                results_dis = np.sum(np.stack([results_dis,distance]),axis=0)
            results_dis = results_dis/len(target_atoms)
            results_dis = np.where(results_dis < reference , 0.0 , results_dis)
            results_dis = results_dis.tolist()
            result.extend(results_dis)
        return pd.concat([cluster_manager.target_combination_df, pd.DataFrame(result, columns=['distance'])], axis=1)
'''


def cal_emd(A:np.array,B:np.array) -> np.array:
    #A,B:原子の座標
    #A,Bは共にm*n*3の配列である必要がある
    if A.shape==B.shape:
        #原点を除く原子の個数
        n = A.shape[1]-1 if (A[0,:]==0).all() else A.shape[1]
        # ユークリッド距離
        d = [cdist(ai, bi) for ai, bi in zip(A, B)]
        # 線形割当問題の解
        assignment = [linear_sum_assignment(di) for di in d]
        # コスト
        distance = np.array([di[assignmenti].sum() / n for di,assignmenti in zip(d, assignment)])
    else:
        distance = np.array([float('nan')]*A.shape[0])
    return distance


def load_files_from_dirpath(dirpath):
    pattern = re.compile(r'[0-9]+\.csv$')
    files = [file for file in glob.glob(os.path.join(dirpath, '*.csv')) if re.search(pattern, file)]
    pattern = re.compile(r'\_([0-9]+)\_')
    return sorted(files, key=lambda x: int(pattern.findall(x)[0]))

def make_combination_indexs_pattern_only_zero(files):
    pattern_list = [i for i, file in enumerate(files) if re.search(r'_0\.csv$', file)]
    return list(itertools.combinations(pattern_list, 2))

def create_combination_indices(files):
    only_zero = make_combination_indexs_pattern_only_zero(files)
    
    result = []
    for target_comb in only_zero:
        match_str = re.sub(r'0\.csv$', '', files[target_comb[1]])#0.csvを削除
        pattern = re.compile(match_str)#正規表現のパターンを作成,0.csvを削除した文字列と一致するファイルを探す
        pattern_list = [i for i, file in enumerate(files) if pattern.match(file)]#一致するファイルのインデックスを取得
        result.extend([(target_comb[0], i) for i in pattern_list])
    return result

def load_files(files):
    return [pd.read_csv(i,index_col=0,usecols=['atom','x','y','z']) for i in files]

def cal_distance(coordinates_a:dict [np.array],coordinates_b:dict [np.arange],target_atoms:list) -> np.array:
    #coordinates_a,coordinates_bは同じサイズのnumpy配列である必要がある。
    distances = list()
    for target_atom in target_atoms:
        distance = cal_emd(coordinates_a[target_atom],coordinates_b[target_atom])
        distances.append(distance)
    distances = np.sum(distances,axis=0)/len(target_atoms)
    distances = np.where(distances < 1e-8 , 0.0 , distances)
    return distances

def format_file_name(path:str) -> dict:
    #dataから、address,cifid,isite,patternを取得する
    address = os.path.dirname(path)
    cifid,isite,pattern = re.search(r'([A-Z]+)_([0-9]+)_([0-9]+)\.csv$',path).groups()
    return {'address':address,'cifid':cifid,'isite':int(isite),'pattern':int(pattern)}

def make_log_format(files,combination_indexs):
    #combination_indexs = [i for j in combination_indexs for i in j]
    all_alist,all_blist = zip(*combination_indexs)
    all_alist,all_blist = list(all_alist),list(all_blist)
    distances_log = pd.DataFrame({'file_i':[files[i] for i in all_alist],'file_j':[files[i] for i in all_blist]})
    distances_log_i = distances_log['file_i'].apply(format_file_name).apply(pd.Series).rename(columns={'address':'address_i','cifid':'cifid_i','isite':'isite_i','pattern':'pattern_i'})
    distances_log_j = distances_log['file_j'].apply(format_file_name).apply(pd.Series).rename(columns={'address':'address_j','cifid':'cifid_j','isite':'isite_j','pattern':'pattern_j'})
    return pd.concat([distances_log_i,distances_log_j],axis=1)


class CalulateSelfDistance:
    def __init__(self,target_atoms=['Si1','O1'],reference=1e-8,chunk:bool | int=30000):
        self.target_atoms = target_atoms
        self.reference = reference
        self.chunk = chunk
    def _calculate_distance(self,coordinates_a:dict [np.array],coordinates_b:dict [np.arange]) -> np.array:
        return cal_distance(coordinates_a,coordinates_b,self.target_atoms)

    def calculate_distance(self,dirpath:str) -> pd.DataFrame:
        #最初に指定したディレクトリに保存されているファイルのパスをリストアップする
        files = load_files_from_dirpath(dirpath)
        #ファイルをロードする
        dfs = load_files(files)
        #必要な座標データを取得する
        coordinates = {atom:np.array([df.query(f'atom=="{atom}"').loc[:,['x','y','z']].to_numpy() for df in dfs]) for atom in self.target_atoms}

        del dfs

        #ファイルの組み合わせを作成する
        combination_indexs = create_combination_indices(files)

        #入力データを作成する
        #チャンクサイズを設定する
        if isinstance(self.chunk,int):
            combination_indexs = [combination_indexs[i:i + self.chunk] for i in range(0, len(combination_indexs), self.chunk)]
        else:
            combination_indexs = [combination_indexs]
        #距離の計算を行う
        distances = list()
        for combination_indexs_i in combination_indexs:
            alist,blist = zip(*combination_indexs_i)
            alist,blist = list(alist),list(blist)
            coordinates_a = {atom:coordinates[atom][alist] for atom in self.target_atoms}
            coordinates_b = {atom:coordinates[atom][blist] for atom in self.target_atoms}
            #距離の計算を行う
            disntance = self._calculate_distance(coordinates_a,coordinates_b)
            distances.extend(disntance)
        #距離の計算結果をデータフレームに変換する
        #return distances
        del coordinates,coordinates_a,coordinates_b,alist,blist
        
        #計算ログのフォーマット作成
        #combnate_indexsを一次元化
        combination_indexs = [i for j in combination_indexs for i in j]
        distances_df = make_log_format(files,combination_indexs)
        distances_df['distance'] = distances
        return distances_df
    
