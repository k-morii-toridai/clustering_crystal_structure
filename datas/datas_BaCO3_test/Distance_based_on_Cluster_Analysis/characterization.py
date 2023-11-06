import polars as pl
import numpy as np
import pandas as pd
import os
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import glob
from itertools import product
import copy
from .distance import cal_emd

#データベースとクラスターとの距離を計算するクラスを作成する
def tilda(u):
    return np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])

def make_inertia_matrix(cluster_coordinate:np.ndarray) -> np.ndarray:
    cluster_coordinate_tilda = np.array([tilda(u) for u in cluster_coordinate])
    inertia_matrixs = -np.array([np.dot(tilda_vec, tilda_vec) for tilda_vec in cluster_coordinate_tilda]) 
    return np.sum(inertia_matrixs,axis=0)

def cal_eigenvalues(cluster_coordinate:np.ndarray) -> list():
    inertia_matrix = make_inertia_matrix(cluster_coordinate)
    eigenvalues, _ = np.linalg.eig(inertia_matrix)
    return sorted(eigenvalues.real)

def _cal_distance(target_coords,database_coords):
    distances = list()
    for target_atom,database_coord in database_coords.items():
        target_cluster_coordinates = np.array([target_coords[target_atom]]*database_coord.shape[0])
        distance = cal_emd(target_cluster_coordinates,database_coord)
        distances.append(distance)
    distances = np.sum(distances,axis=0)/len(target_coords)
    distances = np.where(distances < 1e-8 , 0.0 , distances)
    return distances

class ClusterFeatureCalculator():
    def __init__(self, databasepath,target_atoms=['Si1','O1'],reference=1e-8,sep_value=0.1,offset=5,use_mesh_flag=True):
        self.targets_atoms = target_atoms
        self.sep_value = sep_value
        self.offset = offset
        self.use_mesh_flag = use_mesh_flag

        database_files = glob.glob(os.path.join(databasepath,'*.csv'))
        
        #データベースのファイルパスとそれに対応する固有値を整理し保存する
        database_path_df = pl.LazyFrame({'file_path':database_files})
        database_path_df = database_path_df.with_columns(
            pl.col("file_path").map_elements(os.path.dirname).alias("address_i"),
            pl.col("file_path").map_elements(lambda x: os.path.basename(x).replace('.csv', '').split('_')[0]).alias("cifid_i"),
            pl.col("file_path").map_elements(lambda x: os.path.basename(x).replace('.csv', '').split('_')[1]).alias("isite_i"),
            pl.col("file_path").map_elements(lambda x: os.path.basename(x).replace('.csv', '').split('_')[2]).alias("pattern_i")
        ).select(pl.col(['address_i','cifid_i','isite_i','pattern_i']))

        self.database_path_df = database_path_df.collect().to_pandas()

        #計算用のcluster_coordinateを作成する
        self.reference = reference
        #データベースを数値データに落とし込む
        database_dfs = [pl.scan_csv(file_i) for file_i in database_files]
        database_coordinates = dict()
        for target_atom in target_atoms:
            df_atoms = [df_i.filter(pl.col('atom')==target_atom).select(['x','y','z']) for df_i in database_dfs]
            df_atoms = pl.collect_all(df_atoms)
            database_coordinates[target_atom] = np.array([df_i.to_numpy() for df_i in df_atoms])
        self.database_coordinates = database_coordinates

        #ログフォーマットの作成
        self.log_format = self.database_path_df.copy()

        #データベースのメッシュを作成する
        self.make_mesh_dict()

    def make_mesh_dict(self):
        #データベースのメッシュを作成する
        self.database_mesh_dict ={}
        self.eig_df = {}
        for atom in self.targets_atoms:
            _mesh_dict = [cal_eigenvalues(coordinate) for coordinate in self.database_coordinates[atom]]
            self.eig_df[atom] = pd.DataFrame(_mesh_dict,columns=['eig_1','eig_2','eig_3'])
            _mesh_df = self.eig_df[atom].filter(like='eig').applymap(lambda x: int(x//self.sep_value)).copy()
            _mesh_df.rename(columns={'eig_1':'eig_1_mesh','eig_2':'eig_2_mesh','eig_3':'eig_3_mesh'},inplace=True)
            _mesh_dict = _mesh_df.groupby(['eig_1_mesh','eig_2_mesh','eig_3_mesh']).apply(lambda x:x.index.tolist()).to_dict()
            self.database_mesh_dict[atom] = copy.deepcopy(_mesh_dict)

    def change_offset(self,offset):
        self.offset = offset
    
    def change_sep_value(self,sep_value):
        self.sep_value = sep_value
    
    def change_use_mesh_flag(self,use_mesh_flag):
        self.use_mesh_flag = use_mesh_flag

    def _cal_mesh(self,target_cluster:pd.DataFrame) -> dict():
        #target_clusterのモーメントの固有値を求める
        target_eigs = {}
        for atom in self.targets_atoms:
            target_eig = cal_eigenvalues(target_cluster.query(f'atom=="{atom}"')[['x','y','z']].values)
            target_eigs[atom] = (np.array(target_eig)//self.sep_value).astype(int)
        return target_eigs
    
    def __make_query_keys(self,mesh:dict()) -> list():
        mesh_lenge = [i for i in range(-self.offset,self.offset+1)]
        offsets = product(*[mesh_lenge] * 3)
        target_cells = [tuple(mesh + offset) for offset in offsets]
        return target_cells
    
    def _make_query_keys(self,mesh:dict()) -> list():
        return {atom:self.__make_query_keys(target_eig) for atom,target_eig in mesh.items()}

    def _get_target_databse_indexs(self,querys):
        target_indexs = {}
        for atom,querys_i in querys.items():
            target_indexs[atom] = [i for  query in querys_i if query in self.database_mesh_dict[atom].keys() for i in self.database_mesh_dict[atom][query]]
        
        target_indexs_set = set()
        for val in target_indexs.values():
            if len(target_indexs_set)==0:
                target_indexs_set = copy.deepcopy(set(val))
                continue
            target_indexs_set = target_indexs_set & set(val)
        return list(target_indexs_set)

    def cluster_calculate_features(self,clusterpath):
        #clusterpath:クラスターのファイルパス

        target_cluster = pd.read_csv(clusterpath,index_col=0)
        #target_clusterのモーメントの固有値を求める
        target_mesh = self._cal_mesh(target_cluster)

        #計算対象のメッシュをリストアップ,offsetsは計算対象のメッシュの周囲
        querys = self._make_query_keys(target_mesh)

        #計算対象のメッシュの周囲のクラスターの特徴量を取得する
        if self.use_mesh_flag:
            target_database_indexs = self._get_target_databse_indexs(querys)
        else:
            target_database_indexs = self.database_path_df.index.tolist()

        #計算対象のクラスターの特徴量を取得する
        target_database_coordinates = {key:val[target_database_indexs] for key,val in self.database_coordinates.items()}
        target_coordinates = {key:target_cluster.query(f"atom == @key")[['x','y','z']].to_numpy() for key in self.targets_atoms}
        dis = _cal_distance(target_coordinates,target_database_coordinates)

        if len(dis) == 0:
            #計算対象のクラスターがデータベースに存在しない場合
            self.distances_df = 'no match'
            return np.nan
        self.distances_df = self.log_format.loc[target_database_indexs].copy()
        self.distances_df['distance'] = dis
        return dis.min()

class CrystalFeatureCalculator(ClusterFeatureCalculator):
    def __init__(self, databasepath,method='mean',target_atoms=['Si1','O1'],reference=1e-8,sep_value=0.1,offset=5,use_mesh_flag=True):
        self.method = method
        super().__init__(databasepath,target_atoms,reference,sep_value,offset,use_mesh_flag)
    
    def process_target_cluster(self,target_cluster):
        features = self.cluster_calculate_features(target_cluster)
        return features,self.distances_df

    def calculate_features(self, crystalpath):
        self.target_clusters = glob.glob('{}/*_0.csv'.format(crystalpath))
        result = [self.process_target_cluster(target_cluster) for target_cluster in self.target_clusters]
        self.cluster_features,self.calculate_log = zip(*result)
        self.calculate_log = list(self.calculate_log)
        if self.method == 'mean':
            return np.mean(self.cluster_features, axis=0)
        elif self.method == 'max':
            return np.max(self.cluster_features)