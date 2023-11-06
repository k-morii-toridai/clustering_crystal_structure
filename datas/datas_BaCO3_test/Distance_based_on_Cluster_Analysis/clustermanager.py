
import os
import glob
import re
import pandas as pd
from itertools import combinations
import polars as pl
from .read_info import to_series


class ClusterManager:
    def __init__(self, cluster_list_df):
        self.cluster_list_df = cluster_list_df
        self.target_combination_df = None
    @classmethod
    def from_dirpath(cls, dirpath, dirs=False,ignore_dirs:list=None):
        cluster_list_df = cls.read_cluster_list(dirpath, dirs,ignore_dirs=ignore_dirs)
        return cls(cluster_list_df)
    
    @staticmethod
    def read_cluster_list(dirpath, dirs,ignore_dirs:list=None):
        dirlist = [os.path.join(dirpath,f) for f in os.listdir(dirpath) if os.path.isdir(os.path.join(dirpath, f))] if dirs else [dirpath]
        #ignore_dirsにあるディレクトリをリストから除外する
        if ignore_dirs is not None:
            dirlist = [dirpath for dirpath in dirlist if (dirpath) not in ignore_dirs]
        cluster_list = [path for dirpath in dirlist for path in glob.glob(os.path.join(dirpath,'*.csv'))]
        result_data = []

        for filepath in cluster_list:
            try:
                filename = os.path.basename(filepath)
                dirname = os.path.dirname(filepath)
                cifid, isite, _ = tuple(re.split('_', filename))
                result_data.append((cifid, dirname, int(isite)))
            except:
                #print(f'error:{filepath}')
                continue

        result_df = pd.DataFrame(result_data, columns=['cifid', 'address', 'isite']).drop_duplicates().sort_values(by=['cifid','isite']).reset_index(drop=True)
        return result_df
    
    def __make_target_files(self,target_combination_df):
        target_combination_df = target_combination_df.select(
                                                        pl.concat_str([
                                                            pl.col('address_i'),
                                                            pl.concat_str([
                                                                pl.concat_str([pl.col('cifid_i'),pl.col('isite_i'),pl.col('pattern_i')],separator='_'),    
                                                                pl.lit('csv')],separator='.')]
                                                        ,separator='/').alias('file_i'),
                                                        
                                                        pl.concat_str([
                                                            pl.col('address_j'),
                                                            pl.concat_str([
                                                                pl.concat_str([pl.col('cifid_j'),pl.col('isite_j'),pl.col('pattern_j')],separator='_'),    
                                                                pl.lit('csv')],separator='.')]
                                                        ,separator='/').alias('file_j'),                 
                                                        )
        return target_combination_df

    def __make_target_combination_df(self,df_i,df_j):
        df_i.columns = ['cifid_i', 'address_i', 'isite_i']
        df_j.columns = ['cifid_j', 'address_j', 'isite_j']
        df = pl.concat([df_i, df_j],how='horizontal')
        target_combinations = []
        
        #クラスターのすべての回転パターンを作成する
        dim=len(df_i)
        for i in range(12):
            target_combinations.append(pl.concat([df,pl.DataFrame([i]*dim,schema = ['pattern_j'])],how='horizontal'))
        target_combination_df = pl.concat(target_combinations,how="vertical")
        return pl.concat([target_combination_df,pl.DataFrame([0]*target_combination_df.height,schema = ['pattern_i'])],how='horizontal')

    def calculate_self_distance_file(self):
        combinations_list = list(combinations(range(len(self.cluster_list_df)), 2))
        df_i = pl.from_pandas(self.cluster_list_df.iloc[[i[0] for i in combinations_list]].reset_index(drop=True))
        df_j = pl.from_pandas(self.cluster_list_df.iloc[[i[1] for i in combinations_list]].reset_index(drop=True))
        
        
        target_combination_df = self.__make_target_combination_df(df_i,df_j)

        target_combination_files = self.__make_target_files(target_combination_df)
        self.target_combination_files = target_combination_files.to_pandas()
        self.target_combination_df = target_combination_df.to_pandas()
    
    def to_file_path(self,pattern=0):
        cluster_list_df = pl.from_pandas(self.cluster_list_df)
        cluster_path_list_df = cluster_list_df.select(
                                                pl.concat_str([
                                                    pl.col('address'),
                                                    pl.concat_str([
                                                        pl.concat_str([pl.col('cifid'),pl.col('isite'),pl.lit(f'{pattern}')],separator='_'),    
                                                        pl.lit('csv')],separator='.')]
                                                ,separator='/').alias('file_path')
                                                )
        self.cluster_path_list_df = cluster_path_list_df.to_pandas()
        return

class DCAFormatManager(ClusterManager):
    #DCAの入力フォーマットを管理,現状必要なフォーマットは2つ
    #一つの結晶を構成するクラスター同士の距離を計算するためのフォーマット
    #データベースと一つのクラスターの距離を計算するためのフォーマット
    def __init__(self, cluster_list_df):
        super().__init__(cluster_list_df)
    
    def __make_target_files(self,target_combination_df):
        target_combination_df = target_combination_df.select(
                                                        pl.concat_str([
                                                            pl.col('address_i'),
                                                            pl.concat_str([
                                                                pl.concat_str([pl.col('cifid_i'),pl.col('isite_i'),pl.col('pattern_i')],separator='_'),    
                                                                pl.lit('csv')],separator='.')]
                                                        ,separator='/').alias('file_i'),
                                                        
                                                        pl.concat_str([
                                                            pl.col('address_j'),
                                                            pl.concat_str([
                                                                pl.concat_str([pl.col('cifid_j'),pl.col('isite_j'),pl.col('pattern_j')],separator='_'),    
                                                                pl.lit('csv')],separator='.')]
                                                        ,separator='/').alias('file_j'),                 
                                                        )
        return target_combination_df

    def __make_target_combination_df(self,df_i,df_j):
        #.select(pl.col("Integer").alias("ABC"), pl.col("Integer").alias("DEF"))
        df_i = df_i.select(pl.col('address').alias('address_i'),pl.col('cifid').alias('cifid_i'),pl.col('isite').alias('isite_i'))
        df_j = df_j.select(pl.col('address').alias('address_j'),pl.col('cifid').alias('cifid_j'),pl.col('isite').alias('isite_j'))
        df = pl.concat([df_i, df_j],how='horizontal')
        target_combinations = []
        
        #クラスターのすべての回転パターンを作成する
        dim=len(df_i)
        for i in range(12):
            target_combinations.append(pl.concat([df,pl.DataFrame([i]*dim,schema = ['pattern_i'])],how='horizontal'))
        target_combination_df = pl.concat(target_combinations,how="vertical")
        return pl.concat([target_combination_df,pl.DataFrame([0]*target_combination_df.height,schema = ['pattern_j'])],how='horizontal')

    def calculate_self_distance_file(self):
        combinations_list = list(combinations(range(len(self.cluster_list_df)), 2))
        df_i = pl.from_pandas(self.cluster_list_df.iloc[[i[0] for i in combinations_list]].reset_index(drop=True))
        df_j = pl.from_pandas(self.cluster_list_df.iloc[[i[1] for i in combinations_list]].reset_index(drop=True))
        
        target_combination_df = self.__make_target_combination_df(df_i,df_j)

        target_combination_files = self.__make_target_files(target_combination_df)
        self.target_combination_files = target_combination_files.to_pandas()
        self.target_combination_df = target_combination_df.to_pandas()
    
    def calculate_cluster_featuring(self,clusterpath):
        df_i = pl.from_pandas(self.cluster_list_df)
        if isinstance(clusterpath, str):
            df_j = pl.from_pandas(pd.DataFrame([to_series(clusterpath)]*df_i.height))
        elif isinstance(clusterpath, pd.Series):
            df_j = pl.from_pandas(pd.DataFrame([clusterpath]*df_i.height))
        target_combination_df = self.__make_target_combination_df(df_i,df_j)
        target_combination_files = self.__make_target_files(target_combination_df)
        self.target_combination_files = target_combination_files.to_pandas()
        self.target_combination_df = target_combination_df.to_pandas()