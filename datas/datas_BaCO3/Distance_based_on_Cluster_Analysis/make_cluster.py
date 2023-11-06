from .read_info import Set_Cluster_Info
import pickle,os,re
import pandas as pd
import traceback
def make_cluster_dataset(cifid=None, atom='Si', nn_data_address=None, nn_data=None, adjacent_num=None,
                         isite=None, cluster_address=None, cluster_df=None, rotation=True, outdir=None):
    cwd = os.getcwd()

    try:
        if nn_data_address is not None:
            with open(nn_data_address, 'rb') as f:
                nn_data = pickle.load(f)
            if cifid is None:
                cifid = re.sub('\.pickle', '', os.path.basename(nn_data_address))

        if cluster_address is not None:
            cluster_df = pd.read_csv(cluster_address, index_col=0)
            if cifid is None:
                cifid = re.sub('_[0-9]*_[0-9]*\.csv', '', os.path.basename(cluster_address))

        if cifid is None:
            print('please enter cifid')
            return

        if outdir is not None:
            os.chdir(outdir)

        if isite is not None or cluster_df is not None:
            if cluster_df is not None:
                cluster_info = Set_Cluster_Info(cluster_df=cluster_df)
                isite = cluster_df.iloc[0].isite
            else:
                if adjacent_num is None:
                    print('please enter adjacent_num')
                    return
                cluster_info = Set_Cluster_Info(isite, nn_data, adjacent_num)
            
            for pattern in range(len(cluster_info.shaft_comb)):
                cluster_info.parallel_shift_of_center()
                cluster_info.rotation(pattern=pattern)
                cluster_info.cluster_coords.to_csv('{}_{}_{}.csv'.format(cifid, isite, pattern))
                if not rotation:
                    break
        else:
            if adjacent_num is None:
                print('please enter adjacent_num')
                return
            
            #alllen_ = sum(1 for isite_atom in nn_data.keys() if re.split(r'([a-zA-Z]+)', nn_data[isite_atom][0][0])[1] == atom)
            #cont = 0
            
            for isite in nn_data.keys():
                isite_atom = re.split(r'([a-zA-Z]+)', nn_data[isite][0][0])[1]
                if isite_atom == atom:
                    if adjacent_num is None:
                        print('please enter adjacent_num')
                        return
                    cluster_info = Set_Cluster_Info(isite, nn_data, adjacent_number=adjacent_num)
                    #alllen = alllen_ * len(cluster_info.shaft_comb)
                    
                    for pattern in range(len(cluster_info.shaft_comb)):
                        cluster_info.parallel_shift_of_center()
                        cluster_info.rotation(pattern=pattern)
                        cluster_info.cluster_coords.to_csv('{}_{}_{}.csv'.format(cifid, isite, pattern))
                        if not rotation:
                            break
    except Exception:
        traceback.print_exc()
        
    os.chdir(cwd)
