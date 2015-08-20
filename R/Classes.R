setClass('FuzzyCluster',
         representation(data='matrix',U='matrix',W='matrix',centroids='matrix',phi='numeric',classes='integer',distance='character',
                                       alpha='numeric',`Ue_mean - Ue_req`='numeric',iterations='integer',pred_int='list'))


setClass('FuzzyClusterGroup',representation(clusters='list',confidence='numeric'))