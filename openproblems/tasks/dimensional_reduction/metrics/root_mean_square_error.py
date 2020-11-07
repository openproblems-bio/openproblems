import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import spatial
import pandas as pd
from ....tools.decorators import metric


@metric(metric_name="root mean squared error", maximize=True)
def rmse(adata, method):
    
    """Calculate the root mean squared error (RMSE) between the full (or processed) data matrix and a list of dimensionally-reduced matrices."""
    
    dimensional_reduction_method = "X_" + method
    
    kruskel_matrix = "kruskel_matrix_" + method
    kruskel_score = "kruskel_score_" + method
    RMSE_calculation = "RMSE_" + method

    kruskel_matrices, kruskel_scores, rmse_list = [], [], []
    
    if method == 'pca':
        
        dimensions_to_evaluate = [1, 2, 5, 10, 25, 30, 50]
    
        for i in dimensions_to_evaluate:

            matrix, kruskel_score, rmse_score = get_rmse(
                adata, adata.X, adata.obsm[dimensional_reduction_method][:, 0:i]
            )

            kruskel_matrices.append(matrix)
            kruskel_scores.append(kruskel_score)
            rmse_list.append(rmse_score)

        for i, j in enumerate(dimensions_to_evaluate):

            kruskal_matrix_adata = "kruskal_matrix_" + method + str(j) + "d"
            adata.obsp[kruskal_matrix_adata] = kruskel_matrices[i]
            
            kruskal_score_adata = "kruskal_score_" + method + str(j) + "d"
            adata.uns[kruskal_score_adata] = kruskel_scores[i]
            
            rmse_calculation_adata = "rmse_calculation_" + method + str(j) + "d"
            adata.uns[rmse_calculation_adata] = rmse_list[i]
            
        return adata
    
    if method == "umap" or "tsne":
 
        adata.obsp[kruskel_matrix], adata.uns[kruskel_score], adata.uns[RMSE_calculation] = get_rmse(
                adata, adata.X, adata.obsm[dimensional_reduction_method]
            )

        return adata
    
def get_rmse(adata, high_dimensional_data_matrix, low_dimensional_data_matrix):

    """
    Calculate kruskel's stress. 
    """    

    high_dimensional_data_matrix = pd.DataFrame(high_dimensional_data_matrix)
    low_dimensional_data_matrix = pd.DataFrame(low_dimensional_data_matrix)
    number_of_points = adata.shape[0]

    points = range(number_of_points)

    high_dimensional_distance_matrix = calculate_squareform_pairwise_distance(
        high_dimensional_data_matrix, points
    )
    low_dimensional_distance_matrix = calculate_squareform_pairwise_distance(
        low_dimensional_data_matrix, points
    )

    diff = high_dimensional_distance_matrix - low_dimensional_distance_matrix

    kruskel_matrix = np.sqrt(diff ** 2 / sum(low_dimensional_distance_matrix ** 2))

    kruskel_score = np.sqrt(sum(diff ** 2) / sum(low_dimensional_distance_matrix ** 2))

    y_actual = high_dimensional_distance_matrix
    y_predic = low_dimensional_distance_matrix

    from sklearn.metrics import mean_squared_error
    from math import sqrt

    rms = sqrt(mean_squared_error(y_actual, y_predic))

    return kruskel_matrix, kruskel_score, rms

def calculate_squareform_pairwise_distance(data, points):

    """
    Calculate the pairwise distance between points in a matrix / vector and then format this into a squareform vector.
    """
    
    df = pd.DataFrame(
        sp.spatial.distance.squareform(sp.spatial.distance.pdist(data.loc[points])),
        columns=points,
        index=points,
    )

    return df

