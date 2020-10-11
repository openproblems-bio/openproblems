import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import spatial
import pandas as pd
from ....tools.decorators import metric


@metric(metric_name="compare_rmse", maximize=True)
def compare_RMSE(adata, method):
    
    """Calculate the root mean square error (RMSE) between the full (or processed) data matrix and a list of dimensionally-reduced matrices."""
    
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

def fig_presets_plot(ax, title, x, y, xlab, ylab):

    ax.set_xlabel(xlab, fontsize=15)
    ax.set_ylabel(ylab, fontsize=15)
    ax.set_title(title, fontsize=15)
    ax.plot(
        x, y, lw=2, marker="o", color="mediumpurple",
    )
    ax.spines["left"].set_linewidth(3)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    return ax

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

def plot_rmse(dims, rms_all):
    
    fig = plt.figure()  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    fig_presets_plot(
        ax,
        "RMSE as a function of dimensional representation",
        dims,
        rms_all,
        xlab="Dimensions",
        ylab="RMSE",
    )

    fig.tight_layout(pad=0)
    
def kruskel_plot(data_red, kmat, method):

    import seaborn as sns
    import matplotlib.pyplot as plt

    ax = sns.jointplot(
        x=data_red[0],
        y=data_red[1],
        data=kmat,
        kind="kde",
        color="mediumpurple",
        height=6,
    )
    ax.plot_joint(plt.scatter, c="w", s=2, linewidth=1, alpha=0.85)
    ax.set_axis_labels(str(method) + "-1", str(method) + "-2")
    plt.title(str(method), pad=80)
    
def plot_rmse_comparison(rmse_values, methods_axis):

    # # set the style of the axes and the text color
    plt.rcParams['axes.edgecolor']='#333F4B'
    plt.rcParams['axes.linewidth']=0.8
    plt.rcParams['xtick.color']='#333F4B'
    plt.rcParams['ytick.color']='#333F4B'
    plt.rcParams['text.color']='#333F4B'
    fig, ax = plt.subplots(figsize=(5,3.5))
    plt.plot(rmse_values, methods_axis, "o", markersize=5, color='#007ACC', alpha=0.6)
    plt.hlines(y=methods_axis, xmin=0, xmax=rmse_values, color='#007ACC', alpha=0.2, linewidth=5)

    # set labels
    ax.set_xlabel('RMSE Value', fontsize=15, fontweight='black', color = '#333F4B')
    ax.set_ylabel('')

    # # set axis
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # add an horizonal label for the y axis 
    fig.text(-0.23, 0.96, 'RMSE Comparison', fontsize=15, fontweight='black', color = '#333F4B')

    # change the style of the axis spines
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)

    # set the spines position
    ax.spines['bottom'].set_position(('axes', -0.04))
    ax.spines['left'].set_position(('axes', 0.015))