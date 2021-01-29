from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np
import scprep

_alra = scprep.run.RFunction(
    setup="""
    library(SingleCellExperiment)
    library(Matrix)
    library(rsvd)
    source('https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R')
    
    """,
    args="sce",
    body="""
    assay(sce, "X") <- as(assay(sce, "X"), "CsparseMatrix")
    reducedDim(sce, "train") <- as(reducedDim(sce, "train"), "CsparseMatrix")
    completed <- alra(as.matrix(reducedDim(sce, "train")))
    normCompleted <- completed[[3]]
    reducedDim(sce, "train") <- as.matrix(normCompleted)
    sce
    """,
)

@method(
    method_name="ALRA",
    paper_name="Zero-preserving imputation of scRNA-seq data using low-rank approximation",
    paper_url="https://www.biorxiv.org/content/10.1101/397588v1",
    paper_year=2018,
    code_url="https://github.com/KlugerLab/ALRA",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",
)

def alra(adata):
    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )
    X = scprep.transform.sqrt(X) #note that log transform is used in the examples of alra on github
    adata.obsm['train']=X #we know that the scprep's .run.Rfunction uses anndata.2ri, where reducedDim(d, "PCA")=adata.obsm["X_PCA"], but that these dim-reduction equivalents 
    #do not encompass any given key hold true across all contexts. Thus, we use X_PCA for now to ensure compatability. 
    Y = _alra(X)
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    Y = alra_calc(adata)
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"]=Y
    return adata
