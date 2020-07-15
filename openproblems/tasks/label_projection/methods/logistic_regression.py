import sklearn

def logistic_regression(adata):
    classifier = sklearn.linear_model.LogisticRegression()

    adata_train = adata[adata.obs['is_train']]
    adata_test = adata[~adata.obs['is_train']]

    classifier.fit(adata_train.X, adata_train.obs['labels'])
    adata_test['labels_pred'] = classifier.predict(adata_test.X)

    return adata_test
