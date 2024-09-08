import yaml
from collections import namedtuple


def to_site_donor(data):
    df = data.obs['batch'].copy().to_frame().reset_index()
    df.columns = ['index','batch']
    df['site'] = df['batch'].apply(lambda x: x[:2])
    df['donor'] = df['batch'].apply(lambda x: x[2:]) 
    return df


def split(tr1, tr2, fold):
    df = to_site_donor(tr1) 
    mask = df['site'] == f's{fold+1}'
    maskr = ~mask

    Xt = tr1[mask].layers["normalized"].toarray()
    X = tr1[maskr].layers["normalized"].toarray()

    yt = tr2[mask].layers["normalized"].toarray()
    y = tr2[maskr].layers["normalized"].toarray()

    print(f"{X.shape}, {y.shape}, {Xt.shape}, {yt.shape}")

    return X,y,Xt,yt


def load_yaml(path):
    with open(path) as f:
        x = yaml.safe_load(f)
    res = {}
    for i in x:
        res[i] = x[i]['value']
    config = namedtuple('Config', res.keys())(**res)
    print(config)
    return config
