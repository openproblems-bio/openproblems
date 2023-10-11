from openproblems.datasets.datasets.dataset import Dataset
from openproblems.datasets.functions import SizedSampler, GroupBalancedSampler
from openproblems.datasets.functions import RemoveEmpty
from openproblems.datasets.metadata import Metadata


__all__ = ["Zebrafish"]


# TODO: metaclass will handle custom functions
class Zebrafish(Dataset):
    METADATA = Metadata(
        name="Zebrafish",
        author="Foo",
        year=2016,
        url="https://ndownloader.figshare.com/files/24566651?private_link=e3921450ec1bd0587870"
    )
    _SAMPLER = (
        GroupBalancedSampler(keys="lab", sampler=SizedSampler(500, cells=True)) >>
        RemoveEmpty() >>
        SizedSampler(100, cells=False) >>
        RemoveEmpty()
    )
    # it either inherits parent's processor/verifier/downloader or completely overrides it
