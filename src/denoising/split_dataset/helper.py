# MIT License

# Copyright (c) 2019 Chan Zuckerberg Biohub

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Copied from https://github.com/czbiohub/molecular-cross-validation/blob/master/src/molecular_cross_validation/util.py


from typing import Tuple

import numpy as np

def split_molecules(
    umis: np.ndarray,
    data_split: float,
    overlap_factor: float = 0.0,
    random_state: np.random.RandomState = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Splits molecules into two (potentially overlapping) groups.
    :param umis: Array of molecules to split
    :param data_split: Proportion of molecules to assign to the first group
    :param overlap_factor: Overlap correction factor, if desired
    :param random_state: For reproducible sampling
    :return: umis_X and umis_Y, representing ``split`` and ``~(1 - split)`` counts
             sampled from the input array
    """
    if random_state is None:
        random_state = np.random.RandomState()

    umis_X_disjoint = random_state.binomial(umis, data_split - overlap_factor)
    umis_Y_disjoint = random_state.binomial(
        umis - umis_X_disjoint, (1 - data_split) / (1 - data_split + overlap_factor)
    )
    overlap_factor = umis - umis_X_disjoint - umis_Y_disjoint
    umis_X = umis_X_disjoint + overlap_factor
    umis_Y = umis_Y_disjoint + overlap_factor

    return umis_X, umis_Y