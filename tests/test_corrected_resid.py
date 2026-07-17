import types
import numpy as np
import pytest
from openawsem.openAWSEM import OpenMMAWSEMSystem


def _call(resi, chain_starts, chain_ends, gap):
    mock = types.SimpleNamespace(resi=resi, chain_starts=chain_starts, chain_ends=chain_ends)
    return OpenMMAWSEMSystem.corrected_resid(mock, gap=gap)


@pytest.mark.parametrize("resi, chain_starts, chain_ends, gap, expected", [
    # single chain: no change
    ([0, 0, 1, 1, 2, 2],       [0],    [2],    10, [0,  0,  1,  1,  2,  2 ]),
    # two chains: chain B shifted by gap
    ([0, 0, 1, 1, 2, 2, 3, 3], [0, 2], [1, 3], 10, [0,  0,  1,  1, 12, 12, 13, 13]),
    # three chains: cumulative offsets (+0, +gap, +2*gap)
    ([0, 1, 2, 3, 4, 5],       [0, 2, 4], [1, 3, 5], 10, [0,  1, 12, 13, 24, 25]),
    # non-protein atoms (-1) are preserved
    ([-1, 0, 1, -1, 2, 3],     [0, 2], [1, 3], 10, [-1, 0,  1, -1, 12, 13]),
    # different gap value
    ([0, 1, 2, 3],              [0, 2], [1, 3],  5, [0,  1,  7,  8]),
])
def test_corrected_resid(resi, chain_starts, chain_ends, gap, expected):
    result = _call(resi, chain_starts, chain_ends, gap)
    np.testing.assert_array_equal(result, expected)
