import numpy as np
import pytest

from karstpy.karyotypes import format_karyotypes, parse_karyotypes


def test_integer_round_trip():
    matrix = np.array([[2, 1, 4], [0, 3, 2]], dtype=np.int32)
    np.testing.assert_array_equal(parse_karyotypes(format_karyotypes(matrix)), matrix)
    assert format_karyotypes(parse_karyotypes(["02.01.4"])) == ["2.1.4"]


@pytest.mark.parametrize("labels", [[], ["2..3"], ["2.2."], ["+2.2"], [" 2.2"], ["2.2", "2.02"], ["2", "2.2"]])
def test_bad_labels(labels):
    with pytest.raises(ValueError):
        parse_karyotypes(labels)
