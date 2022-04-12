"""
Basic tests
"""

from longest_gene import find_stop_codon


def test_find_stop_codon() -> None:
    """test find_stop_codon"""
    assert find_stop_codon("TAA") == 0
    assert find_stop_codon("AAATAA") == 3
    assert find_stop_codon("NTAA") is None
    assert find_stop_codon("NNTAA") is None
    assert find_stop_codon("NNNTA") is None
    assert find_stop_codon("NNT") is None
