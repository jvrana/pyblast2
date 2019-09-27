from pyblast.utils import Span


def test_lengths():
    s = Span(8, 10, 10, cyclic=True, allow_wrap=True, index=0)
    assert len(s) == 2
    assert s.b == 10
    # 8, 9

    s = Span(8, 11, 10, cyclic=True, allow_wrap=True, index=1)
    assert len(s) == 3
    assert s.b == 0
    # 8, 9

    s = Span(8, 12, 10, cyclic=True, allow_wrap=True, index=1)
    assert len(s) == 4
    assert s.b == 1
    # 8, 9, 10, 11
