from pyblast.utils import Span


def test_lenghts():
    s = Span(8, 11, 10, cyclic=True, allow_wrap=True, index=0)
    assert len(s) == 3
    assert s.b == 1
    # 8, 9, 0

    s = Span(8, 11, 10, cyclic=True, allow_wrap=True, index=1)
    assert len(s) == 3
    assert s.b == 11
    # 8, 9, 10

    s = Span(8, 12, 10, cyclic=True, allow_wrap=True, index=1)
    assert len(s) == 4
    assert s.b == 2
    # 8, 9, 10, 11