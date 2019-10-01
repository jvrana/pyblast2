from pyblast.utils import Span
import pytest
import numpy as np


class TestInit:
    @pytest.mark.parametrize("a", [-2, -1, 0, 1, 5, 10, 19, 20, 21, 22])
    @pytest.mark.parametrize("b", [-2, -1, 0, 1, 5, 10, 19, 20, 21, 22])
    @pytest.mark.parametrize("index", [-1, 0, 1])
    def test_linear(self, a, b, index):
        """Any linear span outside of bounds or with a > b is invalid"""
        length = 20
        bounds = [index, length + index]
        if a < bounds[0] or a >= bounds[1] or b > bounds[1] or a > b:
            with pytest.raises(IndexError):
                Span(a, b, length, cyclic=False, index=index)
        else:
            s = Span(a, b, length, cyclic=False, index=index)
            assert s.a == a
            assert s.b == b

    def test_cyclic_basic(self):
        s = Span(5, 10, 20, cyclic=True)
        assert s.a == 5
        assert s.b == 10

        s = Span(15, 5, 20, cyclic=True)
        assert s.a == 15
        assert s.b == 5

        s = Span(-1, 5, 20, cyclic=True)
        assert s.a == 19
        assert s.b == 5

    def test_cyclic_edge_case(self):
        s = Span(5, 10, 10, cyclic=True)
        s.a == 5
        s.b == 10

        s = Span(5, 0, 10, cyclic=True)
        assert s.a == 5
        assert s.b == 0

    def test_cyclic_beyond_right_bound(self):
        """Endpoints outside of right bound should be translated."""

        # index 0
        s = Span(5, 10, 10, cyclic=True, index=0)
        assert s.b == 10

        s = Span(5, 11, 10, cyclic=True, index=0)
        assert s.b == 1

        s = Span(5, 12, 10, cyclic=True, index=0)
        assert s.b == 2

        # index 1
        s = Span(5, 10, 10, cyclic=True, index=1)
        assert s.b == 10

        s = Span(5, 11, 10, cyclic=True, index=1)
        assert s.b == 11

        s = Span(5, 12, 10, cyclic=True, index=1)
        assert s.b == 2

        s = Span(5, 13, 10, cyclic=True, index=1)
        assert s.b == 3

        # index 2
        s = Span(5, 12, 10, cyclic=True, index=2)
        assert s.b == 12

        s = Span(5, 13, 10, cyclic=True, index=2)
        assert s.b == 3

    def test_cyclic_beyond_left_bound(self):
        """Start points beyond left bound shoulde be translated"""
        s = Span(5, -1, 10, cyclic=True, index=0)
        assert s.b == 10
        s = Span(5, 0, 10, cyclic=True, index=1)
        assert s.b == 11
        s = Span(5, 1, 10, cyclic=True, index=2)
        assert s.b == 12

        s = Span(-1, 5, 10, cyclic=True, index=0)
        assert s.a == 9
        assert s.b == 5

        s = Span(0, 5, 10, cyclic=True, index=1)
        assert s.a == 10
        assert s.b == 5

        s = Span(1, 5, 10, cyclic=True, index=2)
        assert s.a == 11
        assert s.b == 5

    @pytest.mark.parametrize(
        "abi",
        [
            (5, 20, 0, False),
            (5, 21, 0, True),
            (5, 22, 1, True),
            (5, 21, 1, False),
            (-1, 5, 0, True),
            (0, 5, 1, True),
            (1, 5, 2, True),
        ],
    )
    def test_cyclic_invalid_when_strict(self, abi):
        a, b, index, does_raise = abi
        if does_raise:
            with pytest.raises(IndexError):
                Span(a, b, 20, cyclic=True, index=index, strict=True)
        else:
            Span(a, b, 20, cyclic=True, index=index, strict=True)
        assert Span(a, b, 20, cyclic=True, index=index, strict=False)

    @pytest.mark.parametrize("a", range(0, 10, 3))
    @pytest.mark.parametrize("b", range(0, 10, 3))
    @pytest.mark.parametrize("cyclic", [True, False])
    def test_init(self, a, b, cyclic):
        length = 20
        index = 0
        if not cyclic:
            if a > b or a >= length or b > length or a < index or b < index:
                with pytest.raises(IndexError):
                    Span(a, b, length, cyclic)
        else:
            s = Span(a, b, length, cyclic)

            assert s.a == a
            assert s.b == b
            assert s.context_length == length
            if cyclic and a > b:
                assert len(s) == length - a + b
            else:
                assert len(s) == b - a

    def test_linear_init_should_raise(self):
        """start is beyond bounds, so IndexError should be raised"""
        with pytest.raises(IndexError):
            Span(10, 10, 10, cyclic=False)

    @pytest.mark.parametrize("x", [0, 5, 6, 7, 8, 9])
    @pytest.mark.parametrize("index", [0], ids=["index=0"])
    def test_special_case_empty_cyclic(self, x, index):
        """Special case of empty span should be valid for any endpoints that are the same"""
        s = Span(x, x, 10, cyclic=True, index=index)
        assert len(s) == 0
        assert s.a == x
        assert s.b == x
        assert s.c == x

    @pytest.mark.parametrize("x", [1, 5, 6, 7, 8, 9, 10])
    @pytest.mark.parametrize("index", [1], ids=["index=0"])
    def test_special_case_empty_cyclic_index1(self, x, index):
        """Special case of empty span should be valid for any endpoints that are the same"""
        s = Span(x, x, 10, cyclic=True, index=index)
        assert len(s) == 0
        assert s.a == x
        assert s.b == x
        assert s.c == x

    def test_special_case_empty_cyclic__explicit(self):
        # index=0
        s = Span(0, 0, 10, cyclic=True, index=0)
        assert s.a == 0
        assert s.b == 0
        assert s.c == s.b

        s = Span(10, 10, 10, cyclic=True, index=0)
        assert s.a == 0
        assert s.b == 0
        assert s.c == s.b

        s = Span(20, 20, 10, cyclic=True, index=0)
        assert s.a == 0
        assert s.b == 0
        assert s.c == s.b

        # index=1
        s = Span(10, 10, 10, cyclic=True, index=1)
        assert s.a == 10
        assert s.b == 10
        assert s.c == s.b

        s = Span(11, 11, 10, cyclic=True, index=1)
        assert s.a == 1
        assert s.b == 1
        assert s.c == s.b

        s = Span(21, 21, 10, cyclic=True, index=1)
        assert s.a == 1
        assert s.b == 1
        assert s.c == s.b

        s = Span(0, 0, 10, cyclic=True, index=1)
        assert s.a == 10
        assert s.b == 10
        assert s.c == s.b

        s = Span(10, 10, 10, cyclic=True, index=1)
        assert s.a == 10
        assert s.b == 10
        assert s.c == s.b

    @pytest.mark.parametrize("x", [-1, 0, 1, 5, 10, 11])
    def test_special_case_empty_linear(self, x):
        """Special case of empty span should be valid for any endpoints that are the same"""
        if x < 0 or x >= 10:
            with pytest.raises(IndexError):
                Span(x, x, 10, cyclic=False)
        else:
            s = Span(x, x, 10, cyclic=False)
            assert s.a == x
            assert s.b == x
            assert s.c == x
            assert len(s) == 0

    def test_init_should_raise2(self):
        """Anything outside the bounds while strict should raise exception"""
        with pytest.raises(IndexError):
            Span(9408, 4219, 9408, True, strict=True)

    def test_init_linear(self):
        """Basic constructor for linear span"""
        s = Span(10, 80, 100, False)
        assert s.a == 10
        assert s.b == 80
        assert s.context_length == 100
        assert s.cyclic is False

    def test_init_cyclic(self):
        """Basic constructor for cyclic spans"""
        s = Span(10, 5, 100, True)
        assert s.a == 10
        assert s.b == 5
        assert s.context_length == 100
        assert s.cyclic

    def test_init_linear_raises(self):
        """Linear span raises error when a > b"""
        Span(10, 5, 100, True)
        with pytest.raises(IndexError):
            Span(10, 5, 100, False)

    # TODO: handle wraps better
    def test_invalid_cyclic(self):
        """When strict is True, any indices outside bounds is invalid."""
        print(Span(0, 10, 10, True))
        with pytest.raises(IndexError):
            print(Span(0, 10, 9, True, strict=True))
        print(Span(0, 10, 9, True, ignore_wrap=False))

    class TestInitWraps(object):
        def test_wrapped(self):
            s = Span(5, 9, 20, cyclic=True)
            assert s.a == 5
            assert s.b == 9
            assert s.c == 9
            assert len(s) == 4

            s = Span(5, 29, 20, cyclic=True)
            assert s.a == 5
            assert s.b == 9
            assert s.c == 29
            assert len(s) == 24

            s = Span(5, 49, 20, cyclic=True)
            assert s.a == 5
            assert s.b == 9
            assert s.c == 49
            assert len(s) == 44

            s = Span(25, 49, 20, cyclic=True)
            assert s.a == 5
            assert s.b == 9
            assert s.c == 49 - 20
            assert len(s) == 24

        def test_negative_wrapped(self):
            s = Span(-40, 5, 20, cyclic=True)
            assert len(s) == 45
            assert s.a == 0
            assert s.b == 5

        def test_wrapped_when_allow_wrapping_is_false(self):
            # 5, 6, 7, 8, 9, 0
            s = Span(5, 11, 10, cyclic=True, ignore_wrap=True)
            assert s.a == 5
            assert s.b == 1
            assert s.c == 1
            assert len(s) == 6

            # 5, 6, 7, 8, 9, 0, 1
            s = Span(5, 12, 10, cyclic=True, ignore_wrap=True)
            assert s.a == 5
            assert s.b == 2
            assert s.c == 2
            assert len(s) == 7

        def test_wrapped_reverse(self):
            s = Span(15, 11, 10, cyclic=True)
            assert s.a == 5
            assert s.b == 1
            assert len(s) == 6

        def test_invalid_wrapping(self):
            assert Span(15, 11, 10, cyclic=True, index=0)
            with pytest.raises(IndexError):
                Span(15, 11, 10, cyclic=True, index=1)
            with pytest.raises(IndexError):
                Span(15, 10, 10, cyclic=True, index=0)

        def test_complete_wrap_plus_one(self):
            s = Span(0, 11, 10, cyclic=True)
            assert s.a == 0
            assert s.b == 1
            assert s.c == 11
            assert len(s) == 11


class TestRanges(object):
    def test_linear_range(self):
        s = Span(5, 10, 20)
        assert list(s) == [5, 6, 7, 8, 9]

    def test_empty_range(self):
        s = Span(5, 5, 20)
        assert list(s) == []

    def test_cyclic_range(self):
        s = Span(8, 10, 10, cyclic=True)
        assert list(s) == [8, 9]

        s = Span(8, 0, 10, cyclic=True)
        assert list(s) == [8, 9]

        s = Span(8, 2, 10, cyclic=True)
        assert list(s) == [8, 9, 0, 1]

    def test_cyclic_range_beyond_rbound(self):

        s = Span(8, 11, 10, cyclic=True)
        assert list(s) == [8, 9, 0]

        s = Span(8, 12, 10, cyclic=True)
        assert list(s) == [8, 9, 0, 1]

        s = Span(8, 12, 10, cyclic=True, index=1)
        assert list(s) == [8, 9, 10, 1]

        s = Span(8, 13, 10, cyclic=True, index=2)
        assert list(s) == [8, 9, 10, 11, 2]

    def test_cyclic_range_beyond_lbound(self):

        s = Span(0, 4, 10, cyclic=True)
        assert list(s) == [0, 1, 2, 3]

        s = Span(-1, 4, 10, cyclic=True)
        assert list(s) == [9, 0, 1, 2, 3]

        s = Span(0, 4, 10, cyclic=True, index=1)
        assert list(s) == [10, 1, 2, 3]

        s = Span(0, 4, 10, cyclic=True, index=2)
        assert list(s) == [10, 11, 2, 3]

    def test_complete(self):
        s = Span(0, 10, 10, cyclic=True)
        assert list(s) == list(range(0, 10))

        s = Span(0, 10, 10, cyclic=False)
        assert list(s) == list(range(0, 10))

    def test_complete_wrapped(self):
        s = Span(10, 20, 10, cyclic=True)
        assert list(s) == list(range(0, 10))
        assert len(s) == 10

        s = Span(20, 30, 10, cyclic=True)
        assert list(s) == list(range(0, 10))
        assert len(s) == 10

        s = Span(11, 21, 10, cyclic=True, index=1)
        assert list(s) == list(range(1, 11))
        assert len(s) == 10

        s = Span(21, 31, 10, cyclic=True, index=1)
        assert list(s) == list(range(1, 11))
        assert len(s) == 10

    def test_wrapped_twice(self):

        s = Span(5, 15, 10, cyclic=True)
        assert len(s) == 10
        assert list(s) == list(range(5, 10)) + list(range(5))

        s = Span(5, 25, 10, cyclic=True)
        assert list(s) == list(range(5, 10)) + list(range(10)) + list(range(5))
        assert len(s) == 20
        assert s.a == 5
        assert s.b == 5
        assert s.c == 25

    def test_wrapped_twice_index(self):
        s = Span(5, 15, 10, cyclic=True, index=1)
        assert len(s) == 10
        assert list(s) == list(range(5, 11)) + list(range(1, 5))

        s = Span(5, 25, 10, cyclic=True, index=1)
        assert list(s) == list(range(5, 11)) + list(range(1, 11)) + list(range(1, 5))
        assert len(s) == 20
        assert s.a == 5
        assert s.b == 5
        assert s.c == 25

        s = Span(5, 26, 10, cyclic=True, index=1)
        assert list(s) == list(range(5, 11)) + list(range(1, 11)) + list(range(1, 6))
        assert len(s) == 21
        assert s.a == 5
        assert s.b == 6
        assert s.c == 26

    def test_wrapped_twice_at_index(self):
        s = Span(0, 10, 10, cyclic=True)
        assert list(s) == list(range(0, 10))

        s = Span(0, 11, 10, cyclic=True)
        assert list(s) == list(range(0, 10)) + [0]

        s = Span(0, 20, 10, cyclic=True)
        assert list(s) == list(range(0, 10)) + list(range(0, 10))

        s = Span(0, 21, 10, cyclic=True)
        assert list(s) == list(range(0, 10)) + list(range(0, 10)) + [0]

    def test_cyclic_wrapped_beyond_rbound(self):
        s = Span(9, 15, 10, cyclic=True)
        assert list(s) == [9, 0, 1, 2, 3, 4]

        s = Span(5, 16, 10, cyclic=True)
        assert list(s) == [5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5]

        s = Span(5, 15, 10, cyclic=True)
        assert list(s) == [5, 6, 7, 8, 9, 0, 1, 2, 3, 4]


def test_t():
    """Test translating positions"""
    s = Span(0, 10, 10, cyclic=True)

    assert s.t(0) == 0
    assert s.t(1) == 1
    assert s.t(10) == 0
    assert s.t(-1) == 9


def test_t_index5():
    """Test translating positions"""
    s = Span(0, 10, 10, cyclic=True, index=5)

    assert s.t(0) == 5
    assert s.t(1) == 6
    assert s.t(9) == 14
    assert s.t(10) == 5
    assert s.t(-1) == 14


@pytest.mark.parametrize(
    ("a", "b", "does_span"),
    [(100, 800, False), (0, 100, False), (500, 1000, False), (800, 100, True)],
)
def test_spans_origin(a, b, does_span):
    s = Span(a, b, 1000, True)
    assert s.spans_origin() == does_span


def test_len_linear():
    s = Span(5, 20, 100, False)
    assert len(s) == 15


@pytest.mark.parametrize("b", [5, 9, 10, 11, 12])
@pytest.mark.parametrize("index", [-2, 0, 5])
def test_len_cyclic(b, index):
    s = Span(5, b, 10, cyclic=True, index=index)
    assert len(s) == b - 5


def test_iter_linear():
    s = Span(10, 20, 100, False)
    assert list(s) == list(range(10, 20))


def test_iter_cyclic():
    s = Span(90, 10, 100, True)
    assert list(s) == list(range(90, 100)) + list(range(10))


def test_str():
    print(Span(90, 10, 100, True))


def test_eq():
    """Tests for span equality."""
    s1 = Span(10, 90, 100, True)
    s2 = Span(10, 90, 100, True)
    assert s1 == s2

    s1 = Span(10, 91, 100, True)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2

    s1 = Span(11, 90, 100, True)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2

    s1 = Span(10, 90, 100, False)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2

    s1 = Span(10, 90, 100, True)
    s2 = Span(10, 90, 101, True)
    assert not s1 == s2

    s1 = Span(10, 190, 100, True)
    s2 = Span(10, 90, 100, True)
    assert not s1 == s2
    assert s1.a == s2.a
    assert s2.b == s2.b
    assert s1.c != s2.c


class TestContains:
    @pytest.mark.parametrize("i", range(-20, 20))
    def test_contains_index(self, i):
        s = Span(10, 50, 100, True)
        if 10 <= i < 50:
            assert i in s
        else:
            assert not i in s

    class TestContainsCyclic:
        @pytest.mark.parametrize("i", range(90, 100))
        def test_contains_cyclic_1(self, i):
            s = Span(90, 10, 100, True)
            assert i in s

        @pytest.mark.parametrize("i", range(0, 10))
        def test_contains_cyclic_2(self, i):
            s = Span(90, 10, 100, True)
            assert i in s

        @pytest.mark.parametrize("i", range(85, 90))
        def test_not_contains_cyclic_1(self, i):
            s = Span(90, 10, 100, True)
            assert i not in s

        @pytest.mark.parametrize("i", range(11, 20))
        def test_not_contains_cyclic_2(self, i):
            s = Span(90, 10, 100, True)
            assert i not in s

    class TestContainsSpan:
        def test_contains_region(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(20, 30, 100, True)
            assert s2 in s1
            assert not s1 in s2

        def test_contains_empty(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(30, 30, 100, True)
            assert s2 in s1
            assert not s1 in s2

        def test_contains_empty2(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(10, 10, 100, True)
            assert s2 in s1
            assert not s1 in s2

        def test_contains_self(self):
            s1 = Span(10, 50, 100, True)
            assert s1 in s1

        def test_does_not_contain(self):
            s1 = Span(10, 50, 100, True)
            s2 = Span(11, 50, 100, True)
            assert s1 not in s2
            assert s2 in s1

            s1 = Span(10, 50, 100, True)
            s2 = Span(10, 49, 100, True)
            assert s1 not in s2
            assert s2 in s1

        def test_cyclic_contains(self):
            s1 = Span(80, 20, 100, True)
            s2 = Span(85, 10, 100, True)
            assert s2 in s1
            assert s1 not in s2

        def test_contains_example(self):
            assert Span(5947, 4219, 9408, True)
            with pytest.raises(IndexError):
                Span(9408, 4219, 9408, True, strict=True)

        def test_contains_example2(self):
            assert Span(59, 42, 94, True)
            with pytest.raises(IndexError):
                Span(94, 42, 94, True, strict=True)


class TestIntersection:
    """These test the intersection between two spans, as in below:

    ..code-block::

        |-------|
            |------|
            |---|    << intersection
    """

    @staticmethod
    def x(a1, b1, a2, b2):
        s1 = Span(a1, b1, 100, True)
        s2 = Span(a2, b2, 100, True)
        return s1, s2

    def test_intersection(self):
        s1 = Span(10, 50, 100, True)
        s2 = Span(40, 60, 100, True)

        sliced = s1.intersection(s2)
        assert sliced.a == 40
        assert sliced.b == 50

        sliced = s2.intersection(s1)
        assert sliced.a == 40
        assert sliced.b == 50

    @pytest.mark.parametrize("i", [-1, 0, 1, 2])
    def test_no_intersection(self, i):
        s1, s2 = self.x(10, 50, 50 + i, 80)
        if i >= 0:
            assert not s1.intersection(s2)
            assert not s2.intersection(s1)
        else:
            assert s1.intersection(s2)
            assert s2.intersection(s1)

    def test_encompassed_intersection(self):
        s1, s2 = self.x(10, 50, 20, 30)
        sliced = s1.intersection(s2)
        assert sliced.a == 20
        assert sliced.b == 30

        sliced = s2.intersection(s1)
        assert sliced.a == 20
        assert sliced.b == 30


class TestSlice(object):
    """These test slicing and indexing"""

    @pytest.mark.parametrize("i", list(range(-20, 20)))
    @pytest.mark.parametrize("cyclic", [True, False])
    @pytest.mark.parametrize("index", [0, 1])
    def test_indexing(self, i, cyclic, index):

        s = Span(5, 15, 20, cyclic, index=index)
        assert len(s) == 10
        if i < -len(s) or i >= len(s):
            with pytest.raises(IndexError):
                print(s[i])
        else:
            if i < 0:
                assert s[i] == (i + 15) % 20
            else:
                assert s[i] == (i + 5)

    def test_slice_indexing_basic(self):
        s = Span(5, 9, 10, index=1)
        assert s[0] == 5
        assert s[-1] == 8

        s = Span(5, 12, 10, cyclic=True, index=0)
        print(list(s))
        assert s[0] == 5
        assert s[6] == 1
        assert s[-1] == 1
        assert s[len(s) - 1] == 1

    def test_slice_indexing_basic_index0(self):
        s = Span(5, 11, 10, cyclic=True, index=0)
        assert s.a == 5
        assert s.b == 1
        assert s[0] == 5
        assert s[1] == 6
        assert s[2] == 7
        assert s[3] == 8
        assert s[4] == 9
        assert s[5] == 0
        assert s[-1] == 0

    def test_slice_indexing_basic_index1(self):
        s = Span(6, 12, 10, cyclic=True, index=1)
        a = list(s)
        print(a)
        assert s.a == 6
        assert s.b == 2
        assert s[0] == 6
        assert s[1] == 7
        assert s[2] == 8
        assert s[3] == 9
        assert s[4] == 10
        assert s[5] == 1
        assert s[-1] == 1

    @pytest.mark.parametrize("i", list(range(10, 50, 5)))
    @pytest.mark.parametrize("j", list(range(60, 80, 5)))
    @pytest.mark.parametrize("cyclic", [True, False])
    def test_valid_slices(self, i, j, cyclic):
        s = Span(10, 100, 200, cyclic)
        sliced = s[i:j]
        assert len(sliced) == j - i

    @pytest.mark.parametrize("i", range(10))
    def test_slice_and_invert_slice_should_total_length(self, i):
        s = Span(10, 100, 200, True)
        if i != 4:
            x1 = s[4:i]
            x2 = s[i:4]
            assert len(x1) + len(x2) == 200

    @pytest.mark.parametrize("i", list(range(-20, 20)))
    def test_open_ended_slice_left(self, i):
        s = Span(10, 20, 200, False)
        if i >= 10 or i < -10:
            with pytest.raises(IndexError):
                assert not s[:i]
        else:
            sliced = s[:i]
            assert sliced.a == 10
            if i < 0:
                assert sliced.b == 20 + i
            else:
                assert sliced.b == 10 + i

    @pytest.mark.parametrize("i", list(range(-20, 20)))
    def test_open_ended_slice_right(self, i):
        s = Span(10, 20, 200, False)
        if i >= 10 or i < -10:
            with pytest.raises(IndexError):
                assert not s[i:]
        else:
            sliced = s[i:]
            assert sliced.b == 20
            if i < 0:
                assert sliced.a == 20 + i
            else:
                assert sliced.a == 10 + i

    def test_invalid_slice(self):
        s = Span(90, 10, 100, True)
        sliced = s[:15]
        assert sliced.a == 90
        assert sliced.b == 5

        s = Span(90, 10, 100, True)
        sliced = s[-15:]
        assert sliced.a == 95
        assert sliced.b == 10

    def test_copying(self):
        s = Span(5, 15, 20, True)
        sliced = s[:]
        assert sliced.a == s.a
        assert sliced.b == s.b

    def test_copy_new_ab(self):
        s = Span(5, 15, 20, True)
        copied = s[2, 16]
        assert copied.a == 2
        assert copied.b == 16

    def test_invert_cyclic(self):
        s = Span(5, 15, 20, True)
        i1 = s.invert()[0]

        assert i1.a == 15
        assert i1.b == 5

    def test_invert_linear(self):
        s = Span(5, 15, 20, False)
        s1, s2 = s.invert()
        assert s1.a == 0
        assert s1.b == 5
        assert s2.a == 15
        assert s2.b == 20

    def test_invert_wrapped(self):
        s = Span(5, 15, 10, True)
        s1, s2 = s.invert()
        assert len(s1) == 0
        assert s2 is None

    def test_slice_wrapped_basic(self):

        s = Span(5, 16, 10, True)
        assert len(s) == 11
        assert s.a == 5
        assert s.b == 6
        assert s.c == 16

        s1 = s[:-1]
        assert s1.a == 5
        assert s1.b == 5
        assert s1.c == 15
        assert len(s1) == 10

    def test_slice_wrapped_adv1(self):

        s = Span(5, 36, 10, True)
        assert s.a == 5
        assert s.b == 6
        assert s.c == 36
        assert len(s) == 31

        s1 = s[:-10]
        assert s1.a == 5
        assert s1.b == 6
        assert s1.c == 26
        assert len(s1) == 21

        s1 = s[:-20]
        assert s1.a == 5
        assert s1.b == 6
        assert s1.c == 16
        assert len(s1) == 11

        s1 = s[:-30]
        assert s1.a == 5
        assert s1.b == 6
        assert s1.c == 6
        assert len(s1) == 1

        with pytest.raises(IndexError):
            s1 = s[:-40]

    def test_slice_wrapped_adv2(self):
        s = Span(5, 36, 10, True)
        s1 = s[:4]
        assert s1.a == 5
        assert s1.b == 9
        assert s1.c == 9

        s1 = s[:6]
        assert s1.a == 5
        assert s1.b == 1
        assert s1.c == 11


class TestDifference:
    """These test the 'difference' between spans. The difference of two spans
    is outlined here:

    ..code-block::

        |----------|
           |----|
        |--|    |--|  << difference
    """

    def test_linear_diff(self):

        s1 = Span(20, 80, 100, True)
        s2 = Span(30, 50, 100, True)

        diff = s1.differences(s2)
        assert len(diff) == 2
        assert diff[0].a == 20
        assert diff[0].b == 30
        assert diff[1].a == 50
        assert diff[1].b == 80

    def test_overhang_left_diff(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(10, 30, 100, True)

        diff = s1.differences(s2)
        assert len(diff) == 1
        assert diff[0].a == 30
        assert diff[0].b == 80

        diff = s2.differences(s1)
        print(diff)
        assert len(diff) == 1
        assert diff[0].a == 10
        assert diff[0].b == 20

    def test_overhang_right_diff(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(70, 90, 100, True)
        diff = s1.differences(s2)
        assert len(diff) == 1
        assert diff[0].a == 20
        assert diff[0].b == 70

        diff = s2.differences(s1)
        print(diff)
        assert len(diff) == 1
        assert diff[0].a == 80
        assert diff[0].b == 90

    def test_no_overlap_diff(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(90, 95, 100, True)
        diff = s1.differences(s2)
        assert len(diff) == 1
        assert s1.a == 20
        assert s1.b == 80


class TestConsecutive:
    """These test whether spans are consecutive"""

    def test_consecutive(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(81, 95, 100, True)
        assert not s1.consecutive(s2)
        assert not s2.consecutive(s1)

    def test_not_consecutive(self):
        s1 = Span(20, 80, 100, True)
        s2 = Span(80, 95, 100, True)
        assert s1.consecutive(s2)
        assert not s2.consecutive(s1)

    # TODO: how to represent wrap-around regions?
    # def test_self_consecutive_cyclic(self):
    #     s1 = Span(2, 1, 10, True)
    #     assert s1.consecutive(s1)

    def test_consecutive_over_origin1(self):
        s1 = Span(7, 0, 10, True)
        s2 = Span(0, 2, 10, True)
        assert s1.consecutive(s2)
        assert not s2.consecutive(s1)


class TestConnectingSpans:
    """These tests connecting spans."""

    def test_connecting_span(self):
        s1 = Span(20, 50, 100, False)
        s2 = Span(60, 75, 100, False)
        s3 = s1.connecting_span(s2)
        assert s3.a == 50
        assert s3.b == 60

    def test_connecting_span_cyclic(self):
        s1 = Span(10, 20, 100, True)
        s2 = Span(80, 90, 100, True)
        s3 = s1.connecting_span(s2)
        assert s3.a == 20
        assert s3.b == 80

        s4 = s2.connecting_span(s1)
        assert s4.a == 90
        assert s4.b == 10

    def test_connecting_span_at_origin(self):
        s1 = Span(50, 60, 100, True)
        s2 = Span(0, 30, 100, True)
        s3 = s1.connecting_span(s2)
        assert s3.a == 60
        assert s3.b == 0

    def test_connecting_span_over_origin(self):
        """The connecting span should span the origin"""
        s1 = Span(50, 60, 100, True)
        s2 = Span(5, 30, 100, True)
        s3 = s1.connecting_span(s2)
        assert s3.a == 60
        assert s3.b == 5

    def test_self_connecting_span(self):
        """The connecting span with a cyclic span should be equivalent to its inverse span."""
        s1 = Span(50, 60, 100, True)
        s2 = s1.connecting_span(s1)
        assert s2.a == 60
        assert s2.b == 50

    def test_connecting_span_with_overlap(self):
        """There is no connecting span with overlapping spans."""
        s1 = Span(10, 30, 100, True)
        s2 = Span(20, 50, 100, True)
        assert not s1.connecting_span(s2)

    def test_connecting_span_consecutive_is_empty(self):
        """The connecting span between two consecutive spans is empty."""
        s1 = Span(10, 30, 100, True)
        s2 = Span(30, 50, 100, True)
        assert not s1.overlaps_with(s2)
        s3 = s1.connecting_span(s2)
        assert len(s3) == 0

    def test_connecting_span_linear_no_span(self):
        s1 = Span(10, 99, 100, False)
        s2 = Span(0, 10, 100, False)
        s3 = s1.connecting_span(s2)
        assert s3 is None


class TestEmptySpan:
    @pytest.mark.parametrize("cyclic", [True, False])
    @pytest.mark.parametrize("index", [0, 1])
    @pytest.mark.parametrize("x", [2000, 2999])
    def test_empty_span(self, cyclic, index, x):
        s = Span(x, x, 3000, cyclic=cyclic, index=0)
        assert len(s) == 0
        assert s.a == x
        assert s.b == x

    @pytest.mark.parametrize("cyclic", [True, False])
    @pytest.mark.parametrize("index", [0, 1])
    @pytest.mark.parametrize("x", [2000, 2999])
    def test_empty_span_new(self, cyclic, index, x):
        s = Span(1, 2, 3000, cyclic=cyclic, index=0)
        s = s.new(x, x)
        assert len(s) == 0
        assert s.a == x
        assert s.b == x

    @pytest.mark.parametrize("x", [2000, 2999, 3000, 3001, 0, 1])
    @pytest.mark.parametrize("index", [0, 1])
    def test_empty_span_new_cyclic(self, index, x):
        s = Span(1, 10, 3000, cyclic=True)
        s2 = s.new(x, x, ignore_wrap=False)
        if x >= 3000:
            r = x - 3000
        else:
            r = x
        assert len(s2) == 0
        assert s2.a == r
        assert s2.b == r


class TestAllowWrap:
    @pytest.mark.parametrize("delta", range(0, 21, 5))
    @pytest.mark.parametrize("index", [0, 1, 2])
    def test_allow_wrap(self, delta, index):
        s = Span(
            90 + delta, 98 + delta, 100, cyclic=True, index=index, ignore_wrap=True
        )
        print(list(s.ranges()))
        print(list(s))
        assert len(s) == 8

    @pytest.mark.parametrize("index", [0, 1, 2])
    def test_simple_len(self, index):
        s = Span(90, 10, 100, cyclic=True, index=index)

        assert len(s) == 20
        print(s.t(100))

    @pytest.mark.parametrize("index", [0, 1, 2])
    def test_simple_len_allow_wrap(self, index):
        s = Span(90, 110, 100, cyclic=True, index=index, ignore_wrap=False)
        assert len(s) == 20


class TestSub(object):
    @pytest.mark.parametrize("x", [(100, 1000), (101, 1000), (100, 999)])
    def test_valid_sub1(self, x):
        s = Span(100, 1000, 10000, cyclic=True)
        s2 = s.sub(x[0], x[1])
        assert s2.a == x[0]
        assert s2.b == x[1]

    @pytest.mark.parametrize("x", [(100, 10000), (101, 10000), (100, 9999)])
    def test_valid_sub2(self, x):
        s = Span(100, 10000, 10000, cyclic=True)
        s2 = s.sub(x[0], x[1])
        assert s2.a == x[0]
        assert s2.b == x[1]

    @pytest.mark.parametrize("x", [(0, 9000), (1, 9000), (1, 8999)])
    def test_valid_sub3(self, x):
        s = Span(0, 9000, 10000, cyclic=True)
        s2 = s.sub(x[0], x[1])
        assert s2.a == x[0]
        assert s2.b == x[1]

    @pytest.mark.parametrize("x", [(0, 900, 1000), (900, 100, 1000)])
    @pytest.mark.parametrize("delta", [(0, 0), (1, 0), (0, 1), (1, 1)])
    def test_valid_sub4(self, x, delta):
        s = Span(x[0], x[1], x[2], cyclic=True)
        start = s.a + delta[0]
        end = s.b - delta[1]
        s2 = s.sub(start, end)
        assert s2.a == start
        assert s2.b == end

    @pytest.mark.parametrize("cyclic", [True, False])
    @pytest.mark.parametrize("x", [1432, 1433, 4778, 4779])
    def test_valid_sub_same_indices(self, x, cyclic):
        s = Span(1432, 4779, 4799, cyclic=cyclic)
        s2 = s.sub(x, x)
        assert len(s2) == 0
        assert s2.a == x
        assert s2.b == x

    @pytest.mark.parametrize("cyclic", [True, False])
    @pytest.mark.parametrize("x", [1431, 4780])
    def test_invalid_sub_same_indices(self, x, cyclic):
        s = Span(1432, 4779, 4799, cyclic=cyclic)
        with pytest.raises(IndexError):
            s.sub(x, x)

    @pytest.mark.parametrize("x", [(99, 1000), (100, 1001), (99, 1001)])
    def test_invalid_ranges(self, x):
        s = Span(100, 1000, 10000, cyclic=True)
        with pytest.raises(IndexError):
            s.sub(x[0], x[1])

    @pytest.mark.parametrize(
        "x", [(900, 100), (901, 100), (900, 99), (901, 920), (50, 80)]
    )
    def test_valid_subs_over_origin(self, x):
        s = Span(900, 100, 1000, True)
        s2 = s.sub(x[0], x[1])
        assert s2.a == x[0]
        assert s2.b == x[1]

    @pytest.mark.parametrize("x", [(80, 50), (899, 100), (900, 101), (20, 920)])
    def test_invalid_subs_over_origin(self, x):
        s = Span(900, 100, 1000, True)
        with pytest.raises(IndexError):
            s.sub(x[0], x[1])

    def test_invalid_span(self):
        s = Span(5947, 4219, 10000, cyclic=True, ignore_wrap=False)
        with pytest.raises(IndexError):
            s.sub(28, 5980)

    def test_invalid_span2(self):
        s = Span(5000, 4999, 10000, cyclic=True, ignore_wrap=False)
        with pytest.raises(IndexError):
            s.sub(4998, 5001)

    @pytest.mark.parametrize("x", [(0, 900, 1000), (900, 100, 1000)])
    @pytest.mark.parametrize("delta", [(-11, 0), (0, -1), (-1, -1)])
    def test_invalid_span_over_regions(self, x, delta):
        s = Span(x[0], x[1], x[2], cyclic=True)
        start = s.a + delta[0]
        end = s.b - delta[1]
        with pytest.raises(IndexError):
            s.sub(start, end)

    def test_span_over_region_twice(self):

        s1 = Span(2000, 3000, 3000, cyclic=True, ignore_wrap=False)
        s2 = Span(5000, 6000, 3000, cyclic=True, ignore_wrap=False)
        assert s1.a == 2000
        assert s1.b == 3000
        assert s2.a == 2000
        assert s2.b == 3000

    def test_sub_wrapped(self):

        s = Span(5, 16, 10, cyclic=True)
        s1 = s.sub(5, 6)
        assert len(s1) == 1

        s2 = s.sub(5, 12)
        assert len(s2) == 7


class TestFullWrap(object):

    # def test_full_wrap_same_index(self):
    #     s = Span(0, 0, 1000, cyclic=True, allow_wrap=True)
    #     s.a == 0
    #     s.b == 0
    #     s.c == 0
    #     assert len(s) == 0
    #     # assert not s.spans_origin()
    #
    #     s = Span(0, 1000, 1000, cyclic=True, allow_wrap=True)
    #     s.a == 0
    #     s.b == 1000
    #     s.c == 1000
    #     assert len(s) == 1000
    #     # assert s.spans_origin()
    #
    #     s = Span(0, 2000, 1000, cyclic=True, allow_wrap=True)
    #     s.a == 0
    #     s.b == 1000
    #     s.c == 2000
    #     assert s._nwraps == 1
    #     print(s.ranges())
    #     assert len(s) == 2000

    @pytest.mark.parametrize("start", [0])
    @pytest.mark.parametrize("end", [0, 1])
    @pytest.mark.parametrize("nwraps", [0, 1, 2, 3])
    def test_full_wrap_same_index(self, nwraps, start, end):
        index = 0
        length = 1000
        s = Span(
            start,
            end + nwraps * length,
            length,
            cyclic=True,
            ignore_wrap=False,
            index=index,
        )
        s.a == start
        s.b == end
        s.c == end + nwraps * length
        print(list(s.ranges()))

        assert len(s) == nwraps * length + end - start

    def test_full_wrap_plus_one(self):
        s = Span(0, 1, 1000, cyclic=True, ignore_wrap=False)
        assert s.a == 0
        assert s.b == 1
        assert s.c == 1
        assert len(s) == 1
        assert not s.spans_origin()

        s = Span(0, 1001, 1000, cyclic=True, ignore_wrap=False)
        assert s.a == 0
        assert s.b == 1
        assert s.c == 1001
        assert len(s) == 1001
        assert s.spans_origin()

        s = Span(0, 2001, 1000, cyclic=True, ignore_wrap=False)
        assert s.a == 0
        assert s.b == 1
        assert s.c == 2001
        print(list(s.ranges()))
        assert len(s) == 2001
        assert s.spans_origin()

    def test_full_wrap_sub_region(self):
        s = Span(0, 1000, 1000, cyclic=True)
        assert s.contains_pos(1)
        s2 = s.sub(500, 600)
        assert s2.a == 500
        assert s2.b == 600

    def test_invert(self):
        s = Span(0, 1000, 1000, cyclic=True)
        assert s.contains_pos(1)
        s2 = s.invert()
        assert len(s2[0]) == 0

    def test_copy(self):
        s = Span(0, 1000, 1000, cyclic=True)
        s2 = s[:]
        assert len(s) == 1000
        assert len(s2) == 1000

    def test_copy_wrapped(self):
        s = Span(0, 1002, 1000, cyclic=True)
        s2 = s[:]
        assert len(s) == 1002
        assert len(s2) == 1002
        assert s.a == s2.a
        assert s.b == s2.b
        assert s.c == s2.c


class TestNWraps:
    @pytest.mark.parametrize("i", [0, 1, 2])
    def test_nwraps(self, i):
        l = 1000
        s = Span(100, 100 + i * l, l, cyclic=True, ignore_wrap=False)
        assert len(s) == i * l


class TestSlicingNumpy:
    def test_numpy_slice(self):
        s = Span(50, 120, 100, cyclic=True)
        a = np.arange(100, 200)
        b = a[s]
        print(b)


class TestChangingStartingIndex:
    def test_change_starting_index(self):

        s = Span(1, 10, 10, index=1)
        assert s.a == 1
        assert s.b == 10
        assert s.c == 10
        assert s.index == 1
        assert len(s) == 9

        s2 = s.reindex(0)

        assert s2.a == 0
        assert s2.b == 9
        assert s2.c == 9
        assert s2.index == 0
        assert len(s) == 9

        s2 = s.reindex(-5)
        assert s2.a == 0 - 5
        assert s2.b == 9 - 5
        assert s2.c == 9 - 5
        assert s2.index == -5
        assert len(s) == 9

    def test_change_starting_index_cyclic(self):

        s = Span(90, 110, 100, cyclic=True, index=1)
        s2 = s.reindex(0)
        s2.a == 89
        s2.b == 109
        s2.c == 109
        assert len(s2) == len(s)
        assert s2.index == 0

    def test_change_starting_index_many_wraps(self):

        s = Span(90, 210, 100, cyclic=True, index=1)
        s2 = s.reindex(0)
        s2.a == 89
        s2.b == 109
        s2.c == 209
        assert len(s2) == len(s)
        assert s2.index == 0


def test_ignore_wrap():

    with pytest.raises(IndexError):
        Span(3000, 1, 3000, cyclic=True, ignore_wrap=False)

    s = Span(3000, 1, 3000, cyclic=True, ignore_wrap=True)
    assert s.a == 0
    assert s.b == 1
    assert s.c == 1

    s = Span(3000, 3000 * 10 + 1, 3000, cyclic=True, ignore_wrap=True)
    assert s.a == 0
    assert s.b == 1
    assert s.c == 1


def test_abs_wrap():

    s = Span(6010, 1, 3000, cyclic=True, abs_wrap=True)
    assert s.a == 10
    assert s.b == 1
    assert s.c == 6001

    with pytest.raises(IndexError):
        Span(6010, 1, 3000, cyclic=True, abs_wrap=False)


def test_abs_wrap_doctest():
    # all of the following initializations result in equivalent spans

    # starts at 1, wraps around one full time, ends at 5.
    # Length is 15 - 1
    s1 = Span(1, 15, 10, cyclic=True, abs_wrap=True)
    assert len(s1) == 14

    # starts at 11, wraps around one full time, ends at 15.
    # Then positions are mapped back to context at 1 and 5.
    # Length is still 25 - 11 == 14
    s2 = Span(11, 25, 10, cyclic=True, abs_wrap=True)
    assert len(s2) == 14

    # starts at 11, wraps around one full time, ends at 15.
    # Then positions are mapped back to context at 1 and 5.
    # Length is now 11 + 5 = 14
    s3 = Span(11, 5, 10, cyclic=True, abs_wrap=True)
    assert len(s3) == 14

    # the lengths change with the abs diff in number of times wrapped
    _s = Span(21, 5, 10, cyclic=True, abs_wrap=True)
    assert len(_s) == 24

    _s = Span(21, 15, 10, cyclic=True, abs_wrap=True)
    assert len(_s) == 14

    # all
    assert s1 == s2
    assert s2 == s3


def test_subregion():

    s = Span(5000, 13000, 9000, cyclic=True)

    s1 = s.sub(5000, 4000)
    s2 = s.sub(5000, 4000 + 9000)

    assert s1.a == 5000
    assert s1.b == 4000
    assert s1.c == 4000

    assert s2.a == 5000
    assert s2.b == 4000
    assert s2.c == 13000
