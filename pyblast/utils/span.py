from collections.abc import Container
from collections.abc import Iterable
from collections.abc import Sized
from itertools import chain
from typing import List
from typing import Tuple
from typing import Union


class SpanError(Exception):
    pass


class Span(Container, Iterable, Sized):
    """`Span` maps the provided positions onto a context.

    Spans have no direction and have an underlying context
    that has a certain length, can be linear or cyclic, and has a starting index.

    The following attributes are accessible:

    - `a` - the starting (inclusive) position of the span
    - `b` - the mapped exclusive position of the span
    - `c` - the unmapped exclusive position of the span (used for cyclic spans that
    wrap context multiple times)
    - `cyclic` - whether the underlying context is circular or linear
    - `index` - the first index of the context
    - `context_length` - the length of the underlying context.

    **Linear spans**
    For example, a basic linear span can be represented using the following:

    .. code-block::

        s = Span(0, 10, 20)
        assert s.a == 0
        assert s.b == 10
        # spans are always exclusive at endpoint.
        assert list(s) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    **Indexing and Slicing**

    New indexes can be provided to spans.

    .. code-block::

        s = Span(1, 10, 10, index=1)
        assert s.a == 1
        assert s.b == 10

    Positions can be mapped to the span:

    .. code-block::

        s = Span(5, 8, 10, index=5)
        assert s.a == 5
        assert s.b == 8

        # the 'first' position in the span equivalent to the starting index
        assert s[0] == 5

        # the 'last' position in the span is the last inclusive index
        assert s[-1] == 7

    Positions can be mapped automatically during initialization:

    .. code-block::

        s = Span(-1, 5, 10, cyclic=True, index=0)  # index '-1' is mapped onto last
                                                   # available index on the context,
        '9'
        assert s.a == 9
        assert s.b == 5

        s = Span(6, -1, 10, cyclic=True, index=0)  # index '-1' is mapped onto last
                                                   # available exclusive index on the
         context, '10'
        assert s.a == 6
        assert s.b == 10

    Inclusive positions can be mapped onto the context using `t`:

    .. code-block::

        s = Span(-1, 5, 10, cyclic=True)
        assert s.t(11) == 1
        assert s.t(10) == 0

    **Reindexing**

    Context starting index can be remapped:

    .. code-block::

        s = Span(0, 5, 10)
        s1 = s.reindex(1)
        asset s1.a == 1
        assert s2.b == 6

    **Slicing**

    Spans can be sliced similarly to lists:

    .. code-block::

        s = Span(4, 8, 10)
        s1 = s[1:]
        assert s1.a == 5
        assert s2.b == 8

        s2 = s[:-1]
        assert s2.a == 4
        assert s2.b == 7

        s3 = s[2:-1]
        assert s3.a == 5
        assert s3.b == 7

    **Cyclic spans**
    A cyclic span can be represented in the following way:

    .. code-block::

        s = Span(18, 2, 20, cyclic=True)
        assert s.a == 18
        assert s.b == 2
        assert list(s) == [18, 19, 0, 1]
        assert len(s) == 4

    **Wrapping cyclic spans**
    Spans the wrap around the context multiple times can be represented as well.
    The mapped endpoint
    is found with `span.b` and the non-mapped endpoint is found using `span.c`.

    .. code-block::

        s = Span(9, 21, 10, cyclic=True)
        assert list(s) == [9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        assert len(s) == 12
        assert s.a == 9  # start point
        assert s.b == 1  # mapped endpoint
        assert s.c == 21 # the unmapped endpoint

    The indices are reduced to their lowest 'wrapping' whenever possible. For example,
    the following
    initializations are equivalent:

    .. code-block::

        s1 = Span(1, 5, 10, cyclic=True)
        s2 = Span(11, 15, 10, cyclic=True)
        s3 = Span(21, 25, 10, cyclic=True)
        assert s1 == s2
        assert s2 == s3

    **Iterating through spans**

    Spans can be treated like iterators:

    .. code-block::

        for i in Span(5, 3, 10, cyclic=True):
            print(i)

    **Span inversion**

    Inversion returns two spans that represent everything *except* the span:

    .. code-block::

        s = Span(4, 8, 10)
        s1, s2 = s.invert()  # also s[::-1]
        assert s1.a == 0
        assert s1.b == 4
        assert s2.a == 8
        assert s2.b == 10

        # inverted cyclic span returns one valid span and None
        s = Span(4, 8, 10, cyclic=True)
        s1, s2 = s.invert()
        assert s1.a == 8
        assert s1.b == 4
        assert s2 is None

    **Intersections, differences**

    .. code-block::

        s1 = Span(4, 10, 20)
        s2 = Span(8, 12, 20)
        s3 = s1.intersection(s2)
        assert s3.a == 8
        assert s3.b == 10

    .. code-block::

        s1 = Span(4, 10, 20)
        s2 = Span(8, 12, 20)
        s3, s4 = s1.differences(s2)

        assert s3.a, s3.b == 4, 8
        assert s4,a, s4.b == 10, 12
    """

    __slots__ = [
        "_a",
        "_b",
        "_c",
        "_context_length",
        "_cyclic",
        "_index",
        "_strict",
        "_abs_wrap",
        "_ignore_wrap",
    ]

    def __init__(
        self,
        a: int,
        b: int,
        l: int,
        cyclic=False,
        index=0,
        ignore_wrap=False,
        strict=False,
        abs_wrap=False,
    ):
        """Constructs a new Span. There are several options to customize the
        initialization procedure.

        **strict=True**

        When strict, any index outside the valid bounds of the context raises an IndexError.

        .. code-block::

            Span(1, 10, 10, cyclic=True, strict=True) # no raise
            Span(1, 11, 10, cyclic=True, strict=True) # raises IndexError
            Span(0, 11, 10, index=1, cyclic=True, strict=True) # raises IndexError

        **ignore_wrap=True**

        When wrapping is ignored, indices are simply mapped to the context with no
        consideration
        of the number of times the absolute position would wrap around the context.

        .. code-block::

            # all of the following are equivalent with ignore_wrap == True
            Span(1, 10, 10, cyclic=True, ignore_wrap=True)
            Span(1+10, 10, 10, cyclic=True, ignore_wrap=True)
            Span(1, 10-100, 10, cyclic=True, ignore_wrap=True)
            Span(1+20, 10-100, 10, cyclic=True, ignore_wrap=True)

        **abs_wrap=True**

        When absolute wrapping is used, the absolute difference between starting and
        ending index wrappings
        is calculated, the starting index is to the context while the ending index is
        adjusted such
        that the length will reflect the abs difference between starting and ending
        index wrappings.
        This can be unintuitive
        is best shown with the following example:

        .. code-block::

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

        :param a: start of the span (inclusive)
        :type a: int
        :param b: end of the span (exclusive)
        :type b: int
        :param l: context length of the region
        :type l: int
        :param cyclic: whether the underlying context is cyclic
        :type cyclic: bool
        :param index: the starting index of the region
        :type index: int
        :param strict: if True, positions outside of context bounds are disallowed.
        :type strict: bool
        :param ignore_wrap: if True (default False), initialization indicies that wrap
                            around multiple times will
                            simply be mapped directly to the context (no wrapping used).
        :type ignore_wrap: bool
        :param abs_wrap: if True, the abs difference between start and end wrappings are
                            used. Starting wraps that
                            are greater than ending wraps are valid. If False and the
                            starting wrap is greater than
                            the ending wrap, an IndexError is thrown.
        """

        self._context_length = l
        self._index = index
        self._cyclic = cyclic
        self._strict = strict
        self._abs_wrap = abs_wrap
        self._ignore_wrap = ignore_wrap

        # special empty edge case
        if cyclic and a == b:
            self._set_as_empty(a)
            return

        # check bounds
        if self._strict or not cyclic:
            bounds = self.bounds()
            if not bounds[0] <= a < bounds[1]:
                raise IndexError(
                    "Start {} must be in [{}, {})".format(a, index, index + l)
                )
            if not bounds[0] <= b <= bounds[1]:
                raise IndexError(
                    "End {} must be in [{}, {}]".format(b, index, index + l)
                )

        if self._ignore_wrap:
            start_wrap = 0
            end_wrap = 0
        else:
            start_wrap = int((a - index) / l)
            end_wrap = int((b - index - 1) / (l))

        if self._abs_wrap is False and start_wrap > end_wrap:
            self._a = a
            self._b = b
            self._c = b
            diff = end_wrap - start_wrap
            # return self._set_as_empty(a)
            raise IndexError(
                "Could not interpret span {span}. Starting position wraps around "
                "context {i} times and end position wraps around {j} times."
                " A valid initialization would be Span({a}, {b}, ...)".format(
                    span=self,
                    i=start_wrap,
                    j=end_wrap,
                    a=self._a,
                    b=self._b - diff * self._context_length,
                )
            )

        # set indices
        _a = a - index
        _b = b - index

        if _a >= l or _a < 0:
            self._a = self.t(_a, False)
        else:
            self._a = a
        if _b > l:
            self._b = self.t(_b - 1, False) + 1
        elif _b < 0:
            self._b = self.t(_b, False) + 1
        else:
            self._b = b

        if self._a > self._b and not cyclic:
            raise IndexError(
                "Start {} cannot be greater than end {} for linear spans.".format(
                    self._a, self._b
                )
            )

        # allow wrap mean this will keep track of how many time the span wraps around the context
        if not ignore_wrap and end_wrap - start_wrap:
            if self._abs_wrap:
                _c = self._b + abs(end_wrap - start_wrap) * l
            else:
                _c = self._b + (end_wrap - start_wrap) * l
            self._c = _c
        else:
            self._c = self._b

    @property
    def index(self):
        """Return the starting index of the context."""
        return self._index

    @property
    def a(self):
        """Return the inclusive mapped startpoint."""
        return self._a

    @property
    def b(self):
        """Return the exclusive mapped endpoint."""
        return self._b

    @property
    def c(self):
        """Return the exclusive un-mapped endpoint."""
        return self._c

    @property
    def cyclic(self):
        """Return whether the context is cyclic/circular."""
        return self._cyclic

    @property
    def context_length(self):
        """Return the length of the context."""
        return self._context_length

    def _set_as_empty(self, a: int) -> None:
        _a = self.t(a - self._index, False)
        self._a = self._b = self._c = _a
        return

    @property
    def _nwraps(self):
        return int((self._c - self._index - 1) / (self._context_length))

    def bounds(self) -> tuple:
        """Return the context bounds (end exclusive)"""
        return self._index, self._context_length + self._index

    def t(self, p: int, throw_error=True) -> int:
        """Translates a position 'p' to an index within the context bounds.

        :param p:
        :type p:
        :param throw_error:
        :return:
        :rtype:
        """
        if p >= self.bounds()[1] or p < self.bounds()[0]:
            if throw_error and not self._cyclic:
                raise IndexError(
                    "Position {} outside of linear bounds {}".format(p, self.bounds())
                )
        _x = p % self._context_length
        if _x < 0:
            return self.bounds()[1] + _x
        else:
            return self.bounds()[0] + _x

    @staticmethod
    def _ranges_str(ranges):
        s = ",".join("[{}, {})".format(*r) for r in ranges)
        return "[" + s + "]"

    def ranges(self) -> List[Tuple[int, int]]:
        """Return ranges of valid positions.

        :return:
        :rtype:
        """

        # TODO: replace this with _nwraps > 0
        if self._cyclic and (self._b < self._a or self._nwraps):
            ranges = [(self._a, self.bounds()[1])]
            for _ in range(self._nwraps - 1):
                ranges.append(self.bounds())
            ranges.append((self.bounds()[0], self._b))

            return ranges
        else:
            return [(self._a, self._b)]

    def slices(self):
        """Return list of slices.

        :return:
        :rtype:
        """
        return [slice(*r) for r in self.ranges()]

    def reindex(self, i, strict=None, ignore_wrap=None):
        """Return a new span with positions reindexed.

        :param i: new index
        :param strict: initialize with 'strict'
        :param ignore_wrap: whether to ignore wrapping indices
        :return:
        """
        if ignore_wrap is None:
            ignore_wrap = self._ignore_wrap
        return self.new(None, None, ignore_wrap=ignore_wrap, index=i, strict=strict)

    def new(
        self,
        a: Union[int, None],
        b: Union[int, None],
        ignore_wrap=None,
        index=None,
        strict=None,
        abs_wrap=None,
    ) -> "Span":
        """Create a new span using the same context."""
        if a is None:
            a = self._a
        if b is None:
            b = self._c

        if strict is None:
            strict = self._strict

        if abs_wrap is None:
            abs_wrap = self._abs_wrap

        if ignore_wrap is None:
            ignore_wrap = self._ignore_wrap

        if index is not None:
            d = index - self._index
            a += d
            b += d
        else:
            index = self._index
        return self.__class__(
            a,
            b,
            self._context_length,
            self._cyclic,
            index=index,
            ignore_wrap=ignore_wrap,
            strict=strict,
            abs_wrap=abs_wrap,
        )

    def sub(self, a: int, b: int) -> "Span":
        """Create a sub region starting from a to b.

        :param a: starting pos
        :param b: ending pos
        :return: a new Span.
        """
        # if a == b:
        #     return self.new(a, b)
        if b is not None and a > b and not self._cyclic:
            raise ValueError(
                "Start {} cannot exceed end {} for linear spans".format(a, b)
            )

        valid_ranges = [list(x) for x in self.ranges()]
        valid_ranges[0][0] = max(a, self._a)
        valid_ranges[-1][1] = min(b + 1, self._b + 1)

        if self._nwraps > 0:
            valid_ranges += [(self._a, self._c + 1)]
        # assert len(valid_ranges) <= 2

        def in_range(pos, ranges):
            for i, r in enumerate(ranges):
                if r[0] <= pos < r[1]:
                    return (True, i)
            return (False, None)

        start_in_range, start_range = in_range(a, valid_ranges)

        if not start_in_range:
            raise IndexError(
                "Start {} must be in {}".format(a, self._ranges_str(valid_ranges))
            )

        valid_end_range = valid_ranges[start_range:]
        end_in_range, end_range = in_range(b, valid_end_range)

        if not end_in_range:
            raise IndexError(
                "End {} must be in {}".format(b, self._ranges_str(valid_end_range))
            )
        subregion = self.new(a, b)
        if len(subregion) > len(self):
            raise IndexError("Cannot make subspan.")
        return subregion

    def same_context(self, other: "Span") -> bool:
        """Return if another Span as an equivalent context."""
        return (
            other.context_length == self._context_length
            and self._cyclic == other.cyclic
        )

    def force_context(self, other: "Span") -> None:
        """Raise error if another Span has different context.

        :param other: The other span
        :return: None
        :raises: "Span"Error
        """
        if not self.same_context(other):
            raise SpanError("Cannot compare with different contexts")

    def overlaps_with(self, other: "Span") -> bool:
        """Returns True if other span has an overlap with this span.

        :param other: The other span
        :return: True if the two spans overlap.
        """
        self.force_context(other)
        # if other in self:
        #     return True
        # elif self in other:
        #     return True
        if (
            other.a in self
            or other.b - 1 in self
            or self._a in other
            or self._b - 1 in other
        ):
            return True
        return False

    def differences(self, other: "Span") -> Union[Tuple["Span"], Tuple["Span", "Span"]]:
        """Return a tuple of differences between this span and the other
        span."""
        self.force_context(other)
        if other.a in self and other.b not in self:
            return (self.new(self._a, other.a),)
        elif other.b in self and other.a not in self:
            return (self.new(other.b, self._b),)
        if other in self:
            return self.new(self._a, other.a), self.new(other.b, self._b)
        elif self in other:
            return (self.new(self._a, self._a),)
        else:
            return ((self[:]),)

    def intersection(self, other: "Span") -> "Span":
        """Return the span inersection between this span and the other span."""
        self.force_context(other)
        if other.a in self and other.b not in self:
            return self.sub(other.a, self._b)
        elif other.b in self and other.a not in self:
            return self.sub(self._a, other.b)
        if other in self:
            return other[:]
        elif self in other:
            return self[:]

    def consecutive(self, other: "Span") -> bool:
        """Returns True if other span is immediately consecutive with this
        span."""
        self.force_context(other)
        try:
            return self._b == other.t(other.a)
        except IndexError:
            return False

    def connecting_span(self, other: "Span") -> Union["Span", None]:
        """Return the span that connects the two spans. Returns None.

        :param other:
        :return:
        """
        self.force_context(other)
        if self._cyclic and self._a == other.a and self._b == other.b:
            return self.invert()[0]
        if self.consecutive(other):
            return self[self._b, self._b]
        elif self.overlaps_with(other):
            return None
        else:
            if self._b > other.a and not self._cyclic:
                return None
            return self[self._b, other.a]

    # def union(self, other):
    #     self.force_context(other)
    #     if other in self:
    #         return self[:]
    #     elif self in other:
    #         return other[:]
    #     elif other.a in self and not other.t(other.b - 1) in self:
    #         if self._a == other.b:
    #             return self.new(self._a, None)
    #         return self.new(self._a, other.b)
    #     elif other.t(other.b - 1) in self and other.a not in self:
    #         if other.a == self._b:
    #             return self.new(self._a, None)
    #         return self.new(other.a, self._b)

    # def __ge__(self, other):
    #     self.force_context(other)
    #     return self._a >= other.a

    #     def __lt__(self, other):
    #         return self._a < other.a

    #     def __gt__(self, other):
    #         return self._a > other.a

    #     def __ge__(self, other):
    #         return self._a >= other.a

    # def __invert__(self):
    #     if self._a > self._b:
    #         if self._cyclic:
    #             return self[self._b+1, self._a-1],
    #         else:
    # return

    def invert(self) -> Union[Tuple["Span", "Span"], Tuple["Span", None]]:
        """Invert the region, returning a tuple of the remaining spans from the
        context. If cyclic, a tuple (span, None) tuple is returned. If linear,
        a (span, span) is returned.

        :return: inverted regions
        :rtype: tuple
        """
        if len(self) >= self._context_length:
            return (self[self._a, self._a], None)
        if self._cyclic:
            return (self[self._b, self._a], None)
        else:
            return self[:, self._a], self[self._b, :]

    def __eq__(self, other: "Span") -> bool:
        return (
            self.same_context(other)
            and self._a == other._a
            and self._b == other._b
            and self._c == other._c
        )

    def __ne__(self, other: "Span") -> bool:
        return not (self.__eq__(other))

    @classmethod
    def _pos_in_ranges(cls, pos: int, ranges: List[Tuple[int, int]]):
        for r in ranges:
            if r[0] <= pos < r[1]:
                return True
        return False

    def contains_pos(self, pos: int) -> bool:
        return self._pos_in_ranges(pos, self.ranges())

    def contains_span(self, other: "Span") -> bool:
        if not self.same_context(other):
            return False
        if other.a == other.b == self._a:
            # special case where starting indices and ending indices are the same
            return True
        else:
            if not self.contains_pos(other.a):
                return False
            elif not self.contains_pos(other.t(other.b - 1)):
                return False
            elif not len(other) <= len(self):
                return False
            return True

    def spans_origin(self) -> bool:
        if self._nwraps and self._cyclic:
            return True
        return self._b < self._a and self._cyclic

    def __contains__(self, other: Union["Span", int]):
        if isinstance(other, int):
            return self.contains_pos(other)
        elif issubclass(type(other), Span):
            return self.contains_span(other)

    def __len__(self):
        return sum([r[1] - r[0] for r in self.ranges()])

    def __iter__(self):
        for i in chain(*[range(*x) for x in self.ranges()]):
            yield i

    def __invert__(self):
        return self.invert()

    def _check_index_pos(self, val: int) -> None:
        if val is not None and (val >= len(self) or val < -len(self)):
            raise IndexError(
                "Index '{}' outside of linear span with length of {}".format(
                    val, len(self)
                )
            )

    def __getitem__(self, val):
        if isinstance(val, int):
            self._check_index_pos(val)
            if val < 0:
                return self.t(val + self._c - self._index)
            else:
                return self.t(val + self._a - self._index)
        elif issubclass(type(val), slice):

            self._check_index_pos(val.start)
            self._check_index_pos(val.stop)

            if val.step == -1:
                return self[val.start : val.stop].invert()
            elif val.step is not None and val.step != 1:
                raise ValueError(
                    "{} slicing does not support step {}.".format(
                        self.__class__.__name__, val.step
                    )
                )

            if val.start is None:
                i = self._a
            elif val.start < 0:
                i = self._c + val.start
            else:
                i = self._a + val.start

            if val.stop is None:
                j = self._c
            elif val.stop < 0:
                j = self._c + val.stop
            else:
                j = self._a + val.stop

            return self.new(i, j)
        elif isinstance(val, tuple):
            if len(val) > 2:
                raise ValueError(
                    "{} -- copying only supports (start, stop)".format(val)
                )
            val = list(val)
            for i in [0, 1]:
                if val[i] == slice(None, None, None):
                    val[i] = None
            if val == (None, None):
                val = self.bounds()
            elif val[0] is None:
                val = (self.bounds()[0], val[1])
            elif val[1] is None:
                val = (val[0], self.bounds()[1])
            return self.new(*val)
        else:
            raise ValueError("indexing does not support {}".format(type(val)))

    def __repr__(self):
        return (
            "<{cls} {a} {b} ({c}) len={length} cyclic={cyclic} index={index}, "
            "nwraps={n}>".format(
                cls=self.__class__.__name__,
                a=self._a,
                b=self._b,
                c=self._c,
                length=self._context_length,
                cyclic=self._cyclic,
                index=self._index,
                n=self._nwraps,
            )
        )

    def __str__(self):
        return self.__repr__()


class EmptySpan(Span):
    def ranges(self):
        return [(self.bounds()[0], self.bounds()[0])]
