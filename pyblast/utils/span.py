from collections.abc import Container, Iterable, Sized
from itertools import chain


class SpanError(Exception):
    pass


class Span(Container, Iterable, Sized):
    __slots__ = ["a", "b", "context_length", "cyclic", "index", "_does_wrap_origin"]

    def __init__(
        self, a, b, l, cyclic=False, index=0, allow_wrap=False, does_wrap_origin=False
    ):
        """
        Constructs a new Span.

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
        :param allow_wrap: if True (default False), the region can be initialized with a and b over the origin
        :type allow_wrap: bool
        :param does_wrap_origin: if True (default False), if the span.a == span.b, then this span is specified to
                                 wrap around the origin, fully encompassing the context.
        :type does_wrap_origin: bool
        """
        if a > b and not cyclic:
            raise IndexError(
                "Start {} cannot exceed end {} for linear spans".format(a, b)
            )
        self._does_wrap_origin = False

        self.index = index
        self.context_length = l
        self.cyclic = cyclic

        # check bounds
        if not cyclic or not allow_wrap:
            bounds = self.bounds()
            if not bounds[0] <= a < bounds[1]:
                raise IndexError(
                    "Start {} must be in [{}, {})".format(a, index, index + l)
                )
            if not bounds[0] <= b <= bounds[1]:
                raise IndexError(
                    "End {} must be in [{}, {}]".format(b, index, index + l)
                )

        # empty edge case
        if a == b and a == l + index:
            self.a = self.b = index
            return

        # set indices
        _a = a - index
        _b = b - index
        if _a >= l or _a < 0:
            self.a = self.t(_a)
        else:
            self.a = a
        if _b > l or _b < 0:
            self.b = self.t(_b - 1) + 1
        else:
            self.b = b

        if self.a == self.b:
            self._does_wrap_origin = does_wrap_origin

    def bounds(self):
        """Return the bounds (end exclusive)"""
        return self.index, self.context_length + self.index

    def t(self, p):
        """Translate an index to a valid index"""
        if p >= self.bounds()[1] or p < self.bounds()[0]:
            if not self.cyclic:
                raise IndexError(
                    "Position {} outside of linear bounds {}".format(p, self.bounds())
                )
        return (p % self.context_length) + self.bounds()[0]

    @staticmethod
    def _ranges_str(ranges):
        s = ",".join("[{}, {})".format(*r) for r in ranges)
        return "[" + s + "]"

    def ranges(self):
        """
        Return ranges of valid positions.

        :return:
        :rtype:
        """
        if self.cyclic and (
            (self.a > self.b) or (self.a == self.b and self._does_wrap_origin)
        ):
            return [(self.a, self.bounds()[1]), (self.bounds()[0], self.b)]
        else:
            return [(self.a, self.b)]

    def new(self, a, b, allow_wrap=True, does_wrap_origin=False):
        """Create a new span using the same context."""
        return self.__class__(
            a,
            b,
            self.context_length,
            self.cyclic,
            index=self.index,
            allow_wrap=allow_wrap,
            does_wrap_origin=does_wrap_origin,
        )

    def sub(self, a, b):
        """
        Create a sub region starting from a to b.
        :param a:
        :type a:
        :param b:
        :type b:
        :return:
        :rtype:
        """
        # if a == b:
        #     return self.new(a, b)
        if b is not None and a > b and not self.cyclic:
            raise ValueError(
                "Start {} cannot exceed end {} for linear spans".format(a, b)
            )

        valid_ranges = [list(x) for x in self.ranges()]

        valid_ranges[0][0] = max(a, self.a)
        valid_ranges[-1][1] = min(b + 1, self.b + 1)
        assert len(valid_ranges) <= 2

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
            raise IndexError("Cannot make subspan. S")
        return subregion

    def same_context(self, other):
        return (
            other.context_length == self.context_length and self.cyclic == other.cyclic
        )

    def force_context(self, other):
        if not self.same_context(other):
            raise SpanError("Cannot compare with different contexts")

    def overlaps_with(self, other):
        self.force_context(other)
        # if other in self:
        #     return True
        # elif self in other:
        #     return True
        if (
            other.a in self
            or other.b - 1 in self
            or self.a in other
            or self.b - 1 in other
        ):
            return True
        return False

    def differences(self, other):
        """Return a tuple of differences between this span and the other span."""
        self.force_context(other)
        if other.a in self and other.b not in self:
            return (self.new(self.a, other.a),)
        elif other.b in self and other.a not in self:
            return (self.new(other.b, self.b),)
        if other in self:
            return self.new(self.a, other.a), self.new(other.b, self.b)
        elif self in other:
            return (self.new(self.a, self.a),)
        else:
            return ((self[:]),)

    def intersection(self, other):
        """Return the span inersection between this span and the other span."""
        self.force_context(other)
        if other.a in self and other.b not in self:
            return self.sub(other.a, self.b)
        elif other.b in self and other.a not in self:
            return self.sub(self.a, other.b)
        if other in self:
            return other[:]
        elif self in other:
            return self[:]

    def consecutive(self, other):
        self.force_context(other)
        try:
            return self.b == other.t(other.a)
        except IndexError:
            return False

    def connecting_span(self, other):
        """
        Return the span that connects the two spans. Returns None

        :param other:
        :return:
        """
        self.force_context(other)
        if self.cyclic and self.a == other.a and self.b == other.b:
            return self.invert()[0]
        if self.consecutive(other):
            return self[self.b, self.b]
        elif self.overlaps_with(other):
            return None
        else:
            if self.b > other.a and not self.cyclic:
                return None
            return self[self.b, other.a]

    # def union(self, other):
    #     self.force_context(other)
    #     if other in self:
    #         return self[:]
    #     elif self in other:
    #         return other[:]
    #     elif other.a in self and not other.t(other.b - 1) in self:
    #         if self.a == other.b:
    #             return self.new(self.a, None)
    #         return self.new(self.a, other.b)
    #     elif other.t(other.b - 1) in self and other.a not in self:
    #         if other.a == self.b:
    #             return self.new(self.a, None)
    #         return self.new(other.a, self.b)

    # def __ge__(self, other):
    #     self.force_context(other)
    #     return self.a >= other.a

    #     def __lt__(self, other):
    #         return self.a < other.a

    #     def __gt__(self, other):
    #         return self.a > other.a

    #     def __ge__(self, other):
    #         return self.a >= other.a

    # def __invert__(self):
    #     if self.a > self.b:
    #         if self.cyclic:
    #             return self[self.b+1, self.a-1],
    #         else:
    # return

    def invert(self):
        """
        Invert the region, returning a tuple of the remaining spans from the context.
        If cyclic, a tuple (span, None) tuple is returned. If linear, a (span, span) is returned.

        :return: inverted regions
        :rtype: tuple
        """
        if self.cyclic:
            return (self[self.b, self.a], None)
        else:
            return self[:, self.a], self[self.b, :]

    def __eq__(self, other):
        return self.same_context(other) and self.a == other.a and self.b == other.b

    def __ne__(self, other):
        return not (self.__eq__(other))

    @classmethod
    def _pos_in_ranges(cls, pos, ranges):
        for r in ranges:
            if r[0] <= pos < r[1]:
                return True
        return False

    def contains_pos(self, pos):
        return self._pos_in_ranges(pos, self.ranges())

    def contains_span(self, other):
        if not self.same_context(other):
            return False
        if other.a == other.b == self.a:
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

    def spans_origin(self):
        if self._does_wrap_origin and self.cyclic and self.a == self.b:
            return True
        return self.b < self.a and self.cyclic

    def __contains__(self, other):
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

    def __getitem__(self, val):
        if isinstance(val, int):
            if not self.cyclic and (val >= len(self) or val < -len(self)):
                raise IndexError(
                    "Index '{}' outside of linear span with length of {}".format(
                        val, len(self)
                    )
                )
            else:
                if val < 0:
                    return self.t(val + self.b) - self.index
                else:
                    return self.t(val + self.a) - self.index
        elif issubclass(type(val), slice):
            if val.step:
                raise ValueError(
                    "{} slicing does not support step.".format(self.__class__.__name__)
                )
            if val.start is None:
                i = self.a
            else:
                i = self[val.start]

            if val.stop is None:
                j = self.b
            else:
                j = self[val.stop]
            if i == j and self._does_wrap_origin:
                return self.new(i, j, does_wrap_origin=self._does_wrap_origin)
            return self.new(i, j)
        elif isinstance(val, tuple):
            if len(val) > 2:
                raise ValueError(
                    "{} -- copying does only supports (start, stop)".format(val)
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
        return "<{}={} {} {} start={}>".format(
            self.__class__.__name__,
            id(self),
            (self.a, self.b, self.context_length),
            self.cyclic,
            self.index,
        )

    def __str__(self):
        return "<{} {} {} start={}>".format(
            self.__class__.__name__,
            (self.a, self.b, self.context_length),
            self.cyclic,
            self.index,
        )


class EmptySpan(Span):
    def ranges(self):
        return [(self.bounds()[0], self.bounds()[0])]
