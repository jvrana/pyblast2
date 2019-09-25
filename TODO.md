# alignments with circular templates

1. commit to using 100% JSON for all sequences
2. converting circular sequences for things with gaps doesnt really make much sense
3. refactor span to have nwraps
4. convert span to slices for np indexing
5. make inherited class called WrappingSpan that allows wrapped spans
