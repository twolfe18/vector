
Note:
------
This package now serves as the sparse vector implementation in [prim](https://github.com/mgormley/prim),
which should be preferred over this code, which will not be maintained.


vector
======

a humble vector implementation

This is an implementation of a vector that can be sparse or dense,
and provides some basic stuff you need to work with them. Most of the
operations provided are implemented basically as efficiently as you
could hope to provide them, there shouldn't be any big gotchas. I have
played around with a few implementations before this, and I think this
is the way to go (down to things like not having a separate type/class
for sparse and dense versions, which tends to make things overly complicated).

There are some cool tricks and features provided including:
* "tagged" features indexes (you get two ints to specify a value, helpful for
  things like feature templates and domain adaptation).
* efficient dot products (ok this one is pretty basic...)
* side-effect-free functions to avoid mutation when feasible

