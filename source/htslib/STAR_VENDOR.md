# Bundled HTSlib subset

HTSlib revision: `4a11e1ff234b2d6f2105ddccc4ac4c298f425522` (1.22 development).
HTScodecs revision: `ce66e5f303862aad325e4df908636a31c877b4d8` (1.6.3).
Licenses: `LICENSE` and `htscodecs/LICENSE.md`.

This subset builds STAR's existing BAM I/O, the shared BGZF compression pool,
and indexed VCF/BCF reading for embedded consensus. Runtime C sources and headers
are unmodified. The Makefile retains the upstream static-library build rules and
omits unused statistical, realignment, region-index, VCF utility and external
CRAM accessor objects. Generic HTS file opening requires SAM/CRAM internals and
their codecs even for VCF/BCF consumers. Platform-specific codec sources remain.

Tests, examples, standalone tools, JavaScript codecs, cloud plugins, reference
cache programs, distribution tooling and unrelated documentation are omitted.
STAR uses its existing Makefile and one HTSlib ABI; `htslib_static.mk` supplies
the native link dependencies. This is an embedded library, not a standalone
HTSlib distribution.
