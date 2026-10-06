# Bundled libsais

Revision: `ce90878d784b5ff7d019300535675e4a2e22aae0`.
Source: <https://github.com/IlyaGrebnov/libsais>.
License: `LICENSE` (Apache-2.0).

The imported source and headers are unmodified. STAR builds `libsais.c` and
`libsais64.c` with `CC` and OpenMP; the 64-bit entry points use the 32-bit
implementation for smaller inputs.
