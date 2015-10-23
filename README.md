dgesv
=====

This is a smaller collection of lapack routines, containing mainly:

- dgesv

The collection is based on clapack (3.2.1), but the f2c header has been
modified.

Defines for min/max/abs/dabs have been replaced by fmin/fmax,labs,fabs.
The type definitions use the same data types as OpenModelica (int on
x86_64 instead of long).
