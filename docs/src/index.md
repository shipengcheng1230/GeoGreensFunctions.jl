# Introduction

This package provides commonly used Green's functions for geoscience research.


All functions are translated from their respective original sources, pertaining all original parameter names and coordinate system.
A transform to the ENU (x -> east, y -> north, z -> upward) coordinate is provided if necessary. All methods here are type stable.


While trying my best to test against the original implementations and cross validate
among them, this by no means guarantees them to be error-free. Due to the complexities
in the equations, it is currently difficult to refactor the code myself for better
clarity and performance.

*Some of them hasn't been implemented yet.*

If you encounter any problem, please don't hesitate file an issue with your minimum working example(s).
If you think something can be improved, PR is always appreciated!

## Known Issues

- [issue 15276](https://github.com/JuliaLang/julia/issues/15276) for type stability involving nested function

- segment dislocation sometimes does not match [`dc3d`](@ref), which is reflected by some random tests

- broadcast is mostly not enabled, which significantly increases the compile time

- lots of repeated formula which can be refactored to simplify generated codes
