This repository contains the code associated with the paper RedZeD: Computing persistent homology by Reduction to Zero Differentials by Chris Kapulkin and Nathan Kershaw.

The code is contained in the redzed.jl file and the primary function is redzed(dist, n; threshold). 
The inputs to redzed(dist, n; threshold) are:
    dist: the distance matrix 
    n: the max homology dimension 
    threshold: the radius threshold (default = Inf)
The output is a dictionary of the form: dimension => persistence pairs of that dimension
