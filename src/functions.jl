# Call required packages
using Comodo
# using GLMakie
# using GeometryBasics
# using FileIO
# using Statistics
# using LinearAlgebra
# using Geogram

"""
    comododir()

# Description 

This function simply returns the string for the Comodo path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function simulatedamputationdir()
    joinpath(@__DIR__, "..")
end
