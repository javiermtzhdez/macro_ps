####################################################################
####################################################################
############# Advance Macroeconomics II PS 7
############# Part 1
####################################################################
####################################################################




cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/")


######### Run these commented commands only the first time
using Pkg
Pkg.activate("ps7")
#Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
#Pkg.add("Random")
#Pkg.add("Distributions")
#Pkg.add("GR")
#Pkg.add("Plotly")
#Pkg.add("StatsBase") # https://github.com/JuliaStats/StatsBase.jl

cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps5")

using Random, Distributions
using Parameters
using Statistics
using Plots
using LinearAlgebra
using Latexify
using PrettyTables
using Interpolations
using Dierckx
using ForwardDiff
using Optim
     using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots
include("Scaled_Interpolation_Functions.jl")


# For replicability purposes I will set the seed at 123
Random.seed!(123)

