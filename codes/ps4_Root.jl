####################################################################
####################################################################
############# Advance Macroeconomics II PS 4
####################################################################
####################################################################



cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/")


######### Run these commented commands only the first time
using Pkg
Pkg.activate("ps4")




#Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
#Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
#Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
#Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
#Pkg.add("Latexify")
#Pkg.add("LinearAlgebra")
#Pkg.add("PrettyTables")
#Pkg.add("GR")
#Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
#Pkg.add("Plotly")
#Pkg.add("IntervalArithmetic")
#Pkg.add("IntervalRootFinding")
#Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
#Pkg.add("Random")
#Pkg.add("Interpolations")
#Pkg.add("Distributions")

using Optim 
    using Optim: converged, maximum, maximizer, minimizer, iterations, optimize
using Parameters
using Roots
using IntervalArithmetic, IntervalRootFinding
using Plots
using Dierckx 
using Plots
using LinearAlgebra
using PrettyTables
#using Distributions
#using Random
using Interpolations
using ForwardDiff 

cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps4/")
# Call Scaled Interpolation Functions
include("Scaled_Interpolation_Functions.jl")

#include("functions.jl")

@with_kw struct Par  #  The macro @with_kw which decorates a type definition to allow default values and a keyword constructor, in this case Par
    # This are the structural parameters which wont move throughout of the model
    #z::Float64 = 1 ; # Productivity
    z::Float64 = 1.0 ; # Productivity for the increase for plotting
    α::Float64 = 1/3; # Capital share, to write the alpha you do it as in latex \alp ha
    β::Float64 = 0.98; # Discount factor
    η::Float64 = 1.0; # Parameter for labor
    σ::Float64 = 2.0; # Parameter for the risk aversion
    δ::Float64 = 0.05; # Depreciation rate
    χ:: = (z*(1-α)*(((z*α*β)/(1-β+β*δ))^(1/(1-α)))^α)/(l^(σ+η)*(z*(((z*α*β)/(1-β+β*δ))^(1/(1-α)))^α-δ*(((z*α*β)/(1-β+β*δ))^(1/(1-α)))^σ)
    # VFI values
    max_iter::Int64 = 2000 ; # Maximum number of iterations
    dist_tol::Float64 = 1E-9 ; # Tolerance for the distance (tol for the fixed point)
    # Howard's Policy iterations
    H_tol::Float64 = 1E-9; # Tolerance for policy function iteration.
    N_H::Int64        = 20    ; # Maximum number of policy iterations
    # Minimum consumption and labor for numerical optimization
    c_min::Float64    = 1E-16
    l_min::Float64    = 1E-16
    l_max::Float64    = 1.0
end


p=Par()

# 1. Solve the planner's problem using VFI. Treat all choice variables and states as continuous. 

# a) Direct maximization over a pair of variables: (c,l) or (k',l).
# b) Direct maximization over single variable (k') and the other varibles solved for analytically
# c) Solving the FOC

################################################################################
################################################################################
################################################################################



################################################################################
#                                 Functions
################################################################################

# Steady state values
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β, σ, η, δ = p
    l_ss = ((z*(1-α)*Ψ^α)/(χ*Λ^σ))^(1/(σ+η))
    k_ss = Ψ*l_ss
    c_ss = Λ*l_ss
    y_ss = z*l_ss*Ψ^α
    r_ss = α*z*Ψ^(α-1)
    w_ss = (1-α)*z*Ψ^α
    return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
end

# Utility function
function utility(k,x::Array,p::Par)
#function utility(k,kp,l,p::Par)
    @unpack α,δ,σ,z,η,c_min = p
    kp = x[1]
    l = x[2]
    #u = ((z*(k^α)*l^(1-α)+(1-δ)*k - kp)^(1-σ))/(1-σ) - χ*((l^(1+η))/(1+η))
    u=0
    c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
    if c>c_min
        u_c = (c^(1-σ))/(1-σ)
    else
        u_c = (c_min^(1-σ))/(1-σ)
    end
    # Now the labor part 
    u = u_c - - χ*((l^(1+η))/(1+η))
    return u
end


function Euler_Labor(k,x,p::Par)
    @unpack α,δ,σ,z,η = p
    kp = x[1]
    l = x[2]
    c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
    return -χ*l^η + (1/c^σ)*(1-α)*z*k^α*l^(-α)
end

function Euler_capital(k,x1::Array,x2::Array, p::Par)
    @unpack α,δ,σ,z,η,β = p
    kp = x1[1]
    l = x1[2]
    kpp = x2[1]
    lp = x2[2]
    c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
    cp = z*(kp^α)*lp^(1-α)+(1-δ)*kp - kpp
    return -(1/c^α)+(1/cp^α)*β*(α*z*(lp/kp)^(1-α) + (1-δ))
end


