####################################################################
####################################################################
############# Advance Macroeconomics II PS 6
############# Part 1
####################################################################
####################################################################




cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/")


######### Run these commented commands only the first time
using Pkg
Pkg.activate("ps6")
Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("GR")
#Pkg.add("Plotly")
Pkg.add("StatsBase") # https://github.com/JuliaStats/StatsBase.jl
Pkg.add("Statistics")
Pkg.add("Plots")
Pkg.add("LinearAlgebra")
Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl


cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps6")

using Random, Distributions
using Parameters
using Statistics
using Plots
using LinearAlgebra
#using Latexify
#using PrettyTables
using Interpolations
using Dierckx
using ForwardDiff
using Optim
     using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots
include("Scaled_Interpolation_Functions.jl")


# For replicability purposes I will set the seed at 123
Random.seed!(123)
mkpath("Figures")




# Parameters
    # Generate structure for parameters using Parameters module
    # Set default values for parameters
    @with_kw struct Par
        # Model Parameters
        z_ave::Float64 = 1; # Reference level for productivity ------ it was z_bar
        α::Float64 = 1/3  ; # Production function
        β::Float64 = 0.98 ; # Discount factor
        σ::Float64 = 2 ; # consumption elasticity of subsitution
        η::Float64 = 1 ; # labor/leisure elasticity of substitution
        δ::Float64 = 0.05 ; # Depreciation rate of capital
        ρ::Float64 = 0.90 ; # Persistance of AR(1) productivity process: log(z')=ρlog(z)+η
        σ_η::Float64 = 0.1 ; # Variance of innovations for productivity where η ∼ N(0,σ_η)
        χ::Float64  = 8.217267384235548; # Parameter for labor disutility
        Ψ::Float64  = 10.301099814055622; # Value needed for calculations of χ
        Λ::Float64 = 1.6607895618579493; # Value needed for calculations of χ
        # VFI Paramters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9  ; # Tolerance for distance between current and previous value functions
        # Policy functions
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Minimum consumption for numerical optimization
        c_min::Float64    = 1E-16
    end
# Allocate parameters to object p for future calling
p = Par()


################################################################################
################################################################################
#                                 Functions
################################################################################
################################################################################

######################################################################
############## Steady State
######################################################################
# Steady state values
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z_ave, α, β, σ, η, δ, χ, Λ, Ψ = p
    z = z_ave
    l_ss = ((z*(1-α)*Ψ^α)/(χ*Λ^σ))^(1/(σ+η))
    k_ss = Ψ*l_ss
    c_ss = Λ*l_ss
    y_ss = z*l_ss*Ψ^α
    r_ss = α*z*Ψ^(α-1)
    w_ss = (1-α)*z*Ψ^α
    return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
end

#

# Test steady state function
k_ss,y_ss,c_ss,r_ss,w_ss,l_ss = SS_values(p)
println(" ")
println("------------------------")
println(" Steady State variables")
println("   Quantities: k = $k_ss; y = $y_ss; c = $c_ss;")
println("   Prices:     r = $r_ss; w = $w_ss;")
println("------------------------")
println(" ")





######################################################################
############## K_grid
######################################################################

    # Function to make grid for capital
    function Make_K_Grid(n_k,θ_k,p::Par)
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss,l_ss = SS_values(p)
        # Lower and upper bounds
        lb = 1E-5
        ub = 2*k_ss
        # Get k_grid
        if θ_k≠1
            k_grid = PolyRange(1E-5,2*k_ss;θ=θ_k,N=n_k)
        else
        k_grid = range(1E-5,2*k_ss,length=n_k)
        end
        # Return
        return k_grid
    end

#


######################################################################
############## Rouwenhorst
######################################################################

    function Rouwenhorst95(N,p::Par)
        @unpack ρ,σ_η=p
        # INPUTS:
            # ρ: persistence of unerlying AR(1) process where log(z') = ρlog(z)+η
            # σ_z: Std dev of inovation η in AR(1) process where η∼N(0,σ^2)
            # N: Size of grid for discrete process
        # OUTPUT:
            # z: All possible values of discretized AR(1) process, equally spaced grid of size N
            # Π: Matrix of transition probabilities
            # PDF_z: Stationary PDF of z
        #---------------------------------------------------------------------------
        Π = zeros(N,N)
        Π_Nm = zeros(N-1,N-1)
        P = (1+ρ)/2
        ϕ = σ_η*(sqrt((N-1)/(1-ρ^2)))
        z = range(-ϕ,ϕ;length=N)
        if N==2
            Π = [P 1-P;1-P P]
        else
            Π_Nm = Rouwenhorst95(N-1,p)[2]
            o = zeros(N-1)
            Π = P*[Π_Nm o; o' 0] + (1-P)*[o Π_Nm; 0 o'] + (1-P)*[o' 0; Π_Nm o] + P*[0 o';o Π_Nm]
            Π = Π./repeat(sum(Π,dims=2),1,N)
        end
        PDF_z = pdf.(Binomial(N-1,0.5),(0:N-1))
        return (z,Π,PDF_z)
    end

#