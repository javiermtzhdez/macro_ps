####################################################################
####################################################################
############# Advance Macroeconomics II PS2
####################################################################
####################################################################

# As in Fortran, I think is better to do variable declaration first and then the code...

# IMPLICIT NONE HAHA

######################
# [specification part]
######################

######### This part is just the environment activation
cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/")

######### Run these commented commands only the first time
using Pkg
Pkg.activate("ps2")


#Pkg.add("JuMP")
#Pkg.add("GR")
#Pkg.add("Roots")
#Pkg.add("Plotly")

#Pkg.add("IntervalArithmetic")
#Pkg.add("IntervalRootFinding")

using Parameters
using Plots
using JuMP
using Roots
using IntervalArithmetic, IntervalRootFinding

cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps2/")
#mkpath("Figures") # Creates the folder Figures ensuring the appropriate path exists


@with_kw struct Par  #  The macro @with_kw which decorates a type definition to allow default values and a keyword constructor, in this case Par
    # This are the structural parameters which wont move throughout of the model
    #z::Float64 = 1 ; # Productivity
    z::Float64 = 1.0 ; # Productivity for the increase for plotting
    α::Float64 = 1/3; # Capital share, to write the alpha you do it as in latex \alp ha
    β::Float64 = 0.98; # Discount factor
    η::Float64 = 1.0; # Parameter for labor
    σ::Float64 = 2.0; # Parameter for the risk aversion
    δ::Float64 = 1.0; # Depreciation rate
    l::Float64 = 0.4; # Steady State labor supply
    # VFI values
    max_iter::Int64 = 2000 ; # Maximum number of iterations
    dist_tol::Float64 = 1E-9 ; # Tolerance for the distance (tol for the fixed point)
    # Howard's Policy iterations
    H_tol::Float64 = 1E-9; # Tolerance for policy function iteration.
end

# Allocate parameters to object p for future calling? why not calling them directly? 
# Is this like creating a vector of parameters
p=Par()


####################### Q. 4 Find χ such that l_ss = 0.4
# First we know this values
function psi(p::Par)
    @unpack z, α, β, δ = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    #psi =  ( (1 - β*(1-1) ) / (α*z*β)  )^(1/(α - 1))
    psi =  ( (1 - β*(1-δ) ) / (α*z*β)  )^(1/(α - 1))
    return psi
end
global Ψ =  psi(p)
println("Psi is")
println(Ψ)

# and 
#Λ = z*Ψ^α- δ*Ψ
function lambda(p::Par)
    @unpack z, α, β, δ = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    #psi =  ( (1 - β*(1-1) ) / (α*z*β)  )^(1/(α - 1))
    psi =  ( (1 - β*(1-δ) ) / (α*z*β)  )^(1/(α - 1))
    return z*(psi)^α - δ*psi
end
global Λ =  lambda(p)

println("Lambda is")
println(Λ)

# Since we know what the function for the labor in steady state:

function lss(a,b,c, p::Par) # where a represents Ψ and b represents Λ, and c represents χ
    @unpack z, α, σ, η = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    return ((z*(1-α)*a^α)/(c*b^σ))^(1/(σ + η))
end

lss(Ψ,Λ,χ, p)



# We know that l_ss = 0.4, then we can invert the function to get χ
function f1(a,b,p::Par)# where a represents Ψ and b represents Λ
    @unpack z, α, σ, η = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    return (z*(1-α)*a^α)/(((0.4)^(σ+η))*b^σ)
end

# Then the value of χ that satisfies l_ss = 0.4 is
global χ = f1(Ψ, Λ, p::Par)
println("Chi = ")
println(χ)

# So the way the problem is written, I need some values (i.e. the steady state values for later, I might as well get them now)

# Steady state values
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β, l, σ, η = p
    l_ss = ((z*(1-α)*Ψ^α)/(χ*Λ^σ))^(1/(σ+η))
    k_ss = Ψ*l_ss
    c_ss = Λ*l_ss
    y_ss = z*l_ss*Ψ^α
    r_ss = α*z*Ψ^(α-1)
    w_ss = (1-α)*z*Ψ^α
    return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
end

l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)

# So far I got the steady state values here
# Now I am going to calculate the instantaneous utility function

function utility(k,kp,l, p::Par)
    @unpack z, α, η, σ, δ = p
    u = ((z*(k^α)*l^(1-α)+(1-δ)*k - kp)^(1-σ))/(1-σ) - χ*((l^(1+η))/(1+η))
    if u>0
        return u
    else
        return -Inf
    end
end

# Now, the grid between (10^-5) and 2k_ss
function k_grid_make(lb,ub,n_k)
    grid = range(lb,ub;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
    return grid
end
#k_grid_make(1E-5,2*k_ss,100)
#collect(k_grid_make(1E-5,2*k_ss,100))

# Grid for labor
l_grid = range(0,1,length=11)
#collect(l_grid)






function utility(k,kp,p::Par)
    @unpack z, α = p
    c = z*k.^α  - kp
    if c>0
    return log(c)
    else
    return -Inf
    end
end

























