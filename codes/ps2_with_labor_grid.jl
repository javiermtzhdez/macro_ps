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


# We know that l_ss = 0.4, then we can invert the function to get χ
function f1(a,b,p::Par)# where a represents Ψ and b represents Λ
    @unpack z, α, σ, η = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    return (z*(1-α)*a^α)/(((0.4)^(σ+η))*b^σ)
end

# Then the value of χ that satisfies l_ss = 0.4 is
global χ = f1(Ψ, Λ, p::Par)
println("Chi = ")
println(χ)


lss(Ψ,Λ,χ, p)



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

function consumption(k,kp,l,p::Par)
    @unpack z, α, η, σ, δ = p
    c = ((z*(k^α)*l^(1-α)+(1-δ)*k - kp)^(1-σ))/(1-σ)
    if c>0
        return c
    else
        return -Inf
    end
end

function labor_decision(k,kp,p::Par)
    @unpack z, α, η, σ, δ = p
    inner_k_grid = range(0,1,length=100)
    L_array = [utility(k,kp,inner_k_grid[i],p) for i in 1:100]
    index = findmax(L_array)[2] # Bc findmax get you  the value and the index 
    l = inner_k_grid[index]
    return l
end


function Make_K_Grid(n_k,p::Par)
    # Get SS
    l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    # Get k_grid
    k_grid = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
    # Return
    return k_grid
end



# Define function for Value update and policy functions
# This is the contraction V_n+1(k) = TV_n(k) = max u(c,l) + beta V_n (k')
# So what we are doing here is just generating the new iteration of the value function.
# The inputs are the old value function (V_n(k) in the terminology I am using), the grid for k (because its grid search)
function T_grid_loop(V_old,k_grid,p::Par)
    @unpack z, α, β, δ = p
    n_k  = length(k_grid)
    V    = zeros(n_k)
    G_kp = fill(0,n_k)
    G_l = fill(0.0,n_k)
    for i = 1:n_k   # This calculates the new iteration of the  value function for all points in the grid.
        V_aux = zeros(n_k) ;
        l_aux = zeros(n_k) # Empty vector for auxiliary value of V(i,j)
        for j = 1:n_k
            # Evaluate potential value function for combinations of
            # current capital k_i and future capital k_j
            l_aux[j] = labor_decision(k_grid[i],k_grid[j],p)
            V_aux[j] = utility(k_grid[i],k_grid[j],l_aux[j],p) + β*V_old[j]
            #println(V_aux[j]," ",k_grid[i]," ",k_grid[j]," ",utility(k_grid[i],k_grid[j],z,a,b) + b*V_old[j])
        end
        # Choose maximum value given current capital k_i
        V[i], G_kp[i] = findmax(V_aux)
        G_l[i] = l_aux[G_kp[i]]
    end
    return V, G_kp, G_l
end



# Now, what we are gonna do is iterate the contraction T until I get the fixed point which is the solution of the contraction

# Solve VFI with grid search and loops
function VFI_grid(T::Function,k_grid,p::Par)
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    n_k    = length(k_grid) ; # Number of grid nodes
    V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - Grid Search - n_k=$n_k")
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_l = T(V_old)
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        iter  += 1
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol
            println("VFI - Grid Search - n_k=$n_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            return V_new, G_kp, G_l
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Grid Search - Solution not found")
end



# Solve VFI with grid search and loops
function Solve_VFI_loop(n_k,p::Par)
    println("----------------------------------------------------")
    println("Value Function Iteration - Grid Search with Loops")
    # Get Grid
    k_grid = Make_K_Grid(n_k,p)
    # Solve VFI
    V, G_kp, G_l = VFI_grid(x->T_grid_loop(x,k_grid,p),k_grid,p)
    # Return Solution
    return V,G_kp, G_l, k_grid
end


@time V_20, G_kp_20, G_c_20, k_grid_20 = Solve_VFI_loop(20,p)