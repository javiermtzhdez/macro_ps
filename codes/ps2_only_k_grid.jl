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
#using JuMP
#using Roots
#using IntervalArithmetic, IntervalRootFinding

cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps2/")
#mkpath("Figures") # Creates the folder Figures ensuring the appropriate path exists

#include("functionsps2.jl")
#using .functionsps2

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

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

########################################## FUNCTIONS USED
# Find a way to use them in a module with the parameters... like a global in fortran...


function psi(p::Par)
    @unpack z, α, β, δ = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    #psi =  ( (1 - β*(1-1) ) / (α*z*β)  )^(1/(α - 1))
    psi =  ( (1 - β*(1-δ) ) / (α*z*β)  )^(1/(α - 1))
    return psi
end


function lambda(p::Par)
    @unpack z, α, β, δ = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    #psi =  ( (1 - β*(1-1) ) / (α*z*β)  )^(1/(α - 1))
    psi =  ( (1 - β*(1-δ) ) / (α*z*β)  )^(1/(α - 1))
    return z*(psi)^α - δ*psi
end

function lss(a,b,c, p::Par) # where a represents Ψ and b represents Λ, and c represents χ
    @unpack z, α, σ, η = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    return ((z*(1-α)*a^α)/(c*b^σ))^(1/(σ + η))
end

# We know that l_ss = 0.4, then we can invert the function to get χ
function f1(a,b,p::Par)# where a represents Ψ and b represents Λ
    @unpack z, α, σ, η = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    return (z*(1-α)*a^α)/(((0.4)^(σ+η))*b^σ)
end


# Steady state values
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β, l, σ, η, δ = p
    l_ss = ((z*(1-α)*Ψ^α)/(χ*Λ^σ))^(1/(σ+η))
    k_ss = Ψ*l_ss
    c_ss = Λ*l_ss
    y_ss = z*l_ss*Ψ^α
    r_ss = α*z*Ψ^(α-1)
    w_ss = (1-α)*z*Ψ^α
    return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
end


# So far I got the steady state values here
# Now I am going to calculate the instantaneous utility function

function utility(k,kp,l,p::Par)
    @unpack z, α, η, σ, δ = p
    u = ((z*(k^α)*l^(1-α)+(1-δ)*k - kp)^(1-σ))/(1-σ) - χ*((l^(1+η))/(1+η))
    c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
    #y = (a,b) -> z*(a^α)*b^(1-α)
    #c = (a,b,d) -> y(a,b)-d
    #lb = (l) -> (χ*l.^(η+1)/(η+1))
    #utils = c(k,l,kp).^(1-σ)/(1-σ)-lb(l)
    if c>0
    return u
    else
    return -Inf
    end
end

function labor_allocation(k,kp,p::Par)
    @unpack z, α, η, σ, δ = p
    inner_k_grid = range(0,1,length=100)
    L_array = [utility(k,kp,inner_k_grid[i],p) for i in 1:100]
    index = findmax(L_array)[2] # Bc findmax get you  the value and the index 
    l = inner_k_grid[index]
    return l
end


# Now, the grid between lower bound and upper bound
function k_grid_make(lb,ub,n_k)
    grid = range(lb,ub;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
    return grid
end

function Make_K_Grid(n_k,p::Par)
    # Get SS
    l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    # Get k_grid
    k_grid = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
    # Return
    return k_grid
end

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


# Define function for Value update and policy functions
function T_grid_loop(V_old,k_grid,p::Par)
    @unpack z, α, β, δ = p
    n_k  = length(k_grid)
    V    = zeros(n_k)
    G_kp = fill(0,n_k)
    G_l = fill(0.0,n_k)
    for i = 1:n_k
        V_aux = zeros(n_k) ;
        l_aux = zeros(n_k) # Empty vector for auxiliary value of V(i,j)
        for j = 1:n_k
            # Evaluate potential value function for combinations of
            # current capital k_i and future capital k_j
            l_aux[j] = labor_allocation(k_grid[i],k_grid[j],p)
            V_aux[j] = utility(k_grid[i],k_grid[j],l_aux[j],p) + β*V_old[j]
            #println(V_aux[j]," ",k_grid[i]," ",k_grid[j]," ",utility(k_grid[i],k_grid[j],z,a,b) + b*V_old[j])
        end
        # Choose maximum value given current capital k_i
        V[i], G_kp[i] = findmax(V_aux)
        G_l[i] = l_aux[G_kp[i]]
    end
    return V, G_kp, G_l
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


# Euler Equation
function Euler_Error(k,kp,kpp,l,lp,p::Par)
    # Return percentage error in Euler equation
    @unpack z, α, β, σ, δ = p
    LHS = 1/(z*k^α*l^(1-α)-kp)^σ
    RHS = β*α*z*kp^(α-1)*lp^(1-α)/(z*kp^α*lp^(1-α)-kpp)^σ
    return (RHS/LHS-1)*100
end

# Euler error by number of points in grid
function Euler_grid(k_grid::AbstractArray,G_kp::AbstractArray)
    # Euler error of numerical policy function on grid
    n_k = length(k_grid)
    Euler = zeros(n_k)
    for i=1:n_k
        k   = k_grid[i]
        kp  = k_grid[G_kp[i]]
        kpp = k_grid[G_kp[G_kp[i]]]
        l = labor_allocation(k,kp,p)
        lp = labor_allocation(kp,kpp,p)
        Euler[i] = Euler_Error(k,kp,kpp,l,lp,p)
    end
    # Return
    return Euler
end


#-----------------------------------------------------------
    # VFI - Grid Search - Matrices - Howard's Policy Iteration

    # Define function for Value update and policy functions
    function HT_grid_mat(V_old,U_mat,k_grid,p::Par,n_H)
        @unpack z, α, β, H_tol = p
        # Get Policy Function
        n_k    = length(V_old)
        V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
        V_old  = V
        # "Optimal" U for Howard's iteration
            U_vec = U_mat[G_kp]
        # Howard's policy iteration
        # G_kp is a Cartesian Index
        for i=1:n_H
            V = U_vec .+ β*repeat(V_old',n_k,1)[G_kp]
            if maximum(abs.(V./V_old.-1))<=H_tol
                break
            end
            V_old = V
        end
        # Recover Policy Functions
        G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
        G_l = labor_allocation.(k_grid,k_grid[G_kp],fill(p,n_k))

        # Return output
        return V, G_kp, G_l
    end


# Solve VFI with Howard's policy iteration
    function Solve_VFI_HPI(n_H,n_k,p::Par)
        println("----------------------------------------------------")
        println("Value Function Iteration - Howard's Policy Iteration")
        # Get Grid
        k_grid = Make_K_Grid(n_k,p)
        # Utility matrix
        U_mat = [utility(k_grid[i],k_grid[j],labor_allocation(k_grid[i],k_grid[j],p),p) for i in 1:n_k, j in 1:n_k]
        # Solve VFI
        V, G_kp, G_l = VFI_grid(x->HT_grid_mat(x,U_mat,k_grid,p,n_H),k_grid,p)
        # Return Solution
        return V,G_kp, G_l, k_grid
    end




#-----------------------------------------------------------


# Fixed point with MPB
function VFI_grid_MPB(T::Function,k_grid,p::Par)
    @unpack β, max_iter, dist_tol = p
    # Initialize variables for loop
    n_k    = length(k_grid) ; # Number of grid nodes
    V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
    iter   = 0              ; # Iteration index
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - Grid Search - MPB - n_k=$n_k")
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_l = T(V_old)
        # MPB and Distance
        MPB_l  = β/(1-β)*minimum(V_new-V_old)
        MPB_h  = β/(1-β)*maximum(V_new-V_old)
        V_dist = MPB_h - MPB_l
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check Convergence
        if (V_dist<=dist_tol)
            # Recover value and policy functions
            V = V_old .+ (MPB_l+MPB_h)/2
            # Return
            println("VFI - Grid Search - MPB - n_k=$n_k")
            println("Iterations = $iter and Distance = ",100*V_dist)
            println("------------------------")
            println(" ")
            return V, G_kp, G_l
        end
    end
    # Report error for non-convergence
    error("Error in VFI - Grid Search - MPB - Solution not found")
end

# Define function for Value update and policy functions
    function T_grid_mat(V_old,U_mat,k_grid,p::Par)
        @unpack z, α, β, δ = p
        n_k    = length(V_old)
        V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
        G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
        G_l = labor_allocation.(k_grid,k_grid[G_kp],fill(p,n_k))
        return V, G_kp, G_l
    end


# Solve VFI with MPB
function Solve_VFI_MPB(n_k,p::Par)
    println("----------------------------------------------------")
    println("Value Function Iteration - MacQueen-Porteus Bounds")
    # Get Grid
    k_grid = Make_K_Grid(n_k,p)
    # Utility matrix
    U_mat = [utility(k_grid[i],k_grid[j],labor_allocation(k_grid[i],k_grid[j],p),p) for i in 1:n_k, j in 1:n_k]
    # Solve VFI
    V, G_kp, G_l = VFI_grid_MPB(x->T_grid_mat(x,U_mat,k_grid,p),k_grid,p)
    # Return Solution
    return V,G_kp, G_l, k_grid
end




##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


####################### Q. 4 Find χ such that l_ss = 0.4
# First we know this values
global Ψ =  psi(p)
println("Psi is")
println(Ψ)

# and 
#Λ = z*Ψ^α- δ*Ψ
# Lambda
global Λ =  lambda(p)

println("Lambda is")
println(Λ)

# Then the value of χ that satisfies l_ss = 0.4 is
global χ = f1(Ψ, Λ, p::Par)
println("Chi = ")
println(χ)

# Since we know what the function for the labor in steady state:
lss(Ψ,Λ,χ, p)

# So the way the problem is written, I need some values (i.e. the steady state values for later, I might as well get them now)

l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)


# Grid for capital 
#k_grid_make(1E-5,2*k_ss,100)
#collect(k_grid_make(1E-5,2*k_ss,100))
# Grid for labor
# l_grid = range(0,1,length=11)
#collect(l_grid)


@time V_20, G_kp_20, G_c_20, k_grid_20 = Solve_VFI_loop(20,p)
@time V_50, G_kp_50, G_c_50, k_grid_50 = Solve_VFI_loop(50,p)
@time V_200, G_kp_200, G_c_200, k_grid_200 = Solve_VFI_loop(200,p)
@time V_1000, G_kp_1000, G_c_1000, k_grid_1000 = Solve_VFI_loop(1000,p)