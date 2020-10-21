####################################################################
####################################################################
############# Advance Macroeconomics II PS 4

# I attempt to solve this one with a nested problem

# The logic behind it comes from Adda and Cooper textbook and the updated slides

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
mkpath("Figures")


# Parameters

@with_kw struct Par  #  The macro @with_kw which decorates a type definition to allow default values and a keyword constructor, in this case Par
    # This are the structural parameters which wont move throughout of the model
    #z::Float64 = 1 ; # Productivity
    z::Float64 = 1.0 ; # Productivity for the increase for plotting
    α::Float64 = 1/3; # Capital share, to write the alpha you do it as in latex \alp ha
    β::Float64 = 0.98; # Discount factor
    η::Float64 = 1.0; # Parameter for labor
    σ::Float64 = 2.0; # Parameter for the risk aversion
    δ::Float64 = 0.05; # Depreciation rate
    χ::Float64  = 8.217267384235548;
    Ψ::Float64  = 10.301099814055622;
    Λ::Float64 = 1.6607895618579493;
    # VFI values
    max_iter::Int64 = 2000 ; # Maximum number of iterations
    dist_tol::Float64 = 1E-9 ; # Tolerance for the distance (tol for the fixed point)
    # Howard's Policy iterations
    H_tol::Float64 = 1E-9; # Tolerance for policy function iteration.
    N_H::Int64        = 20    ; # Maximum number of policy iterations
    # Minimum consumption and labor for numerical optimization
    c_min::Float64    = 1E-16
    #l_min::Float64    = 1E-16
    l_max::Float64    = 1.0
end


p=Par()



######################################################################
############## Steady state values
######################################################################

function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β, σ, η, δ, χ, Ψ, Λ = p
    l_ss = ((z*(1-α)*Ψ^α)/(χ*Λ^σ))^(1/(σ+η))
    k_ss = Ψ*l_ss
    c_ss = Λ*l_ss
    y_ss = z*l_ss*Ψ^α
    r_ss = α*z*Ψ^(α-1)
    w_ss = (1-α)*z*Ψ^α
    return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
end

l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)

######################################################################
############## K_grid
######################################################################
function Make_K_Grid(n_k,θ_k,p::Par,scale_type="Poly")
    # Get SS
    k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    # Get k_grid
    if θ_k≠1
        if scale_type=="Poly"
        k_grid = PolyRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
        elseif scale_type=="Exp"
        # k_grid = ExpRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
        else
        error("scale_type must be either Poly or Exp")
        end
    else
    k_grid = range(1E-5,2*k_ss,length=n_k)
    end
    # Return
    return k_grid
end


######################################################################
############## # Generate structure of model objects
######################################################################
@with_kw struct Model
    # Parameters
    p::Par = Par() # Model paramters in their own structure
    # Grids
    θ_k::Float64    = 1     # Curvature of k_grid
    n_k::Int64      = 20    # Size of k_grid
    n_k_fine::Int64 = 1000  # Size of fine grid for interpolation
    k_grid          = Make_K_Grid(n_k,θ_k,p)    # k_grid for model solution
    k_grid_fine     = Make_K_Grid(n_k_fine,1,p) # Fine grid for interpolation
    # Value and policy functions
    V         = Array{Float64}(undef,n_k)       # Value Function
    G_kp      = Array{Float64}(undef,n_k)       # Policy Function
    G_c       = Array{Float64}(undef,n_k)       # Policy Function
    G_l       = Array{Float64}(undef,n_k)       # Policy Function
    V_fine    = Array{Float64}(undef,n_k_fine)  # Value Function on fine grid
    G_kp_fine = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
    G_c_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
    G_l_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
end

M = Model()


######################################################################
############## Consumption Function
######################################################################
function consumption(k,kp,l,p::Par)
    @unpack z, α, β, δ = p
    c = z*k^α*l^(1-α) + (1-δ)*k - kp
    return c
end

######################################################################
############## Euler Labor (scalar)
######################################################################
function Euler_Labor(k, kp, l ,p::Par) # This returns the scalar value of th Euler for Labor
    @unpack α,δ,σ,z,η, χ = p
    c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
    return -χ*l^η + (1/c^σ)*(1-α)*z*k^α*l^(-α)
end

######################################################################
############## Labor lower bound for Bellman Objective 
######################################################################
function l_min(k,kp,c,p::Par)
    @unpack z, α, β, δ, c_min = p
    l_min = ((c_min + kp - (1-δ)*k)/(z*k^α))^(1/(1-α))
    return l_min
end

######################################################################
############## First order condition  -- function for Bellman Objective
############## THIS IS ALGORITHM 7 FROM SERGIO'S HANDOUT
######################################################################
#=
function FOC_l(l, k,kp,p::Par) # l is the only input, the rest are given
    @unpack z, α, β, δ, η, σ, χ = p
    #Compute implied consumption
    #c = l -> z*k^α*l^(1-α) + (1-δ)*k - kp
    # I will return the objective function for minimization
    #foc = l -> -χ*l^η + (1/c^σ)*(1-α)*z*k^α*l^(-α)
    foc = l -> -χ*l^η + (1/(z*k^α*l^(1-α) + (1-δ)*k - kp)^σ)*(1-α)*z*k^α*l^(-α)
    return foc
end
=#

FOC_l = x -> (Euler_Labor(k,kp,x,p))



######################################################################
############## Utility function
######################################################################

function utility(k,kp,l,p::Par)
    #function utility(k,kp,l,p::Par)
        @unpack α,δ,σ,z,η,c_min = p
        #u = ((z*(k^α)*l^(1-α)+(1-δ)*k - kp)^(1-σ))/(1-σ) - χ*((l^(1+η))/(1+η))
        u=0
        c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
        if c>c_min
            u_c = (c^(1-σ))/(1-σ)
        else
            u_c = (c_min^(1-σ))/(1-σ)
        end
        # Now the labor part 
        u = u_c - χ*((l^(1+η))/(1+η))
        return u
    end


######################################################################
############## Bellman Objective Function
######################################################################

function Bellman_Objective(k, kp, p::Par) # kp is 
    # the only input (or the initial guess).
    # This takes as given the k = k_grid[i], the old 
    # iteration of the Value fn = V_old and the Parameters = p
    #@unpack p, n_k, k_grid, V, G_kp, G_c = M
    @unpack z, α, β, δ, c_min, σ = p
    # Check consumption  with full labor and even full depreciation
    # c = f(k,1) - kp
    c = consumption(k,kp,1,p) # This is a fn that comes outside
    if c<c_min 
        return (c_min^(1-σ))/(1-σ) - a*((c-c_min)^2) 
        # !!!!!!!!!! Duda aquí, te regresa una función? 
    end
    # Solve for labor (inside bounds)
    # This is with a root finder for the FOC of labor given the bounds
    # for labor and the k and kp(initial guess) you have
    # l* = Root(FOC_l; l_min, 1, k',k) where l_min = max(0,f^-1(c_min + kp;k))
    # l_min will be a function too
    l_min = max(0,l_min(k,kp,c_min,p)) # This gives you the lower bound for labor
    # Now we use the roots part 
    FOC_l = x -> (Euler_Labor(k,kp,x,p)) # Objective function for the roots finder
    l_star = find_zero(FOC_l,(l_min,1),Roots.Brent())
    # Compute again implied consumption
    c = consumption(k,kp,l_star,p)
    B_obj_fun = -utility(k,kp,l_star,p)
    # B_obj_fun = -utility(k,kp,l_star,p) - β*V_old.(x)
    return B_obj_fun, l_star
end

######################################################################
############## Bellman operator
######################################################################

function T_cts_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c, G_l = M
    @unpack β,α,z,δ,η,σ,c_min, χ = p
    # Define the interpolation of the current value function (V) for values of Vp
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid()))))
    get_kp_max(k,l) = z*(k^α)*(l^(1-α)) + (1-δ)*k - c_min # Function because the max k' depends on the value of capital today k and labor l
    kp_min = k_grid[1]
    for i = 1:n_k
        # Here goes the Bellman Objective function
        kp_max = min(get_kp_max(k_grid[i],1.0),0.9999*k_grid[end])
        Obj_Fun = x->(Bellman_Objective(k_grid[i],x,V,p)[1] - β*Vp.(x))
        min_result = optimize(Obj_Fun,kp_min,kp_max)
        V[i] = -min_result.minimum
        G_kp[i] = min_result.minimizer
        G_l[i] = Bellman_Objective(k_grid[i],G_kp[i],p::Par)[2]
    end
    G_c = z.*(collect(k_grid).^α).*(G_l.^(1-α)) .- G_kp
    # Return results
    return V, G_kp, G_c, G_l
end


        

######################################################################
############## VFI Fixed Point
######################################################################
## VFI fixed point
function VFI_Fixed_Point(T::Function,M::Model)
    # Unpack model structure
    @unpack p, n_k, θ_k, k_grid = M
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    V_old  = zeros(n_k)     ; # Initialize value function, here I just do 0, not the best option
    # V_old  = utility.(collect(k_grid),zeros(n_k),p) ; # Start at utility with zero savings
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - n_k=$n_k - θ_k=$θ_k")
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_c, G_l = T(Model(M,V=copy(V_old)))
            # println("T(V) = $V_new")
            # println("  V  = $V_old")
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol
            println("VFI - n_k=$n_k - θ_k=$θ_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            V_ip = ScaledInterpolations(M.k_grid,V_new, BSpline(Cubic(Line(OnGrid()))))
                V_fine = V_ip.(collect(M.k_grid_fine))
            G_kp_ip = ScaledInterpolations(M.k_grid,G_kp, BSpline(Cubic(Line(OnGrid()))))
                G_kp_fine = G_kp_ip.(collect(M.k_grid_fine))
            G_c_ip = ScaledInterpolations(M.k_grid,G_c, BSpline(Cubic(Line(OnGrid()))))
                G_c_fine = G_c_ip.(collect(M.k_grid_fine))
            G_l_ip = ScaledInterpolations(M.k_grid,G_l, BSpline(Cubic(Line(OnGrid()))))
                G_l_fine = G_l_ip.(collect(M.k_grid_fine))
            # Update model
            M = Model(M; V=V_new,G_kp=G_kp,G_c=G_c,G_l=G_l,V_fine=V_fine,G_kp_fine=G_kp_fine,G_c_fine=G_c_fine,G_l_fine=G_l_fine)
            # Euler error
            Euler = Euler_Error(M.k_grid_fine,G_kp_fine,G_kp_ip.(collect(G_kp_fine)),G_l_fine,G_l_ip.(collect(G_kp_fine)),p)
            # Update model
            M = Model(M; Euler=Euler)
            # Return results
            return M::Model
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Solution not found")
end


@time M_20  = VFI_Fixed_Point(T_cts_max,Model(n_k=20))