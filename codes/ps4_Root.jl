####################################################################
####################################################################
############# Advance Macroeconomics II PS 4
####################################################################
####################################################################


# Adda and Cooper textbook is a good reference for the statement of the problem and how to approach it I think.

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

################################################################################
################################################################################
#                                 Parameters
################################################################################
################################################################################

    @with_kw struct Par  #  The macro @with_kw which decorates a type definition to allow default values and a keyword constructor, in this case Par
        # This are the structural parameters which wont move throughout of the model
        #z::Float64 = 1 ; # Productivity
        z::Float64 = 1.0 ; # Productivity for the increase for plotting
        α::Float64 = 1/3; # Capital share, to write the alpha you do it as in latex \alp ha
        β::Float64 = 0.98; # Discount factor
        η::Float64 = 1.0; # Parameter for labor
        σ::Float64 = 2.0; # Parameter for the risk aversion
        δ::Float64 = 0.05; # Depreciation rate
        χ::Float64  = 8.217267384235548; # Parameter for labor disutility
        Ψ::Float64  = 10.301099814055622; # Value needed for calculations of χ
        Λ::Float64 = 1.6607895618579493; # Value needed for calculations of χ
        a::Float64 = 5.0; # Value for penalty function
        # VFI values
        max_iter::Int64 = 2000 ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9 ; # Tolerance for the distance (tol for the fixed point)
        # Howard's Policy iterations
        H_tol::Float64 = 1E-9; # Tolerance for policy function iteration.
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Minimum consumption 
        c_min::Float64    = 1E-16
    end

p=Par()

################################################################################




################################################################################
################################################################################
#                                 Functions
################################################################################
################################################################################
#1. utility_unr cambió a u_fn
#2. d_utility_kp cambió a euler_kp
#3. d_utility_l cambió a euler_l
#4. get_kp_max(k,l) cambió a kp_max_fn(k,l)
#5. Obj_fn_condl cambió a inner_function_labor


######################################################################
############## Steady State
######################################################################
# Steady state values
    function SS_values(p::Par)
        # This function takes in parameters and provides steady state values
        # Parameters: productivity (z), returns to scale (a) and discount factor (b)
        # Output: values for capital, production, consumption, rental rate, wage
        @unpack z, α, β, σ, η, δ, χ, Λ, Ψ = p
        l_ss = ((z*(1-α)*Ψ^α)/(χ*Λ^σ))^(1/(σ+η))
        k_ss = Ψ*l_ss
        c_ss = Λ*l_ss
        y_ss = z*l_ss*Ψ^α
        r_ss = α*z*Ψ^(α-1)
        w_ss = (1-α)*z*Ψ^α
        return l_ss,k_ss,y_ss,c_ss,r_ss,w_ss
    end

#


######################################################################
############## Euler Error
######################################################################

    function Euler_Error(k,kp,kpp,l,lp,p::Par)
        # Return percentage error in Euler equation
        @unpack z, α, β, σ, δ = p
        LHS = (z.*(k.^α).*(l.^(1-α)).+(1-δ).*k.-kp).^(-σ)
        RHS = β.*(α.*z.*((lp./kp).^(1-α)).+(1-δ)).*((z.*(kp.^α).*(lp.^(1-α)).+(1-δ).*kp.-kpp).^(-σ))
        return (RHS./LHS.-1).*100
    end

#

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

######################################################################
############## Utility function
######################################################################
    function u_fn(k,kp,l,p::Par)
        @unpack α,δ,σ,z,η,c_min, χ = p
        # Utility from consumption
        u_c = 0; 
        # Utility of consumption
        c = z*(k^α)*(l^(1-α))+(1-δ)*k-kp # Consumption from resource constraint
        u_c = (c^(1-σ))/(1-σ)
        # Disutility of labor
        u_l = χ*((l^(1+η))/(1+η))
        return u_c-u_l
    end
#

######################################################################
############## Euler Capital
######################################################################
    function euler_kp(k,kp,l,p::Par)
        @unpack α,δ,σ,z,η,c_min = p
        c = z*(k^α)*(l^(1-α))+(1-δ)*k-kp
        d_u = 0
        if c>c_min
            d_u = c^(-σ)
        else
            d_u = c_min^(-σ)
        end
        return -d_u
    end

#

######################################################################
############## Euler Labor 
######################################################################


    function euler_l(k,kp,l,p::Par)
        @unpack α,δ,σ,z,η, χ, Λ, Ψ,c_min = p
        c = z*(k^α)*l^(1-α)+(1-δ)*k - kp
        d_c = 0
        #return -χ*l^η + (1/c^σ)*(1-α)*z*k^α*l^(-α)
        # Consumption above threshold c_min?
        if c>c_min
            d_c = c^(-σ)
        else
            d_c = c_min^(-σ)
        end 
        return -χ*l^η + d_c*(1-α)*z*(k^α)*(l^(-α))
    end

#

# If Euler labor is written appropriately, by plugging the SS values, I would have it equal to zero
# euler_l(k_ss,k_ss,l_ss,p)
# I get -4.440892098500626e-16
# Which is basically zero...





# Generate structure of model objects
@with_kw struct Model
    # Parameters
    p::Par = Par() # Model paramters
    # Grids
    θ_k::Float64    = 1     # Default Curvature of k_grid
    n_k::Int64      = 20    # Default Size of k_grid
    n_k_fine::Int64 = 1000  # Default Size of fine grid for interpolation
    scale_type::Int64 = 1   # Default grid type (polynomial)
    k_grid          = Make_K_Grid(n_k,θ_k,p)    # k_grid for model solution
    k_grid_fine     = Make_K_Grid(n_k_fine,1,p) # Fine grid for interpolation
    # Value and policy functions
    V         = Array{Float64}(undef,n_k)       # Value Function
    G_kp      = Array{Float64}(undef,n_k)       # Policy Function for capital k'
    G_c       = Array{Float64}(undef,n_k)       # Policy Function for consumption c
    G_l       = Array{Float64}(undef,n_k)       # Policy Function for labor l
    V_fine    = Array{Float64}(undef,n_k_fine)  # Interpolated Value Function
    G_kp_fine = Array{Float64}(undef,n_k_fine)  # Interpolated Policy Function for capial k'
    G_c_fine  = Array{Float64}(undef,n_k_fine)  # Interpolated Policy Function for consumption c
    G_l_fine  = Array{Float64}(undef,n_k_fine)  # Interpolated Policy Function for labor l
    Euler     = Array{Float64}(undef,n_k_fine)  # Errors in Euler equation
end
M = Model()










######################################################################
############## kpmax function
######################################################################
#=
# kp max
function kp_max_fn(k,l,p::Par)
    @unpack β,α,z,δ,η,σ,c_min = p
    return z*(k^α)*(l^(1-α)) + (1-δ)*k - c_min
end
=#

######################################################################
############## Bellman operator
######################################################################

function T_cts_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c, G_l = M
    @unpack β,α,z,δ,η,σ,c_min, χ = p
    # Interpolate current iteration of value function v(k') so I can evaluate it at any k'
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid())))) # Impose monotonicity because I know it is increasing in capital
    # Derivative of value function (of iteration n) wrt capital k'
    dVp(x) = ForwardDiff.derivative(Vp,x)
    # Objective Function
    Obj_fn(k,kp,l,p::Par) = -u_fn(k,kp,l,p) - β*Vp.(kp) # Returns the objective function for a triplet of k,kp and l
    # Euler function for capital k' 
        # [This part combines the contemporary part with the numerical derivative
        # of the Value functionn wrt capital since employing the envelope
        # theorem requires too much]
    d_Obj_fn_kp(k,kp,l,p::Par) = euler_kp(k,kp,l,p) + β*dVp.(kp)
    # Euler function for labor
    d_Obj_fn_l(k,kp,l,p::Par)  = euler_l(k,kp,l,p)
    # Labor min and max
    l_min = 1E-16
    l_max = 1.0
    function kp_max_fn(k,l)
        #unpack β,α,z,δ,η,σ,c_min = p
        return z*(k^α)*(l^(1-α)) + (1-δ)*k - c_min
    end
    function inner_function_labor(k,kp,p::Par)
        Euler_l_min = euler_l(k,kp,l_min,p)
        Euler_l_max = euler_l(k,kp,l_max,p)
        if Euler_l_min <=0 # That is: the agent would like to work less but he cant, we have boundaries, so he is stuck at l_min
            #println("Corner solution - lower bound for labor")
            return -u_fn(k,kp,l_min,p) - β*Vp.(kp), l_min
        elseif Euler_l_max>= 0 # That is: the agent would like to work more but he cant, we have boundaries, so he is stuck at l_max
            return -u_fn(k,kp,l_max,p) - β*Vp.(kp), l_max
            #println("Corner solution - upper bound for labor")
        else # No corner solution, thus return the Bellman Objective  for univariate maximization
            #println("No corner solution - return the Bellman Objective  for univariate maximization")
            min_result = optimize(x->euler_l(k,kp,x,p).^2,l_min,l_max,GoldenSection())
            l_opt = min_result.minimizer
            # Check that consumption is positive with this
            #c_inner = consumption(k,kp,l_opt,p)
            return -u_fn(k,kp,l_opt,p) - β*Vp.(kp), l_opt
            #=
            if c_inner ≥ c_min
                println("Consumption is at least c_min")
                return -ufn(k,kp,l_opt,p) - β*Vp.(kp), l_opt
                println("labor optimal is     ", l_opt)
                #return l_opt
            else
                println(" Something went wrong in the inner function. The consumption at the optimal labor is less than c_min")
            end
            =#
        end
    end
    # Here is where the process for kp begins
    kp_min = 1.001*k_grid[1]
    # Maximization for kp 
    for i in 1:n_k
        kp_max = min(kp_max_fn(k_grid[i],1.0),k_grid[end])
        # Corner solutions for labor and capital with the current boundaries (boundaries for this grid point)
        l_kp_min = inner_function_labor(k_grid[i],kp_min,p)[2]
        l_kp_max = inner_function_labor(k_grid[i],kp_max,p)[2]
        dobj_min = d_Obj_fn_kp(k_grid[i],kp_min,l_kp_min,p)
        dobj_max = d_Obj_fn_kp(k_grid[i],kp_max,l_kp_max,p)
        if dobj_min <= 0.0  # Corner solution for minimum value of capital
            G_kp[i] = kp_min
            G_l[i] = l_kp_min
            V[i] = -Obj_fn(k_grid[i],kp_min,l_kp_min,p)
        elseif dobj_max >= 0.0 # Corner solution for maximum value of capital
            G_kp[i] = kp_max
            G_l[i] = l_kp_max
            V[i] = -Obj_fn(k_grid[i],kp_max,l_kp_max,p)
        else # No corner solution for capital, then find interior solution for kp with maximization
            min_result = optimize(x->inner_function_labor(k_grid[i],x,p)[1],kp_min,kp_max,Brent())
            # Record results
            V[i] = -min_result.minimum
            G_kp[i] = min_result.minimizer
            G_l[i] = inner_function_labor(k_grid[i],G_kp[i],p)[2]
        end
    end
    # Policy for consumption
    G_c = z.*(collect(k_grid).^α).*(G_l.^(1-α)) .- G_kp
    # Return results
    return V, G_kp, G_c, G_l
end





# Graphs
function VFI_Graphs(M::Model,VFI_Type)
gr()
# Value Function
    plot(M.k_grid,M.V,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
    plot!(M.k_grid_fine,M.V_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
    xlabel!("Capital")
    ylabel!("Value")
    savefig("./Figures/VFI_"*VFI_Type*"_V_$(M.n_k)_$(M.θ_k).pdf")
# Capital Policy Function
    plot(M.k_grid_fine,M.k_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(M.k_grid,M.G_kp,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
    plot!(M.k_grid_fine,M.G_kp_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
    xlabel!("Capital")
    ylabel!("Capital")
    savefig("./Figures/VFI_"*VFI_Type*"_G_kp_$(M.n_k)_$(M.θ_k).pdf")
# Labor Policy Function
    plot(M.k_grid_fine,M.k_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(M.k_grid,M.G_l,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
    plot!(M.k_grid_fine,M.G_l_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
    xlabel!("Capital")
    ylabel!("labor")
    savefig("./Figures/VFI_"*VFI_Type*"_G_l_$(M.n_k)_$(M.θ_k).pdf")
# Euler Percentage Error
    plot(M.k_grid_fine,zeros(M.n_k_fine),lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(M.k_grid_fine,M.Euler,linetype=:scatter,marker=(:diamond,2),linecolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
    xlabel!("Capital")
    ylabel!("Percentage Points")
    savefig("./Figures/VFI_"*VFI_Type*"_Euler_$(M.n_k)_$(M.θ_k).pdf")
    println("\n     Graphs Completed for VFI_$VFI_Type - n_k=$(M.n_k) - θ_k=$(M.θ_k)\n")
end

# Function that finds the fixed point of the value function, then interpolates between grid points
function VFI_fixedpoint(T::Function,M::Model)
        ### T : Bellman operator (interior loop) ###
        ### M : Model structure                  ###
# Unpack model structure
@unpack p, n_k, n_k_fine, θ_k, k_grid, k_grid_fine = M
# VFI paramters
@unpack max_iter, dist_tol = p
# Initialize variables for iteration
V_old = zeros(n_k) # Value function
V_dist = 1 # Distance between old and new value function
iter = 1
println(" ")
println("------------------------")
println("VFI - n_k=$n_k - grid curvature θ_k=$θ_k")
# Start VFI
while iter <= max_iter
    # Update value function and policy functions
    V_new, G_kp, G_c, G_l = T(Model(M,V=copy(V_old))) # Call Bellman operator which returns a new value function at each capital grid point
    # Update value function and distance between previous and current iteration
    V_dist = maximum(abs.(V_new./V_old.-1))
    V_old = V_new
    # Report progress
    println("   VFI Loop: iter=$iter, dist=",100*V_dist," %")
    # Report progress every 100 iterations
    #if mod(iter,100)==0
    #    println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
    #end
    # Check if converged
    if V_dist <= dist_tol
        println("VFI - n_k=$n_k - θ_k=$θ_k")
        println("Converged after $iter iterations with a distance of ",100*V_dist," %")
        println("------------------------")
        println(" ")
        # Interpolate to fine grid using natural cubic spline if it converged
        V_ip = ScaledInterpolations(M.k_grid,V_new, FritschButlandMonotonicInterpolation()) # Monotonic spline because I know that the value function is always increasing in capital
            V_fine = V_ip.(collect(M.k_grid_fine))
        G_kp_ip = ScaledInterpolations(M.k_grid,G_kp, BSpline(Cubic(Line(OnGrid()))))
            G_kp_fine = G_kp_ip.(collect(M.k_grid_fine))
        G_c_ip = ScaledInterpolations(M.k_grid,G_c, BSpline(Cubic(Line(OnGrid()))))
            G_c_fine = G_c_ip.(collect(M.k_grid_fine))
        G_l_ip = ScaledInterpolations(M.k_grid,G_l, BSpline(Cubic(Line(OnGrid()))))
            G_l_fine = G_l_ip.(collect(M.k_grid_fine))
        # Percent Euler Error on fine grid
        Euler = Euler_Error(M.k_grid_fine,G_kp_fine,G_kp_ip.(collect(G_kp_fine)),G_l_fine,G_l_ip.(collect(G_kp_fine)),p)
        #Euler = map(x->Euler_Error(x[1],x[2],x[3],x[4],x[5],p),[M.k_grid_fine,G_kp_fine,G_kp_ip.(collect(G_kp_fine)),G_l_fine,G_l_ip.(collect(G_kp_fine))])
        # Update model with solution
        M = Model(M; V=V_new,G_kp=G_kp,G_c=G_c,G_l=G_l,V_fine=V_fine,G_kp_fine=G_kp_fine,G_c_fine=G_c_fine,G_l_fine=G_l_fine,Euler=Euler)
        return M
    end
    # If it didn't converge, go to next iteration
    iter += 1
end
# If it didn't converge, return error
error("Error in VFI - Solution not found")
end






            ### Execute VFI and plot graphs ###
# Curved spaced grid with univariate optimization
@time M_20c  = VFI_fixedpoint(T_cts_max,Model(n_k=20,θ_k=2.6))
     VFI_Graphs(M_20c,"cts_max")


