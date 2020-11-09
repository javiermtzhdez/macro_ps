####################################################################
####################################################################
############# Advance Macroeconomics II PS 5
############# Part B
####################################################################
####################################################################




cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/")


######### Run these commented commands only the first time
using Pkg
Pkg.activate("ps5")
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
############## Utility function
######################################################################
    function u_fn(k,z,kp,l,p::Par)
        @unpack α,δ,σ,η,c_min, χ = p
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
    function euler_kp(k,z,kp,l,p::Par)
        @unpack α,δ,σ,η,c_min = p
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


    function euler_l(k,z,kp,l,p::Par)
        @unpack α,δ,σ,η, χ, Λ, Ψ,c_min = p
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
#=
# Simulation

    function Simulation_MC(Ns,MP::MP,N_0=1000)
        # Compute conditional CDF
        Γ = cumsum(MP.Π,dims=2)
        # Allocate simulated vector
        z_ind    = zeros(Int64,N_0+Ns)
        z_MC     = zeros(N_0+Ns)
        # Starting value for simulation
        z_ind[1] = Int(ceil(length(MP.grid)/2))
        z_MC[1]  = MP.grid[z_ind[1]]
        # Simulate
        for i=2:Ns+N_0
            #= Option 1
            # Draw a uniform random number (r). Compare it with conditional CDF.
            # Conditional CDF given by cumsum of current row of Π, call it Γ
            # Γ[i,j] gives the probability z'<z_j
            # We are looking for j s.t Γ[i,j-1]<r<Γ[i,j]
            # Equivalently, the lowest j s.t. Γ[i,j]-r>0
                z_ind[i] = findmax(sign.(Γ[z_ind[i-1],:] .- rand()))[2]
            =#
            #= Option 2
            # Alternatively we can draw directly from conditional distributional
            # Distributional is categorical with P[z'=z_j] = Π[i,j]
            =#
            z_ind[i] = rand(Categorical(MP.Π[z_ind[i-1],:]))
            z_MC[i]  = MP.grid[z_ind[i]]
        end
        # Throw out first N_0 elements
        z_MC = z_MC[N_0+1:end]
        # Return result
        return z_MC
    end

#

=#
######################################################################
############## Model
######################################################################

# Generate structure of model objects
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model paramters
        # Grids
        θ_k::Float64    = 1     # Default Curvature of k_grid
        n_k::Int64      = 20    # Default Size of k_grid
        n_k_fine::Int64 = 1000  # Default Size of fine grid for interpolation
        n_z::Int64      = 10    # Default size of discretized grid for productivity as a markov process
        scale_type::Int64 = 1   # Default grid type (polynomial)
        k_grid          = Make_K_Grid(n_k,θ_k,p)    # k_grid for model solution
        k_grid_fine     = Make_K_Grid(n_k_fine,1,p) # Fine grid for interpolation
        # Value and policy functions
        V         = Array{Float64}(undef,n_k,n_z)   # Value Function
        G_kp      = Array{Float64}(undef,n_k,n_z)       # Policy Function for capital k'
        G_c       = Array{Float64}(undef,n_k,n_z)       # Policy Function for consumption c
        G_l       = Array{Float64}(undef,n_k,n_z)       # Policy Function for labor l
        V_fine    = Array{Float64}(undef,n_k_fine,n_z)  # Interpolated Value Function
        G_kp_fine = Array{Float64}(undef,n_k_fine,n_z)  # Interpolated Policy Function for capial k'
        G_c_fine  = Array{Float64}(undef,n_k_fine,n_z)  # Interpolated Policy Function for consumption c
        G_l_fine  = Array{Float64}(undef,n_k_fine,n_z)  # Interpolated Policy Function for labor l
        Euler     = Array{Float64}(undef,n_k_fine,n_z)  # Errors in Euler equation
        z_grid    = Array{Float64}(undef,n_z)       # Grid for productivity via Rouwenhorst 
        Π         = Array{Float64}(undef,n_z,n_z)   # Probability transition matrix for productivity
    end
    M = Model()

# Function that finds the fixed point of the value function, then interpolates between grid points
    function VFI_fixedpoint(T::Function,M::Model)
        ### T : Bellman operator (interior loop) ###
        ### M : Model structure                  ###
    # Unpack model structure
    @unpack p,n_k,n_k_fine,θ_k,k_grid,k_grid_fine,n_z,V_fine,G_kp_fine,G_c_fine,G_l_fine,Euler,z_grid,Π = M
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for iteration
    V_old = zeros(n_k,n_z) # Value function
    V_dist = 1 # Distance between old and new value function
    iter = 1
    println(" ")
    println("------------------------")
    println("VFI - n_k=$n_k - grid curvature θ_k=$θ_k - n_z=$n_z")
    # Start VFI
    while iter <= max_iter
    # Update value function and policy functions
    a=1.0
    V_new = a*T(Model(M,V=copy(V_old)))[1] + (1-a)*V_old # Call Bellman operator which returns a new value function at each capital grid point
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
        V_new, G_kp, G_c, G_l = T(Model(M,V=copy(V_new)))
        println("VFI - n_k=$n_k - θ_k=$θ_k")
        println("Converged after $iter iterations with a distance of ",100*V_dist," %")
        println("------------------------")
        println(" ")
        # Interpolate to fine grid on capital using natural cubic spline if it converged
        for i in 1:n_z
            V_ip = ScaledInterpolations(M.k_grid,V_new[:,i], FritschButlandMonotonicInterpolation()) # Monotonic spline because I know that the value function is always increasing in capital
                V_fine[:,i] = V_ip.(collect(M.k_grid_fine))
            G_kp_ip = ScaledInterpolations(M.k_grid,G_kp[:,i], BSpline(Cubic(Line(OnGrid()))))
                G_kp_fine[:,i] = G_kp_ip.(collect(M.k_grid_fine))
            G_c_ip = ScaledInterpolations(M.k_grid,G_c[:,i], BSpline(Cubic(Line(OnGrid()))))
                G_c_fine[:,i] = G_c_ip.(collect(M.k_grid_fine))
            G_l_ip = ScaledInterpolations(M.k_grid,G_l[:,i], BSpline(Cubic(Line(OnGrid()))))
                G_l_fine[:,i] = G_l_ip.(collect(M.k_grid_fine))
            # Percent Euler Error on fine grid
            #Euler[:,i] = Euler_Error(M.k_grid_fine,G_kp_fine,G_kp_ip.(collect(G_kp_fine)),G_l_fine,G_l_ip.(collect(G_kp_fine)),p)
        end
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
#


######################################################################
############## Bellman operator
######################################################################


    function T_cts_max(M::Model)
        @unpack p,n_k,k_grid,n_z,V,G_kp,G_c,G_l,z_grid,Π,z_grid = M
        @unpack β,α,δ,η,σ,c_min = p
        #=
        It is possible to write two fns for boundaries for kp. One for the inner loop, that
        is, when given z, k and kp, we choose l, and one for for the outer
        loop which is when we maximize over kp, but I am going to use a general one =#
        function kp_max_fn(k,l,z)
            #unpack β,α,z,δ,η,σ,c_min = p
            return z*(k^α)*(l^(1-α)) + (1-δ)*k - c_min
        end
        # Define boundaries for k'
        kp_min = 1.001*k_grid[1]
        # Boundaries for labor
        # Labor min and max
        l_min = 1E-16
        l_max = 1.0
        # Interpolation for Value Function
        Vp = [x->ScaledInterpolations(k_grid,x[:,i], BSpline(Cubic(Line(OnGrid())))) for i in 1:n_z]
        #
        #= Now comes the Bellman function (this is a function that returns 
        the objective function (Bellman function) for given z,k,kp,l =#
        function Obj_fn(k,z,kp,l,Π_z::Vector,p::Par)
            # Π_z: Vector of conditional probabilites for productivity next period given z today
            Emax = sum(Π_z[x]*Vp[x](V).(kp) for x in 1:n_z) # sum of z_j probability when starting at z_i : ∑_j π_ij V(kp,z_j)
            return -u_fn(k,z,kp,l,p) - β*Emax
        end
        # Function to get derivative of Emax wrt capital k'
        dVp(x,Π_z::Vector) = sum(Π_z[i]*ForwardDiff.derivative(Vp[i](V),x) for i in 1:n_z)
        # Function that returns derivative of objective function wrt k'
        d_Obj_fn_kp(k,z,kp,l,Π_z::Vector,p::Par) = euler_kp(k,z,kp,l,p) + β*dVp(kp,Π_z)
        # Derivative of objective function wrt labor l
        d_Obj_fn_l(k,z,kp,l,p::Par) = euler_l(k,z,kp,l,p)
        # Define function that finds optimal labor l given (k,z,k') and returns objective function conditional on optimal labor
        function inner_function_labor(k,z,kp,Π_z::Vector,p::Par)
            # Check for corner solutions on labor
            dobj_min = euler_l(k,z,kp,l_min,p)
            dobj_max = euler_l(k,z,kp,l_max,p)
            if dobj_min <= 0
                return Obj_fn(k,z,kp,l_min,Π_z,p),l_min
            elseif dobj_max >= 0
                return Obj_fn(k,z,kp,l_max,Π_z,p),l_max
            else
            # if no corner solutions, find interior solution
                min_result = optimize(x->euler_l(k,z,kp,x,p).^2,l_min,l_max,Brent())
                l = min_result.minimizer
                return Obj_fn(k,z,kp,l,Π_z,p),l
            end
        end
        # Outer loop for all possible values of productivity today
        for j in 1:n_z
            # Inner loop for each capital level in the grid
            for i in 1:n_k
                kp_max = min(kp_max_fn(k_grid[i],1.0,z_grid[j]),k_grid[end])
                # Check for corner solutions on capital
                l_kp_min = inner_function_labor(k_grid[i],z_grid[j],kp_min,Π[j,:],p)[2]
                l_kp_max = inner_function_labor(k_grid[i],z_grid[j],kp_max,Π[j,:],p)[2]
                dobj_min = d_Obj_fn_kp(k_grid[i],z_grid[j],kp_min,l_kp_min,Π[j,:],p)
                dobj_max = d_Obj_fn_kp(k_grid[i],z_grid[j],kp_max,l_kp_max,Π[j,:],p)
                if dobj_min <= 0.0
                    G_kp[i,j] = kp_min
                    G_l[i,j] = l_kp_min
                    V[i,j] = -Obj_fn(k_grid[i],z_grid[j],kp_min,l_kp_min,Π[j,:],p)
                elseif dobj_max >= 0.0
                    G_kp[i,j] = kp_max
                    G_l[i,j] = l_kp_max
                    V[i,j] = -Obj_fn(k_grid[i],z_grid[j],kp_max,l_kp_max,Π[j,:],p)
                else
                # If no corner solution, find interior solution
                    min_result = optimize(x->inner_function_labor(k_grid[i],z_grid[j],x,Π[j,:],p)[1],kp_min,kp_max,Brent())
                    # Check result
                    #converged(min_result) || error("Failed to solve Bellman max for capital =" k_grid[i]" in $(iterations(min_result)) iterations")
                    # Record results
                    V[i,j] = -min_result.minimum
                    G_kp[i,j] = min_result.minimizer
                    G_l[i,j] = inner_function_labor(k_grid[i],z_grid[j],G_kp[i],Π[j,:],p)[2]
                end
            end
            # Fill in policy for consumption
            G_c[:,j] = z_grid[j].*(collect(k_grid).^α).*(G_l[:,j].^(1-α)) .- G_kp[:,j]
        end
        # Return results
        return V, G_kp, G_c, G_l
    end

#

                ### Solve the problem for θ_k=2,n_z=10,n_k=20 ###
# Get discrete grid for productivity and transition matrix
(log_z,Π) = Rouwenhorst95(10,p)[1:2]
z = exp.(log_z)
# Get solution
@time Mc  = VFI_fixedpoint(T_cts_max,Model(n_k=20,θ_k=2,n_z=10,z_grid=z,Π=Π))


# Interpolate the value function along the z dimension for 3d plot
z_grid = range(z[1],z[end],length=size(z)[1])
z_grid_fine = range(z[1],z[end],length=Mc.n_k_fine)
V_fine_3d = zeros(Mc.n_k_fine,Mc.n_k_fine)
V_fine_33d = [ScaledInterpolations(z_grid,Mc.V_fine[i,:], BSpline(Cubic(Line(OnGrid())))).(collect(z_grid_fine)) for i in 1:Mc.n_k_fine]
for i in 1:Mc.n_k_fine
    V_fine_3d[i,:] = V_fine_33d[i]
end
# Surface and contour plot of the value function
gr()
x= Mc.k_grid_fine; y=z_grid_fine; f=V_fine_3d;
plot(x,y,f, st=:surface,xlabel="Capital",ylabel="Productivity") # Surface plot
savefig("./Figures/surface_vf.pdf")
plot(x,y,f, st=:contour,xlabel="Capital",ylabel="Productivity",nlevels=100, width=2, size=[800,480]) # Contour plot
savefig("./Figures/contour_vf.pdf")