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


# Derivative of utility function wrt labor
function d_utility_l(k,x,p::Par)
    @unpack α,δ,σ,z,η,c_min = p
    kp = x[1]
    l = x[2]
    c = z*(k^α)*(l^(1-α))+(1-δ)*k-kp
    d_u = 0
    if c>c_min
        d_u = c^(-σ)
    else
        d_u = c_min^(-σ)
    end
    return d_u*z*(k^α)*(1-α)*(l^(-α))-χ*(l^η)
end

function Euler_Error(k,kp,kpp,l,lp,p::Par)
    # Return percentage error in Euler equation
    @unpack z, α, β, σ, δ = p
    LHS = (z.*(k.^α).*(l.^(1-α)).+(1-δ).*k.-kp).^(-σ)
    RHS = β.*(α.*z.*((lp./kp).^(1-α)).+(1-δ)).*((z.*(kp.^α).*(lp.^(1-α)).+(1-δ).*kp.-kpp).^(-σ))
    return (RHS./LHS.-1).*100
end

# Grid
function Make_K_Grid(n_k,θ_k,p::Par,scale_type)
    # Get SS
    l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    # Get k_grid
    if θ_k≠1
        if scale_type=="Poly"
        k_grid = PolyRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
        elseif scale_type=="Exp"
        k_grid = ExpRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
        else
        error("scale_type must be either Poly or Exp")
        end
    else
    k_grid = range(1E-5,2*k_ss,length=n_k)
    end
    # Return
    return k_grid
end

function Make_L_Grid(n_k,θ_l,p::Par)
    # Get SS
    l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    # Get L_grid
    k_grid = range(0,1,length=n_k)
    # Return
    return k_grid
end


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

# Generate structure of model objects
@with_kw struct Model
    # Parameters
    p::Par = Par() # Model paramters in their own structure
    # Grids
    θ_k::Float64    = 1     # Curvature of k_grid
    θ_l::Float64    = 1     # Curvature of l_grid
    n_k::Int64      = 20    # Size of k_grid
    n_k_fine::Int64 = 1000  # Size of fine grid for interpolation
    n_l_fine::Int64 = 1000  # Size of fine grid for interpolation
    k_grid          = Make_K_Grid(n_k,θ_k,p, "Poly")    # k_grid for model solution
    l_grid          = Make_L_Grid(n_k,θ_l,p)    # l_grid for model solution
    k_grid_fine     = Make_K_Grid(n_k_fine,1,p,"Poly") # Fine grid for interpolation
    l_grid_fine     = Make_L_Grid(n_k_fine,1,p) # Fine grid for interpolation
    # Value and policy functions
    V         = Array{Float64}(undef,n_k)       # Value Function
    G_kp      = Array{Float64}(undef,n_k)       # Policy Function for capital
    G_c       = Array{Float64}(undef,n_k)       # Policy Function for consumption
    G_l       = Array{Float64}(undef,n_k)       # Policy Function for labor  
    V_fine    = Array{Float64}(undef,n_k_fine)  # Value Function on fine grid
    G_kp_fine = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
    G_c_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
    G_l_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
    Euler     = Array{Float64}(undef,n_k_fine)  # Errors in Euler equation
end

M = Model()



#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - Continuous Choice
function T_cts_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c, G_l = M
    @unpack z, α, β, c_min, δ = p
    # Define the interpolation of the current value function (V) for values of Vp
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid()))))
    #dVp(x) = ForwardDiff.derivative(Vp,x[1])
    for i = 1:n_k
        Obj_Fun(k,x,p) = -utility(k_grid[i],x,p) - β*Vp.(x[1])
        # Min and max kp given current k
        #################################################################################
        # Not so sure about these bounds
        #################################################################################
        kp_min = 0.0001 # Suggested by Emmanuel
        # kp_min = k_grid[1]/(1-δ) # The bound I thought of
        kp_max = min(z*k_grid[i]^α  + (1-δ)*k_grid[i] - c_min, k_grid[end])
        #kp_max = 0.0
        # Bounds for l
        l_min = 1E-16
        l_max = 1.0
        # Now check if they bind
        # Lower bound 
        lb = [kp_min,l_min]
        #dObj_min = ForwardDiff.gradient(Obj_Fun,lb)
        dObj_min = -dutility(k_grid[i],lb,p) .- β*dVp(lb[1])
        if dObj_min[1]>0 && dObj_min[2]>0
            V[i]    = -Obj_Fun(lb)
            G_kp[i] = lb[1]
            G_l[i] = lb[2]
        else
        # Upper bound
        ub = [kp_max, l_max]
            dObj_max = -dutility(k_grid[i],ub,p) .- β*dVp(ub[1])
            #dObj_max = ForwardDiff.gradient(Obj_Fun,ub)
        #if dObj_max[1]<0 && dObj_max[2]<0
        #    V[i]    = -bj_Fun(ub)
        #    G_kp[i] = ub[1]
        #    G_l[i] = ub[2]
        #else
        # Now the minimization process
        inner_optimizer = NelderMead()
        l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
        initial_x = [(kp_min+kp_max)/2,l_ss]
        #od = OnceDifferentiable(Obj_Fun, initial_x; autodiff=:forward) # No Funciona bien
        #min_result= optimize(od,lb,ub,initial_x, Fminbox()) # No funciona bien
        #min_result = optimize(x->Obj_fn(k_grid[i],x,p),[kp_min,l_min],[kp_max,l_max],initial_x,Fminbox(inner_optimizer))
        min_result = optimize(x->Obj_fn(k_grid[i],x,p),lb,ub,initial_x,Fminbox(inner_optimizer))
        println(" Cheking it solved the minimization")
        # Check result
        converged(min_result) || error("Failed to solve Bellman max in $(iterations(min_result)) iterations")
        # Record results
        V[i]     = -min_result.minimum
        G_kp[i]  = min_result.minimizer[1]
        G_l[i]  = min_result.minimizer[2]
        #end
        #end
    end
    G_c = z.*(collect(k_grid).^α).*(G_l.^(1-α)) .+(1-δ).*collect(k_grid) .-G_kp
    # Return Results
         println("T(V) = $V")
    return V, G_kp, G_c, G_l

end
#-----------------------------------------------------------


#-----------------------------------------------------------
# VFI Fixed Point
function VFI_Fixed_Point(T::Function,M::Model)
    # Unpack model structure
    @unpack p, n_k, θ_k, k_grid = M
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    V_old  = zeros(n_k)     ; # Initialize value function, here I just do 0, not the best option
    # V_old  = utility.(collect(k_grid),zeros(n_k),p) ; # Start at utility with zero savings
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - n_k=$n_k - θ_k=$θ_k")
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_c, G_l = T(Model(M,V=copy(V_old)))
        # V_new, G_kp, G_c, G_l = T_cts_max(Model(M,V=copy(V_old)))
            # println("T(V) = $V_new")
            # println("  V  = $V_old")
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol #
            #######################################################################################################################################################################
            # OK NOTE HERE THAT ONCE YOU GOT CLOSE ENOUGH USING THE ROUGH GRID, THEN YOU BEGIN THE FINE GRID SEARCH AND INTERPOLATE THE POLICY FUNCTIONS AND THE VALUE FUNCTIONS
            # NOT BEFORE, BECAUSE YOU WOULD WASTE TIME
            #######################################################################################################################################################################
            println("VFI - n_k=$n_k - θ_k=$θ_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            V_ip = ScaledInterpolations(M.k_grid,V_new, BSpline(Cubic(Line(OnGrid()))))
                V_fine = V_ip.(collect(M.k_grid_fine))
            G_kp_ip = ScaledInterpolations(M.k_grid,G_kp, BSpline(Cubic(Line(OnGrid()))))
                G_kp_fine = G_kp_ip.(collect(M.k_grid_fine))
            G_c_ip = ScaledInterpolations(M.k_grid,G_c, BSpline(Cubic(Line(OnGrid()))))
                G_c_fine = G_c_ip.(collect(M.k_grid_fine))
            G_l_ip = ScaledInterpolations(M.k_grid,G_l, BSpline(Cubic(Line(OnGrid()))))
                G_l_fine = G_l_ip.(collect(M.k_grid_fine))
            # Update model
            M = Model(M; V=V_new,G_kp=G_kp,G_c=G_c,V_fine=V_fine,G_kp_fine=G_kp_fine,G_c_fine=G_c_fine, G_l=G_l, G_l_fine=G_l_fine)
            # Return results
            return M
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Solution not found")
end





# Minimum Bracketing Function
    # This function comes from Numerical Recipes, Section 10.1
    # Input: numbers a,b that initiliaze bracket and objective function F
        # Optional: Bounds for the brackeet so that a,b,c ∈ [x_min,x_max]
    # Output: numbers (a,b,c) and images (fa,fb,fc) such that a<b<c and fb<fa,fc
    function mnbrak(a,b,F::Function,x_min=-Inf,x_max=Inf)
        # Auxiliariy variables
        Tiny    = 1E-20
        G_limit = 100
        # Define golden ratio
        γ = 1.618034
        # Evaluate function at end points
        fa, fb = F(a), F(b)
        # Swap a and b so that we can go downhill in the direction from a to b.
            # This way we know for certain that fb<fa, we only need to find c such that fb<fc
        if fb>fa
        a , b = b , a
        fa,fb = fb, fa
        end
        # Guess for value c expanding bracket away from b -> (a,b,c) or (c,b,a)
            # Check bounds
        c  = max( min( b + γ*(b-a) , x_max ) , x_min)
        fc = F(c)
        # While fb>fc we need to keep on bracketing
        while fb>fc
            # Compute u (new candidate) by parabolic extrapolation from a, b, c.
                # TINY is used to prevent any possible division by zero.
            r = (b-a)*(fb-fc)
            q = (b-c)*(fb-fa)
            u = min(max(b-((b-c)*q-(b-a)*r)/(2*sign(q-r)*max(abs(q-r),Tiny)),x_min),x_max)
            u_lim = min(max(b+G_limit*(c-b),x_min),x_max)
            # Test various cases for new candidate
            if ((b-u)*(u-c) > 0) # Parabolic u is between b and c
                fu=F(u)
                if (fu < fc) # Got a minimum between b and c so bounds are (b,u,c)
                    a, fa, b, fb = b, fb, u, fu
                    break
                elseif (fu > fb) # Got a minimum between a and u so bounds are (a,b,u)
                    c, fc = u, fu
                    break
                end
                # If you got here is because candidate u failed
                # Get new candidate by expanding interval with golden ratio
                # Check bounds
                u  = max(min( c+γ*(c-b) , x_max ),x_min)
                fu = F(u)
            elseif ((c-u)*(u-u_lim) > 0.0) # Parabolic u is between c and its limit (ulim)
                fu=F(u)
                if (fu < fc) # Drop c and replace it with u, get new u with golden expansion
                    b, c, fb, fc = c, u, fc, fu
                    u  = max(min( c+γ*(c-b) , x_max ),x_min)
                    fu = F(u)
                end
            elseif ((u-u_lim)*(u_lim-c) >= 0.0) # U is above its limit, reign it in!
                u  = u_lim
                fu = F(u)
            else # Nothing worked, use golden expansion
                u  = max(min( c+γ*(c-b) , x_max ),x_min)
                fu = F(u)
            end
            # If no break then forget the oldest point and go onto next iteration
                # This means assigning b->a, c->b, u-> and forgetting about a
            a, b, c, fa, fb, fc = b, c, u, fb, fc, fu
        end
        # Return solution once out of the loop
            # Swap a and c so that a<b<c
            if a>c
            a , c  = c, a
            fa, fc = fc, fa
            end
        # Minimum bracketed in [a,c] with intermediate point b so that fb<fa,fc
        # println("Minimum bracketed in [$a,$c] with intermediate point b=$b \n Function values: F(a)=$fa, F(b)=$fb, F(c)=$fc")
        return a,c,b,fa,fc,fb
    end
#

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

################################################################################
################################################################################
################################################################################




















# First let's calculate the Steady State values
SS_values(p)

################################################################################
##############          STEADY STATE VARIABLES
################################################################################


# Then the value of χ that satisfies l_ss = 0.4 is
#global χ = f1(Ψ, Λ, p::Par)
println("Chi = ")
println(χ)

# Since we know what the function for the labor in steady state:
lss(Ψ,Λ,χ, p)

# Then the value of χ that guarantees what Sergio asked for is 28.169401999194655


################################################################################
##############          1.a) Direct maximization over a pair of variables: (c,l) or (k',l).
################################################################################


# The value function is 
#\begin{align}
# V(k,l) = \max_{c,l} \left[  \frac{c^{1-\sigma}}{1-\sigma} - \chi \frac{l^{1+\eta}}{1+\eta} + \beta V(k',l')   \right]
#\end{align}
# subject to 
# 0 \leq l \leq 1 
# 0 \leq k' \leq zk^{\alpha}l^{1-\alpha} + (1-\delta)k - k'

# We can rewrite the problem as:
# V(k,l) = \max_{k',l} \left[  \frac{(zk^{\alpha}l^{1-\alpha} + (1-\delta)k - k')^{1-\sigma}}{1-\sigma} - \chi \frac{l^{1+\eta}}{1+\eta} + \beta V(k',l')   \right]
#\end{align}
# subject to 
# 0 \leq l \leq 1 




@time M_20  = VFI_Fixed_Point(T_cts_max,Model(n_k=20,θ_k=3.6))
VFI_Graphs(M_20,"mvariate_max")
























#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - Continuous Choice
function T_cts_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c, G_l = M
    @unpack z, α, β, c_min, δ, η, σ  = p
    # Define the interpolation of the current value function (V) for values of Vp
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid()))))
    dVp(x) = ForwardDiff.derivative(Vp,x[1])
    for i = 1:n_k
        l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values
        # Objective function
        Obj_Fun = x->(-utility(k_grid[i],x,p) - β*Vp.(x[1]))
        # Lower and upper bound for the values
        l_min = 1E-16
        l_max = 1.0
        # We can check the bounds better
        if k_grid[i]<k_ss
            kp_min = k_grid[1]
            kp_max = k_ss
        elseif k_grid[i]>k_ss
            kp_min = k_ss
            kp_max = min(z*k_grid[i]^α  + (1-δ)*k_grid[i] - c_min, k_grid[end])
        #kp_min = k_grid[1]
        #kp_max = min(z*k_grid[i]^α  + (1-δ)*k_grid[i] - c_min, k_grid[end])
        end
        #=
        Now let's be smart as Sergio suggested and check the corner solutions. We know that optimally, we have an Euler equation for labor
        and a Euler equation for capital using the envelope condition
        =#
        # First for labor at the min, this would mean that the agent is working l_min, and still would like to work less
        foc_l_min_bound = Euler_Labor(k_grid[i], [k_grid[1],l_min],p)
        foc_l_max_bound = Euler_Labor(k_grid[i], [k_max,l_max],p)
        foc_k 
        # I have to check the bounds on capital 
        if  foc_l_min_bound<0 # In this case, I would like to work less, but I cant, I am already working the minimum -> Corner Solution
            G_l[i] = l_min
            G_kp[i] = k_min
            V[i] = -Obj_Fun([k_min,l_min])
        elseif foc_l_max_bound>0
            G_l[i] = l_max
            G_kp[i] = k_max
            V[i] = -Obj_Fun([k_max,l_max])
        else
            if 


        # Euler_Labor(k,x,p::Par)
        
        
        
        # Compute the implied consumption
        # I am going to calculate it with two bounds. One is for working full and the other one is working zero
        c_imp = z*(k_grid[i]^α)*(G_l[i]^1-α)+(1-δ)*k_grid[i]-G_kp[i]
        if c_imp<c_min
            V[i] = 
        Obj_Fun(k,x,p) = -utility(k_grid[i],x,p) - β*Vp.(x[1])
        # Min and max kp given current k
        #################################################################################
        # Not so sure about these bounds
        #################################################################################
        kp_min = 0.0001 # Suggested by Emmanuel
        # kp_min = k_grid[1]/(1-δ) # The bound I thought of
        kp_max = min(z*k_grid[i]^α  + (1-δ)*k_grid[i] - c_min, k_grid[end])
        #kp_max = 0.0
        # Bounds for l
        l_min = 1E-16
        l_max = 1.0
        # Now check if they bind
        # Lower bound 
        lb = [kp_min,l_min]
        #dObj_min = ForwardDiff.gradient(Obj_Fun,lb)
        dObj_min = -dutility(k_grid[i],lb,p) .- β*dVp(lb[1])
        if dObj_min[1]>0 && dObj_min[2]>0
            V[i]    = -Obj_Fun(lb)
            G_kp[i] = lb[1]
            G_l[i] = lb[2]
        else
        # Upper bound
        ub = [kp_max, l_max]
            dObj_max = -dutility(k_grid[i],ub,p) .- β*dVp(ub[1])
            #dObj_max = ForwardDiff.gradient(Obj_Fun,ub)
        #if dObj_max[1]<0 && dObj_max[2]<0
        #    V[i]    = -bj_Fun(ub)
        #    G_kp[i] = ub[1]
        #    G_l[i] = ub[2]
        #else
        # Now the minimization process
        inner_optimizer = NelderMead()
        l_ss,k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
        initial_x = [(kp_min+kp_max)/2,l_ss]
        #od = OnceDifferentiable(Obj_Fun, initial_x; autodiff=:forward) # No Funciona bien
        #min_result= optimize(od,lb,ub,initial_x, Fminbox()) # No funciona bien
        #min_result = optimize(x->Obj_fn(k_grid[i],x,p),[kp_min,l_min],[kp_max,l_max],initial_x,Fminbox(inner_optimizer))
        min_result = optimize(x->Obj_fn(k_grid[i],x,p),lb,ub,initial_x,Fminbox(inner_optimizer))
        println(" Cheking it solved the minimization")
        # Check result
        converged(min_result) || error("Failed to solve Bellman max in $(iterations(min_result)) iterations")
        # Record results
        V[i]     = -min_result.minimum
        G_kp[i]  = min_result.minimizer[1]
        G_l[i]  = min_result.minimizer[2]
        #end
        #end
    end
    G_c = z.*(collect(k_grid).^α).*(G_l.^(1-α)) .+(1-δ).*collect(k_grid) .-G_kp
    # Return Results
         println("T(V) = $V")
    return V, G_kp, G_c, G_l

end
#-----------------------------------------------------------