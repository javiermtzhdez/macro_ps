####################################################################
####################################################################
############# Advance Macroeconomics II PS 5
############# Part A
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
mkpath("Figures")

using Random
using Distributions
using Parameters
using Plots
using Plotly
using LinearAlgebra
using Statistics
using StatsBase


######################## PART 1

# Define a markov process struct
    # Generate structure for markov processes using Parameters module
    @with_kw struct MP
        # Model Parameters
        N::Int64 # Number of states
        grid     # Grid of discrete markov process
        Π        # Transition matrix
        PDF      # Stationary distribution
        CDF      # Stationary distribution
    end
#

#= 
Let the AR(1) process be z' = \rho z + \eta 
where z~ Normal(0, sigma_eta)

a) Simulate a markov Chain for z of 10k draws, starting at z_0 = 0 
=#

# For replicability purposes I will set the seed at 123
Random.seed!(123)

# Now create a 10k vector for the z starting at z_0 = 0
z0 = 0
size = 10000
z = zeros(size,2)
z[1,1] = 1
z[1,2] = z0

# The values for sigma_eta and rho are not given in the problem, so I will use the same that are stated for the problem 2
ρ_z = 0.9
d = Normal(0, 0.1)
# pr = rand(d,1)
for i in 2:10000
    z[i,1] = i
    eta_t = rand(d,1)
    z[i,2] = ρ_z*z[i-1,2] + eta_t[1]
end

# Plot
gr()
Plots.plot(1:150, z[1:150,2], label = "zp", linestyle = :dashdotdot, legend = :outerbottomright )
xlabel!("t")
title!("AR (1) process")
Plots.savefig("Figures/original_AR1_z")

# Grid sizes required
N1 = 5
N2 = 15

##################
# Discretization
##################

# Discretization using Tauchen Method

Tauchen_N1 = Tauchen86(0.9,0.1,N1)
Tauchen_N2 = Tauchen86(0.9,0.1,N2)

# Discretization using Rouwenhorst Method

Rouwenhorst_N1 = Rouwenhorst95(0.9,0.1,N1)  # It has a grid size, a grid, transition probability matrix, pdf and cdf
Rouwenhorst_N2 = Rouwenhorst95(0.9,0.1,N2)

##################
### Simulation
##################

# Tauchen
N1_T_sim  = Simulation_MC(size, Tauchen_N1 , 100)
N2_T_sim  = Simulation_MC(size, Tauchen_N2 ,100)

# Rouwenhorst
N1_R_sim  = Simulation_MC(size, Rouwenhorst_N1 ,100)
N2_R_sim  = Simulation_MC(size, Rouwenhorst_N2 ,100)

##################
### Graphs
##################

# Tauchen
gr()
plt = Plots.plot(title="Tauchen Simulation Paths",legend=:outerright,foreground_color_legend = nothing,background_color_legend = nothing)
Plots.plot!(1:100, z[1:100,2], label = "Original", linestyle = :dashdotdot, legend = :outerbottomright )
Plots.plot!(1:100,N1_T_sim[1:100],linewidth=3,label="N=5",linestyle=(:dash))
Plots.plot!(1:100,N2_T_sim[1:100],linewidth=3,label="N=15",linestyle=(:dot))
Plots.xlabel!("Time")
savefig("./Figures/Tauchen_Simulation_Paths.pdf")

# Rouwenhorst
gr()
plt = Plots.plot(title="Rouwenhorst Simulation Paths",legend=:outerright,foreground_color_legend = nothing,background_color_legend = nothing)
Plots.plot!(1:100, z[1:100,2], label = "Original", linestyle = :dashdotdot, legend = :outerbottomright )
Plots.plot!(1:100,N1_R_sim[1:100],linewidth=3,label="N=5",linestyle=(:dash))
Plots.plot!(1:100,N2_R_sim[1:100],linewidth=3,label="N=15",linestyle=(:dot))
Plots.xlabel!("Time")
savefig("./Figures/Rouwenhorst_Simulation_Paths.pdf")

##################
### Moments
##################

# Tauchen
mean_T_n1, std_T_n1, acorr_T_n1 = Moments_MC(N1_T_sim)
mean_T_n2, std_T_n2, acorr_T_n2 = Moments_MC(N2_T_sim)

# Rouwenhorst
mean_R_n1, std_R_n1, acorr_R_n1 = Moments_MC(N1_R_sim)
mean_R_n2, std_R_n2, acorr_R_n2 = Moments_MC(N2_R_sim)

# First four autocorrelations 
# Tauchen
acf_N1_T_sim = StatsBase.autocor(vec(N1_T_sim))
acf_N2_T_sim = StatsBase.autocor(vec(N2_T_sim))

# Rouwenhorst
acf_N1_R_sim = StatsBase.autocor(vec(N1_R_sim))
acf_N2_R_sim = StatsBase.autocor(vec(N2_R_sim))

##################
### Histograms
##################

# Tauchen N=5
plotly()
x = N1_T_sim
Plots.histogram(x, line=(0.3,0.2,:blue), fillcolor=[:blue :black], fillalpha=0.8, title="Tauchen N=5 Histogram")
Plots.xlabel!("grid")
Plots.ylabel!("Count")
savefig("./Figures/N1_hist_T.pdf")

# Tauchen N=15
plotly()
x = N2_T_sim
Plots.histogram(x, line=(0.3,0.2,:blue), fillcolor=[:blue :black], fillalpha=0.8, title="Tauchen N=15 Histogram")
Plots.xlabel!("grid")
Plots.ylabel!("Count")
savefig("./Figures/N2_hist_T.pdf")

# Rouwenhorst N = 5
plotly()
x = N2_R_sim
Plots.histogram(x, line=(0.3,0.2,:blue), fillcolor=[:blue :black], fillalpha=0.8, title="Rouwenhorst N=15 Histogram")
Plots.xlabel!("grid")
Plots.ylabel!("Count")
savefig("./Figures/N1_hist_R.pdf")

# Rouwenhorst N = 15
plotly()
x = N2_R_sim
Plots.histogram(x, line=(0.3,0.2,:blue), fillcolor=[:blue :black], fillalpha=0.8, title="Rouwenhorst N=15 Histogram")
Plots.xlabel!("grid")
Plots.ylabel!("Count")
savefig("./Figures/N2_hist_R.pdf")




# [subprogram part!]

# Tauchen Method

function Tauchen86(ρ,σ,N,Ω::Any=3)
    # Create z grid
        z = range(-Ω*σ,Ω*σ,length=N)
    # Define intermediate step length
        h = (z[2]-z[1])/2
    # Define auxiliary matrices
        z_0 = repeat(z ,1,N) # Matrix of today's z each row is a value, columns are equal
        z_1 = repeat(z',N,1) # Matrix of tomorrow's z each column is a value, rows are equal
    # Define intervals
        z_lim = zeros(N,N,2) # First matrix is lower bounds. Second matrix is uppor bounds.
        z_lim[:,1      ,1] .= -Inf
        z_lim[:,2:end  ,1] .=  ( z_1[:,2:end  ] - ρ*z_0[:,2:end  ] .- h )./σ
        z_lim[:,1:end-1,2] .=  ( z_1[:,1:end-1] - ρ*z_0[:,1:end-1] .+ h )./σ
        z_lim[:,end    ,2] .=  Inf
    # Define reference distribution
        # This line uses "Distributions"
        F(x) = cdf.(Normal(),x)
    # Fill in transition matrix
        Π_z = F.(z_lim[:,:,2]) - F.(z_lim[:,:,1])
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    # Get stationary distribution of markov chain
        PDF_z = real(eigvecs(Π_z')[:,end]); PDF_z = PDF_z/sum(PDF_z) ;
        CDF_z = cumsum(PDF_z)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end

# Rouwenhorst

function Rouwenhorst95(ρ,σ,N)
    # Define paramters for Rouwenhorst's approximation
        p = (1+ρ)/2
        q = p                   # Note: I am leaving q here for comparability with source
        ψ = σ*sqrt((N-1)/(1-ρ^2))
        s = (1-q)/(2-(p+q))     # Note: s=0.5, I leave it for comparability with source
    # Fill in transition matrix
    if N==2
        Π_z = [p 1-p ; 1-q q]
    else
        MP_aux = Rouwenhorst95(ρ,σ,N-1)
        o = zeros(N-1)
        Π_z = p*[MP_aux.Π o ; o' 0] + (1-p)*[o MP_aux.Π ; 0 o'] + (1-q)*[o' 0 ; MP_aux.Π o] + q*[0 o' ; o MP_aux.Π]
        # Adjust scale for double counting
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    end
    # Distribution
        PDF_z = pdf.(Binomial(N-1,1-s),(0:N-1))
        CDF_z = cumsum(PDF_z)
    # Create z grid
        z    = range(-ψ,ψ,length=N)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end

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

# Moments function for sample
function Moments_MC(z_MC)
    mean_MC      = mean(z_MC)
    std_MC       = std(z_MC)
    auto_corr_MC = cor(z_MC[1:end-1],z_MC[2:end])
    return mean_MC, std_MC, auto_corr_MC
end

