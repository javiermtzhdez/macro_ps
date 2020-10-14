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

#Pkg.add("Sobol")

#using Optim 
#    using Optim: converged, maximum, maximizer, minimizer, iterations, optimize
#using Parameters
#using Roots
#using IntervalArithmetic, IntervalRootFinding
#using Plots
#using Dierckx 
#using Plots
#using LinearAlgebra
#using PrettyTables
#using Distributions
#using Random
#using Interpolations
#using ForwardDiff 
using Sobol


cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps4/")

# Rosenbrock Function
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

# The function goes from R^k -> R^l
# Step 0
# Determine bounds for each parameter 
# Generate a sequence of Sobol' points with length N
k=2 # domain
l=1 # range
N = length_seq
s = SobolSeq(lb,ub)
Sseq = zeros(N,2)
Sseq[1,1] = next!(s)[1]
Sseq[1,2] = next!(s)[2]
for i in 2:N
    Sseq[i,1]=next!(s)[1]
    Sseq[i,2]=next!(s)[2]
end
fval = zeros(N,2)
for i in 1:N
    fval[i,1] = f(Sseq[i,:])
    fval[i,2] = i
end

# Now we set N*
N_star = 0.1*N
fval = sortperm(fval[:,1])
fval = fval[1:N_star,:]




function Tiktak(domain,range,length_seq,lb::Array,ub::Array,f::Function)
    k=domain
    l=range
    N = length_seq
    # lb and ub are arrays of length k that gives you the hypercube[lb,up]^k

    



    


end
# Let $p$ be a J-dimensional parameter vector with generic element $p^j$ , with $j=1,...,J$
enerate a sequence of Sobol' points with length N















# Note that this is a function R^2 -> R


