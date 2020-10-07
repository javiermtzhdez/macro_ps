####################################################################
####################################################################
############# Advance Macroeconomics II PS3
####################################################################
####################################################################

######### This part is just the environment activation
cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/")

######### Run these commented commands only the first time
using Pkg
Pkg.activate("ps3")
cd("/Users/cyberdim/Dropbox/github/macro/macro_ps/ps3/")
#Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
#Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
#Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
#Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
#Pkg.add("Latexify")

using Plots
using Dierckx 
using Interpolations 
using ForwardDiff 
using LaTeXStrings
using Latexify

# And using the module I coded
include("interpol.jl")
using .INTERPOL

# The functions to be interpolated are:
# u(c) = log(c)
# u(c) = c^(1/2)
# u(c) = (c^(1-σ))/(1-σ) (for σ=2,5,10)

# We have to interpolate the functions using three different methods.

# For each function and for each method:
# a) Interpolate over the domain [0.05, 2] using grids of different size. Plot them all in the same graph

#mkpath("Figures")
grid_1 = range(0.05,2,length=100)


nom = ["Log","Square","CES σ=2", "CES σ=5", "CES σ=10"]
funs = [lc,sc,crra2, crra5, crra10]
for i in 1:5
    name = nom[i]
    for n = (3,13,15,20)
        fn = funs[i].(grid_1)
        # Grid of nodes for interpolation
        xi = collect(range(0.05,2;length=n)) ; # Collect makes it an array instead of a collection
        yi = funs[i].(xi) # the corresponding y-coordinates
        # Interpolation
        interp=map(z->newton(xi,yi,z),grid_1) # Interpolating poly for the data
        # Plot
        gr()
        plot(title="Interpolation $name n=$n  (Newton Polynomial)")
        plot!(grid_1,fn,linewidth=3,label = "Function: $name",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(grid_1,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("Figures/Newton $name _n_$n")
    end
end



err=[poly_inter_err.(0.05,2,i,funs) for i in [3,13,15]]
Err=[err[1] err[2] err[3], err[4]]

copy_to_clipboard(true)
mdtable(Err, side=nom, head=["n=3","n=13","n=20", "n=30"],fmt = FancyNumberFormatter(4))





nom = ["Log","Square","CES σ=2", "CES σ=5", "CES σ=10"]
funs = [lc,sc,crra2, crra5, crra10]
for i in 1:5
    name = nom[i]
    for n = (3,13,15,20)
        fn = funs[i].(grid_1)
        # Grid of nodes for interpolation
        xi = collect(range(0.05,2;length=n)) ; # Collect makes it an array instead of a collection
        yi = funs[i].(xi) # the corresponding y-coordinates
        # Interpolation
        interp=map(z->CubicNaturalEval(xi,yi,z),grid_1) # Interpolating poly for the data
        # Plot
        gr()
        plot(title="Interpolation $name n=$n (Cubic Spline Natural)")
        plot!(grid_1,fn,linewidth=3,label = "Function: $name",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(grid_1,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("Figures/CSN_$name _n_$n")
    end
end





