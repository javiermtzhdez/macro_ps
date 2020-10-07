####################################################################
####################################################################
############# Advance Macroeconomics II PS1
####################################################################
####################################################################

# As in Fortran, I think is better to do variable declaration first and then the code...

# IMPLICIT NONE HAHA

######################
# [specification part]
######################

#cd("/Users/cyberdim/Dropbox/github/advmacro_ps/ps1/")
#using Pkg
#Pkg.activate("ps1")


#using Parameters
#using Plots
#using Plotly
#using PlotlyBase

#
# Environment activation
#


#
# Creating the folders for our graphs
#
#mkpath("Figures") # Creates the folder Figures ensuring the appropriate path exists


#
# PARAMETERSS
#

@with_kw struct Par  #  The macro @with_kw which decorates a type definition to allow default values and a keyword constructor, in this case Par
    # This are the structural parameters which wont move throughout of the model
    #z::Float64 = 1 ; # Productivity
    z::Float64 = 1.0 ; # Productivity for the increase for plotting
    α::Float64 = 1/3; # Capital share, to write the alpha you do it as in latex \alp ha
    β::Float64 = 0.98; # Discount factor
    # VFI values
    max_iter::Int64 = 2000 ; # Maximum number of iterations
    dist_tol::Float64 = 1E-9 ; # Tolerance for the distance (tol for the fixed point)
    # Howard's Policy iterations
    H_tol::Float64 = 1E-9; # Tolerance for policy function iteration.
end

# Allocate parameters to object p for future calling? why not calling them directly? 
# Is this like creating a vector of parameters
p=Par()


######################
# [execution part]
######################


#
# Utility Function
#
function utility(k,kp,p::Par)
    @unpack z, α = p # The @unpack and @pack! macros work to unpack types, modules, and dictionaries 
    c = z*k.^α - kp # Note here that .^ is a vectorized "dot" operation, automatically defined to 
    # perform ^ element by element in arrays.
    if c>0
        return log(c)
    else
        return -Inf # I dont like this minus infinity stuff. its like fortran, i dont like it personally
    end
end


#
# Steady state values
#

# The steady state values are obtained via a function that takes in parameters and provides ss values
# Parameters: Productivity (z), returns to scale (alpha) and the discount factor (beta)
# It spits ss values for capital, production, consumption, interest rate, and wages

function SS_values(p::Par)
    @unpack z,α, β = p
    k_ss = (β*α*z)^(1/(1-α))
    y_ss = z*k_ss^α
    c_ss = y_ss - k_ss
    r_ss = α*y_ss/k_ss
    w_ss = (1-α)*y_ss
    return k_ss,y_ss,c_ss,r_ss,w_ss
end


# Test steady state function
k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
println(" ")
println("------------------------")
println(" Steady State variables") 
println("   Quantities: k = $k_ss; y = $y_ss; c = $c_ss;") # explanation of this line http://www.jlhub.com/julia/manual/en/function/println
println("   Prices:     r = $r_ss; w = $w_ss;")
println("------------------------")
println(" ")

#Steady state values
# Quantities: k = 0.18670555150547338; y = 0.5715476066494083; c = 0.3848420551439349;
# Prices:     r = 1.0204081632653061; w = 0.38103173776627225;

# So far what we have are the steady state values. 

# Functions to calculate the variables needed. So for any value we enter of capital,x, this gives you back
# policy function (k(t+1), prod, consumption, interest_rate, wages)
function transition_values(x,p::Par)
    @unpack z,α, β = p
    g_t = α*β*(z*1.05)*x^α  
    y_t = (z*1.05)*x^α
    c_t = (z*1.05)*x^α*(1-α*β) # y_t - g_t
    r_t = α*(z*1.05)/(x^(1-α))
    w_t = (1-α)*(z*1.05)*x^α
    return [g_t, y_t, c_t, r_t, w_t]
end




##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# b) Productivity increases 5% permanently
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


# Let start an array (kinda like a numpy array)

m=100
n=7
vals_grid = zeros(m,n)
# The number of columns correspond to time(1), capital(2), policyfn(3), production(4), consumption(5), interest_rate(6), wages(7).
# Note here that Julia starts counting in the naturals, not in zero as Python, (well I guess that depends on your definition of the natural numbers)
# lets numerate the first column
for i in 2:100
    vals_grid[i,1]=i-1
end

vals_grid[1,2]=k_ss

for i in 1:100
    if i<100
        local vals = transition_values(vals_grid[i,2],p::Par)
        vals_grid[i,3] = vals[1]
        vals_grid[i,4] = vals[2]
        vals_grid[i,5] = vals[3]
        vals_grid[i,6] = vals[4]
        vals_grid[i,7] = vals[5]
        vals_grid[i+1,2] = vals[1]
    elseif i==100
        local vals = transition_values(vals_grid[i,2],p::Par)
        vals_grid[i,3] = vals[1]
        vals_grid[i,4] = vals[2]
        vals_grid[i,5] = vals[3]
        vals_grid[i,6] = vals[4]
        vals_grid[i,7] = vals[5]
    end
end

for i in 1:100
    println("    t        ", "        k            ", "       k+1       ", "          y          ", "          c         ", "        r         ", "         w             " )
    println(vals_grid[i,:])
end

# Now comes the plotting
# I am more fond of using Plotly than GR but still, GR is installed already in the environment
# Reaching the steady state is very quick, so I am only gonna use t=15 at most.

# When plotting, in Plots.jl every column is a series, so it is possible to plot multiple lines by plotting a matrix
# of values and each column is interpreted as a separate line

# So I am going to create a grid of size 15, i, which I will have the steady state value of the variable
# and the values for reaching it.

# Consumption
c_grid = zeros(15,3)
for i in 1:15
    c_grid[i,1]=i
    c_grid[i,2]=c_ss
    c_grid[i,3]=vals_grid[i,5]
    println(i, "   C_ss   " , c_ss, "   C_t   ", vals_grid[i,5])
end

c2 = c_grid[1:15,2]
c3 = c_grid[1:15,3]
c4 = fill(c_grid[15,3],15)
gr()
plot(1:15, c2, label = "C_ss original")
plot!(1:15,c3, label = "Consumption path", linestyle = :dashdotdot, legend = :outerbottomright )
plot!(1:15,c4, label = "C_ss new", linestyle = :dashdot, legend = :outerbottomright )
xlabel!("Period")
title!("Consumption")
savefig("Figures/consumption_b")

# Capital
k_grid = zeros(15,3)
for i in 1:15
    k_grid[i,1]=i
    k_grid[i,2]=k_ss
    k_grid[i,3]=vals_grid[i,2]
end
c2 = k_grid[1:15,2]
c3 = k_grid[1:15,3]
c4 = fill(k_grid[15,3],15)
gr()
plot(1:15, c2, label = "K_ss")
plot!(1:15,c3, label = "Capital path", linestyle = :dashdotdot, legend = :outerbottomright )
plot!(1:15,c4, label = "K_ss new", linestyle = :dashdot, legend = :outerbottomright )
xlabel!("Period")
title!("Capital")
savefig("Figures/capital_b")

# Interest rate
r_grid = zeros(15,3)
for i in 1:15
    r_grid[i,1]=i
    r_grid[i,2]=r_ss
    r_grid[i,3]=vals_grid[i,6]
end
c2 = r_grid[1:15,2]
c3 = r_grid[1:15,3]
c4 = fill(r_grid[15,3],15)
gr()
plot(1:15, c2, label = "R_ss")
plot!(1:15,c3, label = "Interest rate path", linestyle = :dashdotdot, legend = :outerbottomright )
plot!(1:15,c4, label = "R_ss new", linestyle = :dashdot, legend = :outerbottomright )
xlabel!("Period")
title!("Interest rate")
savefig("Figures/interest_rate_b")

# Wage
w_grid = zeros(15,3)
for i in 1:15
    w_grid[i,1]=i
    w_grid[i,2]=w_ss
    w_grid[i,3]=vals_grid[i,7]
end

c2 = w_grid[1:15,2]
c3 = w_grid[1:15,3]
c4 = w_grid[15,3]
c4 = fill(w_grid[15,3],15)
gr()
plot(1:15, c2, label = "W_ss")
plot!(1:15,c3, label = "Wage path", linestyle = :dashdotdot, legend = :outerbottomright )
plot!(1:15,c4, label = "W_ss new", linestyle = :dashdot, legend = :outerbottomright )
xlabel!("Period")
title!("Wage")
savefig("Figures/wage_b")

# Production
y_grid = zeros(15,3)
for i in 1:15
    y_grid[i,1]=i
    y_grid[i,2]=y_ss
    y_grid[i,3]=vals_grid[i,4]
end
c2 = y_grid[1:15,2]
c3 = y_grid[1:15,3]
c4 = y_grid[15,3]
c4 = fill(y_grid[15,3],15)
gr()
plot(1:15, c2, label = "Y_ss")
plot!(1:15,c3, label = "Production", linestyle = :dashdotdot, legend = :outerbottomright )
plot!(1:15,c4, label = "Y_ss new", linestyle = :dashdot, legend = :outerbottomright )
xlabel!("Period")
title!("Production")
savefig("Figures/production_b")


# [subprogram part]