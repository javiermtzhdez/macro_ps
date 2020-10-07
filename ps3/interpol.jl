# Interpolation module
module INTERPOL
export newton, poly_inter_err, CubicNatural, CubicNaturalEval

# The Newton Polynomial is taken from http://fsu.digital.flvc.org/islandora/object/fsu%3A657877

function diff(x::Array,y::Array)
    m = length(x) #here m is the number of data points. #the degree of the polynomial is m-1 a=Array{Float64}(undef,m)
    a = [y[i] for i in 1:m]
    # for i in 1:m
    #     a[i]=y[i]
    # end
    #a = [(y[i]-y[i-1])/(x[i]-x[i-j+1]) for j in 2:m, i in reverse(collect(j:m))]
    for j in 2:m
        for i in reverse(collect(j:m))
            a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
        end
    end
    return a
end


function newton(x::Array,y::Array,z)
    m=length(x) #here m is the number of data points, not the degree # of the polynomial
    a=diff(x,y)
    sum=a[1]
    pr=1.0
    for j in 1:(m-1)
        pr=pr*(z-x[j])
        sum=sum+a[j+1]*pr
    end
    return sum
end


function poly_inter_err(a,b,n,f::Function)
    a , b  = min(a,b), max(a,b)
    h = (b-a)/n
    x = [a+(i-1)*h for i in 1:(n+1)]
    Π = h^(n+1)*factorial(n)/4
    ξ = abs(maximum(f.(x)))
    err = Π*ξ/factorial(n+1)
    return err
end

function CubicNatural(x::Array,y::Array)
    m=length(x) # m is the number of data points
    n=m-1
    a=Array{Float64}(undef,m)
    b=Array{Float64}(undef,n)
    c=Array{Float64}(undef,m)
    d=Array{Float64}(undef,n)
    a = y
    # for i in 1:m
    #     a[i]=y[i]
    # end
    h = [x[i+1]-x[i] for i in 1:n]
    # h=Array{Float64}(undef,n)
    # for i in 1:n
    #     h[i]=x[i+1]-x[i]
    # end
    u=Array{Float64}(undef,n)
    u[1]=0
    for i in 2:n
        u[i]=3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1]
    end
    s=Array{Float64}(undef,m)
    z=Array{Float64}(undef,m)
    t=Array{Float64}(undef,n)
    s[1]=1
    z[1]=0
    t[1]=0
    for i in 2:n
        s[i]=2*(x[i+1]-x[i-1])-h[i-1]*t[i-1]
        t[i]=h[i]/s[i]
        z[i]=(u[i]-h[i-1]*z[i-1])/s[i]
    end
    s[m]=1
    z[m]=0
    c[m]=0
    for i in reverse(1:n)
        c[i]=z[i]-t[i]*c[i+1]
        b[i]=(a[i+1]-a[i])/h[i]-h[i]*(c[i+1]+2*c[i])/3
        d[i]=(c[i+1]-c[i])/(3*h[i])
    end
    return a, b, c, d
end

# CubicNaturalEval. The inputs are the value at which
# the spline is evaluated, w, and the x-coordinates of the data. The function first finds the
# interval [xi; xi+1]; i = 1; :::;m-1; w belongs to, and then evaluates the spline at w using the
# corresponding cubic polynomial.

function CubicNaturalEval(x::Array,y::Array,w)
    m=length(x)
    if w<x[1]||w>x[m]
        return print("error: spline evaluated outside its domain")
    end
    n=m-1
    p=1
    for i in 1:n
        if w<=x[i+1]
            break
        else
            p=p+1
        end
    end
    # p is the number of the subinterval w falls into, i.e., p=i means
    # w falls into the ith subinterval $(x_i,x_{i+1}), and therefore
    # the value of the spline at w is
    # a_i+b_i*(w-x_i)+c_i*(w-x_i)^2+d_i*(w-x_i)^3.
    a, b, c, d = CubicNatural(x,y)
    return a[p]+b[p]*(w-x[p])+c[p]*(w-x[p])^2+d[p]*(w-x[p])^3
end

end
