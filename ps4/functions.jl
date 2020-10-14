module functions
export psi, lambda, lss, f1

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























##################################################################
##################################################################
############ Below here is the end for the module

end