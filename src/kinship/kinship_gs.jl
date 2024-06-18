"""

     kinshipGs(climate::Array{Float64,2},ρ::Float64)

Computes a kinship matrix using the Gaussian Kernel.

# Arguments

- `climate` : A matrix of genotype, or climate information data. size(climate)=(r,m), such that `r` genotype markers (or days/years of climate factors,
            i.e. precipitations, temperatures, etc.), and `m` individuals (or environments/sites)
- `ρ` : A free parameter determining the width of the kernel. Could be attained empirically.

# Output

Returns a m x m symmetric (positive definite) matrix containing 1's on the diagonal.

"""
function kinship_gs(climate::Array{Float64,2},ρ::Float64)
 env=axes(climate,2);
 Kc=zeros(env,env);

    @views for c=env, r=c:length(env)
            Kc[r,c]=exp(-mean(abs.(climate[:,c]-climate[:,r]).^2)/ρ)
            Kc[c,r]=Kc[r,c]
    end

    return Kc

end

