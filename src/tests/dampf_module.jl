#   Created by Matteo Melloni
#   dampf_module
#   2 Mar 2024
module DAMPF

using ITensors
using PseudomodesTTEDOPA
using LinearAlgebra
using ProgressMeter

export traceVectorizedMPO, getLocOscGates, getIntGates, getElEvoCoefficients, evolveRho

# Include some variations of the code useful for testing
include("dampf_module/notoptimized.jl")

# General functions
# ===============================

"""
    traceVectorizedMPO(op::MPS, indexes::Vector{Index{Int64}}, localDim::Int64)

Return the trace of the MPO `op` in Generalized Gell-Mann basis with local dimension `localDim`
"""
function traceVectorizedMPO(op::MPS, indexes::Vector{Index{Int64}}, localDim::Int64)
    nindexes = length(indexes)
    el = [localDim^2 for j in 1:nindexes]

    V = ITensor(1.)
    for j in 1:nindexes
        V *= (op[j]*state(indexes[j],el[j]))
    end
    v = scalar(V)

    return v*sqrt(localDim)^nindexes
end

# Evolution gates
# ===============================

# Vibrational gates
# -------------------------------
@doc raw"""
    getLocOscGates(dt::Float64, freqs::Vector{Float64}, 
        temps::Vector{Float64},damps::Vector{Float64}, sys::Vector{Index{Int64}})

Return a vector of Trotter gates that act on each site `sys` of the MPS consisting of all the oscillators.

Each oscillator evolves with a time step `dt` determied by the local Hamiltonian

```math
H_v^m = \omega_m a^\dagger_m a_m
```

and the Lindblad damping

```math
\mathcal{D}_m\rho_E = 2 \gamma_m \left [ (1+n_{T_m}(\omega_m)) (a_m \rho_E a_m^\dagger - \frac{1}{2}\{a_m^\dagger a_m,\rho_E\}) + n_{T_m}(\omega_m) (a_m^\dagger \rho_E a_m - \frac{1}{2}\{a_m a_m^\dagger,\rho_E\})\right ],
```
where ``n_{T_m}(\omega_m) = \frac{1}{e^{\omega_m/T_m}-1}``.
``M`` goes from 1 to `length(sys)` and ``\gamma_m``, ``T_m``, ``\omega_m`` 
are the m-th element of the vectors `damps`, `temps`,  `freqs` respectively.

The time evolution is derived from the Liouville superoperator ``\mathcal{V} = -i[H_v^m,\cdot] + \mathcal{D}_m`` approximated by
the Suzuki Trotter decomposition. Which implied that the local evolution operator is defined as

```math
\rho_E(t+dt) = e^{-i[H_v^m,\rho_E]\textit{dt}} e^{\mathcal{D}_m\textit{dt}}\rho_E(t)
```

This local evolution operator has to be applied to each site of the MPS consisting of the oscillators.
"""
function getLocOscGates(dt::Float64, freqs::Vector{Float64},
    temps::Vector{Float64}, damps::Vector{Float64}, sys::Vector{Index{Int64}})
    ll = length(sys)
    gates = [
       exp(dt * (
       locHamOsc(freqs[j],sys[j]) + 
       2 * damps[j] * locDissipatorOsc(freqs[j],temps[j],sys[j])
       )) 
       for j in 1:ll]
    return gates
 end

#Super-operator for local dissipation:
function locDissipatorOsc(frequency::Float64, temperature::Float64, sysIdx::Index{Int64})
    if(temperature == 0.)
       return op("A⋅ * ⋅Adag",sysIdx)-0.5*op("N⋅",sysIdx) - 0.5*op("⋅N",sysIdx) 
    else
       avg = 1/expm1(frequency/temperature)
    return (1+avg) * (op("A⋅ * ⋅Adag",sysIdx)-0.5*op("N⋅",sysIdx) - 0.5*op("⋅N",sysIdx) ) +
    avg * (op("Adag⋅ * ⋅A",sysIdx)-0.5*op( "A⋅ * Adag⋅",sysIdx) - 0.5*op("⋅Adag * ⋅A",sysIdx) )
    end
 end

#Super-operator for Hamiltonian evolution
function locHamOsc(frequency::Float64, sysIdx::Index{Int64})
    return -1.0im * frequency * (op("N⋅",sysIdx) - op("⋅N",sysIdx))
end

# Interaction gates
# -----------------------------------------
@doc raw"""
    getIntGates(dt::Float64, coups::Matrix{Float64}, Q::Int64, sys::Vector{Index{Int64}}, nsites::Int64)

The system-bath interaction acts locally on each block ``\mathcal{O}^{m,n}`` of the total density matrix ρ.
Where ``m,n`` go from 1 to `nsites`.

`getIntGates` returns a matrix of vectors of Trotter gates for each density matrix block ``\mathcal{O}^{m,n}``.
The Trotter gates of each block act on each site `sys` of the MPS consisting of all the oscillators.

The time evolution with time step `dt` is determined by the interaction Hamiltonian

```math
H_{ev} = \sum_{n=1}^N \sum_{q=1}^Q k_{n,q} \ket{n}\hspace{-4pt}\bra{n} (a_{n,q} + a^\dagger_{n,q}),
```
where ``N`` goes from 1 to `nsites`, `Q` is the total number of oscillators per site,
and ``k_{n,q}`` is the q-th element of the n-th row in `coups`.
`coups` has to be a matrix of this kind [[``1 \; k_{11} ... k_{1Q}``]; ... [``i \; k_{n1} ... k_{nq} ... k_{nQ}``]; ... [``N \; k_{N1} ... k_{NQ}]``],
where ``i`` determines the electronic site coupled to the `Q` oscillators with coupling strength ``k_{nq}``.
"""
function getIntGates(dt::Float64, coups::Matrix{Float64}, Q::Int64, sys::Vector{Index{Int64}}, nsites::Int64)
    
    if size(coups, 2) != Q+1 throw(ArgumentError("Invalid coups")) end

    gateMat = Matrix{Vector{ITensor}}(undef,nsites,nsites)

    for m in 1:nsites
        for n in 1:nsites
            #Every oscillator is to be updated.
            #the update depends on the indices m,n
            appo = Array{ITensor}(undef,length(sys))

            for i in 1:length(sys)
                j = trunc(Int, (i-1)/Q) + 1 # goes from 1 to N
                q = (i-1)%Q + 1 # goes from 1 to Q
                pippo = (Int(coups[j,1]) == m ? coups[j,q+1] : 0.) * op("Asum⋅",sys[i])-
                (Int(coups[j,1]) == n ? coups[j,q+1] : 0.) * op("⋅Asum",sys[i])
                appo[i] = exp(dt * (-1.0im * pippo))
            end
            gateMat[m,n] = appo
        end
    end
    return gateMat
end

# Electron interaction evolution
# ============================================

@doc raw"""
    getElEvoCoefficients(elHam::Matrix{Float64}, dt::Float64, nsites::Int64)

Return the coefficients determined by the time evolution with time step `dt` 
of the density matrix block ``\mathcal{O}^{m,n}``, where ``m,n`` go from 1 to `nsites`.

Applying the time evolution determined by the electron Hamiltonian
```math
H_e = \sum_{n=1}^N \Omega_n \ket{n}\hspace{-4pt}\bra{n} + \sum_{m \neq n } J_{m,n} \ket{m}\hspace{-4pt}\bra{n}
```
to the total density matrix ρ, each density matrix block ``\mathcal{O}^{m,n}`` evolves as
```math
\mathcal{O}^{m,n}(\delta t) = \sum_{i,j=1}^N U_{m,i} \mathcal{O}^{i,j} U_{j,n}^\dagger = \sum_{i,j=1}^N C_{m,i,j,n} \mathcal{O}^{i,j},
```
where ``C_{m,i,j,n} = U_{m,i} * U_{j,n}^\dagger``.
"""
function getElEvoCoefficients(elHam::Matrix{Float64}, dt::Float64, nsites::Int64; cutoff=1e-12)
    evoEl = exp(-1.0im * dt * elHam)
    evoEldaga = evoEl'

    return [
                [
                    [ 
                        [
                            if abs(evoEl[m,i]*evoEldaga[j,l]) < cutoff zero(ComplexF64) else evoEl[m,i]*evoEldaga[j,l] end
                        for l in 1:nsites]
                    for j in 1: nsites]  
                for  i in 1:nsites]   
            for m in 1:nsites]
end

# Total density matrix time evolution
# ============================================

@doc raw"""
    evolveRho(rho::Matrix{MPS}, timestep::Float64, time::Float64, maxBondDim::Int64,
    locOscGates::Vector{ITensor}, intGates::Matrix{Vector{ITensor}},
    elEvoCoefficients::Vector{Vector{Vector{Vector{ComplexF64}}}}, 
    sys::Vector{Index{Int64}}, nsites::Int64, localDim::Int64)

Return the evolved reduced density matrix tracing out the thermal bath, therefore averaging on the states of the oscillators.
The evolution follows the solution, using Suzuki Trotter decomposition, of the master equation

```math
\frac{d\rho}{dt} = - i [H_v+H_{ev}+H_e,\rho] + \mathcal{D}\rho,
```
where ``H_e`` is the electron Hamiltonian, ``H_v`` is the vibrational Hamiltonian of the oscillators,
``H_ev`` is the system-bath interaction Hamiltonian, and ``\mathcal{D}`` is the Lindblad operator acting on each oscillator.

For more details about the parameters `locOscGates`, `intGates`, `elEvoCoefficients`
and the evolution dynamics look at the documentation of the following functions:
- `getLocOscGates` for the evolution determined by ``H_v`` and the Lindblad operator ``\mathcal{D}``
- `getIntGates` for the system-bath interaction evolution determined by ``H_ev``
- `getElEvoCoefficients` for the electronic dynamics determined by the electron Hamiltonian ``H_e``
"""
function evolveRho(rho::Matrix{MPS}, timestep::Float64, time::Float64, maxBondDim::Int64,
    locOscGates::Vector{ITensor}, intGates::Matrix{Vector{ITensor}},
    elEvoCoefficients::Vector{Vector{Vector{Vector{ComplexF64}}}}, 
    sys::Vector{Index{Int64}}, nsites::Int64, localDim::Int64)

    pops = [[[] for _ in 1:nsites] for _ in 1:nsites]

    @showprogress color=:white for _ in 0.:timestep:time
        # Duplicate current rho
        newRho = deepcopy(rho)
        
        # Measure population of every sites
        for m in 1:nsites
            for n in 1:nsites
                push!(pops[m][n], traceVectorizedMPO(rho[m,n], sys, localDim))
            end
        end
        
        # Calculate only the upper rigth part of the density matrix because it is hermitian
        Threads.@threads for m in 1:nsites
            for n in m:nsites
                newRho[m,n] = computePartialRho(rho, m, n, locOscGates, intGates, elEvoCoefficients, nsites, maxBondDim)
            end
        end

        # ρ(m,n) = ρ†(n,m), because the density matrix is hermitian
        for m in 1:nsites
            for n in 1:m-1
                newRho[m,n] = dag(newRho[n,m])
            end
        end
         
        # rho.delete() 
     
        rho = newRho
    end

    return pops
end

# Evolve the density matrix block m,n following the electronic evolution determined by H_e
function applyElEvo(rho::Matrix{MPS}, m::Int64, n::Int64, C::Vector{Vector{Vector{Vector{ComplexF64}}}}, nsites::Int64, maxBondDim::Int64)
    weightOlist = vcat(filter(x->!isnothing(x), [if abs(C[m][i][j][n]) != zero(ComplexF64) C[m][i][j][n] * rho[i,j] end for i in 1:nsites, j in 1:nsites])...)
    return add(weightOlist..., cutoff = 10e-7, maxdim=maxBondDim)
end

# Evolve truncating the resulting MPS at the end of all the sums
function applyElEvoTruncAtTheEnd(rho::Matrix{MPS}, m::Int64, n::Int64, C::Vector{Vector{Vector{Vector{ComplexF64}}}}, nsites::Int64, maxBondDim::Int64)
    weightOlist = vcat(filter(x->!isnothing(x), [if abs(C[m][i][j][n]) != zero(ComplexF64) C[m][i][j][n] * rho[i,j] end for i in 1:nsites, j in 1:nsites])...)
    pureOmn = add(weightOlist..., alg="directsum")
    return truncate(pureOmn, cutoff = 10e-7, maxdim=maxBondDim)
end

# Evolve the density matrix block m,n
function computePartialRho(rho::Matrix{MPS}, m::Int64, n::Int64,
    locOscGates::Vector{ITensor}, intGates::Matrix{Vector{ITensor}},
    C::Vector{Vector{Vector{Vector{ComplexF64}}}}, nsites::Int64, maxBondDim::Int64)
    newRho_mn = rho[m,n]
    
    #apply H_e
    # newRho_mn = applyElEvo(rho, m, n, C, nsites, maxBondDim)
    newRho_mn = applyElEvoTruncAtTheEnd(rho, m, n, C, nsites, maxBondDim)
    
    #apply H_ev
    newRho_mn = apply(intGates[m,n], newRho_mn)
    
    #apply H_v
    newRho_mn = apply(locOscGates, newRho_mn)
    
    return newRho_mn
end


end