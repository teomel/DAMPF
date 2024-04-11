# ==================================================================
#### NOT OPTIMEZED VERSIONS OF THE FUNCTIONS
# ==================================================================

export evolveRhoNOTOPTIMIZED, evolveRhoPARTIALLYOPTIMIZED, getElEvoCoefficientsNOTOPTIMIZED

function evolveRhoNOTOPTIMIZED(rho::Matrix{MPS}, timestep::Float64, time::Float64, maxBondDim::Int64,
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
        for m in 1:nsites
            for n in 1:nsites
                newRho[m,n] = computePartialRho(rho, m, n, locOscGates, intGates, elEvoCoefficients, nsites, maxBondDim)
            end
        end
     
        rho = newRho
    end

    return pops
end

# Without multithreading
function evolveRhoPARTIALLYOPTIMIZED(rho::Matrix{MPS}, timestep::Float64, time::Float64, maxBondDim::Int64,
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
        for m in 1:nsites
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

function getElEvoCoefficientsNOTOPTIMIZED(elHam::Matrix{Float64}, dt::Float64, nsites::Int64)
    evoEl = exp(-1.0im * dt * elHam)
    evoEldaga = evoEl'

    return [
                [
                    [ 
                        [
                            evoEl[m,i]*evoEldaga[j,l]
                        for l in 1:nsites]
                    for j in 1: nsites]  
                for  i in 1:nsites]   
            for m in 1:nsites]
end