#   Created by Matteo Melloni
#   dampf
#   2 Mar 2024

using ITensors
using PseudomodesTTEDOPA
using LinearAlgebra
using Plots
using Revise
using FileIO # to save objects on file
using BenchmarkTools
# includet is in the Revise package, and it's used to not reload dampf_module.jl at every execution on the same REPL
includet("dampf_module.jl")
using .DAMPF

### Setting up system parameters
# ===============================================

include("nsites_params.jl")

println("Testing --- Number of sites: ", nsites)
println("Testing --- Number oscillators per site: ", noscpersite)
println("Testing --- Simulation time: ", time)
println("Testing --- Timestep: ", timestep)
println("Testing --- Maximum bond dimension χ: ", maxBondDim)
println("Testing --- Bath local dimension: ", localDim)
println("Testing --- Total number of oscillators: ", nosc)

### Setting up the initial global density matrix ρ
# =================================================

bath = siteinds("HvOsc", nosc, dim=localDim)
bathMPS = chain(
   [parse_init_state_osc(bath[j], "thermal", frequency=freqs[j], temperature=temps[j]) for j in 1:nosc]...
)

orthogonalize!(bathMPS, 1)

#rho_E is the building block of the matrix ρ...
rho_E = bathMPS;

rho_S = ComplexF64.(zeros(nsites, nsites));
rho_S[1,1] = 1.0;

# Defining the global density matrix ρ
rho = rho_S .* fill(rho_E, nsites, nsites);

### Setting the evolution gates
# ================================================

locOscGates = getLocOscGates(timestep, freqs, temps, damps, bath);
intGates = getIntGates(timestep, coups, noscpersite, bath, nsites);

### Setting the electron evolution coefficients determined by the electron Hamiltonian
# ====================================================================================

elHam = exchange + Diagonal(energies)

CoeffCutoff = 1e-4
elEvoCoefficients = getElEvoCoefficients(elHam, timestep, nsites, cutoff=CoeffCutoff)
# elEvoCoefficientsNOTOPTIMIZED = getElEvoCoefficientsNOTOPTIMIZED(elHam, timestep, nsites)

include("utils.jl") # inlcude the function countRetainedCoefficients() with the variables of the current session

println("Testing --- C coefficients cutoff: ", CoeffCutoff)
println("Testing --- Percentage of retained electronic C coefficients: ", trunc(countRetainedCoefficients()/nsites^4 * 100, digits=1), "%")

println(".")

### Evolve the total density matrix and obtain the reduced density matrix dynamics
# ================================================================================

###  SOME PARTS HAVE TO BE UNCOMMENTED TO RUN DIFFERENT OPTIONS ###

### EvolveRho
maxBondDim = 4
# pops = load("plots/popscompare3.jld2", "pops")
println("χ=$maxBondDim: ")
pops = evolveRho(rho, timestep, time, maxBondDim, locOscGates, intGates, elEvoCoefficients, bath, nsites, localDim)
save("plots/pops$(nsites).jld2", "pops", pops) # save pops with name "pops". To load use: pops = load(filename, "pops")
arrTraceReducedRho = [sum([pops[m][m][j] for m in 1:nsites]) for j in 1:length(pops[1][1])]

# ### EvolveRho to be compared
# maxBondDimCompare = 6
# # popscompare = load("plots/popscompare.jld2", "pops")
# println("χ=$maxBondDimCompare: ")
# popscompare = evolveRho(rho, timestep, time, maxBondDimCompare, locOscGates, intGates, elEvoCoefficients, bath, nsites, localDim)
# save("plots/popscompare$(nsites).jld2", "pops", pops)
# arrTraceReducedRhoCompare = [sum([popscompare[m][m][j] for m in 1:nsites]) for j in 1:length(popscompare[1][1])]

println("")

### Plots of the result
# ================================================================================

#### Plot population
x=range(0,time,length=length(pops[1][1]))
real_plot = [real(pops[m][m])./real(arrTraceReducedRho) for m in 1:nsites]
labels = hcat(["\\rho_{$(m)$(m)}" for m in 1:nsites]...)
display(plot(x, real_plot, title="Populations", label=labels, xlabel="Time [fs]", ylabel="Populations"))
savefig("plots/populations$(nsites).pdf")

# #### Plot population overlapping different bond dimensions
# x=range(0,time,length=length(pops[1][1]))
# real_plot = [real(pops[m][m])./real(arrTraceReducedRho) for m in 1:nsites]
# real_plotcompare = [real(popscompare[m][m])./real(arrTraceReducedRhoCompare) for m in 1:nsites]
# labels = hcat(["\\rho_{$(m)$(m)} (\\chi = $(maxBondDim))" for m in 1:nsites]...)
# labelscompare = hcat(["\\rho_{$(m)$(m)} (\\chi = $(maxBondDimCompare))" for m in 1:nsites]...)
# plot(x, real_plot, title="Populations", label=labels, linewidth=3, xlabel="Time [fs]", ylabel="Populations")
# display(plot!(x, real_plotcompare, label=labelscompare, linewidth=3, linestyle=:dash))
# savefig("plots/populationsCompare$(nsites).pdf")

# Plot coherence
x=range(0,time,length=length(pops[1][1]))
real_plot = [real(pops[m][l])./real(arrTraceReducedRho) for m in 1:nsites for l in (m+1):nsites]
imag_plot = [imag(pops[m][l])./real(arrTraceReducedRho) for m in 1:nsites for l in (m+1):nsites]
# real_labels = hcat(["Re \\rho_{$(m)$(l)}" for m in 1:nsites for l in (m+1):nsites]...)
# imag_labels = hcat(["Im \\rho_{$(m)$(l)}" for m in 1:nsites for l in (m+1):nsites]...)
plot(x, real_plot, title="Coherences", legend=false, xlabel="Time [fs]", ylabel="Coherences")
display(plot!(x, imag_plot, legend=false, linestyle=:dash))
savefig("plots/coherences$(nsites).pdf")
print("")