# ================================================
### FMO complex parameters
# ================================================

nsites = 2 #number of sites
noscpersite = 6 #number of oscillators per site
nosc = nsites * noscpersite # define the total number of oscillator
localDim = 3 #local dimension of the oscillators
maxBondDim = 4 #maximal bond dimension of the MPSs
timestep = 0.5 #integration time-step in fs
time = 300. #simulation time in fs

### Parameters for the system dynamics in cm-1
# ------------------------------------------------

##### Sites
energies = [12430., 12405., 12175., 12315., 12625., 12500., 12450.] #Local energy Ωn of the sites
energies = energies[1:nsites]
exchangepersite = 80. # Exchange energy per every site
##### Oscillators
freqspersite = [247., 763., 1175., 1356., 1521., 160.] #Frequency of the oscillators per site
tempKelvin = 77. #Temperature in Kelvin
huangRhysFactors = [0.056, 0.133, 0.049, 0.019, 0.006, 0.164] #Huang Rhys factors coupling every oscillators per site
dampspersiste = [53., 76., 29., 29., 15., 133.]  #damping rates of the oscillators per site


# ===============================================
### Conversion to fs-1 
### and preparing for the simulation 
### (do not touch this part)
# ===============================================
_OmegaConv = 1.883651567e-4 # Conversion constant from cm-1 to fs-1 (2pi/T)
_Tconv = 0.6950348 # Conversion constant from Kelvin to cm-1

energies = _OmegaConv * energies
exchange = _OmegaConv * hcat([[i == j ? 0 : exchangepersite for i in 1:nsites] for j in 1:nsites]...) 
#Adjacency matrix of the interaction graph (real and symmetric)
freqs = vcat([_OmegaConv * freqspersite for i in 1:nsites]...) #Frequency of all the oscillators
temps = [_OmegaConv *_Tconv * tempKelvin for j in 1:nosc] #Temperature of the oscillators in fs-1
coups = vcat([hcat(i, [w*sqrt(s) for (w,s) in zip(freqs, huangRhysFactors)]...) for i in 1:nsites]...)
#coups = [[1 2. 2. 3.]; [2 2. 2. 3.]; [3 2. 2. 3.]; [4 2. 2. 3.]; [5 2. 2. 3.]] #Coupling of the oscillators to the sites: 
#is a list of couples, indexed by the oscillator number: every couples:
#[site,couplings...] 
damps = vcat([_OmegaConv * dampspersiste for i in 1:nsites]...) #damping rates of the oscillators

