# Testing percentage of retained terms in the electronic sum evolution
function countRetainedCoefficients()
    countnotzeros = 0
    csorted = []
    for i in 1:nsites, j in 1:nsites, n in 1:nsites, m in 1:nsites
       push!(csorted, elEvoCoefficients[i][j][m][n])
       if elEvoCoefficients[i][j][m][n] != zero(ComplexF64) countnotzeros += 1 end
    end
 
    sort!(csorted, rev=true, lt=(x,y)->isless(abs(x),abs(y)))
    csorted = filter(x->!isnothing(x), [if abs(el)<2 el end for el in csorted])
    
    # UNCOMMENT TO PLOT THE DISTRIBUTIONS OF COEFFICIENTS
    # real_plot = abs.(csorted)
    # display(bar(real_plot, linewidth=0, bar_width=2, legend=false, title="(b)", titleloc=:left, yscale=:log10, ylabel=L"|C_{n,i,j,m}|", color=RGB(83/255,146/255,255/255), xlabel="Coefficients number"))
    # savefig("plots/Cdistribution2.pdf")
 
    # for elem in csorted
    #    println(abs(elem))
    # end
 
    return countnotzeros
 end