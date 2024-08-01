/*
function trialmove:
    i = random number from 0 to total_n_atoms
    en_0 = energy(x_i, i) -> energy of atom i
    x_new = x_i + (R-0.5)*step_size
    en_new = energy(x_new, i)

    if R<exp(-beta*(en_new-en_0))
        accept
    else
        reject and revert to x_i

function energy:
    for 1<= j <= total_n_atoms-1
        if i != j
            dr = x_i - x_j
            dr = dr - box*round(dr/box)
            r2 = dr*dr
            if r2 < rc2
                PE formula
    return PE

function widom:
    x_test = randomly generate a coordinate for the test atom
    en_0 = energy(x_test, total_n_atoms+1)
    wtest = wtest + exp(-β * entest)
    μex = − ln(wtest/M)/β

*/