using Plots, LaTeXStrings, Distributions, Arpack, JLD
using QuadGK, SpecialFunctions, HypergeometricFunctions



Ts = [0.2, 0.6, 2.0]
rr = Exponential()

rho(E) = exp(-E);

n = 2^13


for T in Ts
    s_grid = -1.5:0.03:0.0;
    free = similar(s_grid)
    acti = similar(s_grid);
    W = zeros(n,n);
    G(E) = _₂F₁(1, T, 1+T, -exp(E/T))
    
    
    for k in 1:length(s_grid)
        ####Master operator. The energies are sampled directly from the exponential distribution
        s = s_grid[k]
        es = rand(rr, n);
        for i in 1:n
            W[i,i] = -G(es[i])
        end
        for i in 1:n
            for j in 1:i-1
                W[i, j] = exp(-s)*1/(1+ exp(-(es[i] - es[j])/T))*1/(n)
            end
            for j in i+1:n
            W[i, j] = exp(-s)*1/(1+ exp(-(es[i] - es[j])/T))*1/(n)
            end
        end
        ######
        #    Wp = W'
        #top_left = eigs(Wp, nev = 1, which = :LR)
        
        top_right = eigs(W, nev = 1, which = :LR)
        V_0 = real.(top_right[2])[:,1];
        #U_0 = real.(top_left[2])[:,1];
        U_0 = V_0.*exp.(-es/T)
        
        Z = sum(U_0.*V_0)
        
        println(k)
        free[k] = real.(top_right[1][1])
        acti[k] =   sum(exp.(-s)*U_0.*[sum(V_0 ./(1 .+ exp.(-(es[i] .- es)/T))) for i in 1:n])/(Z*n)
    end
    free_mod = deepcopy(free);
    acti_mod = deepcopy(acti);
    
    s_grid_mod = vcat(s_grid, 0.1:0.1:0.5);
    
    acti_mod = vcat(acti_mod, zeros(5));
    free_mod = vcat(free_mod, zeros(5));
    


    save("free_energy_and_activity_T=$(T)_bm.jld", "s_grid", s_grid_mod, "acti", acti_mod, "free", free_mod,
         "n", n)
end
