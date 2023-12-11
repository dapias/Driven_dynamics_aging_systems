using Plots, LaTeXStrings
using QuadGK, SpecialFunctions, HypergeometricFunctions


###Parameters
t0 = 0.0
egrid = collect(0.01:0.02:10.0)
deltat = 0.05;
T = 0.4;
s = 0.25;
tfinal = 10^4;

tstops = round.(exp10.(range(0, stop = log10(tfinal), length = 40))[1:end-1], digits = 1);

###Functions repeatedly called
rho(E) = exp(-E);
G(E) = exp(-E/T)
###Corresponding arrays
rho_grid = rho.(egrid)
escape_grid = G.(egrid);

###Forward master equation with initial condition P(E,0) = ρ(E)
global pf_total = []
global t = t0

pf(E)  = rho(E) + (deltat)*(-G(E)*rho(E) + exp(-s)*rho(E)*quadgk(x-> rho(x)*exp(-(x)/T), 0, Inf)[1])
t = t0 +deltat
global pf_array = pf.(egrid);

while t < tstops[end]
    integral = sum(pf_array.*(exp.(-(egrid)/T)))*(egrid[2]  - egrid[1]) 
    global pf_array += deltat*(-escape_grid.*pf_array .+ exp(-s)*rho_grid.*integral)
    global t += deltat      
    if round(t, digits =2) in tstops
            println(t)
            push!(pf_total, pf_array)
    end
end


###Backward master equation with initial condition q(E, t) = 1
tstops_back = round.(tfinal .- tstops, digits = 2);

global Deltat = 0.0
q(E)  = 1 + (deltat)*(-G(E)*1 + exp(-s)*exp(-E/T))

Deltat += deltat
global q_array = q.(egrid);


global q_deriv = []
global q_total = []


while Deltat < tfinal
    integral = sum(q_array.*rho_grid)*(egrid[2]  - egrid[1])
    q1 = deepcopy(q_array)
    global q_array += deltat*(-escape_grid.*q_array .+ exp(-s)*exp.(-egrid/T)*integral)
    global Deltat += deltat
    if round(Deltat, digits =2) in tstops_back
        println(Deltat)
        push!(q_total, q_array)
        push!(q_deriv, (log.(q_array) .- log.(q1))/deltat)
    end
end

#q_correct = q_total[end:-1:1];
#q_deriv = q_deriv[end:-1:1]


Zs = [sum(q_total[end-i+1].*pf_total[i])*(egrid[2] - egrid[1]) for i in 1:length(tstops_back)];
ptotal = [q_total[end-i+1].*pf_total[i]/Zs[i] for i in 1:length(tstops)];

a_mean = [sum(ptotal[k].*(G.(egrid) .- q_deriv[end-k+1])) *(egrid[2] - egrid[1]) for k in 1:length(tstops)]


plot(tstops, a_mean, xscale = :log10, yscale = :log10, label = "", color = :green)
xx = tstops[20:32];
plot!(xx, xx.^(-(T+1))/(exp(s)-1)^2*exp(s)*T^2*gamma(T), ls = :dash, color = :red, lw = 2, label = "")
plot!(xx, 0.1*xx.^(-(1))/(exp(s)-1)^2*exp(s)*T^2*gamma(T), ls = :dash, color = :blue, lw = 2, label = "")

    

#a_mean = [( sum(G.(egrid).*ptotal[k] ) .- 
#sum( ((log.(q_correct[k+1]) .- log.(q_correct[k]))/(tstops[k+1] - tstops[k])).*ptotal[k]) )*(egrid[2] - egrid[1]) for k in 1:length(tstops)-1];
