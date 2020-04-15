using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;


N=50;
N_steps = 100;
δ_steps = 0.1;

#Build the circulant vector for the a†a part of the Γ with exponential decaying correlations
adaa = zeros(Complex{Float64},N);#
for i=1:div(N,2)
    adaa[i] = exp(-i*0.15)*(rand()+im*rand())
end
adaa[((div(N,2))+1):N]= reverse(adaa[1:div(N,2)]);
#Build the translational invariant a†a part of the Γ
Γ_adaa = Fu.Circulant(adaa);
Γ_adaa = (Γ_adaa+Γ_adaa')/2.

#Build the circulant vector for the aa part of the Γ
aa = zeros(Complex{Float64},N);
aa[2] = rand()+im*rand();
aa[3] = rand()+im*rand();
#Build the translational invariant aa part of the Γ
Γ_aa = Fu.Circulant(aa)
Γ_aa = (Γ_aa-transpose(Γ_aa))/2.;

#Build the translational invariant Γ
Γ= zeros(Complex{Float64}, 2N,2N);
Γ[(1:N),(1:N)]          = Γ_adaa;
Γ[(1:N).+N,(1:N).+N]    = (I-Γ_adaa)';
Γ[(1:N).+N,(1:N)]       = Γ_aa;
Γ[(1:N),(1:N).+N]       = -conj(Γ_aa);
Fu.Print_complex_matrix("Γ",Γ)

H   = Fu.Build_hopping_hamiltonian(N,true);
D,U = Diag_h(H);


Γ_evolved   = copy(Γ);
adaa        = zeros(Complex{Float64}, N_steps)
aa          = similar(adaa);

#Start the time evolution cycle
#at each cycle it saves the value of two correalotors
adaa[1]     = Γ_evolved[1,2];
aa[1]       = Γ_evolved[N+1,2];
for t=2:N_steps
    global Γ_evolved = Fu.Evolve(Γ_evolved,(D,U),δ_steps);
    adaa[t]     = Γ_evolved[1,2];
    aa[t]       = Γ_evolved[N+1,2];
end

figure("Evolutions")
plot(real.(adaa), color="black", label=L"$\mathfrak{R}(\langle a_1^{\dagger}a_2 \rangle)(t)$");
plot(imag.(adaa), color="black",linestyle="--", label=L"$\mathfrak{I}(\langle a_1^{\dagger} a_2 \rangle)(t)$");
plot(real.(aa), color="purple", label=L"$\mathfrak{R}(\langle a_1a_2 \rangle|(t))$");
plot(imag.(aa), color="purple", linestyle="--", label=L"$\mathfrak{I}(\langle a_1 a_2 \rangle)(t)$");
legend()
xlabel("timestep")
