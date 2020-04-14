using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;


N =127;

H   = Fu.Build_hopping_hamiltonian(N,true);

U_ω = Fu.Build_Fourier_matrix(N);
D_ω = U_ω'*H*U_ω;
D,U = Fu.Diag_h(H);


figure("Energies")
plot(diag(real.(D_ω))[(N+1):(2*N)],label="Method Fourier");
plot(real.(diag(D))[(N+1):(2*N)], label="Method Diag_h");
xlabel(L"$k$");
ylabel(L"$\epsilon_k$");
legend();

Γ_ω = Fu.GS_gamma(D_ω,U_ω);
Γ   = Fu.GS_gamma(D,U);
println("")
println("Energy GS Method Fourier:      ",Fu.Energy(Γ_ω,(D_ω,U_ω)))
println("En GS Method Diag_h:           ",Fu.Energy(Γ,(D,U)))

Fu.Print_complex_matrix("Differenza Γ", Γ-Γ_ω)
