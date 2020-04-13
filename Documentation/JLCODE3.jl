using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;

N = 128;

H   = Fu.Build_hopping_hamiltonian(N,true);
U_ω = Fu.Build_Fourier_matrix(N);

D_Fourier = U_ω'*H*U_ω;
D,U = Fu.Diag_h(H);

Γ_GS_Diag_h   = Fu.GS_gamma(D,U);

Γ_GS_Fourier  = zeros(Complex{Float64}, 2*N, 2*N);
for index=1:(2*N)
    if real(D_Fourier[index,index])<=0
        Γ_GS_Fourier[index,index] = 1;
    end
end
Γ_GS_Fourier = U_ω*Γ_GS_Fourier*U_ω';

figure("Energies")
plot(diag(real.(D_Fourier))[(N+1):(2*N)],label="Method Fourier");
plot(real.(diag(D))[(N+1):(2*N)], label="Method Diag_h");
xlabel(L"$k$");
ylabel(L"$\epsilon_k$");
legend()

figure("test", abs.(Γ_GS_Diag_h-Γ_GS_Fourier))
