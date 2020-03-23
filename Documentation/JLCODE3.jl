include("F-utilities.jl")

N = 128;

H   = Build_hopping_hamiltonian(N,true);
U_ω = Build_Fourier_matrix(N);

D_Fourier = U_ω'*H*U_ω;
D,U = Diag_h(H);


figure("Energies")
plot(diag(real.(D_Fourier))[(N+1):(2*N)],label="Method Fourier");
plot(real.(diag(D))[(N+1):(2*N)], label="Method Diag_h");
xlabel(L"$k$");
ylabel(L"$\epsilon_k$");
legend()
