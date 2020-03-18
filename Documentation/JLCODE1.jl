include("F-utilities.jl")

N = 64;

#Generate and diagonalise the hamiltonian
H = Random_NNhamiltonian(N)
H_D, U_D = Diag_h(H)

#Print the energy modes ϵ_k
figure("Energies")
plot(1:N,diag(H_D)[1:N])
xlabel(L"$k$")
ylabel(L"$\epsilon_k$")
