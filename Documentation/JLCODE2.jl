include("F-utilities.jl")

N = 64;
H = Random_NNhamiltonian(N);
H_D, U_D = Diag_h(H);

Γ = GS_Gamma(H_D, U_D);
println("The energy of the ground state is: ", Energy_fermions(Γ,H_D, U_D));

N_A = 32;
#I consider the reduced state over the sites 17,2,...,48
Γ_A = Reduce_gamma(Γ,N_A,17);
#I compute the entangement entropy
S_A = VN_entropy(Γ_A);
#I compute the contour of partition A
c_A = Contour(Γ_A);

lbl_title   = string(L"$S(A) = $", S_A);
lbl_legend  = string(L"$\sum_{i=1}^{N_A} c_A(i) = $", sum(c_A));
figure("Contour of A")
title(lbl_title)
plot(1:N_A, c_A, marker="o", label=lbl_legend);
xlabel("i")
ylabel(L"$c_A(i)$")
legend();
