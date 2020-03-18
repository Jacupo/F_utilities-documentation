include("/Users/jacoposurace/Documents/GitHub/F-utilities/F-utilities.jl")

N_A = 64;
H_A = Random_NNhamiltonian(N_A);
HA_D, UA_D = Diag_h(H_A);


N_B = 32;
H_B = Random_NNhamiltonian(N_A);

Γ_A = GS_Gamma(HA_D, UA_D);
Energy_fermions(Γ_A,HA_D, UA_D)
Γ_B = Thermal_fix_beta(Diag_h(H_B),1)
Γ_B2, β, Δ = Thermal_fix_energy(Diag_h(H_B),Energy(Γ_B,Diag_h(H_B)))
#Diagonalise the hamiltonian
H_D, U_D = Diag_h(H)
