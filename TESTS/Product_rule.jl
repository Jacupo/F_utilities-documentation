include("/Users/jacoposurace/Documents/GitHub/F-utilities/F-utilities.jl")

N_f         = 32;
θ           = π/8;

Jx                  = ones(Float64, N_f);
Jy                  = 0*ones(Float64, N_f);
lambda              = cot(θ)*ones(Float64, N_f);
H_ferm_i            = TFI_Hamiltonian(N_f, Jx, Jy, lambda);
H_ferm_f, U_diag_f  = Diag_h(H_ferm_i);
GS                  = GS_Gamma()
Γ1                  = Thermal_fix_beta((H_ferm_f,U_diag_f), 0.02);
Γ2                  = Thermal_fix_beta((H_ferm_f,U_diag_f), 0.1);

Γ3                  = Product(Γ1,Γ2);

Γ3v                 = Thermal_fix_beta((H_ferm_f,U_diag_f), 0.12)

println("------------------------------");
println("E1:", Energy(Γ1,(H_ferm_f,U_diag_f)));
println("E2:", Energy(Γ2,(H_ferm_f,U_diag_f)));
println("E3:", Energy(Γ3,(H_ferm_f,U_diag_f)));
println("E3v:", Energy(Γ3v,(H_ferm_f,U_diag_f)));
println("------------------------------");



# Print_matrix("Γ1", abs.(Γ1));
# Print_matrix("Γ2", abs.(Γ2));
# Print_matrix("Γ3", abs.(Γ3));
# Print_matrix("Γ3v", abs.(Γ3v));
# Print_matrix("Γ3-Γ3v", abs.(Γ3-Γ3v));
