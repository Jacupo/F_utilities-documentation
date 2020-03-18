D_A = size(Γ_A,1);
D_B = size(Γ_B,1);
D   = D_A+D_B;

Γ_AB = zeros(Complex{Float64}, D,D);
Γ_AB = Inject_gamma(Γ_AB,Γ_A,1);
Γ_AB = Inject_gamma(Γ_AB,Γ_B,D_A+1);
