D_A = size(Gamma_A,1);
D_B = size(Gamma_B,1);
D   = D_A+D_B;

Gamma_AB = zeros(Complex{Float64}, D,D);
Gamma_AB = Inject_gamma(Gamma_AB,Gamma_A,1);
Gamma_AB = Inject_gamma(Gamma_AB,Gamma_B,D_A+1);
