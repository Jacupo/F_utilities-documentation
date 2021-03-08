using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu    = F_utilities;
const LinA  = LinearAlgebra;


N           = 100;
H           = Fu.Build_hopping_hamiltonian(N);
HD, U       = Fu.Diag_h(H);
Gamma       = Fu.GS_gamma(HD,U);
Gamma_RBD   = Fu.RBD(Gamma,5)

prod_modes_border       = zeros(Int64, N);
prod_modes_RBD_border   = zeros(Int64, N);
prod_modes_bulk         = zeros(Int64, div(N,2));
prod_modes_RBD_bulk     = zeros(Int64, div(N,2));
for l=1:N
      DA,UA             = Fu.Diag_gamma(Fu.Reduce_gamma(Gamma,l,1));
      DA_RBD,UA_RBD     = Fu.Diag_gamma(Fu.Reduce_gamma(Gamma_RBD,l,1));
      prod_modes_border[l]     = count(i->i!=0, round.(real.(LinA.diag(DA)[1:l]),digits=14));
      prod_modes_RBD_border[l] = count(i->i!=0, round.(real.(LinA.diag(DA_RBD)[1:l]),digits=14));
end
for l=1:div(N,2)
      DA,UA                   = Fu.Diag_gamma(Fu.Reduce_gamma(Gamma,l, div(N,2)));
      DA_RBD,UA_RBD           = Fu.Diag_gamma(Fu.Reduce_gamma(Gamma_RBD,l, div(N,2)));
      prod_modes_bulk[l]      = count(i->i!=0, round.(real.(LinA.diag(DA)[1:l]),digits=14));
      prod_modes_RBD_bulk[l]  = count(i->i!=0, round.(real.(LinA.diag(DA_RBD)[1:l]),digits=14));
end
figure("prod_eigenvalues");
plot(1:N,prod_modes_border, marker="s", markersize=3, label=L"$\#(eigenval(\Gamma_{1,\dots,\ell}\neq 0,1)$");
plot(1:div(N,2),prod_modes_bulk, marker="s", markersize=3, label=L"$\#(eigenval(\Gamma_{\frac{N}{2},\dots,\frac{N}{2}+\ell})\neq 0,1)$");
plot(1:N,prod_modes_RBD_border, marker="s", markersize=3, label=L"$\#(eigenval(\Gamma(m=5)_{1,\dots,\ell}\neq 0,1)$");
plot(1:div(N,2),prod_modes_RBD_bulk, marker="s", markersize=3, label=L"$\#(eigenval(\Gamma(m=5)_{\frac{N}{2},\dots,\frac{N}{2}+\ell}\neq 0,1)$");
axvline(div(N,2), linestyle="--", linewidth=0.5, color="gray")
xlabel(L"$\ell$")
grid(axis="y", linestyle="--")
legend();
