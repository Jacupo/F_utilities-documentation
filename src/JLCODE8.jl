using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu    = F_utilities;
const LinA  = LinearAlgebra;

function Random_hamiltonian(N)
  A   = rand(N,N)+im*rand(N,N);
  A   = (A+A')/2.;
  bd  = rand(N-1).+im*rand(N-1);
  B   = Tridiagonal(bd, zeros(Complex{Float64}, N), -bd);
  H   = zeros(Complex{Float64}, 2*N, 2*N);
  H[(1:N),(1:N)]          = -conj(A);
  H[(1:N).+N,(1:N)]       = -conj(B);
  H[(1:N),(1:N).+N]       = B;
  H[(1:N).+N,(1:N).+N]    = A;

  return H;
end

N           = 100;
H           = Random_hamiltonian(N);
HD, U       = Fu.Diag_h(H);
Gamma       = Fu.GS_gamma(HD,U);
Gamma_RBD   = Fu.RBD(Gamma,1)
S_modes_border       = zeros(Float64, N);
S_modes_RBD_border   = zeros(Float64, N);
S_modes_bulk         = zeros(Float64, div(N,2));
S_modes_RBD_bulk     = zeros(Float64, div(N,2));
for l=1:N
      S_modes_border[l]     = Fu.VN_entropy(Fu.Reduce_gamma(Gamma,l,1));
      S_modes_RBD_border[l] = Fu.VN_entropy(Fu.Reduce_gamma(Gamma_RBD,l,1));
end
for l=1:div(N,2)
      S_modes_bulk[l]         = Fu.VN_entropy(Fu.Reduce_gamma(Gamma,l, div(N,2)));
      S_modes_RBD_bulk[l]     = Fu.VN_entropy(Fu.Reduce_gamma(Gamma_RBD,l, div(N,2)));
end
figure("Entropies");
plot(1:N,log.(abs.(S_modes_border)), marker="s", markersize=3, label=L"$F=S(\Gamma_{1,\dots,\ell}\neq 0,1)$");
plot(1:div(N,2),log.(abs.(S_modes_bulk)), marker="s", markersize=3, label=L"$F=S(\Gamma_{\frac{N}{2},\dots,\ell+\frac{N}{2}})\neq 0,1)$");
plot(1:N,log.(abs.(S_modes_RBD_border)), marker="s", markersize=3, label=L"$F=S(\Gamma(m=1)_{1,\dots,\ell}\neq 0,1)$");
plot(1:div(N,2),log.(abs.(S_modes_RBD_bulk)), marker="s", markersize=3, label=L"$F=S(\Gamma(m=1)_{\frac{N}{2},\dots,\ell+\frac{N}{2}}\neq 0,1)$");
axvline(div(N,2), linestyle="--", linewidth=0.5, color="gray")
axhline(log(log(2)), linestyle="-.", color="red", label="F=log(D)")
axhline(log(2*log(2)), linestyle="-.", color="red", label="F=2*log(D)")
xlabel(L"$\ell$")
ylabel(L"$\log(F)$")
# grid(axis="y", linestyle="--")
legend();
