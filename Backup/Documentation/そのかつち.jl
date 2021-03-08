using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;

function Mpermute(N,u,d)
    M = diagm(ones(Int64,N));
    M = permute(M,u,d);
    return M;
end
function Mpermute3x3()
    H = diagm(ones(Int64,6));
    H = permutehalfsx(H,1,6);
    H = permutehalfsx(H,2,5);
    H = permutehalfsx(H,4,6);
    H = permutehalfsx(H,3,4);
    # H = permutehalfsx(H,4,5);
    return H;
end


function Mpermute4x4()
    H = diagm(ones(Int64,8));
    H = permutehalfsx(H,1,6);
    H = permutehalfsx(H,2,4);
    H = permutehalfsx(H,1,2);
    H = permutehalfsx(H,5,8);
    # H = permutehalfsx(H,4,5);
    return H;
end

function permute(M,u,d)
    temp = M[u,:];
    M[u,:] = M[d,:];
    M[d,:] = temp;
    temp = M[:,u];
    M[:,u] = M[:,d];
    M[:,d] = temp;
    return M;
end
function permutehalfsx(M,u,d)
    temp = M[u,:];
    M[u,:] = M[d,:];
    M[d,:] = temp;
    # temp = M[:,u];
    # M[:,u] = M[:,d];
    # M[:,d] = temp;
    return M;
end
function permutehalfdx(M,u,d)
    # temp = M[u,:];
    # M[u,:] = M[d,:];
    # M[d,:] = temp;
    temp = M[:,u];
    M[:,u] = M[:,d];
    M[:,d] = temp;
    return M;
end
function permute3x3sx(H)
    H = permutehalfsx(H,1,6);
    H = permutehalfsx(H,2,5);
    H = permutehalfsx(H,4,6);
    H = permutehalfsx(H,3,4);
    return H;
end
function permute3x3dx(H)
    H = permutehalfdx(H,1,6);
    H = permutehalfdx(H,2,5);
    H = permutehalfdx(H,4,6);
    H = permutehalfdx(H,3,4);
    return H;
end

function permute3x3(H)
    H = permute(H,1,6);
    H = permute(H,2,5);
    H = permute(H,4,6);
    H = permute(H,3,4);
    return H;
end

# for N=8:8
# N = 8;
# #Con 4 Fourier è sbagliato, con 16 Diag è sbagliato
#
# H   = Fu.Build_hopping_hamiltonian(N,true);
# U_ω = round.(Fu.Build_Fourier_matrix(N),digits=15);
#
# # U_ω = permute3x3(U_ω);
#
# D_Fourier = U_ω'*H*U_ω;
# D,U = Fu.Diag_h(H);
# # Mp = (Mpermute4x4())';
# # Up = U#*Mp;
# # Dp = Up'*H*Up;
# Γ_GS_Fourier = Fu.GS_gamma(D_Fourier,U_ω)
# Γ_GS_Diag_h   = Fu.GS_gamma(D,U)
#
# for i=1:N
#     println("D:", diag(Γ_GS_Diag_h)[i]+diag(Γ_GS_Diag_h)[i+N], "- F:", diag(Γ_GS_Fourier)[i]+diag(Γ_GS_Fourier)[i+N])
# end
# # if (sum(abs.(Γ_GS_Diag_h-Γ_GS_Fourier))>0.00000001)
# #     println(N, " --> ", sum(abs.(Γ_GS_Diag_h-Γ_GS_Fourier)))
# # end
# println("--------")
# println("A: ",sum(sort(diag(real.(D)))[1:N]))
# println("B: ",sum(sort(diag(real.(D_Fourier)))[1:N]))
# println("En Diag",Fu.Energy(Γ_GS_Diag_h,(D,U)))
# println("En Fu",Fu.Energy(Γ_GS_Fourier,(D_Fourier,U_ω)))
#
# Fu.Print_matrix("D", abs.(U*U'*H*U*U'-U*D*U'))
# # Fu.Print_matrix("Dp",abs.(Dp))
# Fu.Print_matrix("Fou",real.(Γ_GS_Fourier))#Da qua vedo che questo è il sbagliato perchè aad!=1-ada
# Fu.Print_matrix("Diag",real.(Γ_GS_Diag_h))
# Fu.Print_matrix("Diff",real.(Γ_GS_Diag_h-Γ_GS_Fourier))

###################################################
# D = permute(D,1,6);
# D = permute(D,2,4);
# D = permute(D,1,2);
# D = permute(D,5,8);
#
# U = permute(U,1,6);
# U = permute(U,2,4);
# U = permute(U,1,2);
# U = permute(U,5,8);
#######################




# H = permute3x3(H);
# D = permute3x3(D);
# U = permute3x3(U);
# Mp = (Mpermute4x4())';
# Up = U*Mp;
# VV = U_ω*Up';
# D = Up'*H*Up;

# G = Up*U_ω';


#
# Fu.Print_matrix("test 10", real.(Up*Up'));
# Fu.Print_matrix("test 11", real.(U_ω));
# Fu.Print_matrix("test 12", real.(Up-U_ω));
# Fu.Print_matrix("test 13", real.(D-D_Fourier));
# Fu.Print_matrix("test 14", real.(VV*VV'));
#
# # plot(sort(diag(real.(D))))
# # plot(sort(diag(real.(D_Fourier))))
# # Fu.Print_matrix("H",real.(H))
# # Fu.Print_matrix("test back-fourth", real.(U_ω*D_Fourier*U_ω'-H))
# Fu.Print_matrix("test back-fourth", real.(D))
# # Fu.Print_matrix("test back-fourth", real.(Up))
#

#Il problema è che i ground state non sono uguali calcolati coi due modi

# end

# UHU'=VHV'
# quale relazione c'è tra U e V?

# U=V*X con X unitaria misteriosa
# X = U*inv(V)=U*V';

# G = V'U
# G*H*G' = H con G unitaria-> G e H commutano
#
#
# figure("Energies")
# plot(diag(real.(D_Fourier))[(N+1):(2*N)],label="Method Fourier");
# plot(real.(diag(D))[(N+1):(2*N)], label="Method Diag_h");
# xlabel(L"$k$");
# ylabel(L"$\epsilon_k$");
# legend();
#
# Fu.Print_matrix("test Four", real.(Γ_GS_Fourier));
# Fu.Print_matrix("test Diag", real.(Γ_GS_Diag_h));
# Fu.Print_matrix("test Diag", real.(Γ_GS_Diag_h-Γ_GS_Fourier));
#
# println("The energies are the same: ", iszero(round(Fu.Energy(Γ_GS_Diag_h,(D,U))-Fu.Energy(Γ_GS_Fourier,(D_Fourier,U_ω)),digits=10)));
# println("The ground state are the same: ", (round(sum(abs.(Γ_GS_Diag_h-Γ_GS_Fourier)),digits=10)));





N =9;

H   = Fu.Build_hopping_hamiltonian(N,true);

U_ω = Fu.Build_Fourier_matrix(N);
D_ω = U_ω'*H*U_ω;
D,U = Fu.Diag_h(H);
# Fu.Print_Complex_matrix("Dω",D_ω);
# Fu.Print_Complex_matrix("D",D);
Fu.Print_Complex_matrix("Uω",U_ω);
Fu.Print_Complex_matrix("U",U);

# Mp = (Mpermute4x4())';
# Up = U#*Mp;
# Dp = Up'*H*Up;
Γ_ω = Fu.GS_gamma(D_ω,U_ω);
Γ   = Fu.GS_gamma(D,U);
println("Energie----")
println("En ω:      ",Fu.Energy(Γ_ω,(D_ω,U_ω)))
println("En Diag:   ",Fu.Energy(Γ,(D,U)))
Fu.Print_Complex_matrix("Γω",Γ_ω);
Fu.Print_Complex_matrix("Γ",Γ);
Fu.Print_Complex_matrix("Differenza Γ", Γ-Γ_ω)
Fu.Print_Complex_matrix("Differenza Γ UL-BR", Γ[1:N,1:N]+Γ[(1:N).+N,(1:N).+N])
Fu.Print_Complex_matrix("Differenza Γω UL-BR", Γ_ω[1:N,1:N]+Γ_ω[(1:N).+N,(1:N).+N])
fig = plt.figure("Diag Γ",figsize=(5, 5), dpi=80)
plot(real.(diag(Γ)))
fig2 = plt.figure("Diag Γω",figsize=(5, 5), dpi=80)
plot(real.(diag(Γ_ω)))
# println("--------")
# println("A: ",sum(sort(diag(real.(D)))[1:N]))
# println("B: ",sum(sort(diag(real.(D_Fourier)))[1:N]))
# println("En Diag",Fu.Energy(Γ_GS_Diag_h,(D,U)))
# println("En Fu",Fu.Energy(Γ_GS_Fourier,(D_Fourier,U_ω)))
#
# Fu.Print_matrix("D", abs.(U*U'*H*U*U'-U*D*U'))
# # Fu.Print_matrix("Dp",abs.(Dp))
# Fu.Print_matrix("Fou",real.(Γ_GS_Fourier))#Da qua vedo che questo è il sbagliato perchè aad!=1-ada
# Fu.Print_matrix("Diag",real.(Γ_GS_Diag_h))
# Fu.Print_matrix("Diff",real.(Γ_GS_Diag_h-Γ_GS_Fourier))
