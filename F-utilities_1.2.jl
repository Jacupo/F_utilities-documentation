using PyPlot; #Or import?
using LinearAlgebra;

function Print_matrix(title, matrix)
 figure(title)
 pcolormesh(matrix)
 colorbar()
 ylim(size(matrix,1),0)
end

function Build_Omega(N)
#Build the matrix omega of dimension 2*N, that is for N sites.
  eye                 = Matrix{Float64}(I, N, N);
  Ω                   = zeros(Complex{Float64}, 2*N, 2*N);
  Ω[1:N,1:N]          = (1/(sqrt(2)))*eye;
  Ω[1:N,(1:N).+N]      = (1/(sqrt(2)))*eye;
  Ω[(1:N).+N,1:N]      = im*(1/(sqrt(2)))*eye;
  Ω[(1:N).+N,(1:N).+N]  = -im*(1/(sqrt(2)))*eye;
  return Ω;
end

function Build_FxxTxp(N)
 FxxTxp = zeros(Int64, 2*N, 2*N);
 for iiter=1:N
   FxxTxp[2*iiter-1,iiter]   = 1;
   FxxTxp[2*iiter, iiter+N] = 1;
 end
 return FxxTxp;
end

function Build_FxpTxx(N)
 FxpTxx = zeros(Int64, 2*N, 2*N)
 for iiter=1:N
   FxpTxx[iiter,2*iiter-1]    = 1;
   FxpTxx[iiter+N, 2*iiter] = 1;
 end
 return FxpTxx;
end

function Diag_real_skew(M, rand_perturbation::Int64=0)
  #rand_perturbation: 0 do nothing,
  #                   1 add random eps() perturbation to the diagonal,
  #                   2 add random eps() perturbation to the whole matrix.
  #                   Default: default_perturbation=0 -> do nothing.

 N = div(size(M,1),2);

 M = real((M-M')/2.); #Force skew-symmetry
 if (rand_perturbation != 0)
   if (rand_perturbation == 1)
     M += diagm(rand(2*N_f)*eps()) #Se la diagonale è proprio zero ho problemi di convergenza
   end
   if (rand_perturbation == 2)
     random_M = 1*rand(2*N_f,2*N_f)*eps();#Per SETTORI MOMENTO questo bastava 10 e non 100
     random_M = (random_M-random_M')/2.;
     M += random_M;
   end
 end
 Schur_object   = LinearAlgebra.schur(M);
 #Per non avere errori di convergenza
 # Schur_object   = schurfact(M);
 #Forse questo risolve dei problemi numerici
 Schur_ort_i    = Schur_object.vectors;
 Schur_blocks_i = Schur_object.Schur;



 Schur_adjust = zeros(Int64, 2*N, 2*N);
 for iiter=1:N
   if (Schur_blocks_i[2*iiter,2*iiter-1] >= 0.)
     Schur_adjust[2*iiter-1, 2*iiter] = 1;
     Schur_adjust[2*iiter, 2*iiter-1] = 1;
   else
     Schur_adjust[2*iiter-1,2*iiter-1] = 1;
     Schur_adjust[2*iiter, 2*iiter]    = 1;
   end
 end

 M_temp   = Schur_adjust*Schur_blocks_i*Schur_adjust';
 O_temp   = (Schur_ort_i*Schur_adjust);


 not_sorted = true;
 while not_sorted
     not_sorted = false;
     for jiter=2:(N)
       if (M_temp[2*(jiter-1)-1, 2*(jiter-1)] <= M_temp[2*jiter-1, 2*jiter])
         not_sorted = true;
         temp = M_temp[2*(jiter-1)-1, 2*(jiter-1)];
         M_temp[2*(jiter-1)-1, 2*(jiter-1)] = M_temp[2*jiter-1, 2*jiter];
         M_temp[2*(jiter-1), 2*(jiter-1)-1] = -M_temp[2*(jiter-1)-1, 2*(jiter-1)];
         M_temp[2*jiter-1, 2*jiter] = temp;
         M_temp[2*jiter, 2*jiter-1] = -temp;
         temp1 = O_temp[:,2*(jiter-1)];
         O_temp[:,2*(jiter-1)] = O_temp[:,2*jiter];
         O_temp[:,2*jiter] = temp1;
         temp1 = O_temp[:,2*(jiter-1)-1];
         O_temp[:,2*(jiter-1)-1] = O_temp[:,2*jiter-1];
         O_temp[:,2*jiter-1] = temp1;
       end
     end
     N = N-1;
 end

   M_f = M_temp;
   O_f = real.(O_temp);

 return M_f, O_f;
end

function Diag_ferm(M,rand_perturbation::Int64=0)
    N = div(size(M,1),2);

    F_xptxx = Build_FxpTxx(N);

    Ω  = Build_Omega(N);
    M_temp = real(-im*Ω*M*Ω')
    # M_temp += 0.01*rand(2*N_f,2*N_f);
    M_temp = (M_temp-M_temp')/2.;
    #M_temp[1,1] += eps();
    M_temp, O = Diag_real_skew(M_temp,rand_perturbation)
    M_temp = F_xptxx*M_temp*(F_xptxx')
    M_temp = im*Ω'*M_temp*(Ω)

    M_f = M_temp;
    U_f = Ω'*O*(F_xptxx')*Ω;

    return real(M_f), U_f;
end

function Diag_gamma(Γ)
 #Diagonalizzo cosi' invece che con eig cosi' ottengo gli autovettori a coppie
 #ordinate come voglio.
 #La generica unitaria tra Fermioni di Dirac e' fatta come
 #Ω*O*Ω', con O una generica ortogonale NxN specificata da N(N-1)/2 angoli
 Γ = (Γ+Γ')/2.;
 γ, U = Diag_ferm(Γ-0.5*I);

 return U'*Γ*U,U;#real(γ+0.5*eye(size(Γ,1))),U
end

function Evolve_gamma(M,D,U,t)
   N = div(size(M,1),2);

   M = (M+M')/2.;

   M_diag_base = U'*M*U;
   M_diag_base_evolv = exp(im*D*t)*M_diag_base*exp(-im*D*t);
   M_evolv = U*M_diag_base_evolv*(U');

   M_evolv = (M_evolv+M_evolv')/2.

   return M_evolv;
end

function Build_A_NOPBC(Jx, Jy, lambdas)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1   ##SEMBRA CHE NON DEVO
   M_A       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_A[iiter, iiter]   = -2*lambdas[iiter];
       M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);
       M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);
   end
   M_A[1,1]         = -2*lambdas[1];#
   M_A[1,2]         = -(Jx[1]+Jy[1]);#

   # M_A[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#
   # M_A[dimension,1] = (-(-1)^(parity))*-(Jx[1]+Jy[1]);
   M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);#
   M_A[dimension,dimension]    = -2*lambdas[dimension];#

   return M_A
end

function Build_A(Jx, Jy, lambdas)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1   ##SEMBRA CHE NON DEVO
   M_A       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_A[iiter, iiter]   = -2*lambdas[iiter];
       M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);   #Si connettono sempre con il coeff. del più piccolo
       M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);       #Si connettono sempre con il coeff. del più piccolo
   end
   M_A[1,1]         = -2*lambdas[1];#
   M_A[1,2]         = -(Jx[1]+Jy[1]);#

   M_A[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#
   M_A[dimension,1] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#Questo è l'unico caso dove connetto con il coeff del più grande
   M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);   #
   M_A[dimension,dimension]    = -2*lambdas[dimension];#

   return M_A
end


function Build_A_APBC(Jx, Jy, lambdas)
   dimension = size(Jx, 1);
   parity = mod(dimension+1, 2); #!!ATTENZIONE HO CAMBIATO con +1   ##SEMBRA CHE NON DEVO
   M_A       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_A[iiter, iiter]   = -2*lambdas[iiter];
       M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);
       M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);
   end
   M_A[1,1]         = -2*lambdas[1];#
   M_A[1,2]         = -(Jx[1]+Jy[1]);#

   M_A[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#
   M_A[dimension,1] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);
   M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);#
   M_A[dimension,dimension]    = -2*lambdas[dimension];#

   return M_A
end

function Build_B_NOPBC(Jx, Jy)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1 ##SEMBRA CHE NON DEVO
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
       M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
   end
   M_B[1,2]         = Jx[1]-Jy[1];

   # M_B[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]-Jy[dimension]);#
   # M_B[dimension,1] = (-(-1)^(parity))*(Jx[1]-Jy[1]);#

   M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);#

   return M_B
end


function Build_B(Jx, Jy)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1 ##SEMBRA CHE NON DEVO
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
       M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
   end
   M_B[1,2]         = Jx[1]-Jy[1];

   M_B[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]-Jy[dimension]);#
   M_B[dimension,1] = (-(-1)^(parity))*(Jx[dimension]-Jy[dimension]);#

   M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);#

   return M_B
end

function Build_B_APBC(Jx, Jy)
   dimension = size(Jx, 1);
   parity = mod(dimension+1, 2); #!!ATTENZIONE HO CAMBIATO con +1 ##SEMBRA CHE NON DEVO
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
       M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
   end
   M_B[1,2]         = Jx[1]-Jy[1];

   M_B[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]-Jy[dimension]);#
   M_B[dimension,1] = (-(-1)^(parity))*(Jx[dimension]-Jy[dimension]);#

   M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);#

   return M_B
end

function TFI_Hamiltonian(N,Jx,Jy,lambda)
   A = 0.5*Build_A(Jx,Jy,lambda);
   B = 0.5*Build_B(Jx,Jy);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end

function TFI_Hamiltonian_NOPBC(N,Jx,Jy,lambda)
   A = 0.5*Build_A_NOPBC(Jx,Jy,lambda);
   B = 0.5*Build_B_NOPBC(Jx,Jy);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  +B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end




function TFI_Hamiltonian_APBC(N,Jx,Jy,lambda)
   A = 0.5*Build_A_APBC(Jx,Jy,lambda);
   B = 0.5*Build_B_APBC(Jx,Jy);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end

function GS_Gamma(D,U)
   N = div(size(D,1),2);

   Gamma_diag_base = zeros(Complex{Float64}, 2*N, 2*N);
   for iiter=1:N
       Gamma_diag_base[iiter+N, iiter+N] = 1;
   end
   Gamma = U*Gamma_diag_base*U';

   Gamma = (Gamma+(Gamma'))/2.
   return Gamma;
end


function Reduce_gamma(M, N_partition, first_index)
   N_f = div(size(M,1),2);
   first_index = first_index-1;
   periodic_dimension = max((N_partition.+first_index-N_f),0)
   dim_UL = N_partition-periodic_dimension;

   redgamma = zeros(Complex{Float64}, N_partition*2, N_partition*2);
   #Copy the upper left left part of the correlation matrix
   redgamma[1:dim_UL,1:dim_UL] = M[(1:dim_UL).+first_index,(1:dim_UL).+first_index];
   redgamma[(1:dim_UL).+N_partition,1:dim_UL] = M[(1:dim_UL).+N_f.+first_index,(1:dim_UL).+first_index];
   redgamma[1:dim_UL,(1:dim_UL).+N_partition] = M[(1:dim_UL).+first_index,(1:dim_UL).+N_f.+first_index];
   redgamma[(1:dim_UL).+N_partition,(1:dim_UL).+N_partition] = M[(1:dim_UL).+N_f.+first_index,(1:dim_UL).+N_f.+first_index];

   if (periodic_dimension>0)
     redgamma[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension))] = M[1:periodic_dimension,1:periodic_dimension];
     redgamma[1:dim_UL,(dim_UL.+(1:periodic_dimension))] = M[(first_index.+(1:dim_UL)),1:periodic_dimension];
     redgamma[(dim_UL.+(1:periodic_dimension)),1:dim_UL] = M[1:periodic_dimension,(first_index.+(1:dim_UL))];

     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,(dim_UL.+(1:periodic_dimension))] = M[(1:periodic_dimension).+N_f,1:periodic_dimension];
     redgamma[(1:dim_UL).+N_partition,(dim_UL.+(1:periodic_dimension))] = M[(first_index.+(1:dim_UL)).+N_f,1:periodic_dimension];
     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,1:dim_UL] = M[(1:periodic_dimension).+N_f,(first_index.+(1:dim_UL))];

     redgamma[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension)).+N_partition] = M[1:periodic_dimension,(1:periodic_dimension).+N_f];
     redgamma[1:dim_UL,(dim_UL.+(1:periodic_dimension)).+N_partition] = M[(first_index.+(1:dim_UL)),(1:periodic_dimension).+N_f];
     redgamma[(dim_UL.+(1:periodic_dimension)),(1:dim_UL).+N_partition] = M[1:periodic_dimension,(first_index.+(1:dim_UL)).+N_f];

     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,(dim_UL.+(1:periodic_dimension)).+N_partition] = M[(1:periodic_dimension).+N_f,(1:periodic_dimension).+N_f];
     redgamma[(1:dim_UL).+N_partition,(dim_UL.+(1:periodic_dimension)).+N_partition] = M[(first_index.+(1:dim_UL)).+N_f,(1:periodic_dimension).+N_f];
     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,(1:dim_UL).+N_partition] = M[(1:periodic_dimension).+N_f,(first_index.+(1:dim_UL)).+N_f];
   end


   return redgamma
end

function Inject_gamma(gamma, injection, first_index)
 dim_gamma     = div(size(gamma, 1),2);
 dim_injection = div(size(injection, 1), 2);

 first_index = first_index-1;
 periodic_dimension = max((dim_injection+first_index-dim_gamma),0)
 dim_UL = dim_injection-periodic_dimension;

 #Injecto la parte Z nei 4 riquadri
 gamma[(1:dim_UL).+first_index,(1:dim_UL).+first_index]                      = injection[(1:dim_UL),(1:dim_UL)]
 gamma[(1:dim_UL).+(first_index+dim_gamma),(1:dim_UL).+first_index]            = injection[(1:dim_UL).+dim_injection,(1:dim_UL)]
 gamma[(1:dim_UL).+first_index,(1:dim_UL).+first_index.+dim_gamma]            = injection[(1:dim_UL),(1:dim_UL).+dim_injection]
 gamma[(1:dim_UL).+first_index.+dim_gamma,(1:dim_UL).+first_index.+dim_gamma]  = injection[(1:dim_UL).+dim_injection,(1:dim_UL).+dim_injection]


 if (periodic_dimension>0)
   #Injecto A,B,C  per ogni riquadro
   gamma[1:periodic_dimension,1:periodic_dimension] = injection[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension))];
   gamma[(first_index.+(1:dim_UL)),1:periodic_dimension] = injection[1:dim_UL,(dim_UL.+(1:periodic_dimension))];
   gamma[1:periodic_dimension,(first_index.+(1:dim_UL))] = injection[(dim_UL.+(1:periodic_dimension)),1:dim_UL];

   gamma[(1:periodic_dimension).+dim_gamma,1:periodic_dimension] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,(dim_UL.+(1:periodic_dimension))];
   gamma[(first_index.+(1:dim_UL)).+dim_gamma,1:periodic_dimension] = injection[(1:dim_UL).+dim_injection,(dim_UL.+(1:periodic_dimension))];
   gamma[(1:periodic_dimension).+dim_gamma,(first_index.+(1:dim_UL))] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,1:dim_UL];

   gamma[1:periodic_dimension,(1:periodic_dimension).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[(first_index.+(1:dim_UL)),(1:periodic_dimension).+dim_gamma] = injection[1:dim_UL,(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[1:periodic_dimension,(first_index.+(1:dim_UL)).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)),(1:dim_UL).+dim_injection];

   gamma[(1:periodic_dimension).+dim_gamma,(1:periodic_dimension).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[(first_index.+(1:dim_UL)).+dim_gamma,(1:periodic_dimension).+dim_gamma] = injection[(1:dim_UL).+dim_injection,(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[(1:periodic_dimension).+dim_gamma,(first_index.+(1:dim_UL)).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,(1:dim_UL).+dim_injection];
 end

 return gamma
end

function VN_entropy(M)
   N = convert(Int64, size(M,1));

   D,U = LinearAlgebra.eigen((M+M')/2.);   #Se voglio fare con eig
   D = diagm(real(D));

   #Forzo la Hermitianità, se ho risultati strani controlla QUA
   #che M sia abbastanza vicina all'hermitiana prima.
   # M   = (M+M')/2.
   # D,U = Diag_gamma(M);
   S = 0;

   for iiter=1:N
       if (round.(D[iiter,iiter];digits=18)<-0.0000000000001)
        De,Ue = LinerAlgebra.eigen((M+M')/2.);
        for iiter=1:N
         println("DG: ", D[iiter,iiter]);
        end
        for iiter=1:N
         println("DE: ", De[iiter]);
        end
        error = string("Eigenvalue in VE not in [0,1]: ",round(D[iiter,iiter],18))
        throw(ArgumentError(error))
       end
       nu = abs(round.(D[iiter,iiter];digits=18))
       if (nu != 0 && nu != 1)
        #Invece di arrivare fino a N/2 nel ciclo e sommare nu e 1-nu, li passo tutti
        #perchè potrebbe essere che facendo una generica transformazione ortogonale non li abbia
        #ordinati in coppie, anche se Diag_gamma dovrebbe metterli in ordine
           S -= log(nu)*nu;
       end
   end

   return S;
end


function Biggest_eig_H(M,n_eig)
    #La grandezza del vettore restituito scala come 2^(n_gaps-1);
    N = convert(Int64, size(M,1)/2);

    D,U = Diag_ferm(M)

    # d_rev = sort(real(diag(D)),rev=true); #Il primo è il più grande, sono positivi
    d_norev = sort(real(diag(D)));          #Il primo è il più piccolo, sono negativi

    # println(sum(d_norev[1:N]))

    initial_coeff = 0;
    for i=1:(N-n_eig)
        initial_coeff += d_norev[i]
    end
    evor = initial_coeff*ones(Float64,convert(Int64,2^(n_eig-1))+1);

    for i=0:convert(Int64,2^(n_eig-1))
        bin_string = i;
        for j=0:(n_eig-1)
            evor[i+1] += (((-1)^(mod(bin_string,2)))*d_norev[N-j]);
            bin_string = div(bin_string,2);
        end
    end

    return sort(evor,rev=true);
end


function Evolve_gamma_imag(M,D,U,t)
   N_f = convert(Int64, size(M,1)/2.);

   #####AGGIUNTO
   M = (M+M')/2.
   # M[1,1] += eps();
   ####

   M_diag_base = U'*M*U;
   A           = 2*(D*M_diag_base-M_diag_base*D)
   U_t         = exp(A*t)
   M_diag_base_evolv = U_t*M_diag_base*U_t';
   M_evolv             = U*M_diag_base_evolv*U';

   #####AGGIUNTO
   M_evolv = (M_evolv+M_evolv')/2.
   ####

   return M_evolv;
end




function Momentum_sector_H_pi_flux(Lp,Lo,kx)
 H  = zeros(Float64,4*Lo,4*Lo);
 A  = zeros(Float64, 2*Lo, 2*Lo);

 for i=2:(Lo-1)
   A[i,i]       = cos(2*pi*kx/Lp);
   A[i+Lo,i+Lo] = -cos(2*pi*kx/Lp);
   A[i+Lo,i+1]  =  1/2.;
   A[i+Lo,i-1]  =  1/2.;
   A[i,i+1+Lo]  =  1/2.;
   A[i,i-1+Lo]  =  1/2.;
 end
 A[1,1]         = cos(2*pi*kx/Lp);
 A[Lo,Lo]       = cos(2*pi*kx/Lp);
 A[1+Lo,1+Lo]   = -cos(2*pi*kx/Lp);
 A[Lo+Lo,Lo+Lo] = -cos(2*pi*kx/Lp);
 A[1,2+Lo]      = 1/2.;
 A[Lo,Lo-1+Lo]  = 1/2.;
 A[1+Lo,2]      = 1/2.;
 A[Lo+Lo,Lo-1]  = 1/2.;

 H[1:(2*Lo),1:(2*Lo)]                = A;
 H[(2*Lo+1):(4*Lo),(2*Lo+1):(4*Lo)]  = -A;

 return H;
end

function H_fasulla(Lp,Lo,kx)
 H  = zeros(Float64,4*Lo,4*Lo);
 A  = zeros(Float64, 2*Lo, 2*Lo);

 for i=2:(Lo-1)
   A[i,i]       = rand();
   A[i+Lo,i+Lo] = rand();
   A[i+Lo,i+1]  =  mod(i,kx);
   A[i+Lo,i-1]  =  mod(i,kx);
   A[i,i+1+Lo]  =  mod(i,kx);
   A[i,i-1+Lo]  =  mod(i,kx);
 end
 A[1,1]         = rand();
 A[Lo,Lo]       = rand();
 A[1+Lo,1+Lo]   = -rand();
 A[Lo+Lo,Lo+Lo] = -rand();
 A[1,2+Lo]      = mod(Lo,kx);
 A[Lo,Lo-1+Lo]  = mod(Lo,kx);
 A[1+Lo,2]      = mod(Lo,kx);
 A[Lo+Lo,Lo-1]  = mod(Lo,kx);

 H[1:(2*Lo),1:(2*Lo)]                = A;
 H[(2*Lo+1):(4*Lo),(2*Lo+1):(4*Lo)]  = -A;

 return H;
end



function H_fasulla_2(Lp,Lo,kx)
 H  = zeros(Float64,4*Lo,4*Lo);
 A  = zeros(Float64, 2*Lo, 2*Lo);

 for i=2:(Lo-1)
   A[i,i]       = cos(2*pi*kx/Lp);
   A[i+Lo,i+Lo] = -cos(2*pi*kx/Lp);
   A[i+Lo,i+1]  =  1/2.;
   A[i+Lo,i-1]  =  1/2.;
   A[i,i+1+Lo]  =  1/2.;
   A[i,i-1+Lo]  =  1/2.;
 end
 A[1,1]         = cos(2*pi*kx/Lp);
 A[Lo,Lo]       = cos(2*pi*kx/Lp);
 A[1+Lo,1+Lo]   = -cos(2*pi*kx/Lp);
 A[Lo+Lo,Lo+Lo] = -cos(2*pi*kx/Lp);
 A[1,2+Lo]      = 1/3.;
 A[Lo,Lo-1+Lo]  = 1/3.;
 A[1+Lo,2]      = 1/3.;
 A[Lo+Lo,Lo-1]  = 1/3.;

 H[1:(2*Lo),1:(2*Lo)]                = A;
 H[(2*Lo+1):(4*Lo),(2*Lo+1):(4*Lo)]  = -A;

 return H;
end


function build_reorder_sector(L2)
 #Questa O mi trasforma da a†k1,a†k2,..,a†-k1,a†-k2,...,ak1,...,a-k1,....
 #all ordine a†k1,a†-k1,a†k2,a†-k2,.....ak1,a-k1,..

#TODO:
 #Applicare questo alla H e fare già tutto senza questo fin dall'inizio? cioè usare solamente una notazione diversa?
 O = zeros(Int64, 4*L2, 4*L2)

 for iiter=1:(2*L2)
     O[iiter,convert(Int64, ceil(iiter/2)+mod(iiter-1,2)*L2)] = 1;
 end
 for iiter=1:(2*L2)
     O[iiter+(2*L2),convert(Int64, ceil(iiter/2)+mod(iiter-1,2)*L2)+(2*L2)] = 1;
 end

 return O;
end


function diagonalise_block(γi, starting_site, dimension)
 #Questa funzione prende una matrice γ e ritorna
 #la matrice γ_finale con il blocco grande "dimension" che parte dal
 #sito "starting_site" diagonalizzato.
 #Inoltre ritorna U_tot tale che
 #γ =  U_tot*γ_finale*U_tot'

 γ = deepcopy(γi)
 γ = (γ+γ')/2.
 dimension = convert(Int64, dimension)
 if ((starting_site+dimension-1)>size(γ,1)/2)
   println("_diagonalise_block_");
   println("ATTENZIONE IL BLOCCO SFORA IL BORDO");
   println("SS ", starting_site, " d ", dimension, " size(γ,1) ", size(γ,1));
 end
 γ_reduced    = Reduce_gamma(γ,dimension,starting_site);
 D_red, U_red = Diag_gamma(γ_reduced);
 ident        = convert(Array{Complex{Float64},2}, Matrix{Float64}(I, size(γ,1), size(γ,1)));
 U_tot        = Inject_gamma(ident, U_red, starting_site);
 γ_finale     = U_tot'*γ*U_tot;
 γ_finale     = (γ_finale+γ_finale')/2.



 return γ_finale,U_tot
end


function Reduce_bond_dimension(Γ_in,L1,L2,N_sectors,mdL)
 #Qua con m intendo m/L1
 #If the bond dimension is maximal do nothing
 if (mdL>=L2)
  println("L2=", L2, " m/L=", mdL, " impossible to RBD")
  return Γ_in;
 end

 Γ = deepcopy(Γ_in)

 O = build_reorder_sector(L2)
 for kx=1:N_sectors
   Γ[kx,:,:] = O*Γ[kx,:,:]*(O');
 end

 U_Γ  = zeros(Complex{Float64}, N_sectors, 4*L2, 4*L2);
 U    = zeros(Complex{Float64}, N_sectors, 4*L2, 4*L2);
 eye  = Matrix{Float64}(I, 4*L2, 4*L2);
 for kiter=1:N_sectors
  U_Γ[kiter,:,:] = eye;
 end
 v   = ones(Float64, 2, N_sectors);

##################
for L2iter=1:(L2-mdL)
  for k=1:N_sectors
      Γ[k,:,:], U[k,:,:]  = diagonalise_block(Γ[k,:,:],convert(Int64,v[2,k]), 2*(mdL+L2iter)-v[2,k]+1);
      U_Γ[k,:,:]          = U_Γ[k,:,:]*U[k,:,:];
      TBz                 = convert(Int64, v[2,k])
      v[1,k]              = real(Γ[k,TBz,TBz])
  end
  for L1iter=1:L1
    k = sortperm(v[1,:])
    i = 0;
    for j=1:N_sectors
      i=j;
      if (v[2,k[i]]<=(2*(mdL+L2iter)))
        break;
      end
    end
    Γ[k[i],Int(v[2,k[i]]),:] .= 0;
    Γ[k[i],:,Int(v[2,k[i]])] .= 0;
    Γ[k[i],Int(v[2,k[i]])+2*L2,:] .= 0;
    Γ[k[i],:,Int(v[2,k[i]])+2*L2] .= 0;
    Γ[k[i],Int(v[2,k[i]])+2*L2,Int(v[2,k[i]])+2*L2] = 1;

    v[2,k[i]] +=1;
    TBz                 = convert(Int64, v[2,k[i]])
    v[1,k[i]]           = real(Γ[k[i],TBz,TBz])
  end
end
##################

 for kiter=1:N_sectors
  Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_sectors
  Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end

 return Γ;
end


function Initialise_gamma_sector_H_pi_flux(Lp,Lo,kx)
  # M = rand(4*Lo,4*Lo);
  # M = (M+M')/2;




  H         = H_fasulla_2(Lp,Lo,kx);
  H_D, U_D  = Diag_ferm(H);
  M         = GS_Gamma(H_D, U_D);
 return M;
end


function Energy_fermions(Mat,D,U)
   N_f = convert(Int64, size(Mat,1)/2.);

   energy = 0;

   M = deepcopy(Mat)
   #####AGGIUNTO
   M = (M+M')/2.
   # M[1,1] += eps()
   ####

   M_diag_base = real(U'*M*U);
   for iiter=1:(N_f)
       energy += M_diag_base[iiter,iiter]*D[iiter+N_f,iiter+N_f];
       energy += M_diag_base[iiter+N_f,iiter+N_f]*D[iiter,iiter];
   end

   return real(energy);
end
