function [Out_X,U,V,RMSE_1] = HOAT(M_Omega,Omega_array,rak,maxiter,ip)
% This matlab code implements Robust Matrix Completion based on Factorization and Truncated-Quadratic Loss Function.
%
% M_Omega - m x n observed matrix 
% Omega_array is the index of observed entries
% rak is the rank of the object matrix
% maxiter - maximum number of iterations
% ip is the parameter in (35) for our paper 

[m,n] = size(M_Omega);
norm_M_Omega = norm(M_Omega,'fro');

% inital estimation for X, U, V and N
[X,U,V] = LP_F(M_Omega,Omega_array,rak,3);
T = (M_Omega - X).*Omega_array;    
t_m_n = T(find(T));
scale = 10*1.4815*median(abs(t_m_n - median(t_m_n)));
sigma = ones(m,n)*diag(scale);
ONE_1=ones(m,n);
ONE_1(abs(T-median(t_m_n))-sigma<0)=0;
N = T.*ONE_1;

for iter = 1 : maxiter
    D = (M_Omega - N).*Omega_array;
    
    for j = 1:n
        row = find(Omega_array(:,j) == 1);
        U_I =  U(row,:);
        b_I = D(row,j);
        V(:,j) = pinv(U_I,1e-8)* b_I;
    end
    
    for i = 1 : m
        col = find(Omega_array(i,:) == 1);
        V_I = V(:,col);
        b_I = D(i,col);
        U(i,:) = b_I * pinv(V_I,1e-8);
    end
    X = U*V;
    
    
    T = (M_Omega - X).*Omega_array;      
    
    t_m_n = T(find(T));
    scale =min([scale ip*1.4815*median(abs(t_m_n - median(t_m_n)))]);
    
    sigma = ones(m,n)*diag(scale);
    ONE_1=ones(m,n);
    ONE_1(abs(T-median(t_m_n))-sigma<0)=0;
    N = T.*ONE_1;

% stop condition
    RMSE_1= norm((M_Omega-X-N).*Omega_array,'fro')/norm_M_Omega;
    if RMSE_1 < 0.000001
        break;
    end
end
    Out_X = X;
end