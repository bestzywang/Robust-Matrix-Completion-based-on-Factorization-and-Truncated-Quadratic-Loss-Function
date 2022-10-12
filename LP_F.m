function [Out_X,U,V] = LP_F(M_Omega,Omega_array,rak,maxiter)
% This matlab code implements the l2-reg 
% method for Matrix Completion.
%
% M - m x n observed matrix 
% Omega_array is the subset
% rak is the rank of the object matrix
% maxIter - maximum number of iterations


[r,c] = size(M_Omega);
U = randn(r,rak);
V = randn(rak,c);

for iter = 1 : maxiter
    for j = 1:c
        row = find(Omega_array(:,j) == 1);
        U_I =  U(row,:);
        b_I = M_Omega(row,j);
        V(:,j) = pinv(U_I,1e-8)* b_I;
    end
    for i = 1 : r
        col = find(Omega_array(i,:) == 1);
        V_I = V(:,col);
        b_I = M_Omega(i,col);
        U(i,:) = b_I * pinv(V_I,1e-8);
    end
    X = U*V;
end
    Out_X = X;
end