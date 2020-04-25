function [M,lambda,V] = constructMonodromy(stm_T2)
% Function constructs monodromy matrix based on stm input at half period
% Yuri Shimane, 2020/03/26
% ============================================================== %
% FORMAT:
%   [M,eigval,eigvec] = constructMonodromy(stm_T2)
% INPUT: 
%   stm_T2 : 6 by 6 stm mapping from t = 0 to t = P/2
% OUTPUT:
%   M      : 6 by 6 monodromy matrix
%   lambda : 6 by 1 array of eigenvalues
%   V      : 6 by 6 matrix of eigenvectors (each of 6 x 1)
% ============================================================== %

% A-matrix from Howell 1984 eq(7)
A = [1 0 0 0 0 0;
     0 -1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 -1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 -1];
% construct monodromy matrix 
M = A * inv(stm_T2) * A * stm_T2;  % Howell 1984 eq(9)
% eigenvalue analysis
[V,D] = eig(M);   
% >>> returns diagonal matrix D of eigenvalues and matrix V whose columns 
%     are the corresponding right eigenvectors, so that A*V = V*D
lambda = diag(D);


end

