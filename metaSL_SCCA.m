function [u, v, res_iter] = metaSL_SCCA(S_XX, Beta, S_YY, para)
% --------------------------------------------------------------------
% This function used to compute metaSL-SCCA 
% --------------------------------------------------------------------
% Input:
%       - S_XX: genotypic matrix, d x d
%       - Beta: summary statistics matrix from GWAS, d x c
%       - S_YY: phenotypic matrix, c x c
%       - para, parameters: lambda1, lambda2, r1, r2 
% Output:
%       - u: weights for SNPs, d x 1
%       - v: weights for QTs, c x 1

% Initialization
[d, c] = size(Beta);
u0 = ones(d, 1); % d x 1
v0 = ones(c, 1); % c x 1

% scale u and v
scale = sqrt(u0' * S_XX * u0);
u = u0 ./ scale;
scale = sqrt(v0' * S_YY * v0);
v = v0 ./ scale;

% stop criteria
tol = 1e-5;
i = 0;
tol_u = inf;
tol_v = inf;

while (i < 100 && (tol_u > tol || tol_v > tol))% default 100 times of iteration
    i = i + 1;
    
    % update u
    D1 = diag(1 ./ (u + eps));
    u_old = u;
    u = (para.lambda(1) * D1 + para.r(1) * S_XX) \ (Beta * v);
    u = u ./  sqrt(u' * S_XX * u);
    
    % update v
    D2 = diag(1 ./ (v + eps));
    v_old = v;
    v = (para.lambda(2) * D2 + para.r(2) * S_YY) \ (Beta' * u);
    v = v ./  sqrt(v' * S_YY * v);
    
    res_iter(i) = -u' * Beta * v + para.lambda(1) * sum(abs(u)) + para.lambda(2) * sum(abs(v)) +...
                       para.r(1) * (u' * S_XX * u - 1) + para.r(2) * (v' * S_YY * v - 1);
    
     if i == 1
         continue
     else
        tol_u = max(abs(u - u_old));
        tol_v = max(abs(v - v_old));
     end
     
     if tol_u < tol && tol_v < tol
         break
     end
end
 