function [dx,dy,dz,instability] = Newton_backsolve(NS,res_p,res_d,res_mu,pos_vars,free_vars)
% ==================================================================================================================== %
% Newton_backsolve    Solve linear system with factorized matrix, by using backward substitution.
% -------------------------------------------------------------------------------------------------------------------- %
% OUTPUT:
%  [dx,dy,dz,instability] = newtonsolve(NS,res_p,res_d,res_mu,A,A_tr,pos_vars,free_vars)
%  i.e. the Newton direction and a boolean parameter indicating critical ill-conditioning.
%
% Author: Spyridon Pougkakiotis.
% ==================================================================================================================== %
m = size(res_p,1);
n = size(res_d,1);
instability = false;
dx = zeros(n,1);
dz = zeros(n,1);
dy = zeros(m,1);
if (size(pos_vars,1) > 0)
    temp_res = zeros(n,1);
end
% ==================================================================================================================== %
% Solve KKT system with LDL' factors.
% -------------------------------------------------------------------------------------------------------------------- %
if (size(pos_vars,1) > 0)
    temp_res(pos_vars) =  res_mu(pos_vars)./(NS.x(pos_vars));
    rhs = [res_d-temp_res; res_p];
else
    rhs = [res_d; res_p];
end
warn_stat = warning;
warning('off','all');
lhs = NS.L'\(NS.D\(NS.L\rhs(NS.pp)));
if (nnz(isnan(lhs)) > 0 || nnz(isinf(lhs)) > 0)
    instability = true;
    return;
end
warning(warn_stat);
lhs(NS.pp) = lhs;
dx = lhs(1:n,1);
dy = lhs(n+1:n+m,1);
if (size(pos_vars,1) > 0)
    dz(pos_vars) = (res_mu(pos_vars)-NS.z(pos_vars).*dx(pos_vars))./NS.x(pos_vars);
    dz(free_vars) = 0;
end
% ==================================================================================================================== %

% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
end 