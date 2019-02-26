function [ c, A, b, free_variables, objective_const_term ] = LP_Convert_to_Standard_Form( c, A, b, lb, ub, sense )
% ==================================================================================================================== %
% LP_Convert_to_Standard_Form( c, A, b, lb, ub, sense ):
% This function takes as input the data of an LP in the following form:
%                       min    c^T x
%                       s.t.   Ax {<=,=,>=} rhs,
%                              lb <= x <= ub
% and transforms it in a semi-standard form, that is:
%                       min    c_bar^T x
%                       s.t.   A_bar x = b
%                              (x_C)_i >= 0, for i in C, (x_F)_i free for i in F,
% where C in {1,...,n} is a set the indices of constrained variables and F = {1,...,n}\C, is the set of indices
% of free variables, and n is the number of variables of the final model.
%
% Author: Spyridon Pougkakiotis.
% ==================================================================================================================== %

% ==================================================================================================================== %
% Test input data, dimensions, e.t.c.
% -------------------------------------------------------------------------------------------------------------------- %
n = size(A,2);
m = size(A,1);
if (size(lb,2) > 1)
    lb = lb';
elseif (size(ub,2) > 1)
    ub = ub';
elseif (size(sense,2) > 1)
    sense = sense';
elseif (size(c,2) > 1)
    c = c';
elseif (~issparse(A))
    A = sparse(A);
end

if (size(c,1) ~= n || size(sense,1) ~= m || size(lb,1) ~= n || size(ub,1) ~= n)
    error("Incorrect input dimensions")
end
% ==================================================================================================================== %


% ==================================================================================================================== %
% Initialization.
% -------------------------------------------------------------------------------------------------------------------- %
num_of_slacks = 0; % Counter for the slack variables to be added in the inequality constraints.
free_variables = []; % To store the indices of the free variables.
objective_const_term = 0; % To keep constants that need to be added in the objective.
extra_constraints = 0; % Counter for the extra constraints to be added in case of double bounds.
[rows,cols,v] = find(A);
b_new = [];
% ==================================================================================================================== %

% ==================================================================================================================== %
% Make all the constraints to be of equality type (add slack variables)
% -------------------------------------------------------------------------------------------------------------------- %
for i = 1:m    
    if ( sense(i) == '<') % add a slack of the form +x_slack.               
       num_of_slacks = num_of_slacks + 1;
       rows = [rows; i];
       cols = [cols; n + num_of_slacks];
       v = [v; 1]; % assign 1 in the element A(i,n+num_of_slacks).        
    elseif ( sense(i) == '>') % add a slack of the form -x_slack.                 
       num_of_slacks = num_of_slacks + 1;
       rows = [rows; i];
       cols = [cols;  n + num_of_slacks];
       v = [v; -1]; %assign -1 in the element A(i,n+num_of_slacks).    
    end   
end
% ==================================================================================================================== %

% ==================================================================================================================== %
% Add extra constraints to treat the upper and lower bounds of Netlib PBs
% -------------------------------------------------------------------------------------------------------------------- %
for i = 1:n    
    if ((ub(i) == Inf) && (lb(i)> -Inf)) % We only have a lower bound.         
        % In this case we implicitly substitute x_i = w_i + lb(i), w_i >=0
        if (lb(i) ~= 0)
            b(:) = b(:) - A(:,i).*lb(i);  
            objective_const_term = objective_const_term + c(i)*lb(i);
        end
    elseif ((lb(i) == -Inf) && (ub(i) < Inf)) % We only have an upper bound.      
        % In this case we implicitly substitute x_i = ub(i) - w_i, w_i >=0
        k_max = size(cols,1);
        for k = 1:k_max
            if (cols(k) == i)
                v(k) = -v(k);
            end
        end  
        objective_const_term = objective_const_term + c(i)*ub(i);
        c(i) = -c(i); 
        if (ub(i) ~= 0)
            b(:) = b(:) - A(:,i).*ub(i); 
        end
    elseif ((lb(i) > -Inf) && (ub(i) < Inf)) % We have both upper and lower bound.
        % In this case we implicitly substitute x_i = w_i + lb(i)
        if (lb(i) ~= 0)
            b(:) = b(:) - A(:,i).*lb(i);
            objective_const_term = objective_const_term + c(i)*lb(i);
        end
        
        extra_constraints = extra_constraints + 1; %add one constraint, one variable
        num_of_slacks = num_of_slacks + 1;
        b_new = [b_new; ub(i) - lb(i)]; %The RHS of extra constraint w_i + w_i_2 = ub_i - lb_i
        rows = [rows; m + extra_constraints ; m + extra_constraints];
        cols = [cols; i ; n + num_of_slacks];
        v = [v; 1; 1]; % assigns ones in the element A(m+extra_constr,i) and A(m+extra_constr,n+num_of_slacks)  
    else % Otherwise, the variable is free.     
        free_variables = [free_variables; i]; % Simply keep track of them.
    end

end
% ==================================================================================================================== %
b = [b; b_new];
c = [c; zeros(num_of_slacks,1)];
A = sparse(rows,cols,v,m+extra_constraints,n+num_of_slacks);
end

