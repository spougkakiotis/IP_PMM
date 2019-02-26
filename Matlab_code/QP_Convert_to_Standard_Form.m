function [model, b, free_variables, objective_const_term] = QP_Convert_to_Standard_Form(model)
% ==================================================================================================================== %
% QP_Convert_to_Standard_Form(model):
% This function takes as input the data of an LP in the following form:
%                       min    c^T x + (1/2)x^T H x
%                       s.t.   al <= Ax <= au,
%                              lb <= x <= ub
% and transforms it in a semi-standard form, that is:
%                       min    c_bar^T x + (1/2)x^T H_bar x
%                       s.t.   A_bar x = b
%                              (x_C)_i >= 0, for i in C, (x_F)_i free for i in F,
% where C in {1,...,n} is a set the indices of constrained variables and F = {1,...,n}\C, is the set of indices
% of free variables, and n is the number of variables of the final model.
%
% Author: Spyridon Pougkakiotis
% ==================================================================================================================== %

% ==================================================================================================================== %
% Test input data, dimensions, e.t.c.
% -------------------------------------------------------------------------------------------------------------------- %
n = size(model.A,2);
m = size(model.A,1);
if (size(model.xl,2) > 1)
    model.xl = model.xl';
elseif (size(model.xu,2) > 1)
    model.xu = model.xu';
elseif (size(model.au,2) > 1)
    model.au = model.au';
elseif (size(model.al,2) > 1)
    model.al = model.al';
elseif (size(model.g,2) > 1)
    model.g = model.g';
elseif (~issparse(model.A))
    model.A = sparse(model.A);
end

if (size(model.g,1) ~= n || size(model.al,1) ~= m || size(model.au,1) ~= m || size(model.xl,1) ~= n         ...
                         || size(model.xu,1) ~= n || size(model.H,1) ~= size(model.H,2) || size(model.H,1) ~= n)
    error("Incorrect input dimensions")
end
% ==================================================================================================================== %

% ==================================================================================================================== %
% Initialization.
% -------------------------------------------------------------------------------------------------------------------- %
num_of_slacks = 0;        % Counter for the slack variables to be added in the inequality constraints.
free_variables = [];      % To store the indices of the free variables.
objective_const_term = 0; % To keep constants that need to be added in the objective.
extra_constraints = 0;    % Counter for the extra constraints to be added in case of double bounds.
[rows,cols,v] = find(model.A);
if (size(rows,2) > 1)
    rows = rows';
end
if (size(cols,2) > 1)
    cols = cols';
end
if (size(v,2) > 1)
    v = v';
end
b = zeros(m,1);
% ==================================================================================================================== %

% ==================================================================================================================== %
% Make all the constraints to be of equality type (add slack variables)
% -------------------------------------------------------------------------------------------------------------------- %
for i = 1:m   
    if (model.au(i) ~= model.al(i))
        if (model.au(i) == Inf)
            % we only have a lower bound, and hence we should add a slack of the form -x_slack
            num_of_slacks = num_of_slacks + 1;
            rows = [rows; i];
            cols = [cols; n + num_of_slacks];
            v = [v; -1];        % assign -1 in the element A(i,n+num_of_slacks) 
            b(i) = model.al(i); % Fix the RHS        
        elseif (model.al(i) == -Inf)
            % we only have an upper bound, and hence we should add a slack of the form +x_slack
            num_of_slacks = num_of_slacks + 1;
            rows = [rows; i];
            cols = [cols; n + num_of_slacks];
            v = [v; 1];         % assign 1 in the element A(i,n+num_of_slacks)
            b(i) = model.au(i); % Fix the RHS     
        else
            % transform al <=Axi <=au to Axi' = aui, Axi' = ali
            extra_constraints = extra_constraints + 1;
            k_max = size(cols,1);
            for k = 1:k_max
                if (rows(k) == i)
                    cols = [cols; cols(k)];
                    rows = [rows; m + extra_constraints];
                    v = [v; v(k)];
                end
            end
            
            % treat the case of the upper bound
            num_of_slacks = num_of_slacks + 1;
            rows = [rows; i];
            cols = [cols; n + num_of_slacks];
            v = [v; 1];         % assign 1 in the element A(i,n+num_of_slacks)
            b(i) = model.au(i); % Fix the RHS
            
            % Now add a new constraint that will treat the case of the LB
            num_of_slacks = num_of_slacks + 1;
            rows = [rows; m + extra_constraints];
            cols = [cols; n + num_of_slacks];
            v = [v; -1];
            b = [b; model.al(i)]; % The RHS of the extra constraint     
        end      
    else
        b(i) = model.al(i); % Already an equality constraint.
    end
    
end
% ==================================================================================================================== %
model.A = sparse(rows,cols,v,m + extra_constraints, n + num_of_slacks); % Renew the matrix to incude new constraints.
b_new = [];
% ==================================================================================================================== %
% Add extra constraints to treat the upper and lower bounds on the variables.
% -------------------------------------------------------------------------------------------------------------------- %
for i = 1:n % I want the initial n, since only those variables have bounds  
    if ((model.xu(i) == Inf) && (model.xl(i)> -Inf)) % We only have a lower bound 
        % In this case we implicitly substitute x_i = w_i + model.xl(i), w_i >=0
        if (model.xl(i) ~= 0)
            b(:) = b(:) - model.A(:,i).*model.xl(i); 
            objective_const_term = objective_const_term + model.g(i)*model.xl(i);
        end      
    elseif ((model.xl(i) == -Inf) && (model.xu(i) == Inf)) % The variable is free.     
        free_variables = [free_variables; i]; % Simply keep track of them.   
    elseif ((model.xl(i) == -Inf) && (model.xu(i) < Inf)) % We only have an upper bound. 
        % In this case we implicitly substitute x_i = ub(i) - w_i, w_i >=0
        k_max = size(cols,1);
        for k = 1:k_max
            if (cols(k) == i)
                v(k) = -v(k);
            end
        end  
        objective_const_term = objective_const_term + model.g(i)*model.xu(i);
        model.g(i) = -model.g(i); 
        if (model.xu(i) ~= 0)
            b(:) = b(:) - model.A(:,i).*model.xu(i); 
        end
    else % We have both upper and lower bound.
        % In this case we implicitly substitute x_i = w_i + lb(i)
        if (model.xl(i) ~= 0)
            b(:) = b(:) - model.A(:,i).*model.xl(i);
            objective_const_term = objective_const_term + model.g(i)*model.xl(i);
        end
        
        extra_constraints = extra_constraints + 1; %add one constraint, one variable
        num_of_slacks = num_of_slacks + 1;
        b_new = [b_new; model.xu(i) - model.xl(i)]; %The RHS of extra constraint w_i + w_i_2 = ub_i - lb_i
        rows = [rows; m + extra_constraints ; m + extra_constraints];
        cols = [cols; i ; n + num_of_slacks];
        v = [v; 1; 1]; % assigns ones in the element A(m+extra_constr,i) and A(m+extra_constr,n+num_of_slacks)  
    end
end
% ==================================================================================================================== %
b = [b; b_new];
model.g = [model.g; zeros(num_of_slacks,1)];
model.A = sparse(rows,cols,v,size(b,1),size(model.g,1));
end

