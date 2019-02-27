function [D,D_L] = Scale_the_problem(A,scale_option,direction)
% ==================================================================================================================== %
% [D] = Scale_the_problem(A): 
% -------------------------------------------------------------------------------------------------------------------- %
% This function, takes as an input a sparse matrix A, representing the constraints of a quadratic progrmaming problem.
% It checks whether the matrix is well scaled, and if not, it applies some linear transformations to the matrix, in 
% order to improve its numerical properties. The method return a diagonal matrix D, which the user should 
% use in order to recover the solution of the initial problem (after solving the scaled one).
% The optional parameter: scale_option. This parameter can take 3 values:
%   (i)   scale_option = 0, no scaling will be applied.
%   (ii)  scale_option = 1, iterative geometric scaling will be used.
%   (iii) scale_option = 2, equilibrium scaling is employed.
%   (iv)  scale_option = 3, nearest power of 2 scaling is used.
%   (v)   scale_option = 4, mixed strategy, based on the properties of the respective row/column.
% The optional parameter: direction. This parameter can take 3 values:
%   (i)   direction = 'r', right scaling (default).
%   (ii)  direction = 'l', left scaling.
% For more information about these scaling choices, the reader is refered to:
%                                   https://en.wikibooks.org/wiki/GLPK/Scaling
%
% Author: Spyridon Pougkakiotis.
% ==================================================================================================================== %
    if (nargin < 2 || isempty(scale_option))
        scale_option = 1; % Set geometric scaling as the default option.
    end
    if (nargin < 3 || isempty(direction))
        direction = 'r';
    end
    if (direction == 'l')
        A = A';
    end     
    pos_A = abs(A);       % Need it to identify non-zero elements.
    D = zeros(size(A,2),1);
    D_L = [];
    pos_ind = pos_A > 0;   
    % ================================================================================================================ %
    % Based on the input parameters, build the desired scaling factor.
    % ---------------------------------------------------------------------------------------------------------------- %
    if (max(max(pos_A)) <= 10 && min(min(abs(A(pos_ind))))>= 0.1 || scale_option == 0)
        fprintf('No scaling necessary.\n'); % Well scaled or no scale.
        D = ones(size(A,2),1);
    elseif (scale_option == 1) % Geometric scaling (applied on columns for computational efficiency).
        fprintf('The constraint matrix is scaled. Geometric scaling is employed.\n');
        for j = 1:size(A,2)
            rows = pos_A(:,j) > 0; % Find all non-zero elements for this column
            if (any(rows))
                %size(pos_A(rows,j))
                maximum = max(pos_A(rows,j));
                minimum = min(pos_A(rows,j));
                if (maximum*minimum > 10^(-12) && maximum*minimum < 10^(12))
                    D(j) = 1/sqrt(maximum*minimum);
                else
                    D(j) = 1;
                end
            else           
                D(j) = 1; % Extreme case, where one column is all zeros.
            end
        end
    elseif (scale_option == 2) % Equilibrium scaling (applied on columns for efficiency).
        fprintf('The constraint matrix is scaled. Equilibrium scaling is employed.\n');
        for j = 1:size(A,2)
            rows = pos_A(:,j) > 0; % Find all non-zero elements for this column.
            maximum = max(pos_A(rows,j));
            if (maximum > 10^(-6))
                D(j) = 1/maximum;
            else
                D(j) = 1;
            end
        end
    elseif (scale_option == 3) % Nearest power of 2 scaling + geometric scaling (avoiding rounding errors).
        fprintf('The constraint matrix is scaled. Geometric scaling with nearest power of 2 is employed.\n');
        for j = 1:size(A,2)
            rows = pos_A(:,j) > 0; % Find all non-zero elements for this column.
            if (any(rows))
                maximum = max(pos_A(rows,j));
                minimum = min(pos_A(rows,j));
                p = nextpow2(sqrt(maximum*minimum));
                if (maximum*minimum > 10^(-12) && maximum*minimum < 10^(12))
                    D(j) = 1/(2^(p-1));
                else
                    D(j) = 1;
                end
            else
                D(j) = 1; % Extreme case, where one column is all zeros.
            end
        end
    elseif (scale_option == 4)
        fprintf('The constraint matrix is scaled. Mixed scaling is employed.\n');
        for j = 1:size(A,2)
            rows = pos_A(:,j) > 0; % Find all non-zero elements for this column
            if (any(rows))
                %size(pos_A(rows,j))
                maximum = max(pos_A(rows,j));
                minimum = min(pos_A(rows,j));
                if (maximum > 10^3 && minimum < 10^(-3))
                    p = nextpow2(sqrt(maximum*minimum));
                    D(j) = 1/(2^(p-1));
                elseif (1/minimum > maximum && minimum > 10^(-6))
                    p = nextpow2(minimum);
                    D(j) = 1/(2^(p-1));
                elseif (maximum < 10^(6))
                    p = nextpow2(maximum);
                    D(j) = 1/(2^(p-1));
                else
                    D(j) = 1;
                end
            else           
                D(j) = 1; % Extreme case, where one column is all zeros.
            end
        end
    end
    
    % ================================================================================================================ %
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %
