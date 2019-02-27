%This script loads various NETLIB problems and solves them using IP_PMM
clear all;
clc;
%The path on which all the netlib problems lie
Netlib_path = '../NETLIB_PROBLEMS_IN_MATLAB_FORM/netlib'; 
%Finds all the Netlib problems and stores their names in a struct
d = dir(fullfile(Netlib_path,'*.mat')); 


%Open the file to write the results
fileID = fopen('Netlib_tabular_format_final_results.txt','a+');
fields = {'A','obj','sense','rhs','lb','ub','vtype','modelname','varnames','constrnames'};
total_iters = 0;
total_time = 0;
scaling_direction = 'r';
scaling_mode = 3;
pc_mode = true;
tol = 1e-6;
print_mode = 1;
%Each indice k=1..num_of_netlib_files gives the name of each netlib problem through d(i).name
for k = 1:1
    load(fullfile(Netlib_path,d(k).name))
   
    c = model.obj; 
    A = model.A;
    b = model.rhs;

    
    [c,A,b,free_variables,objective_const_term] = LP_Convert_to_Standard_Form(c, A, b, model.lb, model.ub, model.sense);

    n = size(A,2);
    Q = sparse(n,n);
     if (scaling_direction == 'r')
        [D,~] = Scale_the_problem(A,scaling_mode,scaling_direction);
        A = A*spdiags(D,0,n,n); % Apply the right scaling.
        c = c.*D;
    elseif (scaling_direction == 'l')
        [D,~] = Scale_the_problem(A,scaling_mode,scaling_direction);
        A = spdiags(D,0,m,m)*A;  % Apply the left scaling.
        b = b.*D;
    elseif (scaling_direction == 'b')
        [D_R,D_L] = Scale_the_problem(A,scaling_mode,scaling_direction);
        if (size(D_L,1) ~= 0)
            A = (spdiags(D_L.^(1/2),0,m,m)*A)*spdiags(D_R.^(1/2),0,n,n);
            b = b.*D_L;
        else
            A = A*spdiags(D_R,0,n,n); % Apply the right scaling.        
        end
        c = c.*D_R;
    end
    time = 0;
    tic;
    [x,y,z,opt,iter] = IP_PMM(c,A,Q,b,free_variables,tol,200,pc_mode,print_mode); 
    total_iters = total_iters + iter;
    time = time + toc;
    total_time = total_time + time;
    obj_val = c'*x + objective_const_term;
    if (opt == 1)
       fprintf(fileID,'%s & %d & %d & opt  \n',model.modelname, iter, time); 
       fprintf(fileID,'The optimal solution objective is %d.\n',obj_val);
    else
       fprintf(fileID,'%s & %d & %d & non-opt \n',model.modelname, iter, time); 
    end
    
end
fprintf(fileID,'The total iterates were: %d and the total time was %d\n',total_iters,total_time);
fclose(fileID);


