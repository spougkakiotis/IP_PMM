
clc; clear;
fileID1 = fopen('QP_problems_performance_prof_time.txt','r');
fileID2 = fopen('QP_problems_performance_prof_iter.txt','r');

A_t = fscanf(fileID1,'%f');
A_it = fscanf(fileID2,'%f');

num_of_probs = size(A_t,1)/2;
c_t = zeros(num_of_probs,2);
c_it = zeros(num_of_probs,2);

c_t(:,1) = A_t(1:num_of_probs,1);      %2 norm
c_t(:,2) = A_t(num_of_probs+1:end);    %no regularization

c_it(:,1) = A_it(1:num_of_probs,1);                  %2 norm
c_it(:,2) = A_it(num_of_probs+1:2*num_of_probs,1);   %no regularization

c_min_t = zeros(num_of_probs,1);
c_min_it = zeros(num_of_probs,1);

for i = 1:num_of_probs
    c_min_t(i) = min(c_t(i,:));
    c_min_it(i) = min(c_it(i,:));
end

r_t = zeros(num_of_probs,2);
r_it = zeros(num_of_probs,2);

for j = 1:2
    for i = 1:num_of_probs
        r_t(i,j) = c_t(i,j)/c_min_t(i);
        r_it(i,j) = c_it(i,j)/c_min_it(i);
    end
end

t_t = linspace(1,25,200);
t_it = linspace(1,25,200);
perf_func_t = zeros(200,2);
perf_func_it = zeros(200,2);
for k = 1:200
    for j = 1:2
        perf_func_t(k,j) = sum(r_t(:,j) <= t_t(k))/num_of_probs;
         perf_func_it(k,j) = sum(r_it(:,j) <= t_it(k))/num_of_probs;
    end
end
   semilogx(t_t,perf_func_t(:,1),'g^',t_t,perf_func_t(:,2),'b*','MarkerSize',14,'LineWidth',3);
   lgd = legend('IP-PMM','Non-regularized IPM','Location','southeast');
   lgd.FontSize = 25;
   xlabel('Performance ratio (time)','FontSize',18);
   ylabel('Problems solved','FontSize',18);
   set(gca,'FontSize',18)

   figure;
   semilogx(t_it,perf_func_it(:,1),'g^',t_it,perf_func_it(:,2),'b*','MarkerSize',14,'LineWidth',3);
   lgd = legend('IP-PMM','Non-regularized IPM','Location','southeast');
   lgd.FontSize = 25;   
   xlabel('Performance ratio (iterations)','FontSize',18);
   ylabel('Problems solved','FontSize',18);
      set(gca,'FontSize',18)

