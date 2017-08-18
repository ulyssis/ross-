function ilp
% optimize method for clustering

% read data
data_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/';
%  data_dir = './';
file_name = 'numberCC.txt';

score_matrix = dlmread([data_dir file_name ]);
num_node = size(score_matrix, 2);
num_cluster =  size(score_matrix, 1);

cost_matrix = -score_matrix;
min_c = min(cost_matrix(:));

cluster_cost = sum(cost_matrix,2);

cost_matrix(cost_matrix == 0) = 10000;
cost_matrix = cost_matrix + abs(min_c) + 1;
cluster_cost = cluster_cost + abs(min_c)*2 + 2;
%cost_matrix = cost_matrix +3;
%cluster_cost = cluster_cost + 3*2;

% construct function

f_x = reshape(cost_matrix, 1, num_node*num_cluster);
f_w = zeros(1, num_cluster);
% f_w = cluster_cost';
f = [f_x, f_w]';
% f = [f_x]';

dim = length(f);
% construct equality constraint
% equality_1
A_1 = zeros (num_node, dim);

% for every node, it can only belong to one cluster
for i = 1 : num_node
%     A_1(i,  i : num_node: num_node*num_cluster) = 1;
     A_1(i,  (i-1)*num_cluster +1 :(i-1)*num_cluster + num_cluster) = 1;
end
b_1 = ones(num_node,1);

% % equality_2
A_2 = sparse(num_cluster, dim);
for j = 1 : num_cluster
    A_2(j, j: num_cluster: num_node*num_cluster) = 1;
    A_2(j, num_node*num_cluster + j) = 2;
%     A_2(j, ((j-1)*num_node + 1): ((j-1)*num_node + num_node) ) = 1;
%     A_2(j, num_node*num_cluster + j) = 2;
end
b_2 = 2*ones(num_cluster,1);

A = [A_1;A_2];
b = [b_1;b_2];
% lb = zeros(1,dim);
% ub = ones(1,dim);

% call gurobi


[result, fval] = bintprog(f,[],[],A,b);
% result = linprog_gurobi(f',[],[],A,b,lb,ub);
x_matrix = reshape(result,num_cluster,num_node +1);

save('x_matrix.mat','x_matrix');
% convert result
end