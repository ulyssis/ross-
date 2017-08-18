function ilp
% optimize method for clustering

% read data
% data_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/';
data_dir = './numberCC_size_3/';

channel_dir = [data_dir '/numCC_size_is_3/'];
node_per_cluster = 3;
    
    
% output_dir
output_dir = [data_dir 'output/'];
if ~exist(output_dir) 
    mkdir(output_dir);
end

file_list  = dir([data_dir 'number*']);

    
    
% gom = repmat(struct('name',[],'label',[],'channel',[],'score',[]),1,1);
for k = 1 : length(file_list)
    
    file_name = file_list(k).name;
    fprintf(['computing for file ' file_name '\n']);
    score_matrix = dlmread([data_dir file_name ]);
    
    channel = dlmread([channel_dir file_name(1:3) file_name(7:end)]);
    
    num_node = size(score_matrix, 2);
    num_cluster =  size(score_matrix, 1);
    
    cost_matrix = -score_matrix;
    min_c = min(cost_matrix(:));
    
    cluster_cost = sum(cost_matrix,2);
    
    cost_matrix(cost_matrix == 0) = 10000;
    cost_matrix = cost_matrix + abs(min_c) + 1;
    cluster_not_choose_cost = cluster_cost + abs(min_c)*node_per_cluster + node_per_cluster + node_per_cluster;
    %cost_matrix = cost_matrix +3;
    %cluster_cost = cluster_cost + 3*2;
    
    % construct function
    
    f_x = reshape(cost_matrix, 1, num_node*num_cluster);
%     f_w = zeros(1, num_cluster);
    f_w = cluster_not_choose_cost';
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
        A_2(j, num_node*num_cluster + j) = node_per_cluster;
        %     A_2(j, ((j-1)*num_node + 1): ((j-1)*num_node + num_node) ) = 1;
        %     A_2(j, num_node*num_cluster + j) = 2;
    end
    b_2 = node_per_cluster*ones(num_cluster,1);
    
    A = [A_1;A_2];
    b = [b_1;b_2];
    lb = zeros(1,dim);
    ub = ones(1,dim);
    
    % call gurobi
    
    
    % [result, fval] = bintprog(f,[],[],A,b);
%     [result, fval, flag]  = linprog_gurobi(f',[],[],A,b,lb,ub);
    [result, fval, flag]  = linprog_gurobi(f',A_1,b_1,A_2,b_2,lb,ub);
    x_matrix = reshape(result,num_cluster,num_node +1);
    
    
    x_matrix_assignment = x_matrix(:,1:num_node);
    
    note_indicator = find(sum(x_matrix_assignment,1) == 0);
    
    channel_single = channel(note_indicator);
    
    [r c v] = find(x_matrix_assignment);
    score = [];
    for i = 1: length(v)
        score(i) = score_matrix(r(i),c(i));
    end
    
    score(end+1:end+length(channel_single)) = channel_single;
    score_ave = sum(score)/length(score)
    

    save([output_dir file_name '_x_matrix.mat'],'x_matrix');
    save([output_dir file_name '_score_matrix.mat'],'score_matrix');
    save([output_dir file_name 'final.mat'],'score_ave','note_indicator','fval','flag');
    % convert result
end

% save([data_dir 'gom.mat'],'gom');
end