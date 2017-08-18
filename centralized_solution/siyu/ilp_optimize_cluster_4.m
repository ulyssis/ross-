function ilp_optimize_cluster_4
% optimize method for clustering
addpath(genpath('/home/tang/code/gurobi/gurobi562/linux64/'))
% read data

% data_dir = './input_50/';
data_dir = './input_50nodes/';

file_list  = dir([data_dir 'potential_clusters_*']);
% which node number per cluster we want to optimize.

opt_node_number = [4, 3, 2];

para = [0.8];

for p = 1:length(para)
    
    alph = [para(p), 0.15, 0.05];
    
    % output_dir
    output_dir = [data_dir 'output/alph_' num2str(alph(1)) '/'];
    if ~exist(output_dir)
        mkdir(output_dir);
    end
    
    fid = fopen([output_dir 'result.txt'],'w');
    good_cluster_4 = zeros(1,length(file_list));
    good_cluster_3 = zeros(1,length(file_list));
    good_cluster_2 = zeros(1,length(file_list));
    channel_per_cluster_all = 0;
    channel_per_cluster_2_3_all = 0;
    for k = 1 : length(file_list)
        % load score matrix
        file_name = file_list(k).name;
        fprintf(['computing for file ' file_name '\n']);
        score_matrix = dlmread([data_dir file_name ]);
        num_node = size(score_matrix, 2);
        
        % single node channel
        channel_matrix = score_matrix(1:num_node,:);
        channel = zeros(1, num_node);
        [r c v ] = find(channel_matrix);
        channel(r) = v;
        
        % find number of node per cluster and return only such clusters
        score_matrix = score_matrix(num_node+1:end,:);
        [score_matrix, node_per_cluster, cluster_size]= opt_node(score_matrix, opt_node_number);
        num_cluster =  size(score_matrix, 1);
        % possible max cost per cluster
        max_score = max(score_matrix(:));
        max_cost = (max_score + 1) * opt_node_number(1);
        
        non_active_node_cost = max_cost * 1000;
        non_active_cluster_cost = max_cost*200 *ones(1,num_cluster);
        
        
        % compute cost_matrix
        cost_matrix = -score_matrix;
        min_c = min(cost_matrix(:));
        cost_matrix(cost_matrix == 0) = non_active_node_cost;
        cost_matrix = cost_matrix + abs(min_c) + 1;
        % !!!!!!!!!!!!!!!re-weight non_active_cluster_cost !!!!!!!!!!!!!!!!!!
        
        non_active_cluster_cost(find(node_per_cluster == opt_node_number(1))) = non_active_cluster_cost(find(node_per_cluster == opt_node_number(1)))*alph(1);
        non_active_cluster_cost(find(node_per_cluster == opt_node_number(2))) = non_active_cluster_cost(find(node_per_cluster == opt_node_number(2)))*alph(2);
        non_active_cluster_cost(find(node_per_cluster == opt_node_number(3))) = non_active_cluster_cost(find(node_per_cluster == opt_node_number(3)))*alph(3);
        
        % construct function
        f_x = reshape(cost_matrix, 1, num_node*num_cluster);
        %     f_w = zeros(1, num_cluster);
        f_w = non_active_cluster_cost;
        f = [f_x, f_w]';
        %     f = [f_x]';
        
        dim = length(f);
        % construct equality constraint
        % equality_1
        A_1 = zeros (num_node, dim);
        
        % for every node, it can only belong to one cluster
        for i = 1 : num_node
            %             A_1(i,  i : num_node: num_node*num_cluster) = 1;
            A_1(i,  (i-1)*num_cluster +1 :(i-1)*num_cluster + num_cluster) = 1;
        end
        b_1 = ones(num_node,1);
        
        % % equality_2
        
        A_2 = sparse(num_cluster, dim);
        for j = 1 : num_cluster
            fprintf('cluster: %d \n',j);
            A_2(j, j: num_cluster: num_node*num_cluster) = 1;
            A_2(j, num_node*num_cluster + j) = node_per_cluster(j);
            %     A_2(j, ((j-1)*num_node + 1): ((j-1)*num_node + num_node) ) = 1;
            %     A_2(j, num_node*num_cluster + j) = 2;
        end
        b_2 = node_per_cluster;%*ones(num_cluster,1);
        
        %     A = [A_1;A_2];
        %     b = [b_1;b_2];
        lb = zeros(1,dim);
        ub = ones(1,dim);
        
        % call gurobi
        % [result, fval] = bintprog(f,[],[],A,b);
        %     [result, fval, flag]  = linprog_gurobi(f',[],[],A,b,lb,ub);
        
        fprintf('call gurobi\n');
        [result, fval, flag]  = linprog_gurobi(f',A_1,b_1,A_2,b_2,lb,ub);
        
        x_matrix = reshape(result,num_cluster,num_node +1);
        
        
        x_matrix_assignment = x_matrix(:,1:num_node);
        
        single_note_indicator = find(sum(x_matrix_assignment,1) == 0);
        cluster_idx = find(x_matrix(:,end) == 0);
        
        
        if num_node == length(single_note_indicator)
            test = 1;
        end
        
        single_channel = channel(single_note_indicator);
        
        [r c v] = find(x_matrix_assignment);
        score = [];
        for i = 1: length(v)
            score(i) = score_matrix(r(i),c(i));
        end
        
        score(end+1:end+length(single_channel)) = single_channel;
        channel_per_node = sum(score)/length(score);
        
        cluster_channel = [];
        cluster_note_indicator = cell(1,1);

        for c = 1 : length(cluster_idx)
            tmp_node_idx = find(x_matrix_assignment(cluster_idx(c),:));
            cluster_note_indicator{c} = tmp_node_idx;
            cluster_channel(c) = score_matrix(cluster_idx(c),tmp_node_idx(1));
            if length(tmp_node_idx) == opt_node_number(1)
                good_cluster_4(k) = good_cluster_4(k) + 1;
            end
            if length(tmp_node_idx) == opt_node_number(2)
                good_cluster_3(k) = good_cluster_3(k) + 1;
            end
             
            if length(tmp_node_idx) == opt_node_number(3)
                good_cluster_2(k) = good_cluster_2(k) + 1;
            end       
        end
        %     cluster_note_indicator
        %     cluster_channel
        channel_per_cluster = (sum(cluster_channel) + sum (single_channel))/(length(cluster_channel) + length(single_channel));
        channel_per_cluster_2_3 = sum(cluster_channel) /length(cluster_channel);
        
        channel_per_cluster_all = channel_per_cluster_all+channel_per_cluster;
        channel_per_cluster_2_3_all = channel_per_cluster_2_3_all + channel_per_cluster_2_3;
        
        save([output_dir file_name '_' num2str(opt_node_number(1)) '_' num2str(alph(1)) '_x_matrix.mat'],'x_matrix');
        save([output_dir file_name '_' num2str(opt_node_number(1)) '_' num2str(alph(1)) '_score_matrix.mat'],'score_matrix');
        save([output_dir file_name '_' num2str(opt_node_number(1)) '_' num2str(alph(1)) '_final.mat'],'channel_per_node','channel_per_cluster','single_note_indicator','cluster_idx','cluster_note_indicator','single_channel','cluster_channel','fval','flag');
        %
        fprintf(fid, '(%f, %f, %d, %d, %d) \n', channel_per_cluster,channel_per_cluster_2_3, good_cluster_4(k),good_cluster_3(k),good_cluster_2(k) );
        
    end
    
    
    fprintf(fid, '(%f, %f, %f, %f, %f) \n', sum(channel_per_cluster_all)/length(file_list),sum(channel_per_cluster_2_3_all)/length(file_list), sum(good_cluster_4)/length(good_cluster_4), sum(good_cluster_3)/length(good_cluster_3),sum(good_cluster_2)/length(good_cluster_2) );
    fclose(fid);
end
end
function [cost_matrix, cluster_not_choose_cost]=normalize(cost_matrix, num_cluster)
for n = 1 : num_cluster
    [value, ~] = min(cost_matrix(n,:));
    idx = find(cost_matrix(n,:) == value);
    cluster_not_choose_cost(n) = 15 - value;
    for m = 1 : length(idx)
        cost_matrix(n,idx(m)) = cost_matrix(n,idx(m))/length(idx);
    end
end
end
% save([data_dir 'gom.mat'],'gom');
