% obtain average cluster size and CI
function get_clusterSize

input_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/optimazation/output_size20/';
%input_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/output3/';
ouput_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/optimazation/output_size20';

%node_list = [102 90 81 72 60 51];
node_list = [ 20 ];

for n = 1 :  length(node_list)
    node = node_list(n);
    file_list  = dir([input_dir 'potential_clusters_*' num2str(node) '*x_matrix*']);
    
    %fid = fopen([ouput_dir 'number3_node_' num2str(node) '.txt'],'w');
    fid = fopen([ouput_dir '_cluster' '.txt'],'w');
    for i = 1: length(file_list)
        file_name = [input_dir file_list(i).name];
        load(file_name);

        % find out the rows with the element at colum 21 is 1
        TF1= x_matrix(:,21)==1; 
        
        % delete those rows
        x_matrix(TF1,:) =[];
        
        %calculate how many clusters, including the single node clusters
        nc= 20 - sum(sum(x_matrix)) + size(x_matrix, 1);
        
        % average cluster size
        vector_cluster_size = [sum(x_matrix,2)', ones(1, 20 - sum(sum(x_matrix)) )];
        
        size_aver = mean(vector_cluster_size);
        
        fprintf(fid,'%s:', file_list(i).name(1:end-9));
        fprintf(fid,' %f %f', nc, size_aver);
        fprintf(fid, '\n');
    end
    fclose(fid);
end
end

