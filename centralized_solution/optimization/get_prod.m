% obtain average cluster size and CI
function get_prod

input_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/optimazation/output_size20/';
%input_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/output3/';
output_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/optimazation/output_size20';

%node_list = [102 90 81 72 60 51];
node_list = [ 20 ];

for n = 1 :  length(node_list)
    node = node_list(n);
    file_list  = dir([input_dir 'potential_clusters_*' num2str(node) '*final*']);
    
    %fid = fopen([ouput_dir 'number3_node_' num2str(node) '.txt'],'w');
    fid = fopen([output_dir '_pord_CCC' '.txt'],'w');
    for i = 1: length(file_list)
        file_name = [input_dir file_list(i).name];
        load(file_name);


        
        fprintf(fid,'%s:', file_list(i).name(1:end-9));
        fprintf(fid,' %f %f', channel_per_node * 20, channel_per_cluster);
        fprintf(fid, '\n');
    end
    fclose(fid);
end
end

