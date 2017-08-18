function sum_channels_all_nodes

input_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/num_available_channels_on_each_node/';
%input_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/output3/';
ouput_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/num_available_channels_on_each_node/';


    file_list  = dir([input_dir 'numCC*']);
    
    fid = fopen([ouput_dir 'sum_channels_all_nodes.txt'],'w');
    for i = 1: length(file_list)
        file_name = [input_dir file_list(i).name];
        
        channel = dlmread(file_name);
        
        sum_channels_one_node = sum(channel);
        
        fprintf(fid,'%s %d', file_list(i).name, sum_channels_one_node);        
        %fprintf(fid,'%d', sum_channels_one_node);
        fprintf(fid, '\n');
    end
    fclose(fid);
end
