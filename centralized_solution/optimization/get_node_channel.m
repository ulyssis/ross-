function get_node_channel

data_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/optimazation/output_size20_initial/';                                                   
output_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/optimazation/output_size20';                                          
          
 file_list  = dir([data_dir 'potential_clusters_20*']);  
 fid = fopen([output_dir '_pord_singleNode' '.txt'],'w');
 for k = 1 : length(file_list)                                               
                                                                          
     file_name = file_list(k).name;                                          
     fprintf(['computing for file ' file_name '\n']);                        
    score_matrix = dlmread([data_dir file_name ]);                          
                                                                            
     num_node = size(score_matrix, 2);                                       
    channel_matrix = score_matrix(1:num_node,:);  
  % single node channel                                                   
   channel = zeros(1, num_node);                                           
   [r c v ] = find(channel_matrix);                                        
   channel(r) = v;                                                         
   fprintf(fid,' %d', sum(channel));
           fprintf(fid, '\n');
 end
    fclose(fid);
   
end