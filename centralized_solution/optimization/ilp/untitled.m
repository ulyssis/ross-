
data_dir = '/home/di/work/work_on_PC/ROSS_paper/ross_clustering/centralized_solution/';
%  data_dir = './';
file_name = 'numberCC.txt';

score_matrix = dlmread([data_dir file_name ]);

load('x_matrix.mat');
x_matrix = x_matrix(:,1:50);
[r c v] = find(x_matrix);

for i = 1: length(v)
 score(i) = score_matrix(r(i),c(i));
end

all_score = sum(score)/50;
