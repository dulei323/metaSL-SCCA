function [sort_info] = perf_process(u, v, data)
% list the snp_info and QT_info accroding to the sum of weight value(abslute value)
% 
% Input:
%       - u: weight of X, g x 1
%       - v: weight of Y, p x 1
%       
% Output:
%       - sort_info: 6 x g

u = abs(u); % sum of every snp_weight, g x 1
[u_sort, u_sorted_idx] = sort(u, 'descend'); % get the desced list and its line_num index
v = abs(v); % sum of every snp_weight, g x 1
[v_sort, v_sorted_idx] = sort(v, 'descend'); 

snp_id = data.snp_id;
l = 200;
sort_info = cell(6, l);
for i = 1 : l
    sort_info{1, i} = snp_id{u_sorted_idx(i)}; % rsID
    sort_info{2, i} = u_sort(i); 
    sort_info{3, i} = u_sorted_idx(i); 
end     

Y_id = data.Y_id;
for j = 1 : size(v)
    sort_info{4, j} = Y_id{v_sorted_idx(j)}; 
    sort_info{5, j} = v_sort(j); 
    sort_info{6, j} = v_sorted_idx(j); 
end  














