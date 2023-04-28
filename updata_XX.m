function S_XX_new = updata_XX(S_XX, s)
% keep the data in diagonal line and > s
% Input:
%          s: numerical, threshold value 

d = size(S_XX, 1);
S_XX_new = S_XX;
for i = 1: d
    num = find(S_XX_new(i:d, i) < s) + i - 1;  % column data of lower triangle  
    if ~isempty(num) 
        S_XX_new(num(1) : d, i) = 0;
    end
    
    nrow = find(S_XX_new(i, 1: i) < s); % row data of lower triangle
    if ~isempty(nrow)
        S_XX_new(i, 1:nrow(end)) = 0;
    end
    
    num = find(S_XX_new(1:i, i) < s); % column data of up triangle
     if ~isempty(num) 
        S_XX_new(1 : num(end), i) = 0;
     end
    
     nrow = find(S_XX_new(i, i: d) < s) + i - 1; % row data of up triangle
    if ~isempty(nrow)
        S_XX_new(i, nrow(1) : d) = 0;
    end
     
end





