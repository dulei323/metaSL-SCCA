close all; clear; clc;
% -----------------------------------
% Author: Duo Xi, xiduo@mail.nwpu.edu.cn
% Date: 28-Apr-2023
% -----------------------------------

%% Load data
data = load('data_using.mat');

%% metaSL_SCCA
s = 0.2;
S_XX = updata_XX(data.S_XX, s);

paras.lambda = [1, 1];
paras.r = [0.1, 0.1];

[res.u, res.v, res.res_iter] = metaSL_SCCA(S_XX, data.Beta, data.S_YY, paras);
res.sort_info = perf_process(res.u, res.v, data); 
