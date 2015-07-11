%% 22-4-15: Initialize options
function opt = reVERSE_init

opt = struct;

%%% Lambda value (regularization)
opt.lambda = 1;

%%% Number of channels
opt.Nc = 8;

%%% B1 limit for VERSE
opt.b1_limit = 12e-3 * 10;  %<--- limit in Gauss
opt.b1_alpha = 0.9;         %<--- At each iteration of VERSE we limit to b1_limit*b1_alpha
opt.bmax     = opt.b1_limit * opt.b1_alpha;  %<--- overall limit is limit*alpha


%%% Max iterations for reVERSE
opt.Nstop = 10;

%%% Sampling dwell, sec
opt.dt = 6.4e-6;
%%% Max and slew rate limits
opt.Gmax = 30;          % mT/m
opt.Smax = 180;         % T/m/s

%%% Display information
opt.show = false;
end