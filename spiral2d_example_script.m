%%% reVERSE for 2D spiral with and without measured GIRFs

%%% Shaihan Malik July 2015

%% constant definitions

gamma_uT = 267.5221;    % units rad s^-1 uT^-1
gamma_mT = gamma_uT * 1e3;
dt = 6.4e-6;        % sampling dwell time, 6.4us


%% Load in B0 & B1 field maps + FOV info (example data from 7T 8ch head coil)
load b0_b1_fields
% b0 = B0 map, uT
% tx = relative transmit sensitivity map
% X,Y,Z = spatial coordinates

[n1 n2 Nc] = size(tx);
xlist = double(cat(2,X(:),Y(:),Z(:)));

idx = find(m); %<- index of non-masked voxels

%% Define Square target

flip = 60 * pi/180; %<-- target flip angle, rad
P = zeros([n1 n2]);

%%% Offsets
xoff = +6e-3;
yoff = 0e-3;
zoff = 0;
r0 = 15e-3;

%%% Square beam
P = double(((abs(X-xoff)<r0)&(abs(Y-yoff)<r0)));


%%% apodize 
h = fspecial('gaussian',3,0.5);
P = imfilter(P,h);

%%% Now scale to flip
P=P*flip;


%% Define gradient correction model (GIRF)

girf = load('GIRF_3T_London_20140729');

%%% This correction gives G and k after convolution with GIRF
gcor = @(x)(gradient_distort_GIRF(x,girf.ff,girf.Hw,dt,10));

%%% This doesn't apply a girf
gcor0 = @(x)(gradient_distort_GIRF(x));



%% Example design: No GIRF (standard re-VERSE)
    
%%% Set default options
opt = reVERSE_init;
opt.show = 1;

%%% Set up function for system matrix
afun = @(k)(STA_system_matrix(xlist,k,opt.dt*(1:size(k,1)),tx,b0,m,'loopcalc'));

%%% Initial k-space
load K %<--- original demanded k-space (rad/m)

%%% run it
[bb,Gv] = reVERSE_GIRF(P(idx),K,afun,gcor0,opt);

%%% Outputs
rf_out = bb{end}; % mT
G_out  = Gv{end}; % mT/m

%% Example design: Include GIRF

%%% Set default options
opt = reVERSE_init;
opt.show = 1;

%%% Set up function for system matrix
afun = @(k)(STA_system_matrix(xlist,k,opt.dt*(1:size(k,1)),tx,b0,m,'loopcalc'));

%%% Initial k-space
load K %<--- original demanded k-space (rad/m)

%%% run it
[bb,Gv] = reVERSE_GIRF(P(idx),K,afun,gcor,opt);

%%% Outputs
rf_out = bb{end}; % mT
G_out  = Gv{end}; % mT/m