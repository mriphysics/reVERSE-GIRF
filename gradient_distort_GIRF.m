%% 13-2-2014: Gradient distortion with transfer function (GIRF)
% usage: [Gdistorted,ktx,krx] = gradient_distort_GIRF(Ginput,ff,H,dt,Npad)
%                   ff = frequencies for which H is defined
%                   H = transfer function
%                   dt = gradient sample dur (default 6.4us if not set)
%                   Npad = number of padding samples to add at start and
%                   end (default=0)
% Ginput can be MxN - i.e. can accept any number of inputs. H must either
% be Mfx1 in which case it is replicated N times for each gradient input
% (typically these will be the 3 axes) or it must be Mfx3 in which case it
% already includes x,y,z measurements
%
% EDIT 2-3-2015: Remove Padding before compute k-space.
% EDIT 6-3-2015: Explicitly normalize the matrices with df and dt. Note
% that the definition of the forward transform is F = dt*exp(2*pi*1i*ff'*t);
% The supplied GIRF must have been obtained in this manner.
% EDIT 21-4-2015: If only one input, just compute k-space, no distortion
%
% Shaihan Malik 2015  

function [Gdistorted,ktx,krx] = gradient_distort_GIRF(Ginput,ff,H,dt,Npad)

% If only one input argument, skip over girf and just return k-spaces
if nargin==1
    apply_girf = false;
else
    apply_girf = true;
end

if (~exist('dt','var'))||isempty(dt)
    dt=6.4e-6;
end
if ~exist('Npad','var')
    Npad=0;
end
Ngradients = size(Ginput,2);


%%% add some zeros on front and end
G = cat(1,zeros([Npad Ngradients]),Ginput,zeros([Npad Ngradients]));
M = length(G);
t = (0:M-1)*dt;

if apply_girf

    %%% check if frequencies sufficient
    df_req = 1/t(end);
    df_act = ff(2)-ff(1);
    if df_req<df_act
        % New freq range
        ffnew = ff(1):df_req:ff(end);
        ffnew=ffnew-min(abs(ffnew));% centre so we have a zero sample
        % centring may have introduced one frequency out of range
        ffnew(1)=[];
        % resample H
        Hnew = interp1(ff,H,ffnew);
        
        %%% swap
        ff=ffnew;
        H=Hnew;
        if size(H,1)==1
            H=H(:);
        end
        fprintf(1,'GIRF: resampled because frequency resolution was insufficient\n');
        df_act = ff(2)-ff(1);
    end
    
    
    %%% construct Fourier matrix
    F = dt*exp(2*pi*1i*ff'*t);
    % Also define inverse
    Fi = df_act*exp(-2*pi*1i*t'*ff);
    
    %%% Apply transfer function
    Gdistorted = zeros([M Ngradients]);
    %%% check if H has all axes
    if size(H,2)==1
        H=repmat(H,Ngradients);
    end
    if size(H,2)~=Ngradients
        fprintf(1,'Error: GIRF has different number of axes to gradient\n')
        return
    end
    
    %%% Apply H
    Gdistorted = real(Fi*(H.*(F*G)));
   
else
    % If only one input specified, just re-output the gradient
    Gdistorted=G;
end

%%% remove padding
Gdistorted([1:Npad (M-Npad+1):M],:)=[];
M = length(Gdistorted);
t = (0:M-1)*dt;

%%% Now make k-spaces
gamma_mT = 2*pi*42577.46778; % rad s^-1 mT^-1

%%% Tx - integrate to end
ktx = -gamma_mT*flipud(cumtrapz(t,flipud(Gdistorted)));
%%% Rx - simple integration as you go along
krx = gamma_mT*cumtrapz(t,Gdistorted);


end


