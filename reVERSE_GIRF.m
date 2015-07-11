%%% [bb,Gv] = reVERSE_GIRF(mTarget,k0,Afun,kfun,options)
%%% mTarget = target magnetization (STA)
%%% k0      = Initial k-space, mT/m
%%% Afun    = function handle to return A given kspace
%%% kfun    = function handle to return updated G and k given input G (could be GIRF)
%%% options = everything else, reVERSE_init() sets up default options
%%%
%%% Shaihan Malik, 2015

function [bb,Gv] = reVERSE_GIRF(mTarget,k0,Afun,kfun,options)

%%% Translate everything to CGS units for Lustig code
Kcm = k0 / (2*pi*100);% k-space in 1/cm
Gmax = options.Gmax / 10; % G/cm
Smax = options.Smax / 10; % G/cm/ms

%%% Options are specified in CGS units, but RF pulse design uses SI. Keep
%%% some limits here in SI units for plotting
b1_limit_mT = options.b1_limit * 0.1;

%%% First run time optimal code to get starting gradient waveform
% Need g (gradient) and krx: this is the Receive k-space, which we use when
% calling the time optimal code
fprintf(1,'Running initial gradient design ...\n');
[~,~,G,~,krx, ~, ~, ~] = minTimeGradient(Kcm,0,0,0,Gmax,Smax,options.dt*1e3,[],[]);
fprintf(1,'... done. Now being iterations...\n');

%%% Apply distortion model (this needs mT/m to give k-space in correct units)
[~,kedd] = kfun(G*10);

%%% Initialize re-VERSE
Kv = kedd;%<-- distorted version
Mv=length(Kv);
tv=(1:Mv)*options.dt;

% variables
bb = {};
bv = {};
Gv = {};

% 1st gradient is already computed: just convert to mT/m
Gv{1}=G*10;

for ii = 1:options.Nstop
    %%% Regularized STA design
    fprintf(1,'iteration %d: define STA system matrix\n',ii);
    A = Afun(Kv);
  
    fprintf(1,'iteration %d: Solve LLS problem\n',ii);
    
    b = lsqrSOL(size(A,1),size(A,2),A,mTarget,options.lambda,[],[],[],50,0);
    
    % Reshape
    b = reshape(b,[Mv options.Nc]);
  
    % save 
    bb{ii}=b;

    if options.show
        % Display solution
        figure(1);clf
        nr=3;nc=1;
        
        sz=get(0,'screenSize');sz=sz([3 4]);
        set(gcf,'position',[0.2*sz 0.4*sz(2) 0.6*sz(2)])
        subplot(nr,nc,1)
        plot(tv*1e3,max(abs(bb{ii}),[],2),'r');grid on;hold on
        if ii>1
            plot(tv*1e3,max(abs(bv{ii-1}),[],2),'k');grid on;hold on
        end
        xlabel('t/ms');title('RF amplitudes')
        ylabel('B_1/ uT');set(gca,'yticklabel',1e3*get(gca,'ytick'));%<- set axis labels in uT
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        patch([xl(1) xl(1) xl(2) xl(2)],b1_limit_mT*[options.b1_alpha 1 1 options.b1_alpha],[0 1 0],'facealpha',0.1,'edgealpha',0)
        patch([xl(1) xl(1) xl(2) xl(2)],[b1_limit_mT yl(2) yl(2) b1_limit_mT],[1 0 0],'facealpha',0.1,'edgealpha',0)
        ylim(yl);% re-enforce limit
        title(sprintf('Iteration %d',ii))
        
        subplot(nr,nc,2)
        plot((1:length(Gv{1}))*options.dt*1e3,Gv{1},'k');grid on
        xlabel('t/ms');title('Gradients')
        ylabel('G / mT/m')
        if ii>1
            hold on
            plot(tv*1e3,Gv{ii},'r')
            hold off
        end
        
        subplot(nr,nc,3)
        hold on
        for jj=1:ii
            plot(jj,1e3*max(abs(bb{jj}(:))),'*','markersize',15,'linewidth',2)
        end
        grid on
        xlabel('iteration number')
        ylabel('peak B_1 /uT')
        title('Max B_1 per iteration')
        drawnow
    end
    
    % VERSE: check convergence (comparison is in Gauss, hence x10)
    if (ii==options.Nstop)||(max(abs(b(:))*10)<options.b1_limit)
        if ii==options.Nstop
            outstr = 'max iterations were reached';
        else
            outstr = 'specified B1 limit was achieved';
        end
          
        fprintf(1,'reVERSE finished because %s\n',outstr);
        break
    end
    
    fprintf(1,'iteration %d: VERSE\n',ii);
    verse_in = [];
    verse_in.b = b*10;%<-- mT->Gauss
    verse_in.Gmod = sum((Gv{ii}/10).^2,2).^0.5;
    verse_in.bmax = options.b1_limit*options.b1_alpha; % (Gauss)
    verse_in.os = 15; %<--- oversample factor
    
    %%% Call VERSE design:
    % Function returns new gradient, new RF in verse_out and new krx. krx
    % is the receive k-space of the current gradient waveform, and this is
    % what is fed into the code as the 'k-space' curve to follow, not the
    % reverse integral that is the TX k-space
    [~,~,gv,~,krx, ~, ~, ~,~,verse_out] = minTimeGradient_VERSE(krx,0,0,Gmax,Smax,options.dt*1e3,[],[],verse_in);

    %%% get bv from additional outputs
    bv{ii} = verse_out.bt/10;% back to mT

    %%% Get Tx k-space from gcor: TX k-space is used for next pulse design iteration
    Gv{ii+1} = 10*gv; %<- convert back to mT/m
    Mv=length(Gv{ii+1});
    [~,kedd] = kfun(Gv{ii+1});
    tv = (1:Mv)*options.dt;
    Kv = kedd;
   
end



end