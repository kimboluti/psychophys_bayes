% Psychometric Function 1

clear;

sampler = 1; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Data
datarow = 1;
n = dlmread('FP1_data_n.txt','\t'); %numtrials at each x cell
r = dlmread('FP1_data_r.txt','\t'); %count of "longer" resps at each x cell
x = dlmread('FP1_data_x.txt','\t'); %test durations
rprop = dlmread('FP1_data_rprop.txt','\t'); %proportion of "longer" resps at each x level

%odd data shortcut

xmean = 350.*ones(1,18);    %center of oddball presentation levels by participant
nstim = 9.*ones(1,18);  %number of presentation levels by participant
nsubjs = 18;


%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n,'r',r,'xmean',xmean,'nstim',nstim,'nsubjs',nsubjs);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.mua = 0;
    S.mub = 0;
    S.sigmaa = 1;
    S.sigmab = 1;
    S.alpha = -2 + 4.*rand(1,nsubjs);
    S.beta = 0.5.*rand(1,nsubjs);
    init0(i) = S;
end

if ~run_model
    load Psychophysical_1 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'Psychophysical_1.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorparams',{'alpha','beta','mua','mub'},...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'Psychophysical_1.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',{'alpha','beta','mua','mub'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save Psychophysical_1 samples stats
end;

%% Analysis
% Concatenate two chains -- Lee & Waggers
% alpha_all = squeeze([samples.alpha(1,:,:) samples.alpha(2,:,:)]);
% beta_all = squeeze([samples.beta(1,:,:) samples.beta(2,:,:)]);

% Concatenate two chains -- more transparent way
for s = 1:nsubjs
    temp = samples.alpha(:,:,s);
    alphaSample(:,s) = temp(:);
end;

for s = 1:nsubjs
    temp = samples.beta(:,:,s);
    betaSample(:,s) = temp(:);
end;

alpha_avg = stats.mean.alpha;
beta_avg = stats.mean.beta;

% alpha_range = zeros(nsubjs,20);
% beta_range= zeros(nsubjs,20);
% 
% % for each subj, randomly sample 20x from alpha/beta samples for that subj
% for s=1:nsubjs
%     alpha_range(s,:)=randsample(alpha_all(:,s),20);
%     beta_range(s,:)=randsample(beta_all(:,s),20);
% end

%%
% Construct JNDs
for s=1:nsubjs
    JND(:,s) = psychfunc_inv(0.84,xmean(s),alphaSample(:,s),betaSample(:,s)) - psychfunc_inv(0.5,xmean(s),alphaSample(:,s),betaSample(:,s));
end

for s = 1:nsubjs
 plotPosteriorDensity(JND(:,s), 'JND', [2,4,s], 1, []);
end;

% Construct PSEs
for s=1:nsubjs
    PSE(:,s) = psychfunc_inv(0.5,xmean(s),alphaSample(:,s),betaSample(:,s));
end
figure
for s = 1:nsubjs
 plotPosteriorDensity(PSE(:,s), 'PSE', [4,6,s], 1, []);
end;

% %%
% % Plots with psychometric functions are called from psychfunc.m and psychfun_inv.m function-files
% figure(1);clf;hold on;
% set(gcf,'units','norm','position',[.1 .1 .85 .7],'paperpositionmode','auto');
% for s = 1:nsubjs
%     subplot(2,4,s);cla;hold on;
%     plot(x(s,:),rprop(s,:),'s','markeredgecolor','k','markerfacecolor','w','markersize',5);
%     set(gca,'xlim',[190 410],'ylim',[-0.1 1.1],'box','on','xtick',[200:50:400],'ytick',[0 .5 .84 1],'fontsize',18);
%     th = title(['Subject ',num2str(s)]);
%     set(th,'fontsize',18,'vert','mid');
% end;
% ah = gca;
% axes('position',[0,0,1,1],'visible','off');
% th = text(0.5,0.02,'Test Interval (ms)','fontsize',20,'horizontalalignment','center');
% th = text(0.075,0.5,'Proportion of Long Responses','fontsize',20,'rotation',90,'horizontalalignment','center');
% 
% for s = 1:nsubjs
%     subplot(2,4,s);hold on;
%         xscale = x(s,1):0.1:x(s,nstim(s));
% plot(xscale,psychfunc(xscale,xmean(s),alpha_avg(s),beta_avg(s)),'k');
%     line([0 psychfunc_inv(0.5,xmean(s),alpha_avg(s),beta_avg(s))],[0.5 0.5],'color','k','linestyle',':');
%     line([psychfunc_inv(0.5,xmean(s),alpha_avg(s),beta_avg(s)) psychfunc_inv(0.5,xmean(s),alpha_avg(s),beta_avg(s))],[0 0.5],'color','k','linestyle',':');
%     line([0 psychfunc_inv(0.84,xmean(s),alpha_avg(s),beta_avg(s))],[0.84 0.84],'color','k','linestyle',':');
%     line([psychfunc_inv(0.84,xmean(s),alpha_avg(s),beta_avg(s)) psychfunc_inv(0.84,xmean(s),alpha_avg(s),beta_avg(s))],[0 0.84],'color','k','linestyle',':');
%     plot(x(s,:),rprop(s,:),'s','markeredgecolor','k','markerfacecolor','w','markersize',5);
% end;
% ah = gca;
% axes('position',[0,0,1,1],'visible','off');
% th = text(0.5,0.02,'Test Interval (ms)','fontsize',20,'horizontalalignment','center');
% th = text(0.075,0.5,'Proportion of Long Responses','fontsize',20,'rotation',90,'horizontalalignment','center');
% % 
% % 
% print -depsc ../../../Content/CaseStudies/PsychophysicalFunctions/Figures/Psychophysical_1.eps