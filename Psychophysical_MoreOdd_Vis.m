% Psychometric Function Group Odd -- mod for Visual Oddball data

clear;

sampler = 1; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Data

nsubjs = 10;
nCon = 4;

% n = dlmread('FP1_data_n.txt','\t'); %numtrials at each x cell
NN = zeros(nsubjs,9,nCon);
NN(:,:,1) = dlmread('vis_datafiles/dataFP1_n.txt','\t');
NN(:,:,2) = dlmread('vis_datafiles/dataFP2_n.txt','\t');
NN(:,:,3) = dlmread('vis_datafiles/dataFP3_n.txt','\t');
NN(:,:,4) = dlmread('vis_datafiles/dataFP4_n.txt','\t');

% r = dlmread('FP1_data_r.txt','\t'); %count of "longer" resps at each x cell
RR = zeros(nsubjs,9,nCon);
RR(:,:,1) = dlmread('vis_datafiles/dataFP1_r.txt','\t');
RR(:,:,2) = dlmread('vis_datafiles/dataFP2_r.txt','\t');
RR(:,:,3) = dlmread('vis_datafiles/dataFP3_r.txt','\t');
RR(:,:,4) = dlmread('vis_datafiles/dataFP4_r.txt','\t');

% x = dlmread('FP1_data_x.txt','\t'); %test durations
XX = zeros(nsubjs,9,nCon);
XX(:,:,1) = dlmread('vis_datafiles/data_x.txt','\t');
XX(:,:,2) = dlmread('vis_datafiles/data_x.txt','\t');
XX(:,:,3) = dlmread('vis_datafiles/data_x.txt','\t');
XX(:,:,4) = dlmread('vis_datafiles/data_x.txt','\t');

% rprop = dlmread('FP1_data_rprop.txt','\t'); %proportion of "longer" resps at each x level
RPR = zeros(nsubjs,9,nCon);
RPR(:,:,1) = dlmread('vis_datafiles/dataFP1_rprop.txt','\t');
RPR(:,:,2) = dlmread('vis_datafiles/dataFP2_rprop.txt','\t');
RPR(:,:,3) = dlmread('vis_datafiles/dataFP3_rprop.txt','\t');
RPR(:,:,4) = dlmread('vis_datafiles/dataFP4_rprop.txt','\t');



xmean = 350.*ones(1,nsubjs);    %center of oddball presentation levels by participant
nstim = 9.*ones(1,nsubjs);  %number of presentation levels by participant

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',XX,'n',NN,'r',RR,'xmean',xmean,'nstim',nstim,'nsubjs',nsubjs,'nCon',nCon);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.mua = zeros(1,4);
    S.mub = zeros(1,4);
    S.sigmaa = 1.*ones(1,1);
    S.sigmab = 1.*ones(1,1);
    S.alpha = -2 + 4.*rand(nsubjs,4);
    S.beta = 0.5.*rand(nsubjs,4);
    init0(i) = S;
end

if ~run_model
    load Psychophysical_GroupOdd samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'Psychophysical_GroupOdd.txt'), ...
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
            fullfile(pwd, 'Psychophysical_GroupOdd.txt'), ...
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
    save Psychophysical_GroupOdd samples stats
end;

%% Analysis
% Concatenate two chains -- Lee & Waggers
% alpha_all = squeeze([samples.alpha(1,:,:) samples.alpha(2,:,:)]);
% beta_all = squeeze([samples.beta(1,:,:) samples.beta(2,:,:)]);

% Concatenate two chains -- more transparent way
for s = 1:nsubjs
    for c = 1:nCon
        temp = samples.alpha(:,:,s,c);
        alphaSample(:,s,c) = temp(:);
    end;
end;

for s = 1:nsubjs
    for c = 1:nCon
        temp = samples.beta(:,:,s,c);
        betaSample(:,s,c) = temp(:);
    end;
end;

alpha_avg = stats.mean.alpha;
beta_avg = stats.mean.beta;

%%
% Construct JNDs
for s=1:nsubjs
    for c = 1:nCon
        JND(:,s,c) = psychfunc_inv(0.84,xmean(s),alphaSample(:,s,c),betaSample(:,s,c)) - psychfunc_inv(0.5,xmean(s),alphaSample(:,s,c),betaSample(:,s,c));
    end;
end;

%% WORKING HERE
% Construct Condition-level JNDs...but this won't let me do poster
% distribs. something needs to be different...
% for c=1:nCon
%     temp = alphaSample(:,:,c);
%     alphaSample_con(:,c) = temp(:);
% end;
% for c=1:nCon
%     temp = betaSample(:,:,c);
%     betaSample_con(:,c) = temp(:);
% end;
% 
% for c = 1:nCon
%     JND_group(:,c) = psychfunc_inv(0.84,xmean(1),alphaSample_con(:,c),betaSample_con(:,c)) - psychfunc_inv(0.5,xmean(1),alphaSample_con(:,c),betaSample_con(:,c));
% end;


% %% alternative JNDs: 25th to 75th percentiles
% % Construct JNDs
% for s=1:nsubjs
%     for c = 1:nCon
%         JND_alt(:,s,c) = psychfunc_inv(0.75,xmean(s),alphaSample(:,s,c),betaSample(:,s,c)) - psychfunc_inv(0.25,xmean(s),alphaSample(:,s,c),betaSample(:,s,c));
%     end;
% end;
% for JND_alt: plot post dense by subj
% for s = 1:nsubjs
%  plotPosteriorDensity(JND_alt(:,s), 'JND_alt', [4,6,s], 1, []);
% end;

% plot post dense by subj
figure()
for s = 1:nsubjs
 plotPosteriorDensity(JND(:,s), 'JND (ms)', [3,6,s], 1, []);
end;

%% plot post dense by condition
figure()
for c = 1:nCon
    temp = JND(:,:,c);
    plotPosteriorDensity(temp(:),'JND (ms)',[2,2,c],1,[]);
    
end;

%%
% Construct PSEs
for s = 1:nsubjs
    for c = 1:nCon
        PSE(:,s,c) = psychfunc_inv(0.5,xmean(s),alphaSample(:,s,c),betaSample(:,s,c));
    end;
end;

figure()
for s = 1:nsubjs
 plotPosteriorDensity(PSE(:,s), 'PSE', [3,6,s], 1, []);
end;

%% plot post dense by condition
figure()
for c = 1:nCon
    temp = PSE(:,:,c);
    plotPosteriorDensity(temp(:),'PSE (ms)',[2,2,c],1,[]);
    
end;


%% Plots with psychometric functions are called from psychfunc.m and psychfun_inv.m function-files
% this chunk plots rprop 
for c=1:nCon
    
figure(c);clf;hold on;
set(gcf,'units','norm','position',[.1 .1 .85 .7],'paperpositionmode','auto');
for s = 1:nsubjs
    subplot(3,6,s);cla;hold on;
    plot(XX(s,:,c),RPR(s,:,c),'s','markeredgecolor','k','markerfacecolor','w','markersize',5);
    set(gca,'xlim',[270 430],'ylim',[-0.1 1.1],'box','on','xtick',[280:35:420],'ytick',[0 .5 .84 1],'fontsize',10);
    th = title(['Subject ',num2str(s)]);
    set(th,'fontsize',12,'vert','mid');
end;
ah = gca;
axes('position',[0,0,1,1],'visible','off');
th = text(0.5,0.02,'Test Interval (ms)','fontsize',20,'horizontalalignment','center');
th = text(0.075,0.5,'Proportion of Long Responses','fontsize',20,'rotation',90,'horizontalalignment','center');
end

% fills in model curves on those rprop plots
for c=1:nCon
    
figure(c);hold on;
for s = 1:nsubjs
    subplot(3,6,s);hold on;
        xscale = XX(s,1):0.1:XX(s,nstim(s));
plot(xscale,psychfunc(xscale,xmean(s),alpha_avg(s,c),beta_avg(s,c)),'k');
    line([0 psychfunc_inv(0.5,xmean(s),alpha_avg(s,c),beta_avg(s,c))],[0.5 0.5],'color','k','linestyle',':');
    line([psychfunc_inv(0.5,xmean(s),alpha_avg(s,c),beta_avg(s,c)) psychfunc_inv(0.5,xmean(s),alpha_avg(s,c),beta_avg(s,c))],[0 0.5],'color','k','linestyle',':');
    line([0 psychfunc_inv(0.84,xmean(s),alpha_avg(s,c),beta_avg(s,c))],[0.84 0.84],'color','k','linestyle',':');
    line([psychfunc_inv(0.84,xmean(s),alpha_avg(s,c),beta_avg(s,c)) psychfunc_inv(0.84,xmean(s),alpha_avg(s,c),beta_avg(s,c))],[0 0.84],'color','k','linestyle',':');
    plot(XX(s,:,c),RPR(s,:,c),'s','markeredgecolor','k','markerfacecolor','w','markersize',5);
end;
ah = gca;
axes('position',[0,0,1,1],'visible','off');
th = text(0.5,0.02,'Test Interval (ms)','fontsize',20,'horizontalalignment','center');
th = text(0.075,0.5,'Proportion of Long Responses','fontsize',20,'rotation',90,'horizontalalignment','center');
end

% print -depsc ../../../Content/CaseStudies/PsychophysicalFunctions/Figures/Psychophysical_1.eps