% JNDs: plot post dense by subj
for s = 1:nsubjs
 plotPosteriorDensity(JND(:,s), 'JND', [3,6,s], 1, []);
end;

%% JNDs: plot post dense by condition
for c = 1:nCon
    temp = JND(:,:,c);
    plotPosteriorDensity(temp(:),'JND',[2,2,c],1,[]);
    
end;

%% PSEs: plot post dense by subject
for s = 1:nsubjs
 plotPosteriorDensity(PSE(:,s), 'PSE', [3,6,s], 1, []);
end;

%% PSEs: plot post dense by condition
for c = 1:nCon
    temp = PSE(:,:,c);
    plotPosteriorDensity(temp(:),'PSE',[2,2,c],1,[]);
    
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

%% Plot pcurves by condition
% NEED TO MODIFY TO AVERAGE OVER SUBJECTS SOMEHOW!

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


% 
% 
%% plotting from Lee&Wag retention example

% %% Data
% t = x
% % ns = nsubjs;
% k = [18    18    16    13     9     6     4     4     4 nan;
% 	17    13     9     6     4     4     4     4     4 nan;
% 	14    10     6     4     4     4     4     4     4 nan;
% 	nan   nan   nan   nan   nan   nan   nan    nan   nan nan];
% nt = length(t);
% n = 18;

% Draw Posterior Predictive Analysis
figure(6);clf;
sc=20; % Scaling Constant for Drawing Boxes
for i=1:ns
	% Subplots for Subjects
	subplot(2,2,i);hold on;
	% Plot Subject Data
	ph=plot([1:nt],k(i,:),'k-');
	set(ph,'linewidth',2);
	% Plot Posterior Predictive
	for j=1:nt
		count=hist(samples.predk(1,:,i,j),[0:n]);
		count=count/sum(count);
		for x=0:n
			if count(x+1)>0
				ph=plot(j,x,'ks');
				set(ph,'markersize',sc*sqrt(count(x+1)));
				if k(i,j)==x
					set(ph,'markerfacecolor','k');
				end;
			end;
		end;
	end;
	% Set the Axes
	axis([0 nt+1 -1 19]);
	% Title the Subplot
	th=title(['Subject ' int2str(i)]);
	set(th,'fontsize',12,'verticalalignment','mid');
	xlabel('Time Lags','fontsize',12);
	ylabel('Retention Count','fontsize',12);
	% Tidy Up the Subplot
	set(gca,'box','on','xtick',[1:nt],'xticklabel',t,'ytick',[0 18]);
end;

