function p = plotPosteriorDensity(sample, label, theSubPlot, plotDensity, theTitle)

if (~isempty(theSubPlot))
    subplot(theSubPlot(1),theSubPlot(2),theSubPlot(3))
    fs = 12;
else
    figure
    fs = 18;
end;
eps = (max(sample)-min(sample))./100;
bins=[min(sample)-eps:eps:max(sample)-eps];
count=hist(sample,bins);
if(plotDensity)
    density=count/sum(count)/eps;
    y=density;
    x = bins;
    p=plot(bins,density,'k-');
else
   [y,x]= hist(sample,25);
   p = bar(x,y,1);
end;
%axis([0,1.0, 0 max(count)])
%find the max
[theMax theMaxi] = max(count);
maxr = bins(theMaxi);
meanr = mean(sample);
xTheMaxi = find(x<=maxr,1,'last');
text(x(xTheMaxi),y(xTheMaxi),['M = ' num2str(round(100.*meanr)./100)],'BackgroundColor',[1 1 1],'FontSize',fs);
HDIlim = HDIofMCMC(sample, .95 );
hold on
plot(HDIlim,[0,0],'-r','LineWidth',6)
text(HDIlim(1),.5.*y(xTheMaxi),['HDIL: ' num2str(round(100.*HDIlim(1))./100)],'BackgroundColor',[1 1 1],'FontSize',fs);
text(HDIlim(2),.5.*y(xTheMaxi),['HDIU: ' num2str(round(100.*HDIlim(2))./100)],'BackgroundColor',[1 1 1],'FontSize',fs);
title(theTitle)
%find percent below 0
if (~plotDensity)
loc0 = find(bins<=0,1,'last');%find the last point at or below 0
propBelow0 = sum(count(1:loc0))./sum(count);
else
b0 = length(find(sample<=0));
propBelow0 = b0./length(sample);
end;
propAbove0 = 1- propBelow0;
zeroText = [int2str(round(100.*propBelow0)) '%<0<' int2str(round(100.*propAbove0)) '%'];
text(x(xTheMaxi),.75.*y(xTheMaxi),zeroText,'BackgroundColor',[1 1 1],'FontSize',fs);

set(gca,'box','on','fontsize',14);
xlabel(label,'fontsize',16);
if (plotDensity)
    ylabel('Posterior Density','fontsize',16);
else
    ylabel('Frequency','fontsize',16);
end;


