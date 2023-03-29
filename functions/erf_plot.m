function erf_figdata = erf_plot(erpdata_MEG,endIter,timeNew)


%% Grand Average N2pc on each hemisphere: MEG

avgLH_LVF = zeros(endIter,length(timeNew));
avgLH_RVF = zeros(endIter,length(timeNew));
avgRH_LVF = zeros(endIter,length(timeNew));
avgRH_RVF = zeros(endIter,length(timeNew));

for k = 1:endIter
    avgLH_LVF(k,:) = erpdata_MEG(k).avgLH_LVF;
    avgLH_RVF(k,:) = erpdata_MEG(k).avgLH_RVF;
    avgRH_LVF(k,:) = erpdata_MEG(k).avgRH_LVF;
    avgRH_RVF(k,:) = erpdata_MEG(k).avgRH_RVF;    
end


%left hem
GA_avgLH_LVF = mean(avgLH_LVF);
GA_avgLH_RVF = mean(avgLH_RVF);
GA_N2pcLH = GA_avgLH_LVF - GA_avgLH_RVF;

%right hem
GA_avgRH_LVF = mean(avgRH_LVF);
GA_avgRH_RVF = mean(avgRH_RVF);
GA_N2pcRH = GA_avgRH_LVF - GA_avgRH_RVF;

% concatenate
megcontra = mean(cat(1,avgLH_RVF,avgRH_LVF));
megipsi = mean(cat(1,avgRH_RVF,avgLH_LVF));
N2pc_conip = megcontra - megipsi;





%% GA plot MEG
GAleftFig = figure('Name','GA N2pc in Left Hemisphere');
hold on
plot(timeNew,GA_avgLH_LVF,'b');
plot(timeNew,GA_avgLH_RVF,'r--');
plot(timeNew,GA_N2pcLH,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA N2pc in Left Hemisphere')
legend('LVF Target','RVF Target','LVF-RVF')
yL = get(gca,'ylim');
ylim([-1.5e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');
% GAleftFigName = sprintf('GAleftHem');
% saveas(GAleftFig, fullfile(savePath,GAleftFigName),'png');

GArightFig = figure('Name','GA N2pc Right Hemisphere');
hold on
plot(timeNew,GA_avgRH_LVF,'b');
plot(timeNew,GA_avgRH_RVF,'r--');
plot(timeNew,GA_N2pcRH,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA N2pc Right Hemisphere')
legend('LVF Target','RVF Target','LVF-RVF')
yL = get(gca,'ylim');
ylim([-1.5e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');
% line([0.18 0.18],ylim,'Color','k','LineStyle','--');
% line([0.3 0.3],ylim,'Color','k','LineStyle','--');
% GArightFigName = sprintf('GArightHem'); 
% saveas(GArightFig, fullfile(savePath,GArightFigName),'png');

N2pcfig = figure('Name','GA N2pc contra-ipsi');
hold on
plot(timeNew,megcontra,'b');
plot(timeNew,megipsi,'r--');
plot(timeNew,N2pc_conip,'k','LineWidth',2);
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA N2pc contra-ipsi')
legend('contra','ipsi','conta-ipsi')
yL = get(gca,'ylim');
ylim([-1.5e-13 1.5e-13])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');

%% ERF data figure

%left hem
erf_figdata.GA_avgLeftHem_LVF = GA_avgLH_LVF;
erf_figdata.GA_avgLeftHem_RVF = GA_avgLH_RVF;
erf_figdata.GA_N2pcLH = GA_N2pcLH;

%right hem
erf_figdata.GA_avgRightHem_LVF = GA_avgRH_LVF;
erf_figdata.GA_avgRightHem_RVF = GA_avgRH_RVF;
erf_figdata.GA_N2pcRH = GA_N2pcRH;

end