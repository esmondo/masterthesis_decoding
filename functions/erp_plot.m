function erp_figdata = erp_plot(erpdata_EEG,endIter,timeNew)

%% Grand Average N2pc from contralateral and ipsilateral signals: EEG

GA_P3_LVF = zeros(endIter,length(timeNew));
GA_P3_RVF = zeros(endIter,length(timeNew));
GA_P4_LVF = zeros(endIter,length(timeNew));
GA_P4_RVF = zeros(endIter,length(timeNew));
GA_Contra1 = zeros(endIter,length(timeNew));
GA_Ipsi1 = zeros(endIter,length(timeNew));

GA_O9_LVF = zeros(endIter,length(timeNew));
GA_O9_RVF = zeros(endIter,length(timeNew));
GA_O10_LVF = zeros(endIter,length(timeNew));
GA_O10_RVF = zeros(endIter,length(timeNew));
GA_Contra2 = zeros(endIter,length(timeNew));
GA_Ipsi2 = zeros(endIter,length(timeNew));

GA_P7_LVF = zeros(endIter,length(timeNew));
GA_P7_RVF = zeros(endIter,length(timeNew));
GA_P8_LVF = zeros(endIter,length(timeNew));
GA_P8_RVF = zeros(endIter,length(timeNew));
GA_Contra3 = zeros(endIter,length(timeNew));
GA_Ipsi3 = zeros(endIter,length(timeNew));

GA_CP1_LVF = zeros(endIter,length(timeNew));
GA_CP1_RVF = zeros(endIter,length(timeNew));
GA_CP2_LVF = zeros(endIter,length(timeNew));
GA_CP2_RVF = zeros(endIter,length(timeNew));
GA_Contra4 = zeros(endIter,length(timeNew));
GA_Ipsi4 = zeros(endIter,length(timeNew));

GA_PO3_LVF = zeros(endIter,length(timeNew));
GA_PO3_RVF = zeros(endIter,length(timeNew));
GA_PO4_LVF = zeros(endIter,length(timeNew));
GA_PO4_RVF = zeros(endIter,length(timeNew));
GA_Contra5 = zeros(endIter,length(timeNew));
GA_Ipsi5 = zeros(endIter,length(timeNew));

GA_PO7_LVF = zeros(endIter,length(timeNew));
GA_PO7_RVF = zeros(endIter,length(timeNew));
GA_PO8_LVF = zeros(endIter,length(timeNew));
GA_PO8_RVF = zeros(endIter,length(timeNew));
GA_Contra6 = zeros(endIter,length(timeNew));
GA_Ipsi6 = zeros(endIter,length(timeNew));




for k = 1:endIter
    
    GA_P3_LVF(k,:) = erpdata_EEG(k).avgLH_LVF(1,:);
    GA_P3_RVF(k,:) = erpdata_EEG(k).avgLH_RVF(1,:);
    GA_P4_LVF(k,:) = erpdata_EEG(k).avgRH_LVF(1,:);
    GA_P4_RVF(k,:) = erpdata_EEG(k).avgRH_RVF(1,:);
    GA_Contra1(k,:) = erpdata_EEG(k).avgcontra(1,:);
    GA_Ipsi1(k,:) = erpdata_EEG(k).avgipsi(1,:);
     
    GA_O9_LVF(k,:) = erpdata_EEG(k).avgLH_LVF(2,:);
    GA_O9_RVF(k,:) = erpdata_EEG(k).avgLH_RVF(2,:);
    GA_O10_LVF(k,:) = erpdata_EEG(k).avgRH_LVF(2,:);
    GA_O10_RVF(k,:) = erpdata_EEG(k).avgRH_RVF(2,:);
    GA_Contra2(k,:) = erpdata_EEG(k).avgcontra(2,:);
    GA_Ipsi2(k,:) = erpdata_EEG(k).avgipsi(2,:);
    
    GA_P7_LVF(k,:) = erpdata_EEG(k).avgLH_LVF(3,:);
    GA_P7_RVF(k,:) = erpdata_EEG(k).avgLH_RVF(3,:);
    GA_P8_LVF(k,:) = erpdata_EEG(k).avgRH_LVF(3,:);
    GA_P8_RVF(k,:) = erpdata_EEG(k).avgRH_RVF(3,:);
    GA_Contra3(k,:) = erpdata_EEG(k).avgcontra(3,:);
    GA_Ipsi3(k,:) = erpdata_EEG(k).avgipsi(3,:);
    
    GA_CP1_LVF(k,:) = erpdata_EEG(k).avgLH_LVF(4,:);
    GA_CP1_RVF(k,:) = erpdata_EEG(k).avgLH_RVF(4,:);
    GA_CP2_LVF(k,:) = erpdata_EEG(k).avgRH_LVF(4,:);
    GA_CP2_RVF(k,:) = erpdata_EEG(k).avgRH_RVF(4,:);
    GA_Contra4(k,:) = erpdata_EEG(k).avgcontra(4,:);
    GA_Ipsi4(k,:) = erpdata_EEG(k).avgipsi(4,:);
    
    GA_PO3_LVF(k,:) = erpdata_EEG(k).avgLH_LVF(5,:);
    GA_PO3_RVF(k,:) = erpdata_EEG(k).avgLH_RVF(5,:);
    GA_PO4_LVF(k,:) = erpdata_EEG(k).avgRH_LVF(5,:);
    GA_PO4_RVF(k,:) = erpdata_EEG(k).avgRH_RVF(5,:);
    GA_Contra5(k,:) = erpdata_EEG(k).avgcontra(5,:);
    GA_Ipsi5(k,:) = erpdata_EEG(k).avgipsi(5,:);
    
    GA_PO7_LVF(k,:) = erpdata_EEG(k).avgLH_LVF(6,:);
    GA_PO7_RVF(k,:) = erpdata_EEG(k).avgLH_RVF(6,:);
    GA_PO8_LVF(k,:) = erpdata_EEG(k).avgRH_LVF(6,:);
    GA_PO8_RVF(k,:) = erpdata_EEG(k).avgRH_RVF(6,:);
    GA_Contra6(k,:) = erpdata_EEG(k).avgcontra(6,:);
    GA_Ipsi6(k,:) = erpdata_EEG(k).avgipsi(6,:);
    
   
           
end


GA_P3_LVF = mean(GA_P3_LVF);
GA_P3_RVF = mean(GA_P3_RVF);
GA_P4_LVF = mean(GA_P4_LVF);
GA_P4_RVF = mean(GA_P4_RVF);
GA_Contra1 = mean(GA_Contra1);
GA_Ipsi1 = mean(GA_Ipsi1);
GA_ConIpsi1 = GA_Contra1 - GA_Ipsi1;

GA_O9_LVF = mean(GA_O9_LVF);
GA_O9_RVF = mean(GA_O9_RVF);
GA_O10_LVF = mean(GA_O10_LVF);
GA_O10_RVF = mean(GA_O10_RVF);
GA_Contra2 = mean(GA_Contra2);
GA_Ipsi2 = mean(GA_Ipsi2);
GA_ConIpsi2 = GA_Contra2 - GA_Ipsi2;

GA_P7_LVF = mean(GA_P7_LVF);
GA_P7_RVF = mean(GA_P7_RVF);
GA_P8_LVF = mean(GA_P8_LVF);
GA_P8_RVF = mean(GA_P8_RVF);
GA_Contra3 = mean(GA_Contra3);
GA_Ipsi3 = mean(GA_Ipsi3);
GA_ConIpsi3 = GA_Contra3 - GA_Ipsi3;

GA_CP1_LVF = mean(GA_CP1_LVF);
GA_CP1_RVF = mean(GA_CP1_RVF);
GA_CP2_LVF = mean(GA_CP2_LVF);
GA_CP2_RVF = mean(GA_CP2_RVF);
GA_Contra4 = mean(GA_Contra4);
GA_Ipsi4 = mean(GA_Ipsi4);
GA_ConIpsi4 = GA_Contra4 - GA_Ipsi4;

GA_PO3_LVF = mean(GA_PO3_LVF);
GA_PO3_RVF = mean(GA_PO3_RVF);
GA_PO4_LVF = mean(GA_PO4_LVF);
GA_PO4_RVF = mean(GA_PO4_RVF);
GA_Contra5 = mean(GA_Contra5);
GA_Ipsi5 = mean(GA_Ipsi5);
GA_ConIpsi5 = GA_Contra5 - GA_Ipsi5;

GA_PO7_LVF = mean(GA_PO7_LVF);
GA_PO7_RVF = mean(GA_PO7_RVF);
GA_PO8_LVF = mean(GA_PO8_LVF);
GA_PO8_RVF = mean(GA_PO8_RVF);
GA_Contra6 = mean(GA_Contra6);
GA_Ipsi6 = mean(GA_Ipsi6);
GA_ConIpsi6 = GA_Contra6 - GA_Ipsi6;




% combination (only sensors which I think showing high N2pc amplitude)
LH_comb_LVF = mean(cat(1,GA_P7_LVF,GA_PO7_LVF),1);
LH_comb_RVF = mean(cat(1,GA_P7_RVF,GA_PO7_RVF),1);
RH_comb_LVF = mean(cat(1,GA_P8_LVF,GA_PO8_LVF),1);
RH_comb_RVF = mean(cat(1,GA_P8_RVF,GA_PO8_RVF),1);
comb_contra = mean(cat(1,GA_Contra2,GA_Contra4),1);
comb_ipsi = mean(cat(1,GA_Ipsi2,GA_Ipsi4),1);
comb_conipsi = mean(cat(1,GA_ConIpsi2,GA_ConIpsi4),1);

%% GA N2pc 1: P3 & P4

% GA_P3 = figure('Name','GA P3 / Left Hemisphere');
% hold on
% plot(timeNew,GA_P3_LVF)
% plot(timeNew,GA_P3_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA P3 / Left Hemisphere')
% legend('P3 LVF','P3 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% GA_P4 = figure('Name','GA P4 / Right Hemisphere');
% hold on
% plot(timeNew,GA_P4_LVF)
% plot(timeNew,GA_P4_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA P4 / Right Hemisphere')
% legend('P4 LVF','P4 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');

N2pc_1 = figure('Name','GA1 Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_Contra1) % contralateral
plot(timeNew,GA_Ipsi1) % ipsilateral
plot(timeNew,GA_ConIpsi1,'k','LineWidth',2)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA1 Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');


%% GA N2pc 2: O9 & O10

% GA_O9 = figure('Name','GA O9 / Left Hemisphere');
% hold on
% plot(timeNew,GA_O9_LVF)
% plot(timeNew,GA_O9_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA O9 / Left Hemisphere')
% legend('O9 LVF','O9 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% GA_O10 = figure('Name','GA P10 / Right Hemisphere');
% hold on
% plot(timeNew,GA_O10_LVF)
% plot(timeNew,GA_O10_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA O10 / Right Hemisphere')
% legend('O10 LVF','O10 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');

N2pc_2 = figure('Name','GA2 Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_Contra2) % contralateral
plot(timeNew,GA_Ipsi2) % ipsilateral
plot(timeNew,GA_ConIpsi2,'k','LineWidth',2)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA2 Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');


%% GA N2pc 3: P7 & P8

% GA_P7 = figure('Name','GA P7 / Left Hemisphere');
% hold on
% plot(timeNew,GA_P7_LVF)
% plot(timeNew,GA_P7_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA P7 / Left Hemisphere')
% legend('P7 LVF','P7 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% GA_P8 = figure('Name','GA P8 / Right Hemisphere');
% hold on
% plot(timeNew,GA_P8_LVF)
% plot(timeNew,GA_P8_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA P8 / Right Hemisphere')
% legend('P8 LVF','P8 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');

N2pc_3 = figure('Name','GA3 Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_Contra3) % contralateral
plot(timeNew,GA_Ipsi3) % ipsilateral
plot(timeNew,GA_ConIpsi3,'k','LineWidth',2)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA3 Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');


%% GA N2pc 4: CP1 & CP2

% GA_CP1 = figure('Name','GA CP1 / Left Hemisphere');
% hold on
% plot(timeNew,GA_CP1_LVF)
% plot(timeNew,GA_CP1_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA CP1 / Left Hemisphere')
% legend('CP1 LVF','CP1 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% GA_CP2 = figure('Name','GA CP2 / Right Hemisphere');
% hold on
% plot(timeNew,GA_CP2_LVF)
% plot(timeNew,GA_CP2_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA CP2 / Right Hemisphere')
% legend('CP2 LVF','CP2 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');

N2pc_4 = figure('Name','GA4 Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_Contra4) % contralateral
plot(timeNew,GA_Ipsi4) % ipsilateral
plot(timeNew,GA_ConIpsi4,'k','LineWidth',2)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA4 Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');

%% GA N2pc 5: PO3 & PO4

% GA_PO3 = figure('Name','GA PO3 / Left Hemisphere');
% hold on
% plot(timeNew,GA_PO3_LVF)
% plot(timeNew,GA_PO3_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA PO3 / Left Hemisphere')
% legend('PO3 LVF','PO3 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% GA_PO4 = figure('Name','GA OP4 / Right Hemisphere');
% hold on
% plot(timeNew,GA_PO4_LVF)
% plot(timeNew,GA_PO4_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA PO4 / Right Hemisphere')
% legend('PO4 LVF','PO4 RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');

N2pc_5 = figure('Name','GA5 Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_Contra5) % contralateral
plot(timeNew,GA_Ipsi5) % ipsilateral
plot(timeNew,GA_ConIpsi5,'k','LineWidth',2)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA5 Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');

%% GA N2pc 6: PO7 & PO8

GA_PO7 = figure('Name','GA PO7 / Left Hemisphere');
hold on
plot(timeNew,GA_PO7_LVF)
plot(timeNew,GA_PO7_RVF)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA PO7 / Left Hemisphere')
legend('PO7 LVF','PO7 RVF')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');

GA_PO8 = figure('Name','GA OP8 / Right Hemisphere');
hold on
plot(timeNew,GA_PO8_LVF)
plot(timeNew,GA_PO8_RVF)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA PO8 / Right Hemisphere')
legend('PO8 LVF','PO8 RVF')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');

N2pc_6 = figure('Name','GA6 Contralateral and Ipsilateral');
hold on
plot(timeNew,GA_Contra6) % contralateral
plot(timeNew,GA_Ipsi6) % ipsilateral
plot(timeNew,GA_ConIpsi6,'k','LineWidth',2)
hold off
xlabel('time')
ylabel('amplitude')
grid on
title('GA6 Contralateral and Ipsilateral')
legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
yL = get(gca,'ylim');
xlim([-0.1 0.5])
ylim([-4.5e-6 4e-6])
line(xlim,[0 0],'Color','k')
line([0 0],ylim,'Color','k');

%% comb plot

% comb_LH = figure('Name','GA comb LH / Left Hemisphere');
% hold on
% plot(timeNew,LH_comb_LVF)
% plot(timeNew,LH_comb_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA comb LH / Left Hemisphere')
% legend('comb LH LVF','comb LH RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% comb_RH = figure('Name','GA comb RH / Right Hemisphere');
% hold on
% plot(timeNew,RH_comb_LVF)
% plot(timeNew,RH_comb_RVF)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA comb RH / Right Hemisphere')
% legend('comb RH LVF','comb RH RVF')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');
% 
% N2pc_comb = figure('Name','GA comb Contralateral and Ipsilateral');
% hold on
% plot(timeNew,comb_contra) % contralateral
% plot(timeNew,comb_ipsi) % ipsilateral
% plot(timeNew,comb_conipsi,'k','LineWidth',2)
% hold off
% xlabel('time')
% ylabel('amplitude')
% grid on
% title('GA Contralateral and Ipsilateral')
% legend('Contralateral Target','Ipsilateral Target','Contralateral-Ipsilateral')
% yL = get(gca,'ylim');
% ylim([-5.5e-6 4e-6])
% line(xlim,[0 0],'Color','k')
% line([0 0],ylim,'Color','k');

%% figure data

erp_figdata.GA_P3_LVF = GA_P3_LVF;
erp_figdata.GA_P3_RVF = GA_P3_RVF;
erp_figdata.GA_P4_LVF = GA_P4_LVF;
erp_figdata.GA_P4_RVF = GA_P4_RVF;
erp_figdata.GA_Contra1 = GA_Contra1;
erp_figdata.GA_Ipsi1 = GA_Ipsi1;
erp_figdata.GA_ConIpsi1 = GA_ConIpsi1;

erp_figdata.GA_O9_LVF = GA_O9_LVF;
erp_figdata.GA_O9_RVF = GA_O9_RVF;
erp_figdata.GA_O10_LVF = GA_O10_LVF;
erp_figdata.GA_O10_RVF = GA_O10_RVF;
erp_figdata.GA_Contra2 = GA_Contra2;
erp_figdata.GA_Ipsi2 = GA_Ipsi2;
erp_figdata.GA_ConIpsi2 = GA_ConIpsi2;

erp_figdata.GA_P7_LVF = GA_P7_LVF;
erp_figdata.GA_P7_RVF = GA_P7_RVF;
erp_figdata.GA_P8_LVF = GA_P8_LVF;
erp_figdata.GA_P8_RVF = GA_P8_RVF;
erp_figdata.GA_Contra3 = GA_Contra3;
erp_figdata.GA_Ipsi3 = GA_Ipsi3;
erp_figdata.GA_ConIpsi3 = GA_ConIpsi3;

erp_figdata.GA_CP1_LVF = GA_CP1_LVF;
erp_figdata.GA_CP1_RVF = GA_CP1_RVF;
erp_figdata.GA_CP2_LVF = GA_CP2_LVF;
erp_figdata.GA_CP2_RVF = GA_CP2_RVF;
erp_figdata.GA_Contra4 = GA_Contra4;
erp_figdata.GA_Ipsi4 = GA_Ipsi4;
erp_figdata.GA_ConIpsi4 = GA_ConIpsi4;

erp_figdata.GA_PO3_LVF = GA_PO3_LVF;
erp_figdata.GA_PO3_RVF = GA_PO3_RVF;
erp_figdata.GA_PO4_LVF = GA_PO4_LVF;
erp_figdata.GA_PO4_RVF = GA_PO4_RVF;
erp_figdata.GA_Contra5 = GA_Contra5;
erp_figdata.GA_Ipsi5 = GA_Ipsi5;
erp_figdata.GA_ConIpsi5 = GA_ConIpsi5;

erp_figdata.GA_PO7_LVF = GA_PO7_LVF;
erp_figdata.GA_PO7_RVF = GA_PO7_RVF;
erp_figdata.GA_PO8_LVF = GA_PO8_LVF;
erp_figdata.GA_PO8_RVF = GA_PO8_RVF;
erp_figdata.GA_Contra6 = GA_Contra6;
erp_figdata.GA_Ipsi6 = GA_Ipsi6;
erp_figdata.GA_ConIpsi6 = GA_ConIpsi6;



erp_figdata.LH_comb_LVF = LH_comb_LVF;
erp_figdata.LH_comb_RVF = LH_comb_RVF;
erp_figdata.RH_comb_LVF = RH_comb_LVF;
erp_figdata.RH_comb_RVF = RH_comb_RVF;
erp_figdata.comb_contra = comb_contra;
erp_figdata.comb_ipsi = comb_ipsi;
erp_figdata.comb_conipsi = comb_conipsi;




end