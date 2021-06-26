
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2020.1.18
%Plot ROC results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
Pfa = 10.^(-4:0.2:0);

%test-new
load('ofdm-fu12-mp.mat')
ofdm_mp = Pd;
load('ofdm-fu12-doppler.mat')
ofdm_doppler = Pd;
% load('ofdm_mp+doppler_v2.mat')
% ofdm_dop = Pd;
load('ocdm-fu12-mp.mat')
ocdm_mp = Pd;
load('ocdm-fu12-doppler.mat')
ocdm_doppler = Pd;
% load('ocdm_mp+doppler_v2.mat')
% ocdm_dop = Pd;

figure; 
semilogx(Pfa,ocdm_mp,'ro-','Linewidth',1.5,'Markersize',7.0);
hold on; 
semilogx(Pfa,ocdm_doppler,'r+-','Linewidth',1.5,'Markersize',7.0);
% hold on; 
% semilogx(Pfa,ocdm_dop,'rd-','Linewidth',1.5,'Markersize',7.0);
hold on;
semilogx(Pfa,ofdm_mp,'bo:','Linewidth',1.5,'Markersize',7.0);
hold on; 
semilogx(Pfa,ofdm_doppler,'b+:','Linewidth',1.5,'Markersize',7.0);
hold on; 
% semilogx(Pfa,ofdm_dop,'bd:','Linewidth',1.5,'Markersize',7.0);
% hold on;

xlabel('Probability of false alarm (P_{fa})');
ylabel('Probability of detection (P_d)');ylim([0,1]);
% set(gca, 'FontSize', 13)
legend('OCDM-AWGN+Multipath','OCDM-AWGN+Multipath+Doppler',...
'OFDM-AWGN+Multipath','OFDM-AWGN+Multipath+Doppler',...
'Location','southeast');
%,'NumColumns',2
% grid on
% set(gcf, 'OuterPosition' , [0 300 700 600])
set(gcf, 'Color', 'w'); % to set the background color white (gray by default)
% legend('OCDM_AWGN','OFDM_AWGN');