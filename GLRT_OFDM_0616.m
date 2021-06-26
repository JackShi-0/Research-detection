%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20210616
%GLRT - OCDM VS OFDM
%Threshold for OCDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
% parameter_setting1;%参数设置
Nchirp = 512; %information bits of each subcarriers
Npacket = 1;
% Phi_Nchirp = Phi(Nchirp);
load('H_for_Detec.mat')
%Noise parameter
SNR = -12;
%调制方式选择
mod_order = 2; % 调制阶数 QPSK =2, 16QAM =4;
qamMapper = comm.RectangularQAMModulator('ModulationOrder', 2^mod_order,'BitInput',...
                                         true, 'NormalizationMethod', 'Average power'); %modulation method

totalbit = mod_order * Nchirp * Npacket; %the total bits of a packet
for num = 1 : 10
% 产生随机数列
seed = num;
rng(seed);
bit_gen = randi([0 1], totalbit,1); 
%调制
modulation = step(qamMapper,bit_gen); 
%IDFnT
input = modulation;
%input
% ocdm_out_pilot = Phi_Nchirp * input;
ofdm_out_pilot= fft(input);
u = ofdm_out_pilot;


% channel 
P = 128;  %多径数量
uconv = toeplitz([u; zeros(P,1)],[u(1);zeros(P-1,1)]);

[U,S,V] = svd(uconv);
Up = U(:,1:P);
Up2 = Up*Up';

% tmp = uconv * h; % tmp = conv(h,u);
%% 发送端
%mente carlo 
% tic;
times = 1e5;
for time = 1:times
noise = (randn(Nchirp+P,1)+1j*randn(Nchirp+P,1))./2^0.5; 
%D0 = (norm(conj(uconv.') * noise))^2;
% D0 = norm(uconv'*sig_r); %和上一条一个意思
D0 = noise'*Up2*noise;
% noise'*uconv*inv(uconv'*uconv)*uconv'*noise; %the same with the upper
TT(time) = D0;
% DD = sort(TT);
end
% toc;
Pfa =10.^(-4:0.2:0);
threshold = find_threshold_NMF(TT,Pfa,1e3,1e5);

mc=1000;
for time = 1:mc
h = h_even(time,:)';

tmp = uconv * h; % tmp = conv(h,u);
noise = (randn(Nchirp+P,1)+1j*randn(Nchirp+P,1))./2^0.5;
Ps = sum(abs(tmp).^2)/length(tmp);
sigma = sqrt(Ps/(10^(SNR/10)));
sig_r1 = tmp./sigma + noise;
%D1 = (norm(conj(uconv.')*sig_r1))^2;
D1 = sig_r1'*Up2*sig_r1;
D_multipath(time) = D1;

% %doppler+mp %sig_r6
Doppler = -0.003 + (0.003 + 0.003).*rand([1 1]);
s_Dop = resample(u,round((1-Doppler)*6000),6000);
N1 = Nchirp + P;
tmp1 = conv(s_Dop,h);
if (N1 >= length(tmp1))
   s_doppler = [tmp1;zeros(N1-length(tmp1),1)]; 
else
   s_doppler = tmp1(1:N1);
end
sig_r2 = s_doppler./sigma + noise;
% D2 = (norm(conj(uconv.')*sig_r2))^2;
D2 = sig_r2'*Up2*sig_r2;
D_doppler(time) = D2;
end

Pd_mp(num,:) = analyse(D_multipath,threshold);
Pd_doppler(num,:) = analyse(D_doppler,threshold);
end

pd_mp = mean(Pd_mp);
pd_doppler = mean(Pd_doppler);

Pfa = 10.^(-4:0.2:0);
figure();
semilogx(Pfa,pd_mp,'r-o','LineWidth',2);
hold on;
semilogx(Pfa,pd_doppler,'b-o','LineWidth',2);
hold on;
title('ROC曲线','FontSize',20);
xlabel('Pfa','FontSize',20);
ylabel('Pd','FontSize',20);ylim([0,1]);
toc
% save 0616-GLRT-mp-ofdm-fixsymbol pd_mp
% save 0616-GLRT-doppler-ofdm-fixsymbol pd_doppler
save pd_mp_ofdm pd_mp
save pd_doppler_ofdm pd_doppler
save Pd_mp_100_ofdm Pd_mp
save Pd_doppler_100_ofdm Pd_doppler