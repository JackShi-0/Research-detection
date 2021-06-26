%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2021.03.02
%GLRT - OCDM VS OFDM
%Threshold for OFDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
% parameter_setting1;%参数设置
Nchirp = 512; %information bits of each subcarriers
Npacket = 1;
%调制方式选择
mod_order = 2; % 调制阶数 QPSK =2, 16QAM =4;

%Noise parameter
% SNR = -18;

%mente carlo 
times =1e6;
%% 发送端
for time = 1:times
qamMapper = comm.RectangularQAMModulator('ModulationOrder', 2^mod_order,'BitInput',...
                                         true, 'NormalizationMethod', 'Average power'); %modulation method
% 产生随机数列
totalbit = mod_order * Nchirp * Npacket; %the total bits of a packet
bit_gen = randi([0 1], totalbit,1); 

% %信道编码
% Encoded_bit=Func_convoCode( input_pilot ,'encode');%1/2码率
%调制
modulation = step(qamMapper,bit_gen); 
%IDFnT
input = modulation;

%with pilot
% ocdm_out_pilot = Phi_Npilot*input;
ofdm_out_pilot= fft(input);
u = ofdm_out_pilot;

% channel 
P = 128;  %多径数量
% t = 0:1e-3:1e-3*(P-1);
% dB = 40;
% beta = (dB/10)*log(10)/(t(P)-t(1));
% B = exp(-beta*t);
% A_real = raylrnd(B);
% A_img = raylrnd(B);
% A = A_real+1i*A_img;
% A = A/sqrt(sum(A.^2));% 幅值
% h = A.';

uconv = toeplitz([u; zeros(P,1)],[u(1);zeros(P-1,1)]);
% tmp = uconv * h; % tmp = conv(h,u);

%noise
% s = tmp;

noise = (randn(Nchirp+P,1)+1j*randn(Nchirp+P,1))./2^0.5;
    
D0 = (norm(conj(uconv.')*noise))^2;

%NMF
%D = abs(conj(temp1_ocdm.')*sig_r)/(sqrt(sum(abs(temp1_ocdm.').^2)) * sqrt(sum(abs(sig_r).^2)));

TT(time) = D0;

end

%Method2 ROC
% figure();
% ngam =50;
% [Pfa_mf,Pd_mf,gam]=roccurve(TT0,TT1,ngam);
% plot(Pfa_mf,Pd_mf);
% xlabel('Probability of false alarm (Pfa)');
% ylabel('Probability of detection (Pd)');

Pfa =10.^(-4:0.2:0);
threshold = find_threshold_NMF(TT,Pfa,1e3,1e6);

toc