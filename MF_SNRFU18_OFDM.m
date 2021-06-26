%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2021.03.02
% GLRT - OCDM VS OFDM
%Pd vs Pfa for OFDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
Nchirp = 512; %information bits of each subcarriers
Npacket = 1;
mod_order = 2; % 调制阶数 QPSK =2, 16QAM =4;%调制方式选择

%Noise parameter
SNR = -12;

%mente carlo 
times =1e4;
P = 128;  %多径数量

%threshold
load thre-ofdm-norm;
men = threshold;
load HforDetec;
%% 发送端
for time = 1:times

qamMapper = comm.RectangularQAMModulator('ModulationOrder', 2^mod_order,'BitInput',...
                                         true, 'NormalizationMethod', 'Average power'); %modulation method
% 产生随机数列
totalbit = mod_order * Nchirp * Npacket; %the total bits of a packet
bit_gen = randi([0 1], totalbit,1); 

%信道编码
% Encoded_bit=Func_convoCode( input_pilot ,'encode');%1/2码率
%调制
modulation = step(qamMapper,bit_gen); 
%IDFnT
input = modulation;

%with pilot
% ocdm_out_pilot = Phi_Npilot*input;
ofdm_out_pilot= fft(input);
u = ofdm_out_pilot;
uconv = toeplitz([u; zeros(P,1)],[u(1);zeros(P-1,1)]);

% %Multipath
% h = h_even(time,:);
% tmp = uconv * h; % tmp = conv(h,u);
% 
% noise = (randn(Nchirp+P,1)+1j*randn(Nchirp+P,1))./2^0.5;
% Ps = sum(abs(tmp).^2)/length(tmp);
% sigma = sqrt(Ps/(10^(SNR/10)));
% sig_r1 = tmp./sigma + noise;


%doppler+mp %sig_r6
Doppler = -0.003 + (0.003 + 0.003).*rand([1 1]);
s_Dop = resample(ofdm_out_pilot,round((1-Doppler)*6000),6000);
P = 128;  %多径数量
h = h_even(time,:);
N1 = Nchirp + P;
tmp1 = conv(s_Dop,h);
if (N1 >= length(tmp1))
   s_doppler = [tmp1;zeros(N1-length(tmp1),1)]; 
else
   s_doppler = tmp1(1:N1);
end
tmp = s_doppler;
noise = (randn(Nchirp+P,1)+1j*randn(Nchirp+P,1))./2^0.5;
Ps = sum(abs(tmp).^2)/length(tmp);
sigma = sqrt(Ps/(10^(SNR/10)));
sig_r1 = tmp./sigma + noise;


% D1 = abs(conj(s.')*sig_r6);
D1 = (norm(conj(uconv.')*sig_r1))^2;
TT(time) = D1;
% DD = sort(TT);
end

pp = zeros(1,21);
for i = 1:times
    if TT(i)>men(1)
        pp(1) = pp(1)+1;
    end
    if TT(i)>men(2)
        pp(2) = pp(2)+1;
    end
    if TT(i)>men(3)
        pp(3) = pp(3)+1;
    end
    if TT(i)>men(4)
        pp(4) = pp(4)+1;
    end
    if TT(i)>men(5)
        pp(5) = pp(5)+1;
    end
    if TT(i)>men(6)
        pp(6) = pp(6)+1;
    end
    if TT(i)>men(7)
        pp(7) = pp(7)+1;
    end
    if TT(i)>men(8)
        pp(8) = pp(8)+1;
    end
    if TT(i)>men(9)
        pp(9) = pp(9)+1;
    end
    if TT(i)>men(10)
        pp(10) = pp(10)+1;
    end
    if TT(i)>men(11)
        pp(11) = pp(11)+1;
    end
    if TT(i)>men(12)
        pp(12) = pp(12)+1;
    end
    if TT(i)>men(13)
        pp(13) = pp(13)+1;
    end
    if TT(i)>men(14)
        pp(14) = pp(14)+1;
    end
    if TT(i)>men(15)
        pp(15) = pp(15)+1;
    end
    if TT(i)>men(16)
        pp(16) = pp(16)+1;
    end
    if TT(i)>men(17)
        pp(17) = pp(17)+1;
    end
    if TT(i)>men(18)
        pp(18) = pp(18)+1;
    end
    if TT(i)>men(19)
        pp(19) = pp(19)+1;
    end
    if TT(i)>men(20)
        pp(20) = pp(20)+1;
    end
    if TT(i)>0
        pp(21) = pp(21)+1;
    end
end
Pfa = 10.^(-4:0.2:0);
Pd = pp/times;

figure();
semilogx(Pfa,Pd,'-ob','LineWidth',2);
title('ROC曲线','FontSize',20);
xlabel('Pfa','FontSize',20);
ylabel('Pd','FontSize',20);ylim([0,1]);
toc
% filename = 'MF_SNRfu18_ocdm_v1';
% MF_SNRfu18_ocdm = pp/times;
% save(filename,'MF_SNRfu18_ocdm');