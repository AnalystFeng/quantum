%% QPSK
clear all;
clc;
theory = 0;
%SymbolMapping_ = 'Binary';
SymbolMapping_ = 'Gray';                                                  %使用格雷码
qpsk_mo = comm.QPSKModulator('SymbolMapping',SymbolMapping_,...
                                     'BitInput',true ); %Binary | Gray    %比特流输入
qpsk_demo = comm.QPSKDemodulator('SymbolMapping',SymbolMapping_,...
                                         'BitOutput',true);               %比特流输出
snr = -5:1:15;                                                            %仿真信噪比
bit_number = 10^7;                                                        %仿真比特数

for i = 1:length(snr)
    
        if( i == 1)
           tic;
        end
        
        bit_in = randi([0 1],1,bit_number).'; 
        
        bit_channel = qpsk_mo(bit_in);%归一化输出
        bit_channel_awgn = awgn(bit_channel,snr(i));
       %scatterplot(bit_channel_awgn);                                    %查看调制星座图
        bit_out = qpsk_demo(bit_channel_awgn); 
       
        Totalerror_bit = length(find(bit_in ~= bit_out));
        
        if i==1
          p_bit = zeros(1,length(snr));
        end
        
        p_bit(i) = Totalerror_bit/length(bit_in);
        clear bit_in bit_out bit_channel_awgn bit_channel;
        fprintf('%d\n',i);
        toc;
end


Pb_theory = 0.5* erfc(sqrt(10.^(snr/10)/2));
semilogy(snr,p_bit,'ks-',snr,Pb_theory,'k*-');

legend('误比特率','理论误比特率');
title(['AWGN QPSK 硬判决'  ' ' SymbolMapping_  ]);
grid on;
xlabel('SNR');ylabel('BER');




%% 16QAM
clear all;
clc;

qam_16_1 = sqrt(10);                      %能量归一化因子
Constellation_ = [-3-3j -3-1j -3+1j -3+3j -1-3j -1-1j -1+1j -1+3j 1-3j 1-1j 1+1j 1+3j 3-3j 3-1j 3+1j 3+3j]./qam_16_1;
                                          %16QAM星座图
QAMMod = comm.GeneralQAMModulator('Constellation',Constellation_); 
QAMDeMod = comm.GeneralQAMDemodulator('Constellation',Constellation_ ,...
                                              'BitOutput',true,...
                                              'DecisionMethod','Hard decision');
snr = -5:1:15;
bit_number_ = 10^7;
dec_number = bit_number_/4;

for i = 1:length(snr)
    
        if( i == 1)
           tic;
        end
        
        dec_in = randi([0 15],1,dec_number).'; 
        bit_in = dec4_to_bin(dec_in);
        
        dec_channel = QAMMod(dec_in);%归一化输出
       %scatterplot(dec_channel);
        
        bit_channel_awgn = awgn(dec_channel,snr(i));
        
        bit_out = QAMDeMod(bit_channel_awgn).'; 
       
        Totalerror_bit = length(find(bit_in ~= bit_out));
        
        if i==1
          p_bit = zeros(1,length(snr));
        end
        
        p_bit(i) = Totalerror_bit/numel(bit_in);
        clear dec_in bit_in bit_out bit_channel_awgn dec_channel;
        fprintf('%d\n',i);
        toc;
end
semilogy(snr,p_bit,'ks-');
title(['AWGN 16QAM 硬判决']);
grid on;hold on;
xlabel('SNR');ylabel('BER');




%% MISO QPSK
 clc;
 clear all;
%SymbolMapping_ = 'Binary';
SymbolMapping_ = 'Gray';
qpsk_mo = comm.QPSKModulator('SymbolMapping',SymbolMapping_,...
                                 'BitInput',true );
qpsk_demo = comm.QPSKDemodulator('SymbolMapping',SymbolMapping_,...
                                     'BitOutput',true);

SNR_dB=-5:1:15;
bit_number = 10^6;
Frame=bit_number/4;
tic;
for i=1:length(SNR_dB)

ErrorNum=0;

    for j =1:1:Frame
    
    bit_in = randi([0 1],1,4).'; 
    
    bit_channel = qpsk_mo(bit_in);%归一化输出
    x1 = bit_channel(1);
    x2 = bit_channel(2);
   %---------------------------进行空时编码--------------------------------- 
    X_1 = [x1 -conj(x2)]; 
    X_2 = [x2 conj(x1)];
    X_1_channel = awgn(X_1,SNR_dB(i));
    X_2_channel = awgn(X_2,SNR_dB(i));
    R = X_1_channel + X_2_channel;   % R(1)=x1+x2 R(2)=(-conj(x2))+conj(x1)
    %---------------------------进行空时解码--------------------------------
%     A=h1*conj(h1)+h2*conj(h2);
%     X1=(conj(h1)*R(1)+h2*conj(R(2)))/A;
%     X2=(conj(h2)*R(1)-h1*conj(R(2)))/A;
    X1=(R(1)+conj(R(2)))/2;
    X2=(R(1)-conj(R(2)))/2;
    bit_channel_awgn = [X1 X2].';
    
    bit_out = qpsk_demo(bit_channel_awgn); 
    
    ErrorNum=ErrorNum+length(find(bit_in ~= bit_out));
    end
    fprintf('%d\n',i);
    toc;
    P(i)=ErrorNum/(Frame*4);
end

semilogy(SNR_dB,P,'k+-');
title(['AWGN 2*1-MISO-QPSK 硬判决' ' '  'Gray']);
grid on;
xlabel('SNR');ylabel('BER');

%% MIMO QPSK
 clc;
 clear all;
%SymbolMapping_ = 'Binary';
SymbolMapping_ = 'Gray';
qpsk_mo = comm.QPSKModulator('SymbolMapping',SymbolMapping_,...
                                 'BitInput',true );
qpsk_demo = comm.QPSKDemodulator('SymbolMapping',SymbolMapping_,...
                                     'BitOutput',true);

SNR_dB=-5:1:15;
bit_number = 10^6;
Frame=bit_number/4;
tic;
for i=1:length(SNR_dB)

ErrorNum=0;

    for j =1:1:Frame
    
    bit_in = randi([0 1],1,4).'; 
    
    bit_channel = qpsk_mo(bit_in);%归一化输出
    x1 = bit_channel(1);
    x2 = bit_channel(2);
   %---------------------------进行空时编码--------------------------------- 
    X_1 = [x1 -conj(x2)]; 
    X_2 = [x2 conj(x1)];
    X_1_channel = awgn(X_1,SNR_dB(i));
    X_2_channel = awgn(X_2,SNR_dB(i));
    X1_1_channel = awgn(X_1,SNR_dB(i));
    X1_2_channel = awgn(X_2,SNR_dB(i));
    R1 = X_1_channel + X_2_channel;   % R(1)=x1+x2 R(2)=(-conj(x2))+conj(x1)
    R2 = X1_1_channel + X1_2_channel;
    %---------------------------进行空时解码--------------------------------
%     A=h1*conj(h1)+h2*conj(h2);
%     X1=(conj(h1)*R(1)+h2*conj(R(2)))/A;
%     X2=(conj(h2)*R(1)-h1*conj(R(2)))/A;
    X11=(R1(1)+conj(R1(2)))/2;
    X21=(R1(1)-conj(R1(2)))/2;
    bit_channel_awgn1 = [X11 X21].';
    
    X12=(R2(1)+conj(R2(2)))/2;
    X22=(R2(1)-conj(R2(2)))/2;
    bit_channel_awgn2 = [X12 X22].';
    
    bit_channel_awgn = (bit_channel_awgn1 + bit_channel_awgn2)/2;
    bit_out = qpsk_demo(bit_channel_awgn); 
    
    ErrorNum=ErrorNum+length(find(bit_in ~= bit_out));
    end
    fprintf('%d\n',i);
    toc;
    P(i)=ErrorNum/(Frame*4);
end

semilogy(SNR_dB,P,'k+-');
title(['AWGN 2*1-MISO-QPSK 硬判决' ' '  'Gray' ' '  num2str(4*Frame/10^6) 'Mbit']);
grid on;
xlabel('SNR');ylabel('BER');

function bitout = dec4_to_bin(decin)

bitout = zeros(numel(decin),4);

for i = 1:length(decin)
    switch decin(i)
            case 0
            bitout(i,1) = 0;
            bitout(i,2) = 0;
            bitout(i,3) = 0;
            bitout(i,4) = 0;
            
            case 1
            bitout(i,1) = 0;
            bitout(i,2) = 0;
            bitout(i,3) = 0;
            bitout(i,4) = 1;
            
            case 2
            bitout(i,1) = 0;
            bitout(i,2) = 0;
            bitout(i,3) = 1;
            bitout(i,4) = 0;
            
            case 3
            bitout(i,1) = 0;
            bitout(i,2) = 0;
            bitout(i,3) = 1;
            bitout(i,4) = 1;
            
            case 4
            bitout(i,1) = 0;
            bitout(i,2) = 1;
            bitout(i,3) = 0;
            bitout(i,4) = 0;
            
            case 5
            bitout(i,1) = 0;
            bitout(i,2) = 1;
            bitout(i,3) = 0;
            bitout(i,4) = 1;
            
            case 6
            bitout(i,1) = 0;
            bitout(i,2) = 1;
            bitout(i,3) = 1;
            bitout(i,4) = 0;
            
            case 7
            bitout(i,1) = 0;
            bitout(i,2) = 1;
            bitout(i,3) = 1;
            bitout(i,4) = 1;
            
            case 8
            bitout(i,1) = 1;
            bitout(i,2) = 0;
            bitout(i,3) = 0;
            bitout(i,4) = 0;
            
            case 9
            bitout(i,1) = 1;
            bitout(i,2) = 0;
            bitout(i,3) = 0;
            bitout(i,4) = 1;
            
            case 10
            bitout(i,1) = 1;
            bitout(i,2) = 0;
            bitout(i,3) = 1;
            bitout(i,4) = 0;
            
            case 11
            bitout(i,1) = 1;
            bitout(i,2) = 0;
            bitout(i,3) = 1;
            bitout(i,4) = 1;
            
            case 12
            bitout(i,1) = 1;
            bitout(i,2) = 1;
            bitout(i,3) = 0;
            bitout(i,4) = 0;
            
            case 13
            bitout(i,1) = 1;
            bitout(i,2) = 1;
            bitout(i,3) = 0;
            bitout(i,4) = 1;
            
            case 14
            bitout(i,1) = 1;
            bitout(i,2) = 1;
            bitout(i,3) = 1;
            bitout(i,4) = 0;
            
            case 15
            bitout(i,1) = 1;
            bitout(i,2) = 1;
            bitout(i,3) = 1;
            bitout(i,4) = 1;
            
    end
end
    bitout = reshape(bitout.',1,numel(bitout));
end




