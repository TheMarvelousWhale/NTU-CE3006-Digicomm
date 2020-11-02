%--Admin stuff--%
clear all; close all; clc;


%{
        This code aims to encode the original data in two different coding
        methods and compare their coding efficiency with the no encode
        version. 

        Coding method chosen: 
            1> Cyclic
            2> Hamming


%}

%define carrier frequency
fc = 10000; %10kHz
%16 times oversampled -> sample freq = 16 fc
fs = 16 * fc;

%define data rate of 1kbps
dataRate = 1000;
%define number of data bits
nBits = 1024
Enc_nBits = nBits/4*7;   %we doing (7,4) code
%define sampling rate
samplingPeriod = fs / dataRate;

%define Amplitude
Amp = 5;
%define time steps
t = 0: 1/fs : Enc_nBits/dataRate;
t_pure = 0:1/fs : nBits/dataRate;    %pure -- no encode version

%define 6th order LP butterworth filter with 0.2 normalized cutoff frequency
[b_low,a_low] = butter(6, 0.2);
%define 6th order HP butterworth filter with 0.2 normalized cutoff frequency
[b_high,a_high] = butter(6, 0.2, 'high');


%generate carrier frequency
Carrier = Amp .* cos(2*pi*fc*t);
Carrier_pure = Amp.* cos(2*pi*fc*t_pure);

%calculate signal length
SignalLength = fs*Enc_nBits/dataRate + 1;
SignalLength_pure = fs*nBits/dataRate +1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = -30:1:10;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%MODIFY THE VARIABLE BELOW TO CHOOSE AT WHICH SNR VALUE 
%TO PLOT SIGNAL,NOISE and RECEIVE
plot_SNR_dB = 15;



%set run times
Total_Run = 10;

%define placeholder for error calculation
Error_Rate_Hamming = zeros(length(SNR));
Error_Rate_Cyclic = zeros(length(SNR));
Error_Rate_NoEncode = zeros(length(SNR));

%for each SNR value
for i = 1 : length(SNR)
	Avg_Error_Hamming = 0;
    Avg_Error_Cyclic = 0;
    Avg_Error_NoEncode = 0;
    
    %for each SNR value, average the error over %Total_Run times
	for j = 1 : Total_Run
        
        %-----Data generation-----%
        Data = round(rand(1,nBits));
        EncodeHamming = encode(Data, 7, 4, 'hamming/fmt'); 
        EncodeCyclic = encode(Data,7,4,'cyclic');

        %fill the data stream
        DataStream_Hamming = zeros(1, SignalLength);
        DataStream_Cyclic = zeros(1, SignalLength);

        for k = 1: SignalLength - 1
            DataStream_Hamming(k) = EncodeHamming(ceil(k*dataRate/fs));
            DataStream_Cyclic(k) = EncodeCyclic(ceil(k*dataRate/fs));

        end
        DataStream_Hamming(SignalLength) = DataStream_Hamming(SignalLength - 1);
        DataStream_Cyclic(SignalLength) = DataStream_Cyclic(SignalLength-1);

        DataStream_NoEncode = zeros(1, SignalLength_pure);
        for k = 1:SignalLength_pure -1
            DataStream_NoEncode(k)= Data(ceil(k*dataRate/fs));
        end
        DataStream_NoEncode(SignalLength_pure) = DataStream_NoEncode(SignalLength_pure-1);
        
        
        %----- OOK -----%
        resultOOK_Hamming = OOK_transmission(DataStream_Hamming,SNR(i),Carrier,SignalLength,samplingPeriod,Enc_nBits,Amp);
        resultOOK_Cyclic = OOK_transmission(DataStream_Cyclic,SNR(i),Carrier,SignalLength,samplingPeriod,Enc_nBits,Amp);
        resultOOK_NoEncode = OOK_transmission(DataStream_NoEncode,SNR(i),Carrier_pure,SignalLength_pure,samplingPeriod,nBits,Amp);
        
        decodedHamming = decode(resultOOK_Hamming,7,4,'hamming/fmt');
        decodedCyclic = decode(resultOOK_Cyclic,7,4,'cyclic');
        
        %--Calculate Error--%
        ErrorHamming = 0;
        ErrorCyclic = 0;
        ErrorNoEncode = 0;
        for k = 1: nBits
            if(decodedHamming(k) ~= Data(k))
                ErrorHamming = ErrorHamming + 1;
            end
            if(decodedCyclic(k) ~= Data(k))
                ErrorCyclic = ErrorCyclic + 1;
            end
            if (resultOOK_NoEncode(k) ~= Data(k));
                ErrorNoEncode = ErrorNoEncode + 1;
            end
        end
        Avg_Error_Hamming = ErrorHamming + Avg_Error_Hamming;
        Avg_Error_Cyclic = ErrorCyclic + Avg_Error_Cyclic;
        Avg_Error_NoEncode = ErrorNoEncode + Avg_Error_NoEncode;
    end
    Error_Rate_Hamming(i) = Avg_Error_Hamming/Total_Run/nBits;
    Error_Rate_Cyclic(i) = Avg_Error_Cyclic/Total_Run/nBits;
    Error_Rate_NoEncode(i) = Avg_Error_NoEncode/Total_Run/nBits;
end


%Error plot
figure(1);
semilogy (SNR_dB, Error_Rate_Hamming,'r-*');
hold on
semilogy(SNR_dB, Error_Rate_Cyclic, 'b-*');
hold off
hold on
semilogy(SNR_dB, Error_Rate_NoEncode, 'k-*');
hold off
title('Error rate of cyclic and hamming and No Encoding for different SNR');
legend('hamming','cyclic','None');
ylabel('Pe');
xlabel('Eb/No')




%%--HELPER FUNCTION--%%
function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end


%This function simulates the decision device
function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end

%This function is a wrapper for Phase 2 OOK 
function result_OOK = OOK_transmission(DataStream,SNR_val,Carrier,SignalLength,samplingPeriod,Enc_nBits,Amp)
        Signal_OOK = Carrier .* DataStream;
        [b_low,a_low] = butter(6, 0.2);
        %generate noise 
        Signal_Power_OOK = (norm(Signal_OOK)^2)/SignalLength;  %Sum of squared signal amp over signal length
		Noise_Power_OOK = Signal_Power_OOK ./SNR_val;
		NoiseOOK = sqrt(Noise_Power_OOK/2) .*randn(1,SignalLength);
		
        %transmission
		ReceiveOOK = Signal_OOK+NoiseOOK;
        %detection -- square law device
        SquaredOOK = ReceiveOOK .* ReceiveOOK;
        %low pass filter
        FilteredOOK = filtfilt(b_low, a_low, SquaredOOK);
         
        %sample and decision device
        sampledOOK = sample(FilteredOOK, samplingPeriod, Enc_nBits);
        result_OOK = decision_device(sampledOOK,Enc_nBits, Amp/2);  %--OOK threshold is 0.5*(A+0)

end