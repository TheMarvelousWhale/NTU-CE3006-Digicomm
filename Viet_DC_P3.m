%--Admin Stuff--%

clear all; close all; clc;

pskdemod
%define carrier frequency of 10kHz
fc = 10000; 
%define sampling frequency of 16*fc
Fs = 16 * fc;
%define data rate of 1kbps
dataRate = 1000; 


%define number of data bit
nBits = 1024;
%define num bit after encoding
Enc_nBits = 1792;
%Define Amplitude
Amplitude = 5;


%define time steps
t = 0: 1/Fs : Enc_nBits/dataRate;
%define sampling 
samplingPeriod = Fs / dataRate;


%define low pass butterworth filter 6th order with normalized 0.2 cutoff frequency
[b_low, a_low] = butter(6, 0.2, 'low');
%define high pass butterworth filter  6th order with normalized 0.2 cutoff frequency
[b_high, a_high] = butter(6, 0.2, 'high');


%generate carrier frequency
Carrier = Amplitude .* cos(2*pi*fc*t);

%define signal length
SignalLength = Fs*Enc_nBits/dataRate + 1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = -20:1:20;
plot_SNR_dB = -10;
%%SNR = S/N = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

Total_Run = 20;

Error_RateOOK = zeros(length(SNR));
Error_RateBPSK = zeros(length(SNR));

%Different SNR value
for i = 1 : length(SNR)
	Avg_ErrorOOK = 0;
    Avg_ErrorBPSK = 0;
    
	for j = 1 : Total_Run
        
        %------ generate data -----%
        Data = round(rand(1,nBits));
        
        %encode
        EncodeHamming = encode(Data, 7, 4, 'hamming/fmt'); 
        
        DataStream = zeros(1, SignalLength);
        for k = 1: SignalLength - 1
            DataStream(k) = EncodeHamming(ceil(k*dataRate/Fs));
        end
        DataStream(SignalLength) = DataStream(SignalLength - 1);

        %----- OOK -----%
        SignalOOK = Carrier .* DataStream;
        Signal_Power_OOK = (norm(SignalOOK)^2)/SignalLength;
        
        %Generate Noise OOK
		Noise_Power_OOK = Signal_Power_OOK ./SNR(i);
		NoiseOOK = sqrt(Noise_Power_OOK/2) .*randn(1,SignalLength);
		%Received Signal OOK
		ReceiveOOK = SignalOOK+NoiseOOK;
        
        %OOK detection
        SquaredOOK = ReceiveOOK .* ReceiveOOK;
        %low pass filter
        FilteredOOK = filtfilt(b_low, a_low, SquaredOOK);
        
        %sample and decide  
        sampledOOK = sample(FilteredOOK, samplingPeriod, Enc_nBits);
        resultOOK = decision_device(sampledOOK,Enc_nBits, Amplitude/2);  %--OOK threshold is 0.5*(A+0)
        
        %------ BPSK -----%
        DataStreamBPSK = DataStream .* 2 - 1;       %change to -1 and 1 
        SignalBPSK = Carrier .* DataStreamBPSK;
        Signal_Power_BPSK = (norm(SignalBPSK)^2)/SignalLength;
        %Generate Noise BPSK
		Noise_Power_BPSK = Signal_Power_BPSK ./SNR(i);
		NoiseBPSK = sqrt(Noise_Power_BPSK/2) .*randn(1,SignalLength);
		%Received Signal BPSK
		ReceiveBPSK = SignalBPSK+NoiseBPSK;
        
        %non-coherent detection
        SquaredBPSK = ReceiveBPSK .* ReceiveBPSK;
        %high pass filter (supposingly band pass filter)
        FilteredBPSK = filtfilt(b_high, a_high, SquaredBPSK);
        
        %frequency divider
        DividedBPSK = interp(FilteredBPSK, 2);
        DividedBPSK = DividedBPSK(1:length(FilteredBPSK));
        
        %Multiple and Low Pass Filter
        MultipliedBPSK = DividedBPSK .* ReceiveBPSK;
        OutputBPSK = filtfilt(b_low, a_low, MultipliedBPSK);
        
        %demodulate
        %sample and decide
        sampledBPSK = sample(OutputBPSK, samplingPeriod, Enc_nBits);
        resultBPSK = decision_device(sampledBPSK,Enc_nBits,0);           %-- bipolar -- threshold 0        


        %decode
        decodedOOK = decode(resultOOK,7,4,'hamming/fmt');
        decodedBPSK = decode(resultBPSK,7,4,'hamming/fmt');
        
        ErrorOOK = 0;
        ErrorBPSK = 0;
        for k = 1: nBits - 1
            if(decodedOOK(k) ~= Data(k))
                ErrorOOK = ErrorOOK + 1;
            end
            if(decodedBPSK(k) ~= Data(k))
                ErrorBPSK = ErrorBPSK + 1;
            end
        end
        Avg_ErrorOOK = ErrorOOK + Avg_ErrorOOK;
        Avg_ErrorBPSK = ErrorBPSK + Avg_ErrorBPSK;

    end
    
    if (plot_SNR_dB == SNR_dB(i))
            plot_signal = Data;
            plot_mod_OOK = SignalOOK;
            plot_receive_OOK = ReceiveOOK;
            plot_demod_OOK = FilteredOOK;
            plot_decoded_OOK = decodedOOK;
            plot_mod_BPSK = SignalBPSK;
            plot_receive_BPSK = ReceiveBPSK;
            plot_demod_BPSK = FilteredBPSK;
            plot_decoded_BPSK = decodedBPSK;
    end
    
	Error_RateOOK(i) = Avg_ErrorOOK / Total_Run / nBits;
    Error_RateBPSK(i) = Avg_ErrorBPSK / Total_Run/nBits;
end

%Error plot
figure(1);
semilogy (SNR_dB,Error_RateOOK,'k-*');
hold on
semilogy(SNR_dB, Error_RateBPSK, 'c-*');
hold off
title('Error rate of OOK and BPSK for different SNR');
legend('OOK', 'BPSK');
ylabel('Pe');
xlabel('Eb/No')


%OOK plot
figure(2);
subplot(511);title('Generated Data');plot(plot_signal);
subplot(512);title('Modulated OOK');plot(plot_mod_OOK,'k');
subplot(513);title('Received Signal OOK');plot(plot_receive_OOK, 'k')
subplot(514);title('Demodulated OOK');plot(plot_demod_OOK, 'k');
subplot(515);title('Decoded Data');plot(plot_decoded_OOK);

%BPSK plot
figure(3)
subplot(511);title('Generated Data');plot(plot_signal);
subplot(512);title('Modulated BPSK');plot( plot_mod_BPSK,'k');
subplot(513);title('Received Signal BPSK');plot(plot_receive_BPSK, 'k')
subplot(514);title('Demodulated BPSK');plot(plot_demod_BPSK, 'k');
subplot(515);title('Decoded Data');plot(plot_decoded_BPSK);

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