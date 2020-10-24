clear all; close all; clc;


%define carrier frequency
fc = 10000; %10kHz
%Assume carrier signal is 16 times oversampled
%define sampling frequency
fs = 16 * fc;
%define data rate
dataRate = 1000; %10kbps
%Define Signal length
nBits = 1024;
%define num bit after encoding

%Define Amp
Amp = 5;

%low pass butterworth filter
%6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);

%high pass butterworth filter
[d, c] = butter(6, 0.2, 'high');

%time
t = 0: 1/fs : nBits/dataRate;

%Carrier
Carrier = Amp .* cos(2*pi*fc*t);

%signal length
SignalLength = fs*nBits/dataRate + 1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = 0:1:20;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

% Set run times
Total_Run = 20;

Error_RateOOK = zeros(length(SNR));
Error_RateBPSK = zeros(length(SNR));

%Different SNR value
for i = 1 : length(SNR)
	Avg_ErrorOOK = 0;
    Avg_ErrorBPSK = 0;
    
	for j = 1 : Total_Run
        %Generate Data
        Data = round(rand(1,nBits));
        %encode
        
        ContinuousData = zeros(1, SignalLength);
        for k = 1: SignalLength - 1
            ContinuousData(k) = Data(ceil(k*dataRate/fs));
        end
        ContinuousData(SignalLength) = ContinuousData(SignalLength - 1);

        %on-off keying
        SignalOOK = Carrier .* ContinuousData;

        %binary phase shift keying
        ContinuousDataBPSK = ContinuousData .* 2 - 1;
        SignalBPSK = Carrier .* ContinuousDataBPSK;

        Signal_Power_OOK = (norm(SignalOOK)^2)/SignalLength;
        Signal_Power_BPSK = (norm(SignalBPSK)^2)/SignalLength;
        
        %Generate Noise OOK
		Noise_Power_OOK = Signal_Power_OOK ./SNR(i);
		NoiseOOK = sqrt(Noise_Power_OOK/2) .*randn(1,SignalLength);
		%Received Signal OOK
		ReceiveOOK = SignalOOK+NoiseOOK;
        
        %OOK detection
        SquaredOOK = ReceiveOOK .* ReceiveOOK;
        %low pass filter
        FilteredOOK = filtfilt(b, a, SquaredOOK);
        
        %Generate Noise BPSK
		Noise_Power_BPSK = Signal_Power_BPSK ./SNR(i);
		NoiseBPSK = sqrt(Noise_Power_BPSK/2) .*randn(1,SignalLength);
		%Received Signal BPSK
		ReceiveBPSK = SignalBPSK+NoiseBPSK;
        
        %non-coherent detection
        SquaredBPSK = ReceiveBPSK .* ReceiveBPSK;
        %high pass filter (supposingly band pass filter)
        FilteredBPSK = filtfilt(d, c, SquaredBPSK);
        
        %frequency divider
        DividedBPSK = interp(FilteredBPSK, 2);
        DividedBPSK = DividedBPSK(1:length(FilteredBPSK));
        
        %Multiple and Low Pass Filter
        MultipliedBPSK = DividedBPSK .* ReceiveBPSK;
        OutputBPSK = filtfilt(b, a, MultipliedBPSK);
        
        %demodulate
        %sampling AND threshold
        
        
        samplingPeriod = fs / dataRate;
        sampledOOK = sample(FilteredOOK, samplingPeriod, nBits);
        sampledBPSK = sample(OutputBPSK, samplingPeriod, nBits);
        resultOOK = decision_device(sampledOOK,nBits, Amp/2);  %--OOK threshold is 0.5*(A+0)
        resultBPSK = decision_device(sampledBPSK,nBits,0);           %-- bipolar -- threshold 0        
        
		%Calculate the average error for every runtime
		%Avg_ErrorOOK = num_error(resultOOK, EncodeHamming, nBits) + Avg_ErrorOOK;                   
        %Avg_ErrorBPSK = num_error(resultBPSK, EncodeHamming, nBits) + Avg_ErrorBPSK;
        
       
        ErrorOOK = 0;
        ErrorBPSK = 0;
        for k = 1: nBits - 1
            if(resultOOK(k) ~= Data(k))
                ErrorOOK = ErrorOOK + 1;
            end
            if(resultBPSK(k) ~= Data(k))
                ErrorBPSK = ErrorBPSK + 1;
            end
        end
        Avg_ErrorOOK = ErrorOOK + Avg_ErrorOOK;
        Avg_ErrorBPSK = ErrorBPSK + Avg_ErrorBPSK;

	end
	Error_RateOOK(i) = Avg_ErrorOOK / Total_Run;
    Error_RateBPSK(i) = Avg_ErrorBPSK / Total_Run;
end

figure(1);
subplot(221);
semilogy (SNR_dB,Error_RateOOK,'k-*');
hold on
semilogy(SNR_dB, Error_RateBPSK, 'c-*');
%axis([0 20 10^(-5) 1]);
hold off
legend('OOK', 'BPSK');
ylabel('Pe');
xlabel('Eb/No')


 
subplot(222)
title('Received Signal OOK');
plot(ReceiveOOK, 'k')

subplot(223)
title('Squared OOK');
plot(SquaredOOK, 'k');

subplot(224)
title('Filtered Signal OOK');
plot(FilteredOOK, 'k');

figure(2)
title('Captured Data');
plot(sampledOOK);


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