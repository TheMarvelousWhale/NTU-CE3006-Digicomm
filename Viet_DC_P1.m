clear all; 
close all; 
clc;

%---------------Define Declaration---------------%
%Num of bits = 1024
nBits = 1024;

%Signal power = 1                   
Signal_Power = 1; 

%SNR_dB = 10log(Signal_Power/Noise_Power)
%Create matrix to store SNR
SNR_dB = 0:0.5:20;

%SNR = S/N = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%Holder value for plotting later
plot_signal = rand(1,nBits);
plot_noise = rand(1,nBits);
plot_receive = rand(1,nBits);

%MODIFY THE VARIABLE BELOW TO CHOOSE AT WHICH SNR VALUE 
%TO PLOT SIGNAL,NOISE and RECEIVE
plot_SNR_dB = 15;

%Counter for number of run to calculate BER for each SNR value
Total_Run = 20;
%---------------END OF DEFINE--------------------%


%----------------MAIN ROUTINE--------------------%

%Iterate through different SNR value
for i = 1 : length(SNR)
	Sum_Error = 0;
	for j = 1 : Total_Run
        
        
        %-----------TRANSMITTER------------%
        
		%Generate random binary digits(0 or 1), INPUT SIGNAL
		Data = round(rand(1,nBits));
        
		%Convert binary digit to -1 or +1: 
        %if Data is 0 -> signal = 2*-0.5 = -1
		%else if Data is 1 -> signal = 2*0.5 = +1
		Signal = 2 .* (Data - 0.5);
        
		%Get Noise Power from SNR
		Noise_Power = Signal_Power ./SNR(i);
        
        %Get Noise Variance from Noise Power
        %Change the noise variance with respect to SNR (signal to noise ratio) value
        Avg_Noise_Power = Noise_Power/2;
        
        %The randn() function is for normal distribution, generated with equal number of noise samples 
        %(In Matlab, a.*randn(1000,1) + b is a Gaussian Dist, mean b sd a.) 
        Noise = sqrt(Avg_Noise_Power) .*randn(1,nBits);   
        
		%Receiver side Signal
		Receive = Signal+Noise;                      

        
        %-----RECEIVER------%
        
        %Initialize threshold
		Threshold = 0;
        
        %Initialize Error for this run
		Error = 0;
        
        %Fix the threshold value as 0
        %If received signal >= threshold value, threshold = 1
        %If received signal < threshold value, threshold = 0
		for k= 1 : nBits
			if (Receive(k)>= Threshold) && Data(k)==0||(Receive(k)<Threshold && Data(k)==1)
				Error = Error+1;
			end
        end
        
        
        %---------------ERROR STATISTICS CALCULATION--------------------% 
		%BER = TotalError divide by number of bits
		Error = Error ./nBits;  
        
		%Accumulate the error of each run	
		Sum_Error = Error + Sum_Error;
        
        %---------------Choose the value for plotting------------------%
        if (plot_SNR_dB == SNR_dB(i))
            plot_signal = Signal;
            plot_noise = Noise;
            plot_receive = Receive;
        end
        
    end
    
    %Average Error
	Error_Rate(i) = Sum_Error / Total_Run;
end

%Predict BER using 
Pred_BER=(1/2)*erfc(sqrt(SNR)); 


%--------------------------------PLOTTING--------------------------------%
figure("position", [10,100,1400,800]) 

%BER against SNR -- main plot
subplot(221)
semilogy (SNR_dB,Error_Rate,'r*');
xlabel('Normalized SNR')
ylabel('Probability Error');
title('BER against SNR');
hold on
semilogy (SNR_dB,Pred_BER,'m');
legend('Simulation','Prediction');
axis([0 20 10^(-5) 1]);
hold off

%data generation
subplot(222)
plot(plot_signal);
title('Data Generated')

%noise generation
subplot(223)
plot(plot_noise);
title('Noise Generated')

%received data generation
subplot(224)
plot(plot_receive);
title('Received Data')



