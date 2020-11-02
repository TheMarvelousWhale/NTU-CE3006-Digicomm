

fc = 10000; 
fs = 16 * fc;
dataRate = 1000;
nBits = 16;
SymbolRate = 2;
nSymbols = nBits / SymbolRate;
SignalLength = fs*nSymbols/dataRate+1; 
SymbolStream = zeros(1, SignalLength);
Symbol = zeros(1,nSymbols);
Amp = 10;
Data = round(rand(1,nBits));
t = 0: 1/fs : nSymbols/dataRate;
Carrier = Amp .* cos(2*pi*fc*t);

for k = 1: nSymbols
    Symbol(k) = Data(2*k-1)*2+Data(2*k);
end

for k = 1: SignalLength - 1
    
    SymbolStream(k) = Symbol(ceil(k*dataRate/fs));
end
SymbolStream(SignalLength) = SymbolStream(SignalLength - 1);


wave1 = Carrier.* SymbolStream;


%define 6th order LP butterworth filter with 0.2 normalized cutoff frequency
[b,a] = butter(6, 0.2);
[c,d] = butter(6, 0.2, 'high');


%figure(1);
subplot(111);title('wave 1');plot(wave1);

% Display and compare results
%hfvt = fvtool(b,a,'FrequencyScale','log');
