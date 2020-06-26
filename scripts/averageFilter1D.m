function signal = averageFilter1D(signal)

% Usage : signal = averageFilter1D(signal) 
% signal must be column
% simple 5 window averaging filter

newSignal=[];

for i = 1:size(signal,1)
    
    if(i==1 | i==size(signal,1))
        
        newSignal = [newSignal;signal(i,1)]; 
        
    elseif(i==2 | i==size(signal,1)-1)
        
        newSignal = [newSignal;mean([signal(i-1,1); signal(i,1); signal(i+1,1)])];
        
    else
    
        newSignal = [newSignal;mean([signal(i-2,1); signal(i-1,1); signal(i,1); signal(i+1,1); signal(i+2,1)])];
        
    end
    
end

signal = newSignal;
