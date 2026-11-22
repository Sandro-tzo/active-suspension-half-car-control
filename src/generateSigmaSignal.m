function sigma = generateSigmaSignal()
   
    time = 0:0.001:20; 
    
    values = zeros(size(time));

    values(time >=0  & time <= 11) = 2;
    values(time >11  & time <= 13) = 0;
    values(time > 13 & time <= 20) = 1;
    
   
    sigma = timeseries(values, time);
end