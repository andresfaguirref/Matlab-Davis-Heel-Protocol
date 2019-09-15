function data_filtered = Butterworth_LP_filter(data,name,samplig_frequency,LP_cutoff,Order,show)


    %% Filter Design
    Fs = samplig_frequency;      % Sampling Frequency
    Fnyq = Fs/2;    % Nyquist Frequency
    
    % low-pass Filter
    Fco_lp = LP_cutoff;     % low-pass Cutoff Freq in Hz 
    [b_lp,a_lp] = butter(Order,Fco_lp/Fnyq, 'low');     % Butterworth Filter
    
    data_filtered = filtfilt(b_lp, a_lp,data);
    
    if show
        figure()
        
        plot(data_filtered);
        hold on
        plot(data,'r');
        hold off
        title(name);
        legend('Data filtered','Data');
    end
end