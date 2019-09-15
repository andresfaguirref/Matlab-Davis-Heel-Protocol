clear all;
%% options
[fname, pname] = uigetfile({'*.txt','*.mat'},'Select File');
if isequal(fname,0),
    hwarn = warndlg('no File');
end
fh = importdata(fullfile(pname, fname),'\t',5);
data_len = length(fh.data);
%% Filter Design
Fs = 1000;      % Sampling Frequency
Fnyq = Fs/2;    % Nyquist Frequency
% Artifact High-pass Filter
Fco_hp = 10;    % high-pass Cutoff Freq in Hz 
[b_hp,a_hp] = butter(2,Fco_hp/Fnyq, 'high');    % Butterworth Filter
% low-pass Filter
Fco_lp = 300;     % low-pass Cutoff Freq in Hz 
[b_lp,a_lp] = butter(2,Fco_lp/Fnyq, 'low');     % Butterworth Filter

%%
biceps = filtfilt(b_hp, a_hp,fh.data(:,1));
biceps = filtfilt(b_lp, a_lp,fh.data(:,1));
RMS_biceps = RMSmovwin(biceps,0.2*Fs);
t = (0:data_len-1)*(1/Fs);
figure;
[AX,H1,H2] = plotyy(t,[RMS_biceps], t,-fh.data(:,3));
set(get(AX(1),'ylabel'),'String', 'mV');
set(get(AX(2),'ylabel'),'String', 'Angle [°]');
xlabel('time [s]');
title('Right Arm');
legend('BB','ELBOW');
axis tight;
%%vline(1.580,'r');
%%vline(3.450,'b');
%%




