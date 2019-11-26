function filter_signal(x)

Fs = 1;
t = 0:1:length(x)-1;

% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.01, ...
%          'StopbandFrequency',0.02,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
% fvtool(lpFilt)

% y = filter(lpFilt,x);

Hd1 = designfilt('bandpassfir', ...
    'StopbandFrequency1',0.004,'PassbandFrequency1',0.006, ...
    'PassbandFrequency2',0.08 ,'StopbandFrequency2',0.09 , ...
    'StopbandAttenuation1',10,'PassbandRipple',1, ...
    'StopbandAttenuation2',10,'DesignMethod','equiripple','SampleRate',Fs);
fvtool(Hd1)
y = filter(Hd1,x);

% y = filter(b,a,x);

figure(1000)
plot(t,x)
hold on
plot(t,y)
legend('Input Data','Filtered Data')
hold off


% designfilt

end