%% CÁLCULO DAS FREQUÊNCIAS
% Variação de potência no tempo
% Início: 0 W por 10 segundos.
% Fim: 0  W por 100 segundos.
% Demais potências: duração de 900 segundos.
% Período (s) - Potência (W)



% Grupo 1 - Condensador
% T1 - T5
% Grupo 2 - Adiabática
% T6 - T8
% Grupo 3 - Evaporador
% T9 - T13

% Separação por canais.
% T1, T6, T9
% T3, T7, T11
% T5, T8, T13

% Cada patamar de potência tem duração de 900 ms

%% Clean everything.
clear all; clc; close all;

%% Specify the parameters of a signal with a sampling frequency (Fs) of 1 Hz and a signal duration of N seconds.


Fs = 1;            % Sampling frequency
T = 1/Fs;             % Sampling period

%% Import text file.

filename = 'horizontal_circular_26canais_3mm_RP50.txt';
delimiterIn = ' ';
headerlinesIn = 1;
B = importdata(filename,delimiterIn,headerlinesIn);

filename_split = split(filename,'_');
tube_format_cell = filename_split(2);
tube_format = cell2mat(tube_format_cell);


nome = B.colheaders;
tempo = B.data(2:end,1);
t = tempo - tempo(2);
T_amb = B.data(2:end,2);
Pot = B.data(2:end,24);
n_sensores = 13;
L = length(tempo);
Temp = zeros(L,n_sensores);


for i = 1:1:n_sensores
    Temp(:,i) = B.data(2:end,i+2);
end

figure(100)
yyaxis left
plot(tempo,Temp(:,6))
yyaxis right
%hold on
plot(tempo,Pot)
grid on


%% Patamares de T1
j = 1;
m = 1;
k = 0;
flag = 0;

for n = 1:n_sensores
%for n = 1
    
    while (flag == 0)
        k = k+1;
        
        if k > L
            flag = 1;
            break
        end
        
        if Pot(k) > 95
            Patamar(j,m) = Temp(k,n);
            Temp_Pot(k,:) = [Temp(k,n), Pot(k)];
            j = j+1;

            if (Pot(k) ~= Pot(k-1)) 
                if Patamar(5:end,m) == 0
                    Patamar(:,m) = [];
                else
                    m = m + 1;
                    Power(m) = Pot(k+3);
                end
                j = 1;               
            end
        end    
    end

Patamar(:,1) = [];
Power(1) = [];

n_patamares = size(Patamar);
for p = 1:n_patamares(2)
    Xp = Patamar(:,p);    
    Xp = Xp - mean(Xp);
    np = 2^nextpow2(length(Patamar(:,p)));
    f = Fs*(0:(np/2))/np;
    %Y = fft(Xp,np);
    Y = fft(Xp);
    %P = abs((Y.^2)/n);
    P2 = (1/(Fs*np))*abs(Y).^2;
    P1 = P2(1:np/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    
    
         if n == 1
            figure(p)
            subplot(3,2,1)
            plot(f,P1(1:np/2+1),'b')            
            title({['Circular', ' - Power: ', num2str(floor(Power(p))), ' W'];['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            ymax = max(P1(1:np/2+1));
            % ylim([0 100*ceil(ymax)])
            
        elseif n == 6
            figure(p)
            subplot(3,2,3)
            plot(f,P1(1:np/2+1),'b')
            %title(['Sensor T', num2str(n), ' - Power: ', num2str(floor(Power(p))), ' W'])
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 9
            figure(p)
            subplot(3,2,5)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
         end
        
         
         if n == 3
            figure(p+20)
            subplot(3,2,1)
            plot(f,P1(1:np/2+1),'b')
            title({['Circular', ' - Power: ', num2str(floor(Power(p))), ' W'];['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 7
            figure(p+20)
            subplot(3,2,3)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 11
            figure(p+20)
            subplot(3,2,5)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
         end
        
         
         if n == 5
            figure(p+40)
            subplot(3,2,1)
            plot(f,P1(1:np/2+1),'b')
            title({['Circular', ' - Power: ', num2str(floor(Power(p))), ' W'];['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 8
            figure(p+40)
            subplot(3,2,3)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 13
            figure(p+40)
            subplot(3,2,5)
            plot(f,P1(1:np/2+1),'b')
            ylim([0 800])
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        end
    
    
end    
    
j = 1;
m = 1;
k = 0;
flag = 0;
clear Patamar
clear Power
clear Temp_Pot
    
    
    
end


clear all;

%% Import text file.
filename = 'horizontal_ranhurado_26canais_3mm_RP50.txt';
delimiterIn = ' ';
headerlinesIn = 1;
B = importdata(filename,delimiterIn,headerlinesIn);

Fs = 1;            % Sampling frequency
T = 1/Fs;             % Sampling period
nome = B.colheaders;
tempo = B.data(2:end,1);
t = tempo - tempo(2);
T_amb = B.data(2:end,2);
Pot = B.data(2:end,24);
n_sensores = 13;
L = length(tempo);
Temp = zeros(L,n_sensores);


for i = 1:1:n_sensores
    Temp(:,i) = B.data(2:end,i+2);
end

figure(200)
yyaxis left
plot(tempo,Temp(:,6))
yyaxis right
%hold on
plot(tempo,Pot)
grid on


%% Patamares de T1
j = 1;
m = 1;
k = 0;
flag = 0;

for n = 1:n_sensores
%for n = 1
    
    while (flag == 0)
        k = k+1;
        
        if k > L
            flag = 1;
            break
        end
        
        if Pot(k) > 95
            Patamar(j,m) = Temp(k,n);
            Temp_Pot(k,:) = [Temp(k,n), Pot(k)];
            j = j+1;

            if (Pot(k) ~= Pot(k-1)) 
                if Patamar(5:end,m) == 0
                    Patamar(:,m) = [];
                else
                    m = m + 1;
                    Power(m) = Pot(k+3);
                end
                j = 1;               
            end
        end    
    end

Patamar(:,1) = [];
Power(1) = [];

n_patamares = size(Patamar);
for p = 1:n_patamares(2)
    Xp = Patamar(:,p);    
    Xp = Xp - mean(Xp);
    np = 2^nextpow2(length(Patamar(:,p)));
    f = Fs*(0:(np/2))/np;
    %Y = fft(Xp,np);
    Y = fft(Xp);
    %P = abs((Y.^2)/n);
    P2 = (1/(Fs*np))*abs(Y).^2;
    P1 = P2(1:np/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    
    
         if n == 1
            figure(p)
            subplot(3,2,2)
            plot(f,P1(1:np/2+1),'b')
            
            title({['Grooved', ' - Power: ', num2str(floor(Power(p))), ' W'];['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 6
            figure(p)
            subplot(3,2,4)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 9
            figure(p)
            subplot(3,2,6)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
         end
        
         
         if n == 3
            figure(p+20)
            subplot(3,2,2)
            plot(f,P1(1:np/2+1),'b')
            title({['Grooved', ' - Power: ', num2str(floor(Power(p))), ' W'];['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 7
            figure(p+20)
            subplot(3,2,4)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 11
            figure(p+20)
            subplot(3,2,6)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
         end
        
         
         if n == 5
            figure(p+40)
            subplot(3,2,2)
            plot(f,P1(1:np/2+1),'b')
            title({['Grooved', ' - Power: ', num2str(floor(Power(p))), ' W'];['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 8
            figure(p+40)
            subplot(3,2,4)
            plot(f,P1(1:np/2+1),'b')
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        elseif n == 13
            figure(p+40)
            subplot(3,2,6)
            plot(f,P1(1:np/2+1),'b')
            ylim([0 800])
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            grid on
            
        end
    
    
end    
    
j = 1;
m = 1;
k = 0;
flag = 0;
clear Patamar
clear Power
clear Temp_Pot
    
    
    
end
