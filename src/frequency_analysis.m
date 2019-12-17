%% FREQUENCY ANALYSIS OF PULSATING HEAT PIPES

% Sensor placement:
% Group 1 - Condenser
% T1 - T5
% Group 2 - Adiabatic
% T6 - T8
% Group 3 - Evaporator
% T9 - T13

% Channel separation:
% T1, T6, T9
% T3, T7, T11
% T5, T8, T13

% Power profile:
% Start: 0 W for 10 seconds.
% Power step: period 900 seconds.
% End: 0 W for 100 seconds.

% Heat Pipes Laboratory
% Federal University of Santa Catarina
% Florianopolis - Brazil

tic
%% Clear variables; clear command window; close windows.
clearvars; clc; close all;

%% Inputs

% Data to be shown on the left side of the graphic.
filename_l = 'inverted_circular_26canais_3mm_RP50.txt';
% Data to be shown on the right side of the graphic.
filename_r = 'inverted_grooved_26canais_3mm_RP50.txt';
delimiterIn = ' ';
headerlinesIn = 1;

xlim_on = 1; % 0 = Off. 1 = On.
ylim_on = 1; % 0 = Off. 1 = On.
xmax = 0.1; % Limit x axis.

Fs = 1; % Sampling frequency
n_sensors = 13; % Number of sensors

steady_len = 300; % Get the last samples of each power level.

%% Import text files.

% Left
data_l = importdata(filename_l,delimiterIn,headerlinesIn);
% Split filename to get tube format.
filename_l_split = split(filename_l,'_');
tube_format_l_cell = filename_l_split(2);
tube_format_l = cell2mat(tube_format_l_cell);
% Verify if the tube in the left is circular.
tube_l = strcmp(tube_format_l, 'circular');

% Right
data_r = importdata(filename_r,delimiterIn,headerlinesIn);
% Split filename to get tube format.
filename_r_split = split(filename_r,'_');
tube_format_r_cell = filename_r_split(2);
tube_format_r = cell2mat(tube_format_r_cell);
% Verify if the tube in the right is grooved.
tube_r = strcmp(tube_format_r, 'grooved');

% Throw error alert if files do not match and stop script.
if tube_l == 0 || tube_r == 0
   msg = 'Input error. Tube 1 shall be circular and Tube 2 shall be grooved.';
   error(msg)
end

tube_for_l = "Round";
tube_for_r = "Grooved";

% Get filling ratio.
rp_cell = split(filename_l_split(5),'.');
rp = cell2mat(rp_cell(1));

% Create folder to store graphics
position = cell2mat(filename_l_split(1));
folder_name = [position, '_', rp, '_steady'];
[status, msg, msgID] = mkdir(['graphics\',folder_name]);
store_path = ['graphics\',folder_name];

time_v_l = data_l.data(2:end,1);      % Time vector.
power_v_l = data_l.data(2:end,24);    % Power vector.

length_l = length(time_v_l);    % Length of time vector.
temp_l = zeros(length_l,n_sensors);     % Temperature matrix declaration.

time_v_r = data_r.data(2:end,1);      % Time vector.
power_v_r = data_r.data(2:end,24);    % Power vector.

length_r = length(time_v_r);    % Length of time vector.
temp_r = zeros(length_r,n_sensors);     % Temperature matrix declaration.

%% Temperature matrix.
for i = 1:n_sensors
    temp_l(:,i) = data_l.data(2:end,i+2);
    temp_r(:,i) = data_r.data(2:end,i+2);
    
%     f1 = 100 + (2*i - 1);
%     figure(f1)
%     set(f1, 'Position', get(0, 'Screensize'));
%     savename = [tube_format_l, ' sensor T', num2str(i)];
%     title(savename)
%     yyaxis left
%     plot(time_v_l,temp_l(:,i))
%     ylabel('Temperature [\circC]')
%     yyaxis right
%     plot(time_v_l,power_v_l)
%     xlabel('Time [s]')
%     ylabel('Power [W]')
%     grid on
%     set(gca,'FontSize',16)    
%     saveas(f1,fullfile(store_path, savename),'png')
%     close
%     
%     f2 = 100 + 2*i;
%     figure(f2)
%     set(f2, 'Position', get(0, 'Screensize'));
%     savename = [tube_format_r, ' sensor T', num2str(i)];
%     title(savename)
%     yyaxis left
%     plot(time_v_r,temp_r(:,i))
%     ylabel('Temperature [\circC]')
%     yyaxis right
%     plot(time_v_r,power_v_r)
%     xlabel('Time [s]')
%     ylabel('Power [W]')
%     grid on
%     set(gca,'FontSize',16)   
%     saveas(f2,fullfile(store_path, savename),'png')
%     close
    
    %filter_signal(temp_l(:,i))
    %filter_signal(temp_r(:,i))
end

%% Split temperature vectors according to the power level.
j = 1;
m = 1;
k = 0;
flag = 0;

% xi = 0;
% zi = 0;
fi = 0;

for n = 1:n_sensors
    
    while (flag == 0)
        k = k+1;
        
        if k > length_l
            flag = 1;
            break
        end
        
        if power_v_l(k) > 95
            level_l(j,m) = temp_l(k,n);
            level_r(j,m) = temp_r(k,n);
            temp_power_l(k,:) = [temp_l(k,n), power_v_l(k)];
            
            j = j+1;
            
            if (power_v_l(k) ~= power_v_l(k-1))
                if level_l(3:end,m) == 0
                    level_l(:,m) = [];
                    level_r(:,m) = [];
                else
                    m = m + 1;
                    power(m) = power_v_l(k+3);
                end                            
                j = 1;
            end
        end
    end
    
    level_l(:,1) = [];
    power(1) = [];
    level_r(:,1) = [];
%    power_r(1) = [];

    [max_power, max_power_index] = max(power);
    
    n_levels = size(level_l);
    for p = 1:n_levels(2)
        
        if p <= max_power_index
            ascdesc = 'ASC';
        else
            ascdesc = 'DESC';
        end
        
        % Time domain figures
        fi = fi + 200; 
        ftx = 200 + 2*p - 2;
        ftz = 200 + 2*p - 1;
        
        %% Compute Power Spectrum Density 
        
        
        np = 2^nextpow2(length(level_l(steady_len:end-1,p)));
        f = Fs*(0:(np/2))/np;
        
        % Left
        Xp1 = level_l(steady_len:end-1,p);
        Xp = Xp1 - mean(Xp1);               
        Y = fft(Xp);
        P2 = (1/(Fs*np))*abs(Y).^2;
        P1 = P2(1:np/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
                     
        
        %filter_signal(Xp)
        
        % Right
        Zp1 = level_r(steady_len:end-1,p);
        Zp = Zp1 - mean(Zp1);
        W = fft(Zp);
        P4 = (1/(Fs*np))*abs(W).^2;
        P3 = P4(1:np/2+1);
        P3(2:end-1) = 2*P3(2:end-1);
        
        
        
        % Normalize graphics
        [xp1max,xp1max_ind] = max(Xp1);
        [zp1max,zp1max_ind] = max(Zp1);
        if xp1max > zp1max
            yt_max = xp1max;
        else
            yt_max = zp1max;
        end

        yt_max = 5*ceil(yt_max/5);
        
        [xp1min,xp1min_ind] = min(Xp1);
        [zp1min,zp1min_ind] = min(Zp1);
        if xp1min < zp1min
            yt_min = xp1min;
        else
            yt_min = zp1min;
        end

        yt_min = 5*floor(yt_min/5);
   
        figure(ftx)
        set(ftx, 'Position', get(0, 'Screensize'));
        savename = ['Round - Sensor T', num2str(n), ' - Power ', num2str(floor(power(p))), ' W - ', ascdesc];      
        title_str = ['Round - Sensor T', num2str(n), ' - Power ', num2str(floor(power(p))), ' W'];
        plot_time_domain_seg(Xp1, ftx, title_str, savename, store_path, yt_min, yt_max)
        close 
        
        figure(ftz)
        set(ftz, 'Position', get(0, 'Screensize'));
        savename = ['Grooved - Sensor T', num2str(n), ' - Power ', num2str(floor(power(p))), ' W - ', ascdesc];  
        title_str = ['Grooved - Sensor T', num2str(n), ' - Power ', num2str(floor(power(p))), ' W']; 
        plot_time_domain_seg(Zp1, ftz, title_str, savename, store_path, yt_min, yt_max)
        close
        
        %filter_signal(Zp)
        
        %% Plot PSD
        
        % Normalize graphics
        [P1max,P1max_ind] = max(P1(1:np/2+1));
        [P3max,P3max_ind] = max(P3(1:np/2+1));
        if P1max > P3max
            Pmax = P1max;
        else
            Pmax = P3max;
        end
%         ymax = 100*ceil(Pmax/100);
        ymax = 10*ceil(Pmax/10);
        
        
        % Compute local maximum points
        TF1 = islocalmax(P1);
        TF3 = islocalmax(P3);
        
        if n == 1
            h1 = p;
            figure(h1)
                       
            subplot(3,2,1)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
           
            subplot(3,2,2)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)         
            
        elseif n == 6
            h1 = p;
            figure(h1)
            
            subplot(3,2,3)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,4)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
        elseif n == 9
            h1 = p;
            figure(h1)
            
            subplot(3,2,5)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,6)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
                               
            savename = ['Group 169 - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
%             saveas(h1,fullfile(store_path, savename),'fig')
            saveas(h1,fullfile(store_path, savename),'png')
            close
        end

        
        if n == 3
            h2 = p+20;
            figure(h2)
            
            subplot(3,2,1)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,2)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
        elseif n == 7
            h2 = p+20;
            figure(h2)
            
            subplot(3,2,3)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,4)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
        elseif n == 11
            h2 = p+20;
            figure(h2)

            subplot(3,2,5)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,6)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
            savename = ['Group 3711 - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
%             saveas(h2,fullfile(store_path, savename),'fig')
            saveas(h2,fullfile(store_path, savename),'png')
            close
        end

        
        if n == 5
            h3 = p+40;
            figure(h3)
            
            subplot(3,2,1)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,2)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
        elseif n == 8
            h3 = p+40;
            figure(h3)
            
            subplot(3,2,3)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,4)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
        elseif n == 13
            h3 = p+40;
            figure(h3)
            
            subplot(3,2,5)
            plot_psd(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, TF1, tube_for_l)
            
            subplot(3,2,6)
            plot_psd(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, TF3, tube_for_r)
            
            savename = ['Group 5813 - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
%             saveas(h3,fullfile(store_path, savename),'fig')
            saveas(h3,fullfile(store_path, savename),'png')
            close
            
        end

        
    end
    
    j = 1;
    m = 1;
    k = 0;
    flag = 0;
    clear level_l
    clear level_r 
    clear power
    clear temp_power
           
end
toc