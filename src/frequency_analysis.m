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
%% Clear workspace; clear command window; close windows.
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
xmax = 0.1;

Fs = 1;     % Sampling frequency
n_sensors = 13;     % Number of sensors

%% Import text files.

% Left
data_l = importdata(filename_l,delimiterIn,headerlinesIn);
% Split filename to get tube format.
filename_l_split = split(filename_l,'_');
tube_format_l_cell = filename_l_split(2);
tube_format_l = cell2mat(tube_format_l_cell);

% Right
data_r = importdata(filename_r,delimiterIn,headerlinesIn);
% Split filename to get tube format.
filename_r_split = split(filename_r,'_');
tube_format_r_cell = filename_r_split(2);
tube_format_r = cell2mat(tube_format_r_cell);

%if tube_format_l ~= 'circular' && tube_format_r ~= 'ranhurado'
%    msg = 'File error. Tube 1 shall be circular and Tube 2 shall be grooved.';
%    error(msg)
%end

rp_cell = split(filename_l_split(5),'.');
rp = cell2mat(rp_cell(1));

% Create folder to store graphics
position = cell2mat(filename_l_split(1));
folder_name = [position, '_', rp, '_png_steady'];
[status, msg, msgID] = mkdir(['graphics\',folder_name]);
store_path = ['graphics\',folder_name];

%nome = data_l.colheaders;
%t = tempo - tempo(2);
%T_amb = data_l.data(2:end,2);

time_v_l = data_l.data(2:end,1);      % Time vector.
power_v_l = data_l.data(2:end,24);    % Power vector.

length_l = length(time_v_l);    % Length of time vector.
temp_l = zeros(length_l,n_sensors);     % Temperature matrix declaration.

time_v_r = data_r.data(2:end,1);      % Time vector.
power_v_r = data_r.data(2:end,24);    % Power vector.

length_r = length(time_v_r);    % Length of time vector.
temp_r = zeros(length_r,n_sensors);     % Temperature matrix declaration.

% Temperature matrix.
for i = 1:n_sensors
    temp_l(:,i) = data_l.data(2:end,i+2);
    temp_r(:,i) = data_r.data(2:end,i+2);
    
    f1 = 100 + (2*i - 1);
    figure(f1)
    set(f1, 'Position', get(0, 'Screensize'));
    savename = [tube_format_l, ' sensor T', num2str(i)];
    title(savename)
    yyaxis left
    plot(time_v_l,temp_l(:,i))
    ylabel('Temperature [\circC]')
    yyaxis right
    plot(time_v_l,power_v_l)
    xlabel('Time [s]')
    ylabel('Power [W]')
    grid on
    set(gca,'FontSize',16)
    
    saveas(f1,fullfile(store_path, savename),'png')
    close
    
    f2 = 100 + 2*i;
    figure(f2)
    set(f2, 'Position', get(0, 'Screensize'));
    savename = [tube_format_r, ' sensor T', num2str(i)];
    title(savename)
    yyaxis left
    plot(time_v_r,temp_r(:,i))
    ylabel('Temperature [\circC]')
    yyaxis right
    plot(time_v_r,power_v_r)
    xlabel('Time [s]')
    ylabel('Power [W]')
    grid on
    set(gca,'FontSize',16)   
    saveas(f2,fullfile(store_path, savename),'png')
    close
end


%% Split temperature vectors according to the power level.
j = 1;
m = 1;
k = 0;
flag = 0;

xi = 0;
zi = 0;

f = 1;
g = 1;
h = 0;
flag_r = 0;

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
                
%                 if level_r(5:end,m) == 0
%                     level_r(:,m) = [];
%                 else
%                     m = m + 1;
%                     power_r(m) = power_v_r(k+3);
%                 end
                
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
        
        % Power Spectrum Density 
        % Left
        Xp = level_l(200:end-1,p);
        Xp = Xp - mean(Xp);
        np = 2^nextpow2(length(level_l(:,p)));
        f = Fs*(0:(np/2))/np;
        Y = fft(Xp);
        P2 = (1/(Fs*np))*abs(Y).^2;
        P1 = P2(1:np/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        xi = xi + 1;
        xp1 = 200 + xi;
        figure(xp1)
        set(xp1, 'Position', get(0, 'Screensize'));
        savename = [tube_format_l, ' sensor T', num2str(n), ' - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
        %title([tube_format_l,' Sensor ', num2str(n), ' - Power: ', num2str(floor(power(p))), ' W'])       
        plot(Xp)
        title(savename)
        xlabel('Time [s]')
        ylabel('Temperature [\circC]')
        grid on
        set(gca,'FontSize',16)       
        saveas(xp1,fullfile(store_path, savename),'png')
        close
        
        % Right
        Zp = level_r(200:end-1,p);
        Zp = Zp - mean(Zp);
        %np = 2^nextpow2(length(level_l(:,p)));
        %f = Fs*(0:(np/2))/np;
        W = fft(Zp);
        P4 = (1/(Fs*np))*abs(W).^2;
        P3 = P4(1:np/2+1);
        P3(2:end-1) = 2*P3(2:end-1);
        
        zi = zi + 1;
        zp1 = 200 + zi;
        figure(zp1)
        set(zp1, 'Position', get(0, 'Screensize'));
        savename = [tube_format_r, ' sensor T', num2str(n), ' - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
        %title([tube_format_r,' Sensor ', num2str(n), ' - Power: ', num2str(floor(power(p))), ' W'])       
        plot(Zp)
        title(savename)
        xlabel('Time [s]')
        ylabel('Temperature [\circC]')
        grid on
        set(gca,'FontSize',16)       
        saveas(zp1,fullfile(store_path, savename),'png')
        close
        
        % Normalize graphic
        [P1max,P1max_ind] = max(P1(1:np/2+1));
        [P3max,P3max_ind] = max(P3(1:np/2+1));
        if P1max > P3max
            Pmax = P1max;
        else
            Pmax = P3max;
        end
        ymax = 100*ceil(Pmax/100);
        ystep = ymax/5;
        
        
        % Plot
        if n == 1
            h1 = p;
            figure(h1)
            %set(gcf, 'Position', get(0, 'Screensize'),'Visible','off');
            set(gcf, 'Position', get(0, 'Screensize'));
            
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p)
            subplot(3,2,1)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            %title({['Circular', ' - Power: ', num2str(floor(power(p))), ' W'];['Sensor T', num2str(n)]}, 'Fontsize',20)
            title({'Circular';['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            %xticklabels({})
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            hold on
            plot(f(P1max_ind),P1max,'r*')
            set(gca,'FontSize', 16)
            
            subplot(3,2,2)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)            
            title({'Grooved';['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            hold on
            plot(f(P3max_ind),P3max,'r*')
            set(gca,'FontSize',16)            
            
        elseif n == 6
            h1 = p;
            figure(h1)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p)
            subplot(3,2,3)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,4)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end          
            set(gca,'FontSize',16)
            
        elseif n == 9
            h1 = p;
            figure(h1)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p)
            subplot(3,2,5)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,6)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
                               
            savename = ['Group 169 - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
            saveas(h1,fullfile(store_path, savename),'png')
            close
        end
        %suptitle(['Power: ', num2str(floor(power(p))), ' W'])
        
        if n == 3
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p+20)
            h2 = p+20;
            figure(h2)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            subplot(3,2,1)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title({'Circular';['Sensor T', num2str(n)]})         
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,2)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title({'Grooved';['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
        elseif n == 7
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p+20)
%             figure(p+20)
            h2 = p+20;
            figure(h2)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            subplot(3,2,3)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,4)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
        elseif n == 11
%             figure(p+20)
            h2 = p+20;
            figure(h2)
            set(gcf, 'Position', get(0, 'Screensize'));
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p+20)
            subplot(3,2,5)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,6)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            savename = ['Group 3711 - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
            saveas(h2,fullfile(store_path, savename),'png')
            close
        end
        %suptitle(['Power: ', num2str(floor(power(p))), ' W'])
        
        if n == 5
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p+40)
            h3 = p+40;
            figure(h3)
%             figure(p+40)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            subplot(3,2,1)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title({'Circular';['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,2)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title({'Grooved';['Sensor T', num2str(n)]})
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
        elseif n == 8
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p+40)
            h3 = p+40;
            figure(h3)
%             figure(p+40)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            subplot(3,2,3)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,4)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
        elseif n == 13
            %figure('Name',[folder_name,' - Power: ', num2str(floor(power(p))), ' W'],'NumberTitle',p+40)
            h3 = p+40;
            figure(h3)
%             figure(p+40)
            set(gcf, 'Position', get(0, 'Screensize'));
            
            subplot(3,2,5)
            plot(f,P1(1:np/2+1),'b','LineWidth',1)
            ylim([0 800])
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            subplot(3,2,6)
            plot(f,P3(1:np/2+1),'b','LineWidth',1)
            ylim([0 800])
            title(['Sensor T', num2str(n)])
            xlabel('Frequency [Hz]')
            ylabel('Magnitude')
            xticks(0:xmax/10:xmax)
            yticks(0:ystep:ymax)
            grid on
            if xlim_on == 1
                xlim([0 xmax])
            end
            if ylim_on == 1
                ylim([0 ymax])
            end
            set(gca,'FontSize',16)
            
            savename = ['Group 5813 - Power ', num2str(floor(power(p))), ' W - ', ascdesc];
            saveas(h3,fullfile(store_path, savename),'png')
            close
            
        end
        %suptitle(['Power: ', num2str(floor(power(p))), ' W'])
        
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