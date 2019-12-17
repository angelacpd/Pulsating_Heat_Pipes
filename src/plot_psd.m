function plot_psd(f, P, n, np, xlim_on, ylim_on, xmax, ymax, TF, tube_format)

ystep = ymax/5;
set(gcf, 'Position', get(0, 'Screensize'));

plot(f,P(1:np/2+1),'b','LineWidth',1)
title({tube_format;['Sensor T', num2str(n)]})
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

plot_points = 0;
if plot_points == 1
    % Get vector of local maximum points.
    a = f(TF); % Frequency
    b = P(TF); % Power
    l = min([20, length(a)]);
    format_a = zeros(l,1);
    format_b = zeros(l,1);
    for k = 1:l
        % Plot local maximum points
        plot(a(k),b(k),'r*')
        format_a(k) = round(a(k)*10000)/10000;
        format_b(k) = round(b(k));
        txt = ['\leftarrow (' num2str(format_a(k)) ', ' num2str(format_b(k)) ')'];
        % Write the coordinates of local maximum points on graphic.
        text(a(k),b(k),txt,'FontSize',14)
    end
end

% Set font size of graphic axes.
set(gca,'FontSize', 16)

end