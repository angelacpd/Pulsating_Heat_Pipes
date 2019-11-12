function plot_psd_l(f, P1, n, np, xlim_on, ylim_on, xmax, ymax, ystep, TF1)

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
a = f(TF1);
b = P1(TF1);
hold on
% plot(f(P1max_ind),P1max,'r*')
% plot(f(TF1),P1(TF1),'r*')
% plot(a,b,'r*')
l = min([6, length(a)]);
format_a = zeros(l,1);
format_b = zeros(l,1);
for k = 1:l
    plot(a(k),b(k),'r*')
    format_a(k) = round(a(k)*10000)/10000;
    format_b(k) = round(b(k));
    txt = ['\leftarrow (' num2str(format_a(k)) ', ' num2str(format_b(k)) ')'];    
    text(a(k),b(k),txt,'FontSize',14)
end
% legend(num2str(a(1)))
set(gca,'FontSize', 16)

end