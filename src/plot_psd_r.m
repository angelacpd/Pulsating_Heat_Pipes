function plot_psd_r(f, P3, n, np, xlim_on, ylim_on, xmax, ymax, ystep, TF3)

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
a = f(TF3);
b = P3(TF3);
hold on
% plot(f(P3max_ind),P3max,'r*')
% plot(f(TF3),P3(TF3),'r*')

l = min([6, length(a)]);
format_a = zeros(l,1);
format_b = zeros(l,1);
for k = 1:l
    plot(a(k),b(k),'r*')
    format_a(k) = round(a(k)*10000)/10000;
    format_b(k) = round(b(k));
    txt = ['\leftarrow (' num2str(format_a(k)) ', ' num2str(format_b(k)) ')'];    
    text(a(k),b(k),txt,'FontSize',14)
    %text(a(k),b(k),txt,'FontSize',14)
end
set(gca,'FontSize',16)

end