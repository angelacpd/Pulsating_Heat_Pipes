function plot_time_domain_seg(X, fi, title_str, savename, store_path, ymin, ymax)

ystep = (ymax-ymin)/10;
xmax = length(X);
plot(X,'b','LineWidth',1)
title(title_str)
xlabel('Time [s]')
ylabel('Temperature [\circC]')
ylim([ymin ymax])
xticks(0:50:xmax)
yticks(ymin:ystep:ymax)
grid on
set(gca,'FontSize',16)
saveas(fi, fullfile(store_path, savename),'png')

end