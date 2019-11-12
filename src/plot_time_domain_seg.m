function plot_time_domain_seg(X, fi, savename, store_path)

plot(X)
title(savename)
xlabel('Time [s]')
ylabel('Temperature [\circC]')
grid on
set(gca,'FontSize',16)
saveas(fi, fullfile(store_path, savename),'png')
%close

end