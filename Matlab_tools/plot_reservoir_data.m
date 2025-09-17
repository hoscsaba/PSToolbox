function plot_reservoir_data(fname)

%addpath('/Users/hoscsaba/Documents/GitHub/PSToolbox/Matlab_tools')

d=load_reservoir_data(fname)

figure
subplot(2,1,1), plot(d.t,d.p)
ylabel('p, bar'), grid on
subplot(2,1,2), plot(d.t,d.mp_in,'r-',d.t,d.mp_out,'k')
xlabel('t, s'), ylabel('mp, kg/s'), grid on
legend('inlet','outlet')
end