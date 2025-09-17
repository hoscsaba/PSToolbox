function plot_valve_data(fname)

%addpath('/Users/hoscsaba/Documents/GitHub/PSToolbox/Matlab_tools')

d=load_valve_data(fname)

figure
subplot(3,1,1), plot(d.t,d.x)
ylabel('x, m'), grid on
subplot(3,1,2), plot(d.t,d.v)
ylabel('v, m/s'), grid on
subplot(3,1,3), plot(d.t,d.mp)

ylabel('m dot, kg,s'), grid on, xlabel('t, s')
end