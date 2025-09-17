function plot_LWP_data(fname)

%addpath('/Users/hoscsaba/Documents/GitHub/PSToolbox/Matlab_tools')

d=load_LWP_data(fname);

figure
subplot(2,2,1), plot(d.t,d.p_head,'b',d.t,d.p_tail,'r')
ylabel('p, bar'), grid on, legend('head','tail')
subplot(2,2,2), plot(d.t,d.v_head,'b',d.t,d.v_tail,'r')
ylabel('v, m/s'), grid on
subplot(2,2,3), plot(d.t,d.T_head,'b',d.t,d.T_tail,'r')
ylabel('T, C'), grid on, xlabel('t, s')
subplot(2,2,4), plot(d.t,d.mp_head,'b',d.t,d.mp_tail,'r')
ylabel('m dot, kg,s'), grid on, xlabel('t, s')
end