d1 = load_LWP_data('../examples/p1.dat');
d2 = load_LWP_data('../examples/p2.dat');
d3 = load_LWP_data('../examples/p3.dat');

figure(1)
subplot(2,1,1)
plot(d1.t,d1.p_head,'r'), hold on
plot(d1.t,d1.p_tail,'r--'), hold on
plot(d2.t,d2.p_head,'k-'), hold on
plot(d2.t,d2.p_tail,'k--'), hold on
plot(d3.t,d3.p_head,'b-'), hold on
plot(d3.t,d3.p_tail,'b--'), hold off
xlabel('t, s'), ylabel('p, bar') 

subplot(2,1,2)
plot(d1.t,d1.v_head,'r'), hold on
plot(d1.t,d1.v_tail,'r--'), hold on
plot(d2.t,d2.v_head,'k-'), hold on
plot(d2.t,d2.v_tail,'k--'), hold on
plot(d3.t,d3.v_head,'b-'), hold on
plot(d3.t,d3.v_tail,'b--'), hold off
xlabel('t, s'), ylabel('v, m/s') 