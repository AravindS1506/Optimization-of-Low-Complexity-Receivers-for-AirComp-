O=2:2:10;
time_zf=[5.968543783	5.849813437	5.878325996	6.098526723	5.834788009];
time_zfas=[10.06909235	8.831538783	7.713535383	7.12003781	5.540798308];
figure;
plot(O,time_zfas); hold on;
plot(O,time_zf);
xlabel("Number of RF chains O");
ylabel("Time (s)");
legend("ZF+AS algorithm","ZF+Our Algorithm");

