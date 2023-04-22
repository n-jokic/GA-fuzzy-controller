clear all;
close all; 
clc;

Num = [2];
Den = [1 12 24];


sim_file_name = 'PID_simulink';

open_system(sim_file_name);

y = [];

t_start = 0;
t_end  = 10;




Kp = 50;

Ki =70;

Kd = 3;
sim(sim_file_name);

steady_state = y(end);

k = find(y >= 0.9*steady_state);
Tr = time(k(1));
k = find(y >= 0.98*steady_state);
for i = 0 : length(k)-1
    if k(length(k) - i) ~= k(length(k) -i -1)
        Ts = time(k(i+1));
        break
    end
end
P = (max(y) - steady_state)/steady_state;
Peak = max(y);

f = figure();
plot(time, y); xlabel('t[s]'); ylabel('y(t)'); title('Odskocni odziv PID'); 
grid on;

figure(f);
print -deps PID;

savefig('PID.fig');

save_system(sim_file_name);
close_system(sim_file_name);

disp('Ref. Time-Span No.Osc. Amplitude     Tr         P[%]         Ts        Peak');

disp(sprintf('                              %3.4f     %10.4f      %10.4f      %10.4f',Tr, P*100, Ts, Peak));