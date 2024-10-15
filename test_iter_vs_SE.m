clc;                                                                                                                                                                                                                                                                  clc;
clear;
close all;

iteration  = 5;
SE_arr = zeros(iteration,200);


% Progressbar
% ppm = ParforProgressbar(iteration);
% ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);

for times = 1:iteration
    temp_SE = zeros(iteration,200);
    [SE] = pgam_main_iter_vs_SE;
    temp_SE (times,:) = SE;

    SE_arr = SE_arr + temp_SE;
    fprintf('%d complete\n',times);
    %Progressbar
    % pause(100/iteration);
    % ppm.increment();

end

%Delete Progressbar
% delete(ppm);

x_1 = find(SE_arr(1,:));
x_2 = find(SE_arr(2,:));
x_3 = find(SE_arr(3,:));
x_4 = find(SE_arr(4,:));
x_5 = find(SE_arr(5,:));

max_1 = max(x_1);
max_2 = max(x_2);
max_3 = max(x_3);
max_4 = max(x_4);
max_5 = max(x_5);

% SE_arr = SE_arr ./ iteration;
% 
figure(1)
plot(x_1,SE_arr(1,(1:max_1)),'-o',x_2,SE_arr(2,(1:max_2)),'-x',x_3,SE_arr(3,(1:max_3)),'-^',x_4,SE_arr(4,(1:max_4)),'-d',x_5,SE_arr(5,(1:max_5)),'-o');
grid on;
xlim([0 60])
colororder([1 0 0;0 0 0;0 0 1;0 0 0;0 0 0]);
xlabel("Number of iterations")
ylabel("Downlink Achievable Sum SE [bit/s/Hz]")