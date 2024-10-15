clc;                                                                                                                                                                                                                                                                  clc;
clear;
close all;


N_arr = [2,3,4,5,6,7,8,9,10].^2;
iteration  = 30;
SE_arr = zeros(1,length(N_arr));
SE6_arr = zeros(1,length(N_arr));
SE_MS_arr = zeros(1,length(N_arr));
SE6_MS_arr = zeros(1,length(N_arr));

% Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);

parfor times = 1:iteration
    temp_SE = zeros(1,length(N_arr));
    temp_SE_6 = zeros(1,length(N_arr));
    temp_MS_SE = zeros(1,length(N_arr));
    temp_MS_SE_6 = zeros(1,length(N_arr));
    for i = 1:length(N_arr)
        [SE] = pgam_ES_N_vs_SE(N_arr(i),1/4);
        temp_SE (i) = SE;
        [SE] = pgam_ES_N_vs_SE(N_arr(i),1/6);
        temp_SE_6 (i) = SE;
        [SE] = pgam_MS_N_vs_SE(N_arr(i),1/4);
        temp_MS_SE (i) = SE;
        [SE] = pgam_MS_N_vs_SE(N_arr(i),1/6);
        temp_MS_SE_6 (i) = SE;
    end

    SE_arr = SE_arr + temp_SE;
    SE6_arr = SE6_arr + temp_SE_6;
    SE_MS_arr = SE_MS_arr + temp_MS_SE;
    SE6_MS_arr = SE6_MS_arr + temp_MS_SE_6;
    fprintf('%d complete\n',times);
    %Progressbar
    pause(100/iteration);
    ppm.increment();

end

%Delete Progressbar
delete(ppm);


SE_arr = SE_arr ./ iteration;
SE6_arr = SE6_arr ./ iteration;
SE_MS_arr = SE_MS_arr ./ iteration;
SE6_MS_arr = SE6_MS_arr ./ iteration;

figure(1)
plot(N_arr,SE_arr,'--o',N_arr,SE6_arr,'--o',N_arr,SE_MS_arr,'-o',N_arr,SE6_MS_arr,'-o');
grid on;
colororder([1 0 0;0 0 0;1 0 0;0 0 0]);
xlabel("Number of RIS elements (N)")
ylabel("Downlink Achievable Sum SE [bit/s/Hz]")
legend('Random phase shifts (ES protocol) (N=64, lambda/4)','Random phase shifts (ES protocol) (N=64, lambda/6)', ...
    'Random phase shifts (MS protocol) (N=64, lambda/4)','Random phase shifts (MS protocol) (N=64, lambda/6)',"Location","Best")
