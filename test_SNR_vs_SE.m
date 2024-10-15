clc;                                                                                                                                                                                                                                                                  clc;
clear;
close all;


SNR_arr = [-10,0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100];

ES_N64_arr = zeros(1,length(SNR_arr));
ES_N100_arr = zeros(1,length(SNR_arr));
MS_N64_arr = zeros(1,length(SNR_arr));
MS_N100_arr = zeros(1,length(SNR_arr));

iteration  = 20;
% Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);

parfor times = 1:iteration
    temp_ES_64 = zeros(1,length(SNR_arr));
    temp_ES_100 = zeros(1,length(SNR_arr));
    temp_MS_64 = zeros(1,length(SNR_arr));
    temp_MS_100 = zeros(1,length(SNR_arr));
    for i = 1:length(SNR_arr)
        [SE] = pgam_ES_SNR_vs_SE(SNR_arr(i),64,(1/4));
        temp_ES_64 (i) =   temp_ES_64 (i) + SE;
        [SE] = pgam_ES_SNR_vs_SE(SNR_arr(i),100,(1/4));
        temp_ES_100 (i) =   temp_ES_100 (i) + SE;
        [SE] = pgam_MS_SNR_vs_SE(SNR_arr(i),64,(1/4));
        temp_MS_64 (i) =   temp_MS_64 (i) + SE;
        [SE] = pgam_MS_SNR_vs_SE(SNR_arr(i),100,(1/4));
        temp_MS_100 (i) =   temp_MS_100 (i) + SE;
    end
    ES_N64_arr = ES_N64_arr + temp_ES_64;
    ES_N100_arr = ES_N100_arr + temp_ES_100;
    MS_N64_arr = MS_N64_arr + temp_MS_64;
    MS_N100_arr = MS_N100_arr + temp_MS_100;
    fprintf('%d complete\n',times);
    %Progressbar
    pause(100/iteration);
    ppm.increment();

end

%Delete Progressbar
delete(ppm);

ES_N64_arr =  ES_N64_arr ./ iteration;
ES_N100_arr =  ES_N100_arr ./ iteration;
MS_N64_arr =  MS_N64_arr ./ iteration;
MS_N100_arr =  MS_N100_arr ./ iteration;

figure(1)
plot(SNR_arr,ES_N64_arr,'--o',SNR_arr,ES_N100_arr,'--o',SNR_arr,MS_N64_arr,'-*',SNR_arr,MS_N100_arr,'-*');
grid on;
colororder([1 0 0;0 0 0;1 0 0;0 0 0]);
xlabel("SNR (dB)")
ylabel("Downlink Achievable Sum SE [bit/s/Hz]")
legend('Random phase shifts (ES protocol, N=64)','Random phase shifts (ES protocol, N=100)', ...
    'Random phase shifts (MS protocol, N=64)','Random phase shifts (MS protocol, N=100)',"Location","Best")