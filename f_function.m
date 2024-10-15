function [signal,interference,gamma,sum_SE]=f_function(Rk,Psik,K,noise_power,signal_power)


    sum_SE = 0;
    sum_inter2 = 0;
    signal = zeros(K,1);
    interference = zeros(K,1);
    gamma = zeros(K,1);
    for i = 1:K
        sum_inter2 = sum_inter2 + trace(Psik(:,:,i));
    end

    for k = 1:K
        signal(k) = (trace(Psik(:,:,k)))^2;
        sum_inter1 = 0;
        for i = 1:K           
            sum_inter1 = sum_inter1 + trace(Rk(:,:,k)*Psik(:,:,i));
        end
        interference(k) = sum_inter1 - trace((Psik(:,:,k))^2) + ((K*noise_power)/(signal_power)) * sum_inter2;
        gamma(k) = abs(signal(k)/interference(k));
        sum_SE = sum_SE+log2(1+gamma(k));
    end

end