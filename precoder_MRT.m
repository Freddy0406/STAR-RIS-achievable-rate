function [precoder_matrix] = precoder_MRT(estimated_channel,bessel_coef,tao,times,users,BS_antennas)

    % Let f_k = h_hat_k
    precoder_matrix = zeros(BS_antennas,users,times);
    for n = tao+1 : times
        for k = 1:users
            if(k<3)
                precoder_matrix(:,k,n) = estimated_channel(:,k);
            else
                precoder_matrix(:,k,n) = bessel_coef(n-users) * estimated_channel(:,k);
            end
        end
    end
end