function [estimated_channel,Qk] = channel_ES(bessel_coef,agg_channel,tao,user_num,BS_antenna,times,noise_var,R_channel,pilot_P)


    % Channel estimation 只量測tao+1個時間點
    tao = tao+1;
    A = randn(tao, user_num);
    pilot_sq = orth(A);
    pilot_sq = pilot_sq*sqrt(tao*pilot_P);

    %% Recieved at BS

    Y = zeros(BS_antenna,tao,tao);

    for n = 1:tao
        for i = 1:user_num       
            Z = complex_gaussian_gen(BS_antenna,tao,zeros(BS_antenna,1),noise_var*eye(BS_antenna)); 
            Y(:,:,n) = Y(:,:,n) + agg_channel(:,i,n)*pilot_sq(:,i)'+Z; 
        end
    end



    %% Correlate with pilot sequence
    
    rk = zeros(BS_antenna,user_num,times);

    for n = 1:tao
        for i = 1:user_num
            rk(:,i,n) = (1/sqrt(tao*pilot_P)) * Y(:,:,n) * pilot_sq(:,i);
        end
    end
    %% LMMSE estimation

    estimated_channel = zeros(size(agg_channel));
    Qk = zeros(BS_antenna,BS_antenna,user_num);

    for k = 1:user_num
        Qk(:,:,k) = pinv(R_channel(:,:,k)+(noise_var/(tao*pilot_P)*eye(BS_antenna)));
        for n = 1:tao
            if(k<3)
                estimated_channel(:,k,n) = R_channel(:,:,k)*Qk(:,:,k)*rk(:,k,n);
                % estimated_channel(:,k,n) = eye(length(R_channel(:,:,k)))*rk(:,k,n);
            else
                estimated_channel(:,k,n) = bessel_coef(tao+1-n)*R_channel(:,:,k)*Qk(:,:,k)*rk(:,k,n);
                % estimated_channel(:,k,n) = bessel_coef(tao+1-n)*eye(length(R_channel(:,:,k)))*rk(:,k,n);
            end 
        end
    end


end