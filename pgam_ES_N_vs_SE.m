function [final_SE] = pgam_ES_N_vs_SE(RIS_element,wavelen_coef)
    %% Configuration
    
    d_0 = 20;
    Kt = 2;     % T zone user
    Kr = 2;     % R zone user
    K = Kt + Kr;
    SNR = 70;
    carrier_freq = 3*10^8;
    wave_length = (3*10^8)/carrier_freq;
    
    % BS, RIS configuration
    [loc_RIS,~] = configRIS;
    [loc_BS,M,co_bandwidth,co_time] = configBS;
    [loc_TUE1,loc_TUE2,loc_RUE1,loc_RUE2] = configUsers(loc_RIS,d_0);
    tao_c = co_bandwidth * co_time;         % tao_c
    tao = K;                                % tao (Must >= K)
    dH = wave_length * wavelen_coef;
    element_size = dH^2;
    N = RIS_element;
    
    pilot_P_dBm = 20;
    pilot_P_linear = 10^(pilot_P_dBm/10) * power(10,-3);
    % Noise
    noise_dBm = -174+10*log10(co_bandwidth);        % noise power in dBm  
    noise_var = 10^((noise_dBm)/10)*(10^-3);        % noise power in Watt
    
    rho_dBm = noise_dBm+SNR;
    rho_linear = 10^((rho_dBm)/10)*(10^-3);
    
    % Aquire correlation matrix at the BS and RIS
    [R_BS] = corr_matrix_BS(M);
    [R_RIS] = corr_matrix_RIS(dH,sqrt(N));
    %% Channel pathloss
    
    % Path loss exponent
    PLE = -2;
    
    % BS to RIS pathloss
    [pathloss_BS2RIS] = pathloss_cal(loc_BS,loc_RIS,element_size,PLE);
    
    % BS to UE pathloss
    pathloss_BS2UE = zeros(1,K);
    [pathloss_BS2UE(1)] = pathloss_cal(loc_BS,loc_TUE1,element_size,PLE);
    [pathloss_BS2UE(2)] = pathloss_cal(loc_BS,loc_TUE2,element_size,PLE);
    [pathloss_BS2UE(3)] = pathloss_cal(loc_BS,loc_RUE1,element_size,PLE);
    [pathloss_BS2UE(4)] = pathloss_cal(loc_BS,loc_RUE2,element_size,PLE);
    
    % penetration loss (NLoS)
    penetration_loss_dB = 15;
    penetration_loss_linear = power(10,-(penetration_loss_dB/10));
    pathloss_BS2UE = pathloss_BS2UE.*penetration_loss_linear;
    
    % RIS to UE pathloss
    pathloss_RIS2UE = zeros(1,K);
    [pathloss_RIS2UE(1)] = pathloss_cal(loc_TUE1,loc_RIS,element_size,PLE);
    [pathloss_RIS2UE(2)] = pathloss_cal(loc_TUE2,loc_RIS,element_size,PLE);
    [pathloss_RIS2UE(3)] = pathloss_cal(loc_RUE1,loc_RIS,element_size,PLE);
    [pathloss_RIS2UE(4)] = pathloss_cal(loc_RUE2,loc_RIS,element_size,PLE);
    
    % Aggreate non-direct channel pathloss into pathloss_hat
    pathloss_hat = zeros(1,4);
    for i = 1:K
        pathloss_hat(i) = pathloss_BS2RIS*pathloss_RIS2UE(i);
    end
    
    %% Initial random phase matrix 
    min_phase = 0;
    max_phase = 2*pi;
    theta_t = exp(1i*(min_phase+rand(N,1)*(max_phase-min_phase)));
    theta_r = exp(1i*(min_phase+rand(N,1)*(max_phase-min_phase)));
    
    %% Initial random amplitude matrix
    beta_t = ones(N,1).*sqrt(0.5);
    beta_r = ones(N,1).*sqrt(0.5);
    
    bPhi_t = diag(beta_t.*theta_t);
    bPhi_r = diag(beta_r.*theta_r);
    
    %% Covariance matrix
    Rk = zeros(M,M,K);
    Psi = zeros(M,M,K);         % Covariance matrix for estimated channel
    Qk = zeros(M,M,K);
    
    
    for i = 1:K
        if(i<=Kt)
            Rk(:,:,i) =pathloss_BS2UE(i)*R_BS + pathloss_hat(i)*...
            trace(R_RIS*bPhi_t*R_RIS*bPhi_t')*R_BS;
        else
            Rk(:,:,i) = pathloss_BS2UE(i)*R_BS + pathloss_hat(i)*...
            trace(R_RIS*bPhi_r*R_RIS*bPhi_r')*R_BS;
        end   
    end
    
    for k = 1:K
        Qk(:,:,k) = pinv(Rk(:,:,k)+(noise_var/(tao*pilot_P_linear))*eye(M));
        Psi(:,:,k) = Rk(:,:,k)*Qk(:,:,k)*Rk(:,:,k);
    end
   
   
    
    
    %% PGAM initial parameters
    init_steps = 1000;            % Initial step
    kappa = 0.95;                % Adjust step
    epsilon = 10^-5;            % Stop criterion
    max_iter = 200;
    nabla_sum_SE_coef = (tao_c-tao)/(tao_c*log10(2));
    sum_SE_coef = (tao_c-tao)/tao_c;
    iter_n = 1;                 % Iteration times
    step_len = init_steps;
    sum_SE = zeros(1,max_iter);
    
    prev = struct('theta_t',theta_t,'theta_r',theta_r,'theta',[theta_t',theta_r']',...
                  'beta_t',beta_t,'beta_r',beta_r,'beta',[beta_t',beta_r']','Psi',Psi,...
                  'bPhi_t',bPhi_t,'bPhi_r',bPhi_r,'Rk',Rk,'Qk',Qk);


    while true  % Outer loop
        %% Initialize
        miu = zeros(1,K);
        nabla_theta_t_sk = zeros(N,1,K);
        nabla_theta_r_sk = zeros(N,1,K);
        nabla_theta_t_Ik = zeros(N,1,K);
        nabla_theta_r_Ik = zeros(N,1,K);
        nabla_theta_t_f = zeros(N,1);
        nabla_theta_r_f = zeros(N,1);
    
        nabla_beta_t_sk = zeros(N,1,K);
        nabla_beta_r_sk = zeros(N,1,K);
        nabla_beta_t_Ik = zeros(N,1,K);
        nabla_beta_r_Ik = zeros(N,1,K);
        nabla_beta_t_f = zeros(N,1);
        nabla_beta_r_f = zeros(N,1);
    
        sum_Psi = 0;
        R_bar = zeros(M,M,K);
        Psi_hat = zeros(M,M,K);
        A_tilda_t = zeros(N,N,K);
        A_tilda_r = zeros(N,N,K);
    
        new_Rk = zeros(M,M,K);
        new_Qk = zeros(M,M,K);
        new_Psik = zeros(M,M,K);
    
        %% ============= gradient_theta ============= %%
        %% Calculate for nabla_theta_Sk
        A_t = R_RIS*prev.bPhi_t*R_RIS;
        A_r = R_RIS*prev.bPhi_r*R_RIS;
        for k = 1:K
            miu(k) = 2*pathloss_hat(k)*trace(prev.Psi(:,:,k))*trace((prev.Qk(:,:,k)*prev.Rk(:,:,k)+...
                prev.Rk(:,:,k)*prev.Qk(:,:,k)-prev.Qk(:,:,k)*(prev.Rk(:,:,k))^2*prev.Qk(:,:,k))*R_BS);
            if(k<=Kt)    
                nabla_theta_t_sk(:,:,k)= miu(k)*diag(A_t*diag(prev.beta_t));
            else
                nabla_theta_r_sk(:,:,k)= miu(k)*diag(A_r*diag(prev.beta_r));
            end
        end
    
        %% Calculate for nabla_theta_Ik
        for k = 1:K
            sum_Psi = sum_Psi + prev.Psi(:,:,k);
            R_bar(:,:,k) = prev.Rk(:,:,k)+(((K*noise_var)/rho_linear)*eye(M));
        end
    
        for k = 1:K
            Psi_hat(:,:,k) = sum_Psi-2*(prev.Qk(:,:,k)*prev.Rk(:,:,k)*prev.Psi(:,:,k)+prev.Psi(:,:,k)*prev.Rk(:,:,k)*prev.Qk(:,:,k)...
                             -prev.Qk(:,:,k)*prev.Rk(:,:,k)*prev.Psi(:,:,k)*prev.Rk(:,:,k)*prev.Qk(:,:,k));
        end
    
        %% Calculate for A_tilda
        % Build A_tilda_t
        for k = 1:K
            if(k<=Kt)
                nu_k_bar = pathloss_hat(k)*trace(Psi_hat(:,:,k)*R_BS);
                nu_ki = 0;
                for i = 1:Kt
                    R_ki = prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k) - prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i) + ...
                        R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i);
                    nu_ki = nu_ki + pathloss_hat(k)*trace(R_ki * R_BS);
                end
                A_tilda_t(:,:,k) = (nu_k_bar + nu_ki)*A_t;
            else
                nu_ki = 0;
                for i = 1:Kt
                    R_ki = prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k) - prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i) + ...
                        R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i);
                    nu_ki = nu_ki + pathloss_hat(k)*trace(R_ki * R_BS);
                end
                 A_tilda_t(:,:,k) = nu_ki * A_t;
            end
        end
    
        % Build A_tilda_r
        for k = 1:K
            if(k<=Kt)
                nu_ki = 0;
                for i = Kt+1:K
                    R_ki = prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k) - prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i) + ...
                        R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i);
                    nu_ki = nu_ki + pathloss_hat(k)*trace(R_ki * R_BS);
                end
                 A_tilda_r(:,:,k) = nu_ki * A_r;
            else
                nu_k_bar = pathloss_hat(k)*trace(Psi_hat(:,:,k)*R_BS);
                nu_ki = 0;
                for i = Kt+1:K
                    R_ki = prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k) - prev.Qk(:,:,i)*prev.Rk(:,:,i)*R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i) + ...
                        R_bar(:,:,k)*prev.Rk(:,:,i)*prev.Qk(:,:,i);
                    nu_ki = nu_ki + pathloss_hat(k)*trace(R_ki * R_BS);
                end
                A_tilda_r(:,:,k) = (nu_k_bar + nu_ki)*A_r;
            end
        end
    
        for k = 1:K 
            nabla_theta_t_Ik(:,:,k) = diag(A_tilda_t(:,:,k)*diag(prev.beta_t));
            nabla_theta_r_Ik(:,:,k) = diag(A_tilda_r(:,:,k)*diag(prev.beta_r));
        end
    
    
        %% 計算此次(第n次)之sum_SE
        [signal,interference,gamma_k,sum_SE(iter_n)] = f_function(prev.Rk,prev.Psi,K,noise_var,rho_linear);
        fprintf('Current SE: %d\n',sum_SE(iter_n))
    
    
        for k = 1:K
            nabla_theta_t_f = nabla_theta_t_f+nabla_sum_SE_coef*(interference(k)*nabla_theta_t_sk(:,:,k)-signal(k)*nabla_theta_t_Ik(:,:,k))/((1+gamma_k(k))*(interference(k)^2));
            nabla_theta_r_f = nabla_theta_r_f+nabla_sum_SE_coef*(interference(k)*nabla_theta_r_sk(:,:,k)-signal(k)*nabla_theta_r_Ik(:,:,k))/((1+gamma_k(k))*(interference(k)^2));
        end
        nabla_theta_f = [nabla_theta_t_f', nabla_theta_r_f']';
        %% ============gradient_beta============== %%
    
        for k = 1:K
            % Calculate for nabla_beta_Sk
            if(k<=Kt)
                nabla_beta_t_sk(:,:,k)= 2*miu(k)*real(diag(A_t'*diag(prev.theta_t)));
            else
                nabla_beta_r_sk(:,:,k)= 2*miu(k)*real(diag(A_r'*diag(prev.theta_r)));
            end
        end
    
        for k = 1:K 
            nabla_beta_t_Ik(:,:,k) = 2*real(diag(A_tilda_t(:,:,k)'*diag(prev.theta_t)));
            nabla_beta_r_Ik(:,:,k) = 2*real(diag(A_tilda_r(:,:,k)'*diag(prev.theta_r)));
        end
    
        for k = 1:K
            nabla_beta_t_f = nabla_beta_t_f+nabla_sum_SE_coef*(interference(k)*nabla_beta_t_sk(:,:,k)-signal(k)*nabla_beta_t_Ik(:,:,k))/((1+gamma_k(k))*(interference(k)^2));
            nabla_beta_r_f = nabla_beta_r_f+nabla_sum_SE_coef*(interference(k)*nabla_beta_r_sk(:,:,k)-signal(k)*nabla_beta_r_Ik(:,:,k))/((1+gamma_k(k))*(interference(k)^2));
        end
        nabla_beta_f = [nabla_beta_t_f', nabla_beta_r_f']';
    
        %% PGAM
        % inner loop start
        while true  
            new_theta = prev.theta + step_len * nabla_theta_f;
            [proj_theta] = projection_theta(new_theta);
            new_beta = prev.beta + step_len * nabla_beta_f;
            [proj_beta] = projection_beta(new_beta);
    
            new_Phi_t = diag(proj_beta(1:N) .* proj_theta(1:N));
            new_Phi_r = diag(proj_beta(N+1:2*N) .* proj_theta(N+1:2*N));
            for k = 1:K
                if(k<=Kt)
                    new_Rk(:,:,k) = pathloss_BS2UE(k)*R_BS + pathloss_hat(k)*...
                                trace(R_RIS*new_Phi_t*R_RIS*new_Phi_t')*R_BS;
                else
                    new_Rk(:,:,k) = pathloss_BS2UE(k)*R_BS + pathloss_hat(k)*...
                                trace(R_RIS*new_Phi_r*R_RIS*new_Phi_r')*R_BS;
                end
                new_Qk(:,:,k) = pinv(new_Rk(:,:,k)+(noise_var/(tao*pilot_P_linear))*eye(M));
                new_Psik(:,:,k) = new_Rk(:,:,k) * new_Qk(:,:,k) * new_Rk(:,:,k) ;
            end
            
            Q_func = sum_SE(iter_n) + 2*real(nabla_theta_f'*(proj_theta-prev.theta))-(1/step_len)*norm(proj_theta-prev.theta)^2+...
                (nabla_beta_f'*(proj_beta-prev.beta))-(1/step_len)*norm(proj_beta-prev.beta)^2;
            [~,~,~,new_f_func] = f_function(new_Rk,new_Psik,K,noise_var,rho_linear);
    
            if(new_f_func<=Q_func)
                step_len = step_len*kappa;
                %fprintf('change step into: %d\n',step_len);
            else
                prev.theta_t = proj_theta(1:N);
                prev.theta_r = proj_theta(N+1:2*N);
                prev.theta = [prev.theta_t',prev.theta_r']';
                prev.beta_t = real(proj_beta(1:N));
                prev.beta_r = real(proj_beta(N+1:2*N));
                prev.beta = [prev.beta_t',prev.beta_r']';
                prev.Psi = new_Psik;
                prev.bPhi_t = new_Phi_t;
                prev.bPhi_r = new_Phi_r;
                prev.Rk = new_Rk;
                prev.Qk = new_Qk;
                break;
            end
        end
        % inner loop end
        if(abs(sum_SE(iter_n)-new_f_func)<=epsilon)
            sum_SE(iter_n+1) = new_f_func;
            break;
        elseif(iter_n>=max_iter)
            break;
        else
            step_len = init_steps;
            iter_n = iter_n + 1;
        end
    end
    final_SE = sum_SE_coef * max(sum_SE);
end


