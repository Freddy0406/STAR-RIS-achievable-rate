function [projected_beta] = projection_beta(beta)
    len = length(beta);
    projected_beta = zeros(len,1);
    N = len/2;
    for i = 1:len
        if(i<=N)
            projected_beta(i) = beta(i)/sqrt(beta(i)^2+beta(i+N)^2);
        else
            projected_beta(i) = beta(i)/sqrt(beta(i-N)^2+beta(i)^2);
        end
    end
end