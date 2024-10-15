function [projected_theta] = projection_theta(theta)
    len = length(theta);
    projected_theta = zeros(len,1);
    for i = 1:len
        projected_theta(i) = theta(i)/abs(theta(i));
    end
end