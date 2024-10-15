function [R]=corr_matrix_BS(M)

    P = M/2;
    r = zeros(M);
    omega = 0.3;

    for a = 1:M 
        for b = 1:P
            phi = (-1*pi/2)+(b-1)*pi/P;
            r(a,b) = (1/sqrt(P))*exp(-2i*pi*omega*(a-1)*sin(phi));
        end
    end

    R = r*r';

end