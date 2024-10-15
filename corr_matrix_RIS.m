function [R]=corr_matrix_RIS(lambda,N_HV)

    %Extract the element size    
    d = lambda/4;
    
    %Generate a grid for the elements
    gridPoints = (0:N_HV-1)*d;
    
    [X,Y] = meshgrid(gridPoints,gridPoints);
    
    locations = X(:)+1i*Y(:);
    
    
    %Total number of elements
    N = length(locations);
    
    
    %Compute the spatial correlation matrix
    R = zeros(N,N);
    
    for m = 1:N
        for l = 1:N
            
            R(m,l) = sinc(2*abs(locations(m)-locations(l))/lambda);
            
        end
    end
end