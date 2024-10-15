function [Z] = standard_complex_gaussian_gen(row,column,variance)

    z_real = sqrt(variance) * randn(row, column);
    z_imag = sqrt(variance) * randn(row, column);
    Z = (z_real + 1i * z_imag)./sqrt(2);

end