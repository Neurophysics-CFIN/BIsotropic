function Y = sph_harm(l, theta, phi)
% Compute the associated Legendre polynomial with Condon-Shotley phase
P_lm = legendre(l, cos(theta)); % m = 0...l rows , theta columns
if l~=0 %add m< 0
    m = (l:-1:1)';
    %c_m = (-1).^m .*(factorial(l - m)./factorial(l + m));
    % P_lm = cat(1,c_m'.*P_lm(l+1:-1:2,:),P_lm);
    P_lm = cat(1,(-1).^m.*P_lm(l+1:-1:2,:),P_lm);
end
 % Compute the normalization constant
m = (-l:l)';
K_lm = sqrt((2*l+1)/(4*pi) * factorial(l-abs(m))./factorial(l+abs(m)));
% Compute the spherical harmonic
Y = K_lm .* P_lm .* exp(1i * m .* phi');
end