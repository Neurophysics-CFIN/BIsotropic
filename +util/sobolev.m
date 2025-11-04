function Vdeg = sobolev(kappa,s,Lmax)
ells = 0:2:Lmax; % Generate a vector of even spherical harmonic degrees from 0 to Lmax
Vdeg     = (1 + (ells.*(ells+1))/kappa^2).^(-s);
end