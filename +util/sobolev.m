function Vdeg = sobolev(kappa,s,Lmax)
%   SOBOLEV  Sobolev-type spectral weights as a function of spherical-harmonic degree.
%
%   Vdeg = SOBOLEV(kappa, s, Lmax)
%
%   Generates a vector of Sobolevspectral weights indexed
%   by spherical-harmonic degree ℓ = 0..Lmax. These weights are used
%   to define band weightings or smoothness penalties in spectral energies on
%    SO(3).
%
%   INPUTS
%     kappa : Positive scale parameter controlling the spectral roll-off.
%             Larger kappa shifts weight toward higher degrees.
%     s     : Smoothness/exponent parameter. Larger s produces faster decay
%             with ℓ (stronger smoothing).
%     Lmax  : Maximum degree ℓ.
%
%   OUTPUT
%     Vdeg  : (Lmax+1)-by-1 vector of weights, where
%               Vdeg(ell+1) corresponds to degree ℓ.
%
%   DEFINITION
%     For ℓ = 0..Lmax,
%       Vdeg(ell+1) = ( 1 + ℓ(ℓ+1) / kappa^2 )^(−s).

%   EXAMPLE
%     % Sobolev weights up to Lmax = 8
%     Vdeg = sobolev(3.0, 2.5, 8);

ells = 0:Lmax; % Generate a vector of even spherical harmonic degrees from 0 to Lmax
Vdeg     = (1 + (ells.*(ells+1))/kappa^2).^(-s);
end