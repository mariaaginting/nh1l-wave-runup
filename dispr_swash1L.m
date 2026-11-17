function [omega,norm_omega,k_in,Vg] = dispr_swash1L(k,depth,g,om_in)
% DISPR_SWASH1L  Dispersion relation of the 1-layer non-hydrostatic model.
%   omega = sqrt(g*h)*k/sqrt(1+(k*h/2)^2)
%
%   [omega,norm_omega,k_in,Vg] = dispr_swash1L(k,depth,g,om_in)
%     k      wavenumbers
%     depth  depth profile (vector)
%     g      gravitational acceleration
%     om_in  target frequencies
%   Returns omega(k), norm_omega = omega^2/(g*h*k^2), the wavenumbers k_in
%   matching om_in, and the group velocity Vg (length(k_in) x length(depth)).

omega       = sqrt(g*depth).*k./sqrt(1+(k.*depth./2).^2);
norm_omega  = 1./(1+(k.*depth./2).^2);
k_in        = interp1(omega,k,om_in,'linear');
Vg          = zeros(length(k_in),length(depth));

for j = 1:length(depth)
    for i = 1:length(k_in)
        if k_in(i) == 0
            Vg(i,j) = sqrt(g*depth(j));
        else
            Vg(i,j) = sqrt(g*depth(j))/ ...
                ((1+(k_in(i)*depth(j)/2)^2)).^(3/2);
        end
    end
end
Vg(isnan(Vg)) = 0;

end
