function [inv_x,inv_y] = inv_label(p,Nx)
% INV_LABEL  Convert linear index p (x-fastest ordering) to grid indices.
%   [i,j] = inv_label(p,Nx) returns the x-index i and y-index j of cell p
%   on a grid with Nx cells in the x-direction, so that p = i + (j-1)*Nx.

inv_y = ceil(p/Nx);
inv_x = p - (inv_y-1)*Nx;
end
