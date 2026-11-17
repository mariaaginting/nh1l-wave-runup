function [pd,t_fact,t_solve] = Pnonhydros1Nwd_new(hh,U1,V1,W1,Wb1,Wb0,options)
% PNONHYDROS1NWD_NEW  Non-hydrostatic pressure (1 layer, 2D, wet-dry).
%
%   [pd,t_fact,t_solve] = Pnonhydros1Nwd_new(hh,U1,V1,W1,Wb1,Wb0,options)
%
%   Assembles and solves the pentadiagonal pressure system A*pd = b on the
%   Nx-by-Ny cell grid (cells numbered x-fastest, see inv_label.m).
%
%   options.wet_wet_coupling (default false):
%     false - original unsymmetric stencil (cases 1-4)
%     true  - SPD wet-wet graph Laplacian + Helmholtz. Dry cells are p=0,
%             shoreline/domain edges are Neumann. Use this on a moving
%             shoreline; the original stencil is a singular Neumann problem
%             there (RCOND ~ 1e-17) and the pressure blows up.
%
%   The linear solve is sparse mldivide (UMFPACK / CHOLMOD), not
%   decomposition/SparseLDL, whose rcond check crashes on wet-dry matrices.

dx    = options.dx;
dy    = options.dy;
dt    = options.dt;
thres = options.thres;
wet_wet = isfield(options,'wet_wet_coupling') && options.wet_wet_coupling;

[Nx,Ny] = size(Wb1);
Nxy = Nx*Ny;
rx = dy*dt/dx;
ry = dx*dt/dy;

if wet_wet
    [A,b] = assemble_spd_wet(hh,U1,V1,W1,Wb1,Wb0,Nx,Ny,Nxy,dx,dy,dt,rx,ry,thres);
else
    [A,b] = assemble_original(hh,U1,V1,W1,Wb1,Wb0,Nx,Ny,Nxy,dx,dy,dt,rx,ry,thres);
end

t0 = tic;
pd = A \ b;
t_fact = toc(t0);
t_solve = 0;

pd(~isfinite(pd)) = 0;
if isfield(options,'g'), gcap = options.g; else, gcap = 9.81; end
p_cap = 50 * gcap * max(max(hh(:)), thres);
pd = max(min(pd, p_cap), -p_cap);

end

%% ------------------------------------------------------------------------
function [A,b] = assemble_spd_wet(hh,U1,V1,W1,Wb1,Wb0,Nx,Ny,Nxy,dx,dy,dt,rx,ry,thres)
% Symmetric positive definite Poisson on wet cells:
%   a_ij = -c_ij,  a_ii = sum_j c_ij + 2*dt*dx*dy/h
%   c_ij = 0.5*(h_i+h_j)*rx/2  (x-face), similarly for y.
% Missing neighbours (dry cell or domain edge) are Neumann. Dry rows are p=0.

wet = hh > thres;
b  = zeros(Nxy,1);
D  = zeros(1,Nxy);
aa = zeros(1,Nxy-1);
cc = zeros(1,Nxy-1);
ab = zeros(1,Nxy-Nx);
cb = zeros(1,Nxy-Nx);

for m = 1:Nxy
    [j,k] = inv_label(m,Nx);
    if wet(j,k)
        b(m) = dx*dy*(2*Wb1(j,k)-Wb0(j,k)-W1(j,k)) ...
            - hh(j,k)*dy*(U1(j+1,k)-U1(j,k)) - hh(j,k)*dx*(V1(j,k+1)-V1(j,k));
        if ~isfinite(b(m)), b(m) = 0; end
        D(m) = 2*dt*dx*dy/max(hh(j,k),thres);
    else
        D(m) = 1;
    end
end

for k = 1:Ny
    for j = 1:Nx-1
        m = j + (k-1)*Nx;
        if wet(j,k) && wet(j+1,k)
            cij = 0.5*(hh(j,k)+hh(j+1,k))*rx/2;
            cc(m) = -cij;
            aa(m) = -cij;
            D(m)   = D(m)   + cij;
            D(m+1) = D(m+1) + cij;
        end
    end
end

for k = 1:Ny-1
    for j = 1:Nx
        m = j + (k-1)*Nx;
        if wet(j,k) && wet(j,k+1)
            cij = 0.5*(hh(j,k)+hh(j,k+1))*ry/2;
            cb(m) = -cij;
            ab(m) = -cij;
            D(m)    = D(m)    + cij;
            D(m+Nx) = D(m+Nx) + cij;
        end
    end
end

D(D <= 0 | ~isfinite(D)) = 1;

A = sparse(1:Nxy,1:Nxy,D,Nxy,Nxy) ...
  + sparse(1:Nxy-1,2:Nxy,cc,Nxy,Nxy) ...
  + sparse(2:Nxy,1:Nxy-1,aa,Nxy,Nxy) ...
  + sparse(1:Nxy-Nx,1+Nx:Nxy,cb,Nxy,Nxy) ...
  + sparse(1+Nx:Nxy,1:Nxy-Nx,ab,Nxy,Nxy);
A = (A + A.')/2;
end

%% ------------------------------------------------------------------------
function [A,b] = assemble_original(hh,U1,V1,W1,Wb1,Wb0,Nx,Ny,Nxy,dx,dy,dt,rx,ry,thres)
% Original unsymmetric pentadiagonal stencil (cases 1-4).

b  = zeros(Nxy,1);
D  = zeros(1,Nxy);
aa = zeros(1,Nxy-1);
cc = zeros(1,Nxy-1);
ab = zeros(1,Nxy-Nx);
cb = zeros(1,Nxy-Nx);

for m = 1:Nxy
    [j,k] = inv_label(m,Nx);
    if hh(j,k) > thres
        b(m) = dx*dy*(2*Wb1(j,k)-Wb0(j,k)-W1(j,k)) ...
            - hh(j,k)*dy*(U1(j+1,k)-U1(j,k)) - hh(j,k)*dx*(V1(j,k+1)-V1(j,k));
        D(m) = hh(j,k)*(rx+ry) + 2*dt*dx*dy/max(hh(j,k),thres);
        if m < Nxy
            cc(m) = -hh(j,k)*rx/2;
            aa(m) = -hh(j,k)*rx/2;
        end
    else
        b(m) = 0;
        D(m) = 1;
    end
end

for m = 1:Nxy-Nx
    [j,k] = inv_label(m,Nx);
    if hh(j,k) > thres
        ab(m) = -hh(j,k)*ry/2;
        cb(m) = -hh(j,k)*ry/2;
    end
end

for m = 1:Nx
    [j,k] = inv_label(m,Nx);
    if hh(j,k) > thres
        cb(m) = -hh(j,k)*ry;
    else
        cb(m) = 0;
    end
end

for m = (Ny-2)*Nx+1:Nxy-Nx
    [j,k] = inv_label(m,Nx);
    if hh(j,k) > thres
        ab(m) = -hh(j,k)*ry;
    else
        ab(m) = 0;
    end
end

for s = 1:Ny-1
    cc(s*Nx) = 0;
    aa(s*Nx) = 0;

    [j,k] = inv_label(s*Nx-1,Nx);
    if hh(j,k) > thres
        aa(s*Nx-1) = -hh(j,k)*rx;
    else
        aa(s*Nx-1) = 0;
    end

    [j,k] = inv_label(s*Nx+1,Nx);
    if hh(j,k) > thres
        cc(s*Nx+1) = -hh(j,k)*rx;
    else
        cc(s*Nx+1) = 0;
    end
end

[j,k] = inv_label(1,Nx);
if hh(j,k) > thres
    cc(1) = -hh(j,k)*rx;
else
    cc(1) = 0;
end
[j,k] = inv_label(Nxy-1,Nx);
if hh(j,k) > thres
    aa(Nxy-1) = -hh(j,k)*rx;
else
    aa(Nxy-1) = 0;
end

A = sparse(1:Nxy,1:Nxy,D,Nxy,Nxy) ...
  + sparse(1:Nxy-1,2:Nxy,cc,Nxy,Nxy) ...
  + sparse(2:Nxy,1:Nxy-1,aa,Nxy,Nxy) ...
  + sparse(1:Nxy-Nx,1+Nx:Nxy,cb,Nxy,Nxy) ...
  + sparse(1+Nx:Nxy,1:Nxy-Nx,ab,Nxy,Nxy);
end
