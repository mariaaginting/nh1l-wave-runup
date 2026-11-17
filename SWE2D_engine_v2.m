function [V1,U1] = SWE2D_engine_v2(E1,U0,V0,hh,qu,qv,options)
% SWE2D_ENGINE_V2  Hydrostatic momentum step of the 2D shallow-water solver.
%   [V1,U1] = SWE2D_engine_v2(E1,U0,V0,hh,qu,qv,options)
%
%   Updates the face velocities U (Nx+1 x Ny) and V (Nx x Ny+1) from the new
%   surface elevation E1, with upwind advection and a wet/dry mask. Written
%   with matrix operations (no loops); improved from the code of Bu Cisca.
%   options: struct with g, dx, dy, dt, thres.
%
%   Written by: N. Subasita (29-July-2019)

[Nx,Ny] = size(E1);
g       = options.g;
dx      = options.dx;
dy      = options.dy;
dt      = options.dt;
thres   = options.thres;

U1 = zeros(Nx+1,Ny);
V1 = zeros(Nx,Ny+1);

%% Momentum equation -> update U
WL   = (hh(2:Nx,1:Ny)>thres | hh(1:Nx-1,1:Ny)>thres);
advx = 0;
ufd  = WL.*(qu(1:Nx-1,1:Ny)+qu(2:Nx,1:Ny))./max(thres,hh(2:Nx,1:Ny)+hh(1:Nx-1,1:Ny));     % uu_x downstream
ufu  = WL.*(qu(3:Nx+1,1:Ny)+qu(2:Nx,1:Ny))./max(thres,hh(2:Nx,1:Ny)+hh(1:Nx-1,1:Ny));     % uu_x upstream
advx = advx + WL.*((ufu<0).*ufu.*(U0(3:Nx+1,1:Ny)-U0(2:Nx,1:Ny))+(ufd>0).*ufd.*(U0(2:Nx,1:Ny)-U0(1:Nx-1,1:Ny))); % momentum conservative
advy = zeros(size(U0));
% 1 < j < Ny
ufd  = WL(:,2:Ny-1).*(qv(2:Nx,2:Ny-1)+qv(2:Nx,1:Ny-2))./max(thres,hh(2:Nx,2:Ny-1)+hh(1:Nx-1,2:Ny-1)); % vu_y downstream
ufu  = WL(:,2:Ny-1).*(qv(2:Nx,2:Ny-1)+qv(2:Nx,3:Ny))./max(thres,hh(2:Nx,2:Ny-1)+hh(1:Nx-1,2:Ny-1));   % vu_y upstream
advy(2:Nx,2:Ny-1) = advy(2:Nx,2:Ny-1)+...
    WL(:,2:Ny-1).*((ufu<0).*ufu.*(U0(2:Nx,3:Ny)  -U0(2:Nx,2:Ny-1)))+...
    WL(:,2:Ny-1).*((ufd>0).*ufd.*(U0(2:Nx,2:Ny-1)-U0(2:Nx,1:Ny-2)));
% j = 1
ufd  = WL(:,1).*qv(2:Nx,1)./max(thres,hh(2:Nx,1));
ufu  = WL(:,1).*(qv(2:Nx,1)+qv(2:Nx,2))./max(thres,hh(2:Nx,1)+hh(1:Nx-1,1));
advy(2:Nx,1) = advy(2:Nx,1)+WL(:,1).*((ufu<0).*ufu.*(U0(2:Nx,2)-U0(2:Nx,1))+(ufd>0).*ufd.*(U0(2:Nx,1)));
% j = Ny
ufd  = WL(:,Ny).*(qv(2:Nx,Ny)+qv(2:Nx,Ny-1))./max(thres,hh(2:Nx,Ny)+hh(1:Nx-1,Ny));
ufu  = WL(:,Ny).*qv(2:Nx,Ny)./max(thres,hh(2:Nx,Ny));
advy(2:Nx,Ny) = advy(2:Nx,Ny)+WL(:,Ny).*((ufu<0).*ufu.*(-U0(2:Nx,Ny))+(ufd>0).*ufd.*(U0(2:Nx,Ny)-U0(2:Nx,Ny-1)));
% U at the next time step
U1(2:Nx,1:Ny) = WL.*(U0(2:Nx,1:Ny)-dt*(g/dx*(E1(2:Nx,1:Ny)-E1(1:Nx-1,1:Ny))+advx/dx+advy(2:Nx,1:Ny)/dy));

%% Momentum equation -> update V
WL   = (hh(1:Nx,2:Ny)>thres | hh(1:Nx,1:Ny-1)>thres);
advy = 0;
ufd  = WL.*(qv(1:Nx,1:Ny-1)+qv(1:Nx,2:Ny))./max(thres,hh(1:Nx,2:Ny)+hh(1:Nx,1:Ny-1)); % vv_y downstream
ufu  = WL.*(qv(1:Nx,3:Ny+1)+qv(1:Nx,2:Ny))./max(thres,hh(1:Nx,2:Ny)+hh(1:Nx,1:Ny-1)); % vv_y upstream
advy = advy+WL.*((ufu<0).*ufu.*(V0(1:Nx,3:Ny+1)-V0(1:Nx,2:Ny))+(ufd>0).*ufd.*(V0(1:Nx,2:Ny)-V0(1:Nx,1:Ny-1))); % momentum conservative
advx = zeros(size(U0));
% 1 < i < Nx
ufd  = WL(2:Nx-1,:).*(qu(1:Nx-2,2:Ny)+qu(2:Nx-1,2:Ny))./max(thres,hh(2:Nx-1,2:Ny)+hh(2:Nx-1,1:Ny-1)); % uv_x downstream
ufu  = WL(2:Nx-1,:).*(qu(3:Nx,2:Ny)+qu(2:Nx-1,2:Ny))./max(thres,hh(2:Nx-1,2:Ny)+hh(2:Nx-1,1:Ny-1));   % uv_x upstream
advx(2:Nx-1,2:Ny) = advx(2:Nx-1,2:Ny)+ ...
    WL(2:Nx-1,:).*((ufu<0).*ufu.*(V0(3:Nx,2:Ny)-V0(2:Nx-1,2:Ny)))+...
    WL(2:Nx-1,:).*((ufd>0).*ufd.*(V0(2:Nx-1,2:Ny)-V0(1:Nx-2,2:Ny)));
% i = Nx
ufu  = WL(Nx,:).*qu(Nx,2:Ny)./max(thres,hh(Nx,2:Ny));
ufd  = WL(Nx,:).*(qu(Nx-1,2:Ny)+qu(Nx,2:Ny))./max(thres,hh(Nx,2:Ny)+hh(Nx,1:Ny-1));
advx(Nx,2:Ny) = advx(Nx,2:Ny)+WL(Nx,:).*((ufu<0).*ufu.*(-V0(Nx,2:Ny))+(ufd>0).*ufd.*(V0(Nx,2:Ny)-V0(Nx-1,2:Ny)));
% i = 1
ufu = WL(1,:).*(qu(2,2:Ny)+qu(1,2:Ny))./max(thres,hh(1,2:Ny)+hh(1,1:Ny-1));
ufd = WL(1,:).*qu(1,2:Ny)./max(thres,hh(1,2:Ny));
advx(1,2:Ny)  = advx(1,2:Ny)+WL(1,:).*((ufu<0).*ufu.*(V0(2,2:Ny)-V0(1,2:Ny))+(ufd>0).*ufd.*(V0(1,2:Ny)));
% V at the next time step
V1(1:Nx,2:Ny) = WL.*(V0(1:Nx,2:Ny)-dt*(g/dy*(E1(1:Nx,2:Ny)-E1(1:Nx,1:Ny-1))+advy/dy+advx(1:Nx,2:Ny)/dx));

end
