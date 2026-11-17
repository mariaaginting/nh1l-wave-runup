% CASE 2 -- Propagation of a solitary wave over a flat bottom (2DH, 1 layer)
%
% A sech^2 hump of height 2*amp at rest splits into two solitary waves of
% amplitude amp. The right-going wave is compared with the analytical
% solitary wave at t = 32 s.
%
% Requires: SWE2D_engine_v2.m, Pnonhydros1Nwd_new.m, closest.m, inv_label.m

clear; close all;

%% Parameters
Lx = 600;
Ly = 10;
h  = 10;                       % still-water depth
g  = 9.81;
dx = 1; dy = 1;
C  = 0.5;                      % Courant number
dt = C*dx/sqrt(2*g*h)
thres = 0;                     % wet/dry threshold

options.g     = g;
options.dx    = dx;
options.dy    = dy;
options.dt    = dt;
options.thres = thres;

x = -Lx:dx:Lx;  y = 0:dy:Ly;
[xx,yy] = meshgrid(x,y);
Nx = length(x); Ny = length(y); Nxy = Nx*Ny;
T  = 32+dt; Nt = floor(T/dt);
t  = 0:dt:T;

%% Bottom
d = h*ones(Nx,Ny);

figure
surf(xx,yy,d','edgecolor','none')
title('Depth')

figure
surf(xx,yy,-d','edgecolor','none')
title('Bottom Topography')

%% Preallocation
E0 = zeros(Nx,Ny);   E1 = zeros(Nx,Ny);   hh = zeros(Nx,Ny);
U0 = zeros(Nx+1,Ny); U1 = zeros(Nx+1,Ny);
V0 = zeros(Nx,Ny+1); V1 = zeros(Nx,Ny+1);
qu = zeros(Nx+1,Ny); qv = zeros(Nx,Ny+1);
W0 = zeros(Nx,Ny);   W1 = zeros(Nx,Ny);    % vertical velocity at the surface
Wb0 = zeros(Nx,Ny);  Wb1 = zeros(Nx,Ny);   % vertical velocity at the bed
pd = zeros(Nxy,1);   P = zeros(Nx,Ny);

%% Initial condition: solitary hump at rest (splits into two waves)
amp   = 2;
x0    = 80+Lx;
gamma = sqrt(3*amp/(4*h^2*(h+amp)));
for k = 1:Ny
    for j = 1:Nx
        E0(j,k) = max(2*amp*(sech(gamma*(j*dx-x0)))^2,-d(j,k)+thres);
        hh(j,k) = E0(j,k)+d(j,k);
    end
end

figure
surf(xx,yy,E0','edgecolor','none')
title('Initial eta \eta_0')
view(3,25)
xlabel('x'); ylabel('y'); zlabel('z')
figure
surf(xx,yy,E0','edgecolor','none')
title('Initial eta \eta_0')
view(0,1)
xlabel('x'); zlabel('z')

figure
plot(x,E0(:,closest(y,5)))
grid on
xlabel('x[m]'); ylabel('\eta[m]')
axis([-Lx Lx 0 2*amp])

%% Time loop: finite volume with flooding & drying
colormap('winter'); caxis([-2, 2]);
fig = figure;
set(fig,'Position',[600 200 1200 400],'color',[1 1 1])
for n = 1:Nt

    % Face depths and fluxes
    for j = 2:Nx
        for k = 1:Ny
            if U0(j,k) > 0
                hh(j,k) = d(j-1,k)+E0(j-1,k);       % upwind: take depth from the left
            elseif U0(j,k) < 0
                hh(j,k) = d(j,k)+E0(j,k);           % upwind: take depth from the right
            else
                hh(j,k) = min(d(j-1,k),d(j,k))+max(E0(j-1,k),E0(j,k));
            end
            qu(j,k) = U0(j,k)*hh(j,k);
        end
    end
    qu(1,1:Ny)    = hh(1,1:Ny).*U0(1,1:Ny);
    qu(Nx+1,1:Ny) = hh(Nx,1:Ny).*U0(Nx+1,1:Ny);
    qv(1:Nx,1)    = hh(1:Nx,1).*V0(1:Nx,1);
    qv(1:Nx,Ny+1) = hh(1:Nx,Ny).*V0(1:Nx,Ny+1);
    for k = 2:Ny
        for j = 1:Nx
            if V0(j,k) > 0
                hh(j,k) = d(j,k-1)+E0(j,k-1);
            elseif V0(j,k) < 0
                hh(j,k) = d(j,k)+E0(j,k);
            else
                hh(j,k) = min(d(j,k-1),d(j,k))+max(E0(j,k-1),E0(j,k));
            end
            qv(j,k) = V0(j,k)*hh(j,k);
        end
    end

    % Surface elevation (mass equation)
    for k = 1:Ny
        for j = 1:Nx
            E1(j,k) = E0(j,k)-dt/dx*(qu(j+1,k)-qu(j,k))-dt/dy*(qv(j,k+1)-qv(j,k));
            hh(j,k) = d(j,k)+E1(j,k);
        end
    end

    % Hydrostatic momentum step and vertical velocity at the bed
    [V1,U1] = SWE2D_engine_v2(E1,U0,V0,hh,qu,qv,options);
    for k = 2:Ny-1
        for j = 2:Nx-1
            if hh(j,k) > thres || hh(j-1,k) > thres
                Wb1(j,k) = -max(0,U1(j,k))*(d(j,k)-d(j-1,k))/dx-min(0,U1(j+1,k))*(d(j+1,k)-d(j,k))/dx ...
                    -max(0,V1(j,k))*(d(j,k)-d(j,k-1))/dy-min(0,V1(j,k+1))*(d(j,k+1)-d(j,k))/dy;
                W1(j,k) = W0(j,k)+(Wb1(j,k)-Wb0(j,k));   % W_t = -Wb_t
            else
                Wb1(j,k) = 0; W1(j,k) = 0;
            end
        end
    end

    % Non-hydrostatic pressure
    pd = Pnonhydros1Nwd_new(hh,U1,V1,W1,Wb1,Wb0,options);
    for m = 1:Nxy
        [i,j] = inv_label(m,Nx);
        if hh(i,j) > thres
            P(i,j) = pd(m);
        else
            P(i,j) = 0;
        end
    end

    % Pressure correction of U, V, W; update of Wb and W
    for k = 1:Ny
        for j = 2:Nx
            U1(j,k) = U1(j,k)-dt/dx/2*(P(j,k)-P(j-1,k));
        end
    end
    for k = 2:Ny
        for j = 1:Nx
            V1(j,k) = V1(j,k)-dt/dy/2*(P(j,k)-P(j,k-1));
        end
    end
    for k = 1:Ny
        for j = 1:Nx
            W1(j,k) = W1(j,k)+2*dt/hh(j,k)*P(j,k);
        end
    end
    for k = 2:Ny-1
        for j = 2:Nx-1
            Wb1(j,k) = -max(0,U1(j,k))*(d(j,k)-d(j-1,k))/dx-min(0,U1(j+1,k))*(d(j+1,k)-d(j,k))/dx ...
                -max(0,V1(j,k))*(d(j,k)-d(j,k-1))/dy-min(0,V1(j,k+1))*(d(j,k+1)-d(j,k))/dy;
            W1(j,k) = W1(j,k)-dt*(U1(j,k)+U1(j+1,k))/2*(W1(j,k)-W1(j-1,k))/dx ...
                -dt*(V1(j,k)+V1(j,k+1))/2*(W1(j,k)-W1(j,k-1))/dy;   % advection -u*w_x - v*w_y
        end
    end
    Wb1(1,1:Ny) = Wb1(2,1:Ny); Wb1(Nx,1:Ny) = Wb1(Nx-1,1:Ny);

    % Plotting
    if mod(n,100) == 1
        figure(5)
        colormap('winter')
        F1 = surf(xx,yy,-d'); hold on;     % bottom
        F2 = surf(xx,yy,E1'); hold off;    % surface elevation
        set(F1,'FaceColor','None','EdgeColor','k');
        set(F2,'FaceColor','interp');
        view(10,15)
        zlim([0,2*amp]);
        title(sprintf('surface profile at t= %1.4s',num2str(n*dt)))
        pause(0.01);
    end

    figure(6)
    plot(x,E1(:,closest(y,5)))
    grid on
    xlabel('x[m]'); ylabel('\eta[m]')
    title([' \eta [m] at t = ', num2str(t(n)) ,' s'])
    axis([-Lx Lx 0 2*amp])
    pause(0.01)

    % Advance in time
    E0 = E1;
    U0 = U1; V0 = V1;
    W0 = W1; Wb0 = Wb1;
end

%% Analytical solitary wave at t = 32 s
c    = sqrt(g*(h+amp));
time = 32;
ii   = 1:Nx;
eta_ex  = amp*(sech(gamma*(ii*dx-x0-c*time))).^2;
eta_ex0 = amp*(sech(gamma*(ii*dx-x0))).^2;

figure
plot(x,eta_ex0,'-.b',x,eta_ex,'r')
xlabel('x[m]'); ylabel('\eta[m]')
title(' \eta [m] at t = 32 s')
axis([0 Lx 0 amp+0.5])
grid on
