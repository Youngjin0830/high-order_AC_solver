clc; clear; close all;

% Define grid and spacing in x and y directions
Nx = 128; Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ny = 256; Ly = 0; Ry = 2;

% Generate cell-centered uniform grids
x = linspace(Lx+0.5*h,Rx-0.5*h,Nx);
y = linspace(Ly+0.5*h,Ry-0.5*h,Ny);

% Set polynomial order alpha 
alp = 1; % alp = 5

% Interface thickness parameter
if alp == 1
ep = 0.019;
else
ep = 0.0051;
end

% Compute maximum allowable time step from monotonicity condition of $g(\phi) with $\alpha=5$$
alp5 = 5; eps5 = 0.0051;
psi = ((alp5-1)./(4*alp5-1)).^(1./(2*alp5));
max_dt = alp5*eps5^2./(psi.^(2*alp5-2)*(2*alp5-1-(4*alp5-1)*psi.^(2*alp5)));

% Set time step size and number of iterations
dt = 0.99*max_dt;
Nt = 550;

% Precompute $g(\phi)$ on uniformly spaced values for interpolation
psi = linspace(-1,1,301);
g = dt/(alp*ep^2)*(psi.^(4*alp-1)-psi.^(2*alp-1))+psi;

% Initial condition
phi = -ones(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if x(i) >= 0.1 && x(i) <= 0.45 && y(j) >= 0.1 && y(j) <= 0.65
            phi(i,j) = 1;
        elseif x(i) >= 0.55 && x(i) <= 0.9 && y(j) >= 0.1 && y(j) <= 1.25
            phi(i,j) = 1;
        elseif x(i) >= 0.55 && x(i) <= 0.9 && y(j) >= 1.35 && y(j) <= 1.9
            phi(i,j) = 1;
        elseif x(i) >= 0.1 && x(i) <= 0.45 && y(j) >= 0.75 && y(j) <= 1.9
            phi(i,j) = 1;
        end
    end
end

% Visualize the initial condition
figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
surf(x,y,phi','EdgeColor','none')
shading interp
colormap jet
clim([-1,1])
view(2)
axis([Lx Rx Ly Ry -1 1])
pbaspect([1 2 0.8])
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.739299610894942 -0.124513618677043])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.16147859922179 1.79182879377432])
drawnow

% Set up tridiagonal matrix coefficients for Thomas algorithm
for i=1:Nx
    ax(i) = -1/h^2;
    cx(i) = -1/h^2;
end
for j=1:Ny
    ay(j) = -1/h^2;
    cy(j) = -1/h^2;
end
bx(1)=1/dt+1.0/h^2;
for i=2:Nx-1
    bx(i) = 1/dt+2.0/h^2;
end
bx(Nx)=1/dt+1.0/h^2;
by(1) = 1/dt+1.0/h^2;
for j=2:Ny-1
    by(j) = 1/dt+2.0/h^2;
end
by(Ny) = 1/dt+1.0/h^2;

% Compute numerical solution using the operator splitting method combined with an interpolation approach
nphi = phi;
for it = 1:Nt
    it
    % Step 1: Solve Eq. (17) using the interpolation technique
    for i = 1:Nx
        phi(i,:) = interp1([-1-10*eps g 1+10*eps],[-1 psi 1],phi(i,:));
    end
    % Step 2: Solve Eq. (18) using the Thomas algorithm
    for j=1:Ny
        dx = phi(:,j)/dt;
        nphi(1:Nx,j) = thomas(ax, bx, cx, dx);
    end
    % Step 3: Solve Eq. (19) using the Thomas algorithm
    for i=1:Nx
        dy = nphi(i,:)/dt;
        phi(i,1:Ny) = thomas(ay, by, cy, dy);
    end

    % Plot numerical solution
    if mod(it,50) == 0
        figure(1); clf; hold on; box on;
        set(gca,'fontsize',21);
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gcf,'PaperPositionMode','auto')
        surf(x,y,phi','EdgeColor','none')
        shading interp
        colormap jet
        clim([-1,1])
        view(2)
        axis([Lx Rx Ly Ry -1 1])
        pbaspect([1 2 0.8])
        text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.739299610894942 -0.124513618677043])
        text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.16147859922179 1.79182879377432])
        drawnow
    end
end

%% Thomas algorithm for tridiagonal systems
function x = thomas(alpha,beta,gamma,f)
n = length(f);
for  i = 2:n
    mult = alpha(i)/beta(i-1);
    beta(i) = beta(i)-mult*gamma(i-1);
    f(i) = f(i)-mult*f(i-1);
end
x(n) = f(n)/beta(n);
for i = n-1:-1:1
    x(i) = (f(i)-gamma(i)*x(i+1))/beta(i);
end
end