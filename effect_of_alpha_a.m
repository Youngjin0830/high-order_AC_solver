clc; clear; close all;

% Fix the random seed for reproducibility
rand('seed',30);

% Define grid and spacing in x and y directions
Nx = 128; Lx = 0; Rx = 2; h = (Rx-Lx)/Nx;
Ny = 64; Ly = 0; Ry = 1;

% Generate cell-centered uniform grids
x = linspace(Lx+0.5*h,Rx-0.5*h,Nx);
y = linspace(Ly+0.5*h,Ry-0.5*h,Ny);

% Set polynomial order alpha 
alp = 1; % alp = 3 or alp = 5

% Interface thickness parameter
ep0 = 3*h/(2*sqrt(2)*atanh(0.9));
ep = ep0/sqrt(alp);

% Compute maximum allowable time step from monotonicity condition of $g(\phi) with $\alpha=5$$
alp5 = 5; eps5 = ep0/sqrt(alp5);
psi = ((alp5-1)./(4*alp5-1)).^(1./(2*alp5));
max_dt = alp5*eps5^2./(psi.^(2*alp5-2)*(2*alp5-1-(4*alp5-1)*psi.^(2*alp5)));

% Set time step size and number of iterations
dt = 0.99*max_dt;
Nt = 350;

% Precompute $g(\phi)$ on uniformly spaced values for interpolation
psi = linspace(-1,1,301);
g = dt/(alp*ep^2)*(psi.^(4*alp-1)-psi.^(2*alp-1))+psi;

% Initial condition
phi = 1*(2*rand(Nx,Ny)-1);
for i = 1:Nx
    for j = 1:Ny
        if x(i) >= 0.2 && x(i) <= 0.4 && y(j) >= 0.2 && y(j) <= 0.4
            phi(i,j) = -1;
        elseif x(i) >= 1.6 && x(i) <= 1.8 && y(j) >= 0.6 && y(j) <= 0.8
            phi(i,j) = 1;
        end
    end
end

% Visualize the initial condition
figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
surf(x,y,phi','EdgeColor','k','FaceColor','w')
view([-60 35])
axis([Lx Rx Ly Ry -1 1])
pbaspect([2 1 0.8])
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.590096037906264 -0.628970736110547 0.701766193301204])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.518990665212574 0.661326914441162 -0.424144402544677])
text('Interpreter','latex','String','$\phi$','FontSize',24,'Position',[-0.270860285091494 1.03504304117355 1.01023483223916])
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
        phi(i,:) = interp1([-1-eps g 1+eps],[-1 psi 1],phi(i,:));
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
        surf(x,y,phi','EdgeColor','k','FaceColor','w')
        view([-60 35])
        axis([Lx Rx Ly Ry -1 1])
        pbaspect([2 1 0.8])
        text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.590096037906264 -0.628970736110547 0.701766193301204])
        text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.518990665212574 0.661326914441162 -0.424144402544677])
        text('Interpreter','latex','String','$\phi$','FontSize',24,'Position',[-0.270860285091494 1.03504304117355 1.01023483223916])
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