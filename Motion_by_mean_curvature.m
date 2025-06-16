clc; clear; close all;

% Define grid and spacing in x and y directions
Nx = 128; Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ny = 128; Ly = 0; Ry = 1;

% Generate cell-centered uniform grids
x = linspace(Lx+0.5*h,Rx-0.5*h,Nx);
y = linspace(Ly+0.5*h,Ry-0.5*h,Ny);

% Set polynomial order alpha 
alp = 4;

% Interface thickness parameter
ep = h;

% Compute maximum allowable time step from monotonicity condition of $g(\phi)$
psi = ((alp-1)./(4*alp-1)).^(1./(2*alp));
max_dt = alp*ep^2./(psi.^(2*alp-2)*(2*alp-1-(4*alp-1)*psi.^(2*alp)));

% Set time step size and number of iterations
dt = 0.1*max_dt;
Nt = 2000;

% Precompute $g(\phi)$ on uniformly spaced values for interpolation
psi = linspace(-1,1,101);
g = dt/(alp*ep^2)*(psi.^(4*alp-1)-psi.^(2*alp-1))+psi;

% Initialize $\phi$ with small random values in [-0.2, 0.2]
R0 = 0.3;
phi = tanh((R0-sqrt((x'-0.5).^2+(y-0.5).^2))/(sqrt(2)*ep));

% Visualize the initial condition
figure(1); clf; hold on; box on;
contour(x,y,phi,[0 0],'EdgeColor','k','LineWidth',1)
view(2)
axis image
axis([Lx Rx Ly Ry -1 1])
drawnow

% Initial radii of analytic and numerical solutions
it = 0;
analytic_R(1) = sqrt(R0^2-2*it*dt);
c = contour(x,y,phi,[0 0],'EdgeColor','none');
numerical_R(1) = mean(sqrt(sum((c(:,2:end)-[0.5; 0.5]).^2,1)));

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
    if it == 500 || it == 1000 || it == 1500 || it == 2000
        figure(1); hold on; box on;
        contour(x,y,phi,[0 0],'EdgeColor','k','LineWidth',1)
        view(2)
        axis image
        axis([Lx Rx Ly Ry -1 1])
        drawnow
    end

    analytic_R(it+1) = sqrt(R0^2-2*it*dt);
    c = contour(x,y,phi,[0 0],'EdgeColor','none');
    numerical_R(it+1) = mean(sqrt(sum((c(:,2:end)-[0.5; 0.5]).^2,1)));
end
figure(1);
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
text('Interpreter','latex','String','$x$','FontSize',23,'Position',[0.747446675565907 -0.0622809508932377])
text('Interpreter','latex','String','$y$','FontSize',21,'Position',[-0.0671547047629869 0.913862932772871])
axis image
axis([Lx Rx Ly Ry -1 1])
annotation('arrow',[0.314880952380952 0.472619047619047],[0.296825396825397 0.466666666666668]);

% Plot the evolution of analytic and numerical radii over time
ms = 10;
figure(2); clf; hold on; grid on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto','Position',[500 500 800 400])
text('Interpreter','latex','String','$t$','FontSize',21,'Position',[0.0348258353283817 0.0619631901840491])
text('Interpreter','latex','String','$R(t)$','FontSize',21,'Position',[-0.00346722242084707 0.276441717791411])
t = [0:it]*dt;
plot(t,analytic_R,'k-*','markersize',ms,'LineWidth',1,'MarkerIndices',[100:200:2000])
plot(t,numerical_R,'ro','markersize',ms,'LineWidth',1,'MarkerIndices',[1:200:2000])
axis([0 Nt*dt 0.08 0.32])
drawnow
leg = legend('Analytic','Numerical');
set(leg,'Interpreter','latex','Location','best')
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