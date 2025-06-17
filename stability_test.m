clc; clear; close all;

% Fix the random seed for reproducibility
rand('seed',96);

% Define grid and spacing in x and y directions
Nx = 128; Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ny = 128; Ly = 0; Ry = 1;

% Generate cell-centered uniform grids
x = linspace(Lx+0.5*h,Rx-0.5*h,Nx);
y = linspace(Ly+0.5*h,Ry-0.5*h,Ny);

% Set polynomial order alpha 
alp = 3;

% Interface thickness parameter
ep0 = 2*h/(2*sqrt(2)*atanh(0.9));
ep = ep0/sqrt(alp);

% Compute maximum allowable time step from monotonicity condition of $g(\phi)$
psi = ((alp-1)./(4*alp-1)).^(1./(2*alp));
max_dt = alp*ep^2./(psi.^(2*alp-2)*(2*alp-1-(4*alp-1)*psi.^(2*alp)));

% Set time step size and number of iterations
dt = max_dt;
Nt = 550;

% Precompute $g(\phi)$ on uniformly spaced values for interpolation
psi = linspace(-1,1,101);
g = dt/(alp*ep^2)*(psi.^(4*alp-1)-psi.^(2*alp-1))+psi;

% Initialize $\phi$ with small random values in [-0.2, 0.2]
phi = 0.2*(2*rand(Nx,Ny)-1);

% Add 100 random circular perturbations
nind = 100;
randind = ceil(Nx*rand(nind,2));
randradius = 0.1*rand(nind,1);
for k = 1:nind
    val = 0.25*(2*rand(1)-1);
    for i = 1:Nx
        for j= 1:Ny
            if sqrt((x(i)-x(randind(k,1)))^2+(y(j)-y(randind(k,2)))^2) <= randradius(k)
                phi(i,j) = val+0.05*(2*rand(1)-1);
            end
        end
    end
end

% Initial maximum and minimum values of $\phi^0$
it = 0;
maxi(it+1) = max(phi,[],'all');
mini(it+1) = min(phi,[],'all');

% Visualize the initial condition
figure(1); clf; hold on; box on;
set(gca,'fontsize',21);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
surf(x,y,phi','EdgeColor','k','FaceColor','w')
view([-35 35])
axis([Lx Rx Ly Ry -1 1])
pbaspect([1 1 0.8])
text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.36277315404724 -0.697962928056722 0.11834709868257])
text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.310900714385596 0.549302625060356 -0.55962592008332])
text('Interpreter','latex','String','$\phi$','FontSize',24,'Position',[-0.180986474838287 0.947117932288588 0.917082453628964])
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
    % Maximum and minimum values of the numerical solution $\phi^n$
    maxi(it+1) = max(phi,[],'all');
    mini(it+1) = min(phi,[],'all');

    % Plot numerical solution
    if mod(it,50) == 0
        figure(2); clf; hold on; box on;
        set(gca,'fontsize',21);
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gcf,'PaperPositionMode','auto')
        surf(x,y,phi','EdgeColor','k','FaceColor','w')
        view([-35 35])
        axis([Lx Rx Ly Ry -1 1])
        pbaspect([1 1 0.8])
        text('Interpreter','latex','String','$x$','FontSize',24,'Position',[0.36277315404724 -0.697962928056722 0.11834709868257])
        text('Interpreter','latex','String','$y$','FontSize',24,'Position',[-0.310900714385596 0.549302625060356 -0.55962592008332])
        text('Interpreter','latex','String','$\phi$','FontSize',24,'Position',[-0.180986474838287 0.947117932288588 0.917082453628964])
        drawnow
    end
end

%% Plot the temporal evolution of the maximum and minimum values of $\phi^n$
lw = 1; ms = 8;
t = [0:it]*dt;
figure(3); clf; box on; hold on; grid on;
set(gca,'fontsize',16);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto')
plot(t,maxi,'b^-','linewidth',lw,'markersize',ms,'MarkerIndices',[1:50:Nt+1],'MarkerFaceColor','b');
plot(t,mini,'rv-','linewidth',lw,'markersize',ms,'MarkerIndices',[1:50:Nt+1],'MarkerFaceColor','r');
axis([t(1) t(end) -1.1 1.1])
xticks([0:2:8]*10^(-3));
xticklabels({'$0$','$0.002$','$0.004$','$0.006$','$0.008$'})
text('Interpreter','latex','String','$t$','FontSize',18,'Position',[0.00699087507725343 -1.27985196975726])
leg = legend('$Max(\phi^n)$','$Min(\phi^n)$');
set(leg,'interpreter','latex','Location','best')

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