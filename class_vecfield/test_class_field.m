set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultFunctionLineLineWidth',2)
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex');

%% Testing funcitonality of field class
%

%--------------------------------------------------
% Initialize a coordinate system
%

X = field( 10, 11, 'coords' );

dx = X.x(1,2) - X.x(1,1);
dy = X.y(2,1) - X.y(1,1);

%--------------------------------------------------
% Make Vortex using field init
%

vort_x = @(x,y,z) -y;
vort_y = @(x,y,z)  x;
vort_z = @(x,y,z) 0*x;

U = field( X, vort_x, vort_y, vort_z );

%--------------------------------------------------
% Plotting
%
figure(1)
plotField( X, U );

%-------------------------------------------------
% Test Algebur
%
V =  - U;
V = U + X;
V = U - X;
V = U .* X;
V = U ./ X;
V = scalarMult( U, X.x );
V = dot(U,X);
V = cross(U,X);
V = mag(U);


% %
% figure()
% plotField( X, U );



%--------------------------------------------------
%% Divergence
%

% Make a denser coordinate sys
%
X = field( 12, 11, 'coords' ) + field( 12, 11, 'coords' );

dx = X.x(1,2) - X.x(1,1);
dy = X.y(2,1) - X.y(1,1);


% Init the field
%
f = @(a,b) atan( a./sqrt(1+b.^2) ) ./(sqrt(1+b.^2));
fx = @(x,y,z) f( x, y);
fy = @(x,y,z) f( y, x);
fz = @(x,y,z) 0*x;

fdiv = @(x,y,z) 2./( 1 + x.^2 + y.^2 );

DivU_known = fdiv( X.x, X.y, X.z ); %<-- known val of divg

U = field( X, fx, fy, fz );

% Get Divergence
%
DivU = divg( U, dx, dy, 'neumann' );

figure(2)
plotScalarField( X, DivU);
shading interp
clim([0,2])
colorbar
hold on
plotField( X, U );
hold off
grid on
title('Field and divergence')
% 
% figure(3)
% plotScalarField( X, DivU_known);
% shading interp
% clim([0,2])
% colorbar
% hold on
% plotField( X, U );
% hold off
% title('Field and analytical divergence')

figure(3) 
plotScalarField( X, DivU - DivU_known );
shading interp
colorbar
title('pointwise errors in divergence')

%% Try Laplacian
%

% Make a denser coordinate sys
%
X = scalarMult( field( 2*24,2*24, 'coords' ), 2);

dx = X.x(1,2) - X.x(1,1);
dy = X.y(2,1) - X.y(1,1);

% Init the field
%
f = @(x,y) sin(pi*x).*cos(pi*y);
fx = @(x,y,z) f( x, y);
fy = @(x,y,z) f( y, x);
fz = @(x,y,z) 0*x;

fdiv = @(x,y,z) 2./( 1 + x.^2 + y.^2 );

U = field( X, fx, fy, fz );

D2U_known = scalarMult( U, (-2*pi^2) ); %<-- known val of lapl

% Get vectir laplacian
%
D2U = laplacian( U, dx, dy );

% Plots
%
figure(4)
plotField( X,D2U );
xlim([-2,2])
ylim([-2,2])
title('vec laplacian')
xlabel('x')
ylabel('y')

figure(5)
plotField( X,D2U_known );
xlim([-2,2])
ylim([-2,2])
title('vec laplacian (analytical)')
xlabel('x')
ylabel('y')


figure(6)
%plotField( X, ( D2U-D2U_known ) );
plotScalarField( X, (mag( D2U-D2U_known )) );
shading interp
clim([0,1])
colorbar
xlim([-2,2])
ylim([-2,2])
title('pointwise errors in vector laplacian')
xlabel('x')
ylabel('y')


%% Try curl
%
%
%--------------------------------------------------
% Initialize a coordinate system
%

X = field( 15, 16, 'coords' );
X = scalarMult( X, 2);

dx = X.x(1,2) - X.x(1,1);
dy = X.y(2,1) - X.y(1,1);

%--------------------------------------------------
% Make Vortex using field init
%

vort_x = @(x,y,z) -y;
vort_y = @(x,y,z)  x;
vort_z = @(x,y,z) 0*x;

U = field( X, vort_x, vort_y, vort_z );

CU = curl(U, dx, dy, 'neumann'); % Should be cnst 2*k

CU_known = field( X, vort_z, vort_z, @(x,y,z)2 );

figure(7)
plotField( X,CU );
xlim([-2,2])
ylim([-2,2])
title('curl ')
xlabel('x')
ylabel('y')

figure(8)
plotScalarField( X, mag(CU-CU_known) );
xlim([-2,2])
ylim([-2,2])
clim([0,1])
shading interp
colorbar
title('pntwise error in curl')
xlabel('x')
ylabel('y')

%% Try gradient

%--------------------------------------------------
% Initialize a coordinate system
%
m = 5;
X = field( m*8, m*8+1, 'coords' );
X = scalarMult( X, 2);

dx = X.x(1,2) - X.x(1,1);
dy = X.y(2,1) - X.y(1,1);

%---------------------------------------------------
% Make a Scalar field
%
phi = 1./( 1 + dot(X,X) );

grad_phi_known = scalarMult(X, -2./( 1+dot(X,X) ).^2 );

%----------------------------------------------------
% Take Gradient
%
grad_phi = grad( phi, dx, dy, 'neumann');

figure(7)
plotScalarField( X, phi );
shading interp
%clim([-2,0])
hold on
plotField( X, grad_phi );
hold off
xlim([-2,2])
ylim([-2,2])
title('gradient of scalar field')
xlabel('x')
ylabel('y')

figure(8)
plotScalarField( X, mag( grad_phi - grad_phi_known ) );
xlim([-2,2])
ylim([-2,2])
%clim([0,1])
shading interp
colorbar
% hold on
% plotField( X, grad_phi - grad_phi_known );
% hold off
title('pntwise error in gradient')
xlabel('x')
ylabel('y')