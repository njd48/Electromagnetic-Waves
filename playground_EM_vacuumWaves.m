

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultFunctionLineLineWidth',2)
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex');

%%

clear all, close all;

addpath('class_vecfield');
addpath('Maxwell');

colmap = flipud(brewermap([],'RdBu'));

%%

% ----------------------------------------- %
% --- Coordinates ------------------------- %
%
    N_x = 2^7 * 3; 
    X   = field( N_x, N_x, 'coords' );
    dx  = X.x(1,2) - X.x(1,1);
    Xb  = X.addcnst( (1/2)*dx );

% ----------------------------------------- %
% --- Example E-M wave -------------------- %
%
    % % rhat = X.scalarMult(1./ mag(X));
    % % 
    % % vort_x = @(x,y,z) -y;
    % % vort_y = @(x,y,z)  x;
    % % vort_z = @(x,y,z) 0*x;
    % % 
    % % that = field( X, vort_x, vort_y, vort_z );
    % % that = that.scalarMult( 1./ mag(that) );
    % % 
    % % zhat = field( X, @(x,y,z) 0*x, @(x,y,z) 0*x, @(x,y,z) 1 + 0*x );
    % % 

    envelope = @( r ) exp( -90*r.^2 ); 
     
    sincblock = @(x,y,z) sinc( 40*sqrt( x.^2 + y.^2 ) ); 

    sinblock = @(x,y,z) sin( 30*(x + y/3) );

    % % coscblock = @(x,y,z) cosc( 12*sqrt( x.^2 + y.^2 ) ); 
    % % 

    fE   = @(x,y,z) sinblock( (x+0.5), y, z ) .*envelope( (y) ).*envelope( (x+0.5) );

    % % 
    % % E = that.scalarMult( fE( X.x, X.y, X.z ) );
    % % B = zhat.scalarMult( fE( Xb.x, Xb.y, Xb.z ) );
    % % 
    E = field( X,   @(x,y,z) 0.*x,   @(x,y,z) 0.*x,   @(x,y,z) fE(x,y,z) );
    B = field( Xb,  @(x,y,z) 0.*x,   @(x,y,z) 0.*x,   @(x,y,z) 0.*x );

   
    figure(1)
    plotScalarField( X, mag(E) );
        clim([0,0.2])
        colormap( colmap )
        shading interp
        colorbar
    hold on
    plotField( X, E , 'interp', 40, 'r');    
    plotField( Xb, B, 'interp', 40 , 'b' );
    hold off
    xlim([-1,1])
    ylim([-1,1])
    zlim([-1,1])
    xlabel('$x$')
    ylabel('$y$')

    %%
% ----------------------------------------- %
% --- Try Propagating  -------------------- %
%

    t  = 0;
    dt = dx/6;

    T_final = 0.4;

    P       = round(T_final/dt);
    P_write = 24;

    wtime = 0;

    for p = 1:P

        tic(); % Time performance

        [E,B] = vaccum_evolve_once( X, E, B, dx, dx, dt, t );

        t = t + dt;

        wtime = wtime + toc();  % Time performance
    
        if 0==mod(p,P_write)
            figure(2)
            plotScalarField( X, E.mag() );
            clim([0,0.1])
            colormap( colmap )
            shading interp
            colorbar
            hold on
            plotField( X, E, 'interp', 40, 'r' );
            plotField( Xb, B, 'interp', 40 , 'b' );
            %plotField( X, E.cross(B), 'interp', 45 , 'm' );
            hold off
            xlim([-1,1])
            ylim([-1,1])
            zlim([-1,1])
            xlabel('$x$')
            ylabel('$y$')
            title( sprintf('time $t = %2.4f $', p*dt) )
            drawnow;
        end

    end
%%
    wtime_per_dt = wtime/p;
    fprintf( 'first test of EM vacuum waves...\n')
    fprintf( 'simulation Box: %5i x %5i \n',    N_x, N_x);
    fprintf( 'walltime per timestep: %1.8g (s)\n', wtime_per_dt );


    %% Fnc


function [y] = sinc(x)

    y = 0*x;

    for i = 1:numel(x)
        if x(i) == 0
            y(i) = 1;
        else 
            y(i) = sin(pi*x(i))/(pi*x(i));            
        end
    end

end

function [y] = cosc(x)

    y = 0*x;

    for i = 1:numel(x)
        if x(i) == 0
            y(i) = 0;
        else 
            y(i) = cos(pi*x(i))/(pi*x(i));            
        end
    end

end