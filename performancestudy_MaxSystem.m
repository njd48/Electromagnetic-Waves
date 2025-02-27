

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

% Time study points

T_final = 0.6;

Nvals = [  2^7, 3*2^6, 2^8, 3*2^7, 2^9, 3*2^8, 2^10 ]';

wtime        = 0*Nvals;
wtime_per_dt = wtime;

Mtests = numel(Nvals);

for m = 1:Mtests

    [wtime(m), wtime_per_dt(m)] = time_study( Nvals(m), T_final );

end

table( Nvals, wtime, wtime_per_dt)

%%

    % Fit last 4 data
    p = polyfit( log(Nvals(end-3:end)), log(wtime_per_dt(end-3:end)), 1  );
    fline = @(x) exp(p(2))*x.^p(1);
    eqnstr  = ['fitline: $Ax^\alpha$, $\alpha = ', num2str(p(1),2), '$'];

figure(1)
    loglog( Nvals, wtime_per_dt, '--sk' )
    hold on
    loglog( Nvals(end-4:end), fline( Nvals(end-4:end) ) , '--m' )
    hold off
    text( Nvals(end-3)*1.15, wtime_per_dt(end-3), eqnstr )
    xlabel('simulation box width: $N$')
    ylabel('wall-time per evolution step: $t_w$ (s)')
    title('Performance study: Maxwells Eqs in a Plane')
    grid on

%%
function [wtime, wtime_per_dt] = time_study( N_x, T_final )

% ----------------------------------------- %
% --- Coordinates ------------------------- %
%
    %  N_x = 2^7 * 3; %---
    M = N_x;
    N = N_x;
    H = 1;
    L = 1;

    dt = 1/6 * ( max([L,H])/min([N,M]) );

% ----------------------------------------- %
% --- Example E-M wave -------------------- %
% --- Set Initial Condition --------------- %
%

    envelope  = @(  r  ) exp( -90*r.^2 );      
    sincblock = @(x,y,z) sinc( 40*sqrt( x.^2 + y.^2 ) ); 
    sinblock  = @(x,y,z) sin( 30*(x + y/3) );
    fE        = @(x,y,z) sinblock( (x+0.5), y, z ) .*envelope( (y) ).*envelope( (x+0.5) );

    %--------------------------------------------------------------
    sys = MaxSystem( M, N, H, L, ...
         @(x,y,z) 0.*x,   @(x,y,z) 0.*x,   @(x,y,z) fE(x,y,z), ...
         @(x,y,z) 0.*x,   @(x,y,z) 0.*x,   @(x,y,z) 0.*x       );
    %--------------------------------------------------------------

% ----------------------------------------- %
% --- Try Propagating  -------------------- %
%
    % T_final = 0.4; ---------------------
    t  = 0;
    

    P       = round(T_final/dt);
    % P_write = 4;

    wtime = 0;

    tic(); 

    for p = 1:P    

        sys = sys.eulerEvolve(dt);
        t   = t + dt;
        
         % if 0==mod(p,P_write)
         %     figure(2)
         %     sys.plotsys();
         %     title( sprintf('time $t = %2.4f $', t) );
         %     drawnow;
         % end
    end

    wtime = wtime + toc();  % Time performance

    wtime_per_dt = wtime/p;

end

    % fprintf( 'EM vaccum waves with sim environment as class...\n')
    % fprintf( 'simulation Box: %5i by %5i \n',    N_x, N_x);
    % fprintf( 'walltime per timestep: %1.8g (s)\n', wtime_per_dt );


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