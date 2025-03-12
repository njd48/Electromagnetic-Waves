

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

save_gif_anim = 1;

%%

% ----------------------------------------- %
% --- Coordinates ------------------------- %
%
    N_x = 2^7 * 4; 
    M = N_x;
    N = N_x;
    H = 1;
    L = 1;

    % Note: edit so dt is not a property of the system
    %       but provided externally for alternate TS routines
    dt = 1/6 * ( max([L,H])/min([N,M]) );

% ----------------------------------------- %
% --- Example E-M wave -------------------- %
%

    envelope = @( r ) exp( -90*r.^2 ); 
     
    sincblock = @(x,y,z) sinc( 40*sqrt( x.^2 + y.^2 ) ); 

    sinblock = @(x,y,z) sin( 30*(x + y/3) );

    fE   = @(x,y,z) sinblock( (x+0.1), y, z ) .*envelope( (y) ).*envelope( (x+0.1) );

    %--------------------------------------------------------------
    sys = MaxSystemDielectric( M, N, H, L, ...
         @(x,y,z) 0.*x,   @(x,y,z) 0.*x,   @(x,y,z) 0.*x , ...
         @(x,y,z) 0.*x,   @(x,y,z) 0.*x,   @(x,y,z) fE(x,y,z)      );
    %--------------------------------------------------------------
    % System at this point is a vaccum

    % Append shapes with more permeability 
    shape1 = shape( 'square', [1,0], 1,  5 );
    shape2 = shape( 'circle', [-1,1], 1, 2 );
    sys = sys.appendShape(shape1).appendShape(shape2);

    % Set
    sys = sys.setPermeabilityDomain();

    % Plot surface of the system permeability
    figure(1)
    plotScalarField( sys.Xe, sys.ep )
    colormap(colmap)
    shading interp
    colorbar

    % Plot initial condition (E,B)
    figure(2)
    sys.plotsys('withshapes')
    % colormap(colmap)
    % shading interp
    % colorbar   
    

    %%
% ----------------------------------------- %
% --- Try Propagating  -------------------- %
%



    t  = 0;
    %dt = dx/6;

    T_final = 0.8;

    P       = round(T_final/dt);
    P_write = 16;

    wtime = 0;

    if save_gif_anim
        fname   = 'EM_wave_dielec.gif';
        nframes = floor(P/P_write);

        %Frames(nframes) = getframe();
        q = 1;
    end

    for p = 1:P    

        tic(); 

        sys = sys.eulerEvolve(dt);
        t   = t + dt;

        wtime = wtime+toc();  % Time performance


        if 0==mod(p,P_write)
            figure(2)
            sys.plotsys('withshapes');
            legend('$|\vec{E}|$','$\vec{P}$','$\vec{B}$' )            
            title( sprintf('time $t = %2.4f $', t) );
            %annotations
            text( -0.65, -0.85,'vacuum, $\varepsilon = 1$', 'Color','w' )
            text( -0.8, 0.7, ['$\varepsilon = ', num2str(shape2.ep), '$'], 'Color','w' )
            text( 0.46, -0.55, ['$\varepsilon = ', num2str(shape1.ep), '$'], 'Color','w' )
            drawnow;

            if save_gif_anim
                %Frames(q) = getframe(gca);
                exportgraphics(gcf, fname,'Append',true);
                %q = q+1;
            end
        end
    end

    if save_gif_anim
    end


%%
    wtime_per_dt = wtime/p;
    fprintf( 'EM vaccum waves with sim environment as class...\n')
    fprintf( 'simulation Box: %5i by %5i \n',    N_x, N_x);
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