
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultFunctionLineLineWidth',2)
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex');

%% Test class shape


%----------------------------------------------------
% Initialize shapes 
%
%        shape( type, [c_x, c_y], radius, permeability_property );
%
shape1 = shape( 'square', [1,0], 1, 1 );
shape2 = shape( 'circle', [-1,1], 1, 1 );

%Array of shapes
shapes = [shape1, shape2];

%----------------------------------------------------
% Test plots
%
figure(1)
for i = 1:numel(shapes)
    shapes(i).plot()
    hold on
end
hold off
axis([ -3, 3, -3, 3 ])
grid on
title('shapes')

%----------------------------------------------------
% test locating points
%
    n_test_pnts  = 4;

    pnts = -1 + 3*rand(n_test_pnts, 2);

    pnt_class= cell(1,n_test_pnts);

     % pnts        = [0.75498    -0.25251];
     % n_test_pnts = 1;

    % pnts        = [0.9330    0.0522];
    % n_test_pnts = 1;

    for n = 1:n_test_pnts

        in_shape = zeros( 1, numel(shapes) );

        tail = ' ';

        for i = 1:numel(shapes)

            in_shape(i) = shapes(i).isPntInside( pnts(n,:) );
            
            if in_shape(i)
                tail = [tail, num2str(i), ', ' ];
            end
        end
        
        if sum(in_shape) == 0
            pnt_class{n} = 'outside all shapes';
        else
            pnt_class{n} = ['inside shape ', tail(1:end-2)];
        end

    end

    % plot these
        horiz_lag = 0.175;
        vert_lag  = -0.1;
    figure(1)
    for n = 1:n_test_pnts

        hold on
        plot( pnts(n,1), pnts(n,2), 'o' )
        text( pnts(n,1)+horiz_lag, pnts(n,2)+vert_lag, pnt_class{n} )

    end
    hold off

%----------------------------------------------------
% test compatibility with field class
%

M = 40;  N = 41;

X =  field( 40,41, 'coords' );
X =  X.scalarMult(3);

indicator = zeros( M, N);

for i = 1:numel(shapes)
    indicator = indicator + shapes(i).isinside(X);
end

% Draw a heatmap of the points which are inside shapes
figure(1)
    plotScalarField(X, indicator )
    shading interp
    hold on
    for i = 1:numel(shapes)
        shapes(i).plot()
    end
    hold off
    axis([ -3, 3, -3, 3 ])
    grid on
    title('shapes')