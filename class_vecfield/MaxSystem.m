classdef MaxSystem
    %MAXSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % sim parameters
        M
        N
        dx
        dy
        dt

        % vector fields
        Xe
        Xb
        E
        B        
    end
    
    methods
        function obj = MaxSystem(  M, N, H, L, f1, f2, f3, g1, g2, g3 )
          
            obj.M  = M;
            obj.N  = N;           

            obj.Xe   = field( M, N, 'coords' );
            obj.Xe.x = L*obj.Xe.x;
            obj.Xe.y = H*obj.Xe.y;
            obj.dx   = obj.Xe.x(1,2) - obj.Xe.x(1,1);
            obj.dy   = obj.Xe.y(2,1) - obj.Xe.y(1,1);
            obj.Xb   = obj.Xe;
            obj.Xb.x = obj.Xb.x + (1/2)*obj.dx ;
            obj.Xb.y = obj.Xb.y + (1/2)*obj.dy ;

            obj.E = field( obj.Xe,  f1,  f2,   f3 );
            obj.B = field( obj.Xb,  g1,  g2,   g3 );

        end
        
        %------------------------------------------------------------------ 
        %

        function obj = eulerEvolve(obj, dt)
        % Time step routine
        % euler method
        % EM field in a vaccum
        % Peroidic BC
        %
            obj.B = obj.B.minus( scalarMult( curl( obj.E, obj.dx, obj.dy, 'periodic' ) , dt ) );
            obj.E = obj.E.plus(  scalarMult( curl( obj.B, obj.dx, obj.dy, 'periodic' ) , dt ) );
        end
        %------------------------------------------------------------------

        function [H] = plotsys(obj)

            colmap = flipud(brewermap([],'RdBu'));

            %figure(1)
            H = plotScalarField( obj.Xe, mag(obj.E) );
                clim([0,0.2])
                colormap( colmap )
                shading interp
                colorbar
            hold on
            plotField( obj.Xe, obj.E , 'interp', 40, 'r');    
            plotField( obj.Xb, obj.B , 'interp', 40, 'b' );
            hold off
            xlim([-1,1])
            ylim([-1,1])
            zlim([-1,1])
            xlabel('$x$')
            ylabel('$y$')
        end

        function [H] = plotsysZ(obj)

            colmap = flipud(brewermap([],'RdBu'));

            %figure(1)
            H = plotScalarField( obj.Xe, obj.E.z );
                clim([-0.3,0.3])
                colormap( colmap )
                shading interp
                colorbar
            hold on
            plotField( obj.Xe, obj.E , 'interp', 40, 'r');    
            plotField( obj.Xb, obj.B , 'interp', 40, 'b' );
            hold off
            xlim([-1,1])
            ylim([-1,1])
            zlim([-1,1])
            xlabel('$x$')
            ylabel('$y$')
        end
    end
end

