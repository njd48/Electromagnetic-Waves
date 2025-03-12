classdef MaxSystemDielectric
    %MAXSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % sim parameters
        M
        N
        dx
        dy
        dt

        % Dielectric Properties
        % mu let mu_0 = 1
        % ep % const scalar field
        % reduce number of divisions by storing inverse
        invep

        % array of class shape
        shapes

        % vector fields
        Xe
        Xb
          
        D
        B
    end
    
    methods

        function obj = MaxSystemDielectric(  M, N, H, L, f1, f2, f3, g1, g2, g3 )
          
            obj.M  = M;
            obj.N  = N;           

            % Set up coordinate system
            obj.Xe   = field( M, N, 'coords' );
            obj.Xe.x = L*obj.Xe.x;
            obj.Xe.y = H*obj.Xe.y;
            obj.dx   = obj.Xe.x(1,2) - obj.Xe.x(1,1);
            obj.dy   = obj.Xe.y(2,1) - obj.Xe.y(1,1);
            obj.Xb   = obj.Xe;
            obj.Xb.x = obj.Xb.x + (1/2)*obj.dx ;
            obj.Xb.y = obj.Xb.y + (1/2)*obj.dy ;

           
            % dielectric permeability is automatically set
            % to 1 in this function
            obj.invep = ones(M,N);   %<- we store its inverse to avoid division

            % Set initial condition for E and B fields
            obj.D = field( obj.Xe,  f1,  f2,   f3 );    %< E field temp goes here
            obj.B = field( obj.Xb,  g1,  g2,   g3 );

            % Start with empty array of shape
            obj.shapes = shape('empty',0,0,0);
            obj.shapes(1) = [];

        end

        function obj = appendShape(  obj, shape )
            %Appends a shape to the shapes list
            obj.shapes(end+1) = shape; 
        end

        function obj = setPermeabilityDomain( obj )

            domain_ep = obj.ep;

            for n = 1:numel(obj.shapes)

                %Set values for scalar field ep
                domain_ep = domain_ep + (obj.shapes(n).ep-1) *obj.shapes(n).isinside( obj.Xe );
            end

            % Store it in the simulation environment
            obj.invep = 1./domain_ep; 

            % Rescale E field to the D field
            obj.D = obj.D.scalarMult( domain_ep );
        end

        %------------------------------------------------------------------ 
        %
        %  Methods for auxilliary properties, to reduce storage

        function Epermeability = ep( obj )
            Epermeability = 1./(obj.invep) ;
        end

        function Efield = E( obj )
            Efield = obj.D.scalarMult( obj.invep );
        end

        function Pfield = P( obj )
            Pfield = obj.D.scalarMult( ( obj.ep-1 ).*obj.invep );
        end
        
        
        %------------------------------------------------------------------ 
        %

        function obj = eulerEvolve(obj, dt)
        % Time step routine
        % euler method
        % EM field in a dielectric
        % Peroidic BC
        %
                                                  % vv-- Need E = 1/ep D --vv
            obj.B = obj.B.minus( scalarMult( curl( obj.D.scalarMult(obj.invep),...
                            obj.dx, obj.dy, 'periodic' ) , dt ) );

            obj.D = obj.D.plus(  scalarMult( curl( obj.B,                      ...
                            obj.dx, obj.dy, 'periodic' ) , dt ) );
        end
        %------------------------------------------------------------------

        function [H] = plotsys(obj, varargin)

            colmap = flipud(brewermap([],'RdBu'));

            %figure(1)
            H = plotScalarField( obj.Xe, mag(obj.E) );
                clim([0,0.2])
                colormap( colmap )
                shading interp
                colorbar
            hold on
            plotField( obj.Xe, obj.P , 'interp', 40, 'r');    
            plotField( obj.Xb, obj.B , 'interp', 40, 'b' );
            hold off
            xlim([-1,1])
            ylim([-1,1])
            zlim([-1,1])
            xlabel('$x$')
            ylabel('$y$')

            for k = 1:numel(varargin)

                switch varargin{k}

                    case 'withshapes'
                        hold on
                        for i = 1:numel(obj.shapes)
                            obj.shapes(i).plot()
                        end
                        hold off

                    otherwise
                        error(['option, ', varargin{k}, 'not recognized\n']);
                end
            end

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
            plotField( obj.Xe, obj.P , 'interp', 40, 'r');    
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

