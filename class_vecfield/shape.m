classdef shape
    %SHAPE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        center
        radius
        ep = 1;
    end
    
    methods
        %------------------------------------------------------------------
        % Initializations
        %
        function obj = shape( shapetype, center, radius, Epermeability )
                % Presently have methods for 
                %   shapetype = 
                %       'circle'
                %       'square'
            obj.type   = shapetype;
            obj.center = center;     %(2-tuple)
            obj.radius = radius;
            obj.ep     = Epermeability;
        end
        
        %------------------------------------------------------------------
        % locating methods
        %
        function inside = isPntInside( obj, p )
                % p is a 2-tuple
                switch obj.type

                    case 'circle' %-------------------------------------
                        if norm( p - obj.center , 2 ) <= obj.radius
                            inside = 1;
                        else 
                            inside = 0;
                        end

                    case 'square' %-------------------------------------
                        a = 1/sqrt(2)*obj.radius;
                        if abs( p(1) - obj.center(1) ) <= a % within x bounds
                            if abs( p(2) - obj.center(2) ) <= a % within y bounds
                                inside = 1;
                            else 
                                inside = 0;
                            end
                        else
                            inside = 0;
                        end

                    otherwise
                        error(['shape has unknown type \"', obj.type, '\" \n']);
                end
            


        end
        function inside = isinside( obj, coordfield )

            inside = zeros( coordfield.M, coordfield.N );

            for m = 1:coordfield.M
                for n = 1:coordfield.N
                    inside(m,n) = obj.isPntInside( [coordfield.x(m,n), coordfield.y(m,n)] );
                end
            end
        end

        %------------------------------------------------------------------
        % drawing methods
        %
        function [x,y] = boundary(obj)      

            switch obj.type
                case 'square' % draw the corners
                    x = obj.center(1) + obj.radius/sqrt(2)*[1, -1, -1,  1, 1];
                    y = obj.center(2) + obj.radius/sqrt(2)*[1,  1, -1, -1, 1];
                case 'circle' % draw a curve 
                    t = linspace( 0, 2*pi, 50 );
                    x = obj.center(1) + obj.radius * cos(t);
                    y = obj.center(2) + obj.radius * sin(t);
                otherwise
                    error(['shape has unknown type \"', obj.type, '\" \n']);
            end
        end

        function [H] = plot(obj)

            [x,y] = obj.boundary();

            plot( x, y, '-k')

        end
    end
end

