classdef field
    %FIELD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M %(y)
        N %(x)
        x
        y
        z
    end
    
    methods
        
        %------------------------------------------------------------------
        % Routines for initializing a field
        % 
        function obj = field( varargin )

            if nargin <= 2 

                obj.M = varargin{1};
                obj.N = varargin{2};

                obj.x = zeros(obj.M,obj.N);
                obj.y = zeros(obj.M,obj.N);
                obj.z = zeros(obj.M,obj.N);

            elseif nargin == 3
                obj.M = varargin{1};
                obj.N = varargin{2};

                switch varargin{3}
                    case 'coords'
                        X = linspace(-1,1,obj.N);
                        Y = linspace(-1,1,obj.M);
                        [X,Y] = meshgrid(X,Y);
                        obj.x = X;
                        obj.y = Y;
                        obj.z = 0;  % <-- careful, scalar cnst zero to save mem
                    otherwise 
                        error('this case of field initialization is not programmed')
                end


            elseif nargin == 4
                X = varargin{1}; % coordinates, also a field
                f = varargin{2}; % function handles
                g = varargin{3};
                h = varargin{4};

                obj = field( X.M, X.N );

                obj.x = f( X.x, X.y, X.z );
                obj.y = g( X.x, X.y, X.z );
                obj.z = h( X.x, X.y, X.z );

            else
                error('unable to parse narg when initializing field');
            end
            
        end
        %------------------------------------------------------------------
        


        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % Algebra
        % 
        function c = uminus(a)
            c = a;
            c.x = -a.x;
            c.y = -a.y;
            c.z = -a.z;
        end
        function c = plus(a,b)
            c   = a;
            c.x = a.x + b.x;
            c.y = a.y + b.y;
            c.z = a.z + b.z;
        end
        function a = addcnst(a, S)
            a.x = a.x + S;
            a.y = a.y + S;
            a.z = a.z + S;
        end
        function c = minus(a,b)
            c = plus(a,-b);
        end

        % Multiplications
        function c = times(a,b)
            c   = a;
            c.x = a.x .* b.x;
            c.y = a.y .* b.y;
            c.z = a.z .* b.z;
        end
        function c = rdivide(a,b)
            c   = a;
            c.x = a.x ./ b.x;
            c.y = a.y ./ b.y;
            c.z = a.z ./ b.z;
        end
        function a = scalarMult(a, S)
            a.x = a.x .*S;
            a.y = a.y .*S;
            a.z = a.z .*S;
        end        
        function c = power(a,b)
            c = a;
            c.x = a.x .^ b.x;
            c.y = a.y .^ b.y;
            c.z = a.z .^ b.z;
        end
        

        % vec products
        function c = dot(u,v)
            c = (u.x .* v.x) + (u.y .* v.y) + (u.z .* v.z) ;
        end
        function c = cross(u,v)
            c = u;
            c.x = (u.y.*v.z) - (u.z.*v.y);
            c.y = (u.z.*v.x) - (u.x.*v.z);
            c.z = (u.x.*v.y) - (u.y.*v.x);
        end
        function c = mag( u )
            c   = sqrt( dot(u,u) );
        end
        %------------------------------------------------------------------


        %------------------------------------------------------------------
        % Grad, Divergence, curl, and laplacian
        %
        % (grad.m and scalar_laplacian.m are functions of scalar fields
        %  which output vector fields. programmed external to this
        %  field.m file )
        
        function DivU = divg( U, dx, dy, bc )

            DivU = zeros( U.M, U.N );
            ep_x  = (1/2/dx);
            ep_y  = (1/2/dy);

            switch bc
                case 'neumann'
                    
                    DivU(:,2:end-1) = ep_x*( U.x(:,3:end) - U.x(:,1:end-2) );
                    DivU(2:end-1,:) = DivU(2:end-1,:) ...
                                    + ep_y*( U.y(3:end,:) - U.y(1:end-2,:) );

                otherwise 
                    error('cannot take div with provided b.c.');
            end
        end
        %------------------------------------------------------------------

        function CU = curl( U, dx, dy, bc )

            CU   = field( U.M, U.N );
            ep_x  = (1/2/dx);
            ep_y  = (1/2/dy);

            switch bc
                case 'neumann'
                    
                    % x
                    CU.x(2:end-1,:) = ep_y*( U.z(3:end,:) - U.z(1:end-2,:) );

                    % y
                    CU.y(:,2:end-1) = -ep_x*( U.z(:,3:end) - U.z(:,1:end-2) );

                    % z                    
                    CU.z(:,2:end-1) = ep_x*( U.y(:,3:end) - U.y(:,1:end-2) );
                    CU.z(2:end-1,:) = CU.z(2:end-1,:) ...
                                    - ep_y*( U.x(3:end,:) - U.x(1:end-2,:) );

                case 'periodic'
                    % x
                    CU.x(2:end-1,:) = ep_y*( U.z(3:end,:) - U.z(1:end-2,:) );

                    CU.x(1,:)   = ep_y*( U.z(2,:) - U.z(end,:) );
                    CU.x(end,:) = ep_y*( U.z(1,:) - U.z(end-1,:) );

                    % y
                    CU.y(:,2:end-1) = -ep_x*( U.z(:,3:end) - U.z(:,1:end-2) );

                    CU.y(:,1)   = -ep_x*( U.z(:,2) - U.z(:,end) );
                    CU.y(:,end) = -ep_x*( U.z(:,1) - U.z(:,end-1) );

                    % z                    
                    CU.z(2:end-1,2:end-1) = ep_x*( U.y(2:end-1,3:end) - U.y(2:end-1,1:end-2) );
                    CU.z(2:end-1,2:end-1) = CU.z(2:end-1,2:end-1) ...
                                    - ep_y*( U.x(3:end,2:end-1) - U.x(1:end-2,2:end-1) );

                    CU.z(:, 1)   = CU.z(:, 1) ...
                                    + ep_x*( U.y(:,2) - U.y(:,end) );
                    CU.z(:, end) = CU.z(:, end) ...
                                    + ep_x*( U.y(:,1) - U.y(:,end-1) );

                    CU.z( 1, :)   = CU.z( 1, :) ...
                                    - ep_y*( U.x(2,:) - U.x(end,:) );
                    CU.z( end, :) = CU.z( end, :) ...
                                    - ep_y*( U.x(1,:) - U.x(end-1,:) );


                otherwise 
                    error('cannot take curl with provided b.c.');
            end
        end

        %------------------------------------------------------------------
        function [D2U] = laplacian( U, dx, dy )
            
            D2U   = field( U.M, U.N );
            D2U.x = scalar_laplacian( U.x, dx, dy );
            D2U.y = scalar_laplacian( U.y, dx, dy );
            D2U.z = scalar_laplacian( U.z, dx, dy );

        end

        

        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % Plotting function
        %
        function H = plotField(X, U, varargin)

            if nargin == 2
                H = quiver3( X.x, X.y, 0*X.x, U.x, U.y, U.z );
            else 
                switch varargin{1}
                    case 'interp'
                        n_r = varargin{2};
                        k_x = log(n_r/X.N)/log(2); 
                        k_y = log(n_r/X.M)/log(2);

                        xm   = interp2( X.x, k_x );
                        ym   = interp2( X.y, k_y );
                        um   = interp2( X.x, X.y, U.x, xm, ym );
                        vm   = interp2( X.x, X.y, U.y, xm, ym );
                        wm   = interp2( X.x, X.y, U.z, xm, ym );

                        if nargin == 4
                            H = quiver3( xm, ym, 0*xm, um, vm, wm );
                        else 
                            H = quiver3( xm, ym, 0*xm, um, vm, wm, ...
                                varargin{3:end} );
                        end
                    otherwise
                        H = quiver3( X.x, X.y, 0*X.x, U.x, U.y, U.z, ...
                            varargin{:} );
                end
            end
            
        end

        function H = plotScalarField(X, u, varargin )

            if nargin == 2
                H = pcolor( X.x, X.y, u );
            else 
                H = pcolor( X.x, X.y, u, varargin{:} );
            end
            
        end
        %------------------------------------------------------------------
    end
end

