
function [D2u] = scalar_laplacian( u, dx, dy )

    ep_x  = (1/dx^2);
    ep_y  = (1/dy^2);

    D2u = 0*u;

    D2u(2:end-1,2:end-1) = ...
          ep_x*( u(2:end-1,3:end) + u(2:end-1,1:end-2) - 2*u(2:end-1,2:end-1) )  ...
        + ep_y*( u(3:end,2:end-1) + u(1:end-2,2:end-1) - 2*u(2:end-1,2:end-1) );
end
