function [ J ] = kcurrent( X, dx, dy, I, t )

    % For a signal use fnc handle
    % @(t)I(t)

    J = field( X, @(x,y,z) 0*x, @(x,y,z) 0*x, @(x,y,z) I(t)*gau(x,y,dx,dy) );


end

function [z] = gau(x,y,dx,dy)

    z = 4/sqrt(pi)/(dx*dy) * exp( -(x/dx/2).^2 -(y/dy/2).^2 )  ;

end
