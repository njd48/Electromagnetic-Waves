function [ E2, B2 ] = vaccum_evolve_once( X, E, B, dx, dy, dt, t )

    % % Current source?
    % tint = 0.002;
    % H    = 1;
    % I = @(t) H*( tint^2 - (t-tint).^2 ).*(t<=tint) ; %+ 0*(t>dt)


    %sys
   
    B2 =  B - scalarMult( curl(  E, dx, dy, 'periodic' ) , dt );

    E2 =  E + scalarMult( curl( B2, dx, dy, 'periodic' )  , dt);

    

end

