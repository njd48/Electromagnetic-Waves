function Du = grad( u, dx, dy, bc)

    [M,N]  = size(u);

    Du   = field(M,N);

    ep_x  = (1/2/dx);
    ep_y  = (1/2/dy);

    switch bc
        case 'neumann'
            
            %x
            Du.x(:,2:end-1) = ep_x*( u(:,3:end) - u(:,1:end-2) );

            %y
            Du.y(2:end-1,:) = ep_y*( u(3:end,:) - u(1:end-2,:) );

            %z                    
            %zero

        otherwise 
            error('cannot take grad with provided b.c.');
    end
end
