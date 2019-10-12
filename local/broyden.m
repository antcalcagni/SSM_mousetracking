function x = broyden(x0,bounds,a,b,mu,sigma,sigmax,z)

    

    opt.maxiter=750;
    opt.tolfun=1e-8;
    opt.tolx=1e-6;
    
    x      = x0(:);
    it     = 0;
    %F      = feval(f, x)
    F = eq1(a,b,mu,sigma,sigmax,x,z);
        
    
    normf =  norm(F);
    J      = jacobi(F,x,a,b,mu,sigma,sigmax,z);  % Intial Jacobian matrix
    
%     if nargout > 1
%         ithist.x = [x(:)';zeros(opt.maxiter,length(x))];
%         ithist.f = [F(:)';zeros(opt.maxiter,length(x))];
%         ithist.normf = [normf;zeros(opt.maxiter,1)];
%     end
    
    normdx  = 2*opt.tolx;
    while( it < opt.maxiter+1 && normdx > opt.tolx && normf > opt.tolfun)           
        if rcond(J) < 1e-15
           error('Singular jacobian at iteration %d\n',it)
        end
       dx     = -J\F;
       normdx = norm(dx);

       if length(bounds)>1
           % make sure x stays within bounds
           for j = 1:20
               jl = find(x+dx<bounds(:,1));
               dx(jl) = dx(jl)/2;
               ju = find(x+dx>bounds(:,2));
               dx(ju) = dx(ju)/2;
               if isempty(jl) && isempty(ju)
                   break
               end
           end
       end
       

       x  = x+dx;
       it = it+1;
       %F  = feval(f, x);
       F = eq1(a,b,mu,sigma,sigmax,x,z);
       normf = norm(F);
       J  = J + F*dx'/(dx'*dx);  
       

       if nargout > 1
           ithist.x(it+1,:) = x(:)';
           ithist.f(it+1,:) = F(:)';
           ithist.normf(it+1,:) = normf;
       end
    end

    % Check if the iterations converged and issue warning if needed
    if it >= opt.maxiter && norm(F) > opt.tolfun
        warning('No convergence in %d iterations.\n',it+1)
%     elseif normf>opt.tolfun
%         %warning('Newton step < %g, but function norm > %g\n',...
%             opt.tolx,opt.tolfun)
%     elseif normdx>opt.tolx
%         %warning('Function norm < %g, but newton step norm > %g\n',...
%             opt.tolfun,opt.tolx)        
    end
%     if nargout > 1
%         ithist.x(it+2:end,:) = [];
%         ithist.f(it+2:end,:) = [];
%         ithist.normf(it+2:end) = [];
%     end
end

function J = jacobi(y0,x,a,b,mu,sigma,sigmax,z)
% Quick and dirty numerical Jacobian for function f at x
% y0: f(x);

    delta = 1e-6*(max(1,sqrt(norm(x))));
    n = length(y0);
    m = length(x);
    J = zeros(n,m);
    for i = 1:m
        dx = zeros(m,1);
        dx(i) = delta/2;
        %J(:,i) = (feval(f,x+dx)-feval(f,x-dx))/delta;
        f1 = eq1(a,b,mu,sigma,sigmax,x+dx,z);
        f2 = eq1(a,b,mu,sigma,sigmax,x-dx,z);
        J(:,i) = (f1-f2)/delta;
    end
end
