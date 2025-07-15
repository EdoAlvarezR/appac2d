function [wakes,gamma,iter,E] = solveWake(foils,Ainv,RHS,CT,opts)
% SOLVEWAKE  Solve for the global circulation solution and wake shape using
% Shollenberger's algorithm:
%   1. Guess initial wake shape and bound circulation
%   2. Calculate induced velocities on the airfoil surfaces due to the wake
%   3. Solve circulation strengths on the airfoil surfaces such that there is
%      no flow through the solid boundaries
%   4. Calculate induced velocities on the wake
%   5. Move the wake panels to satisfy flow tangency on the wake boundaries
%   6. Calculate bound circulation on the wake from the newly calculated
%      average velocity across the wake boundary
%   7. If the new wake circulation is close to the previous solution, exit the
%      program, else return to step 2

% Attach wake to the trailing edge of each propulsive element %%%%%%%%%%%%%%%%
N = opts.NumPanels + 1; % add far-field panel to the panel count
wakes.m = [N N];
[xinit,yinit,traj] = initWake(foils,Ainv*RHS); % obtain starting wake shape
for i = 2:-1:1
    k = (i-1)*N+(1:N);
    x = xinit{i}; y = yinit{i}; n = numel(x);
    stencil = linspace(0,opts.WakeLengthChords,N-n+1).';
    % if strcmpi(opts.NodeSpacing,'cosine')
    %     wakes.xo(k,:) = foils.xo(1+(i-1)*foils.m(1)) + ...
    %         opts.WakeLengthChords*(1-cos(linspace(0,pi/2,N))).';
    % elseif strcmpi(opts.NodeSpacing,'uniform')
    %     wakes.xo(k,:) = foils.xo(1+(i-1)*foils.m(1)) + ...
    %         linspace(0,opts.WakeLengthChords,N).';
    % end
    wakes.xo(k,1) = [x(1:n-1); x(n)+traj(1)*stencil];
    wakes.yo(k,1) = [y(1:n-1); y(n)+traj(2)*stencil];
    wakes.dx(k,1) = [diff(wakes.xo(k)); 1e3*traj(1)];
    wakes.dy(k,1) = [diff(wakes.yo(k)); 1e3*traj(2)];
end
wakes.theta = atan2(wakes.dy,wakes.dx);
wakes.co = [wakes.xo+wakes.dx/2 wakes.yo+wakes.dy/2];
wakes.ds = sqrt(wakes.dx.^2 + wakes.dy.^2);
wakes.so = [0; cumsum(wakes.ds(1:N-1)); 0; cumsum(wakes.ds(N+1:2*N-1))];

% Create initial guess for wake circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaInf = sqrt(2*CT + 1) - 1;
wakes.gamma = repelem([gammaInf;-gammaInf],N+1);
gnew = wakes.gamma;

if strcmpi(opts.Display,'iter') || strcmpi(opts.Display,'final')
    figure;
    hold on; axis image; k = 0;
    for i = 1:numel(foils.m)
        plot(foils.xo(k+[1:foils.m(i) 1]),foils.yo(k+[1:foils.m(i) 1]),'k-');
        k = k + foils.m(i);
    end
    h(1) = plot(wakes.xo(1:N),wakes.yo(1:N),'b-');
    h(2) = plot(wakes.xo(N+1:2*N),wakes.yo(N+1:2*N),'r-');
end

iter = 0;
E = 1;
while (E > opts.FunctionTolerance) && (iter < opts.MaxIterations)
    iter = iter + 1;

    if strcmpi(opts.Display,'iter')
        set(h(1),'XData',wakes.xo(1:N),'YData',wakes.yo(1:N));
        set(h(2),'XData',wakes.xo(N+1:2*N),'YData',wakes.yo(N+1:2*N));
        drawnow;
    end

    % Solve airfoil circulation distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,V] = influence(foils.co,wakes,1);
    A = -U.*sin(foils.theta) + V.*cos(foils.theta);
    gamma = Ainv*(RHS - [A*wakes.gamma; ...
                         -wakes.gamma(1); ...
                         -wakes.gamma(N+2); ...
                         zeros(numel(foils.m)-2,1)]);

    % Calculate induced velocities on wake boundaries %%%%%%%%%%%%%%%%%%%%%%%%
    [U,V] = influence(wakes.co,foils,1);
    u = U*gamma + 1;
    v = V*gamma;
    [U,V] = influence(wakes.co,wakes,0);
    u = u + U*wakes.gamma;
    v = v + V*wakes.gamma;
    Vbar = sqrt(u.*u + v.*v);

    % Update wake shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wakes.dx = u./Vbar.*wakes.ds;
    wakes.dy = v./Vbar.*wakes.ds;
    % wakes.dy([N 2*N]) = 0; % far-field panels remain flat
    wakes.xo(  2:N  ) = wakes.xo(1)   + cumsum(wakes.dx(  1:N-1  ));
    wakes.xo(N+2:2*N) = wakes.xo(N+1) + cumsum(wakes.dx(N+1:2*N-1));
    wakes.yo(  2:N  ) = wakes.yo(1)   + cumsum(wakes.dy(  1:N-1  ));
    wakes.yo(N+2:2*N) = wakes.yo(N+1) + cumsum(wakes.dy(N+1:2*N-1));
    % % Handle wakes crossing
    % if wakes.yo(N) >= wakes.yo(2*N)
    %     yd = wakes.yo(2*N) - wakes.yo(N) - 0.025;
    %     % Rebuild YO of the lower wake directly rather than through CUMSUM
    %     wakes.yo(2:N) = wakes.yo(2:N) + yd* ...
    %         (wakes.xo(2:N)-wakes.xo(1))/opts.WakeLengthChords;
    %     wakes.dy(1:N-1) = diff(wakes.yo(1:N));
    % end
    wakes.theta = atan2(wakes.dy,wakes.dx);
    wakes.co(:,1) = wakes.xo + wakes.dx/2;
    wakes.co(:,2) = wakes.yo + wakes.dy/2;

    % Update wake circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = CT./Vbar;
    g(N+1:2*N) = -g(N+1:2*N);
    % Interpolate circulation values from panel midpoints to endpoints
    gnew(  2:N    ) = g(  1:N-1  ) + diff(g(  1:N  )).* ...
        wakes.ds(  1:N-1  )./(wakes.ds(  1:N-1  )+wakes.ds(  2:N  ));
    gnew(N+3:2*N+1) = g(N+1:2*N-1) + diff(g(N+1:2*N)).* ...
        wakes.ds(N+1:2*N-1)./(wakes.ds(N+1:2*N-1)+wakes.ds(N+2:2*N));
    gnew(1)   = 2*g(1)   - gnew(2);
    gnew(N+2) = 2*g(N+1) - gnew(N+3);
    % Decay circulation to the far-field value
    k = find(wakes.so(  1:N  ) >= wakes.so(N)   - opts.GammaDecayChords, 1);
    gnew(    k:N    ) =  gammaInf + ( gammaInf-gnew(k)    )* ...
        (wakes.so(  k:N  )-wakes.so(N)  )/(wakes.so(N)  -wakes.so(k)  +eps);
    k = find(wakes.so(N+1:2*N) >= wakes.so(2*N) - opts.GammaDecayChords, 1);
    gnew(N+1+k:2*N+1) = -gammaInf + (-gammaInf-gnew(k+N+1))* ...
        (wakes.so(N+k:2*N)-wakes.so(2*N))/(wakes.so(2*N)-wakes.so(k+N)+eps);

    % Calculate residual between current and previous iteration %%%%%%%%%%%%%%
    E = sum(abs(gnew - wakes.gamma))/(2*N*gammaInf);

    % Adjust wake circulation by a relaxation factor %%%%%%%%%%%%%%%%%%%%%%%%%
    wakes.gamma = opts.RelaxationFactor*gnew + ...
        (1 - opts.RelaxationFactor)*wakes.gamma;
end

if strcmpi(opts.Display,'iter') || strcmpi(opts.Display,'final')
    set(h(1),'XData',wakes.xo(1:N),'YData',wakes.yo(1:N));
    set(h(2),'XData',wakes.xo(N+1:2*N),'YData',wakes.yo(N+1:2*N));
    drawnow;
end
end

function [xout,yout,traj] = initWake(foils,gamma)
    idx = [0 0];
    for i = 1:2
        idx = idx(2) + [1 foils.m(i)];
        r = [foils.dx(idx) foils.dy(idx)];
        R = sqrt(r(:,1).^2 + r(:,2).^2);
        c = 0.5*(-r(1,:)/R(1) + r(2,:)/R(2)); % camber direction vector
        y0([i i+2]) = [foils.xo(idx(1)) foils.yo(idx(1))] + ...
            c/sqrt(c(1)*c(1) + c(2)*c(2)) * 0.25*min(R);
    end
    [t,y] = ode45(@objfun,[0 1],y0);
    N = size(y,1);
    % Anchor first point to exactly the trailing edge
    y(1,1:2) = foils.xo([1 foils.m(1)+1]);
    y(1,3:4) = foils.yo([1 foils.m(1)+1]);

    % Trim wakes to have a clean trajectory
    k1 = N; k2 = N;
    prp = [y(k2,2)-y(k1,1) y(k2,4)-y(k1,3)];
    par = [y(k1-1,1)-y(k1,1) y(k1-1,3)-y(k1,3)];
    E = dot(prp,par)/norm(prp)/norm(par);
    while abs(E) > 0.05
        if E > 0; k1 = k1 - 1; else; k2 = k2 - 1; end
        prp = [y(k2,2)-y(k1,1) y(k2,4)-y(k1,3)];
        par = [y(k1-1,1)-y(k1,1) y(k1-1,3)-y(k1,3)];
        E = dot(prp,par)/norm(prp)/norm(par);
    end

    xout = {y(1:k1,1), y(1:k2,2)};
    yout = {y(1:k1,3), y(1:k2,4)};
    traj = [prp(2) -prp(1)]/norm(prp);

    function dydt = objfun(t,y)
        dydt = zeros(4,1);
        [U,V] = influence(reshape(y,[2 2]),foils,1);
        dydt(1:2) = U*gamma + 1;
        dydt(3:4) = V*gamma;
    end
end
