function visout = flowVis(options,foils,varargin)
% FLOWVIS  Visualize flow around surfaces with bound vorticity.
%   FLOWVIS(OPTIONS,FOILS,WAKES,K1,K2) creates a contour plot with overlaid
%   streamlines using the style specified in OPTIONS. At minimum, the unified
%   struct FOILS containing the solid surfaces must be given. If there is a
%   powered wake, all of WAKES, K1, and K2 must be given where WAKES is the
%   unified struct for the wake boundaries and K1 and K2 are indices specify-
%   ing the nodes (one on each surface) where the actuator disk is located.

% Meshing settings (not user defined)
opts.kind = 'delfront';
opts.rho2 = 1;

bb = [min(foils.co),max(foils.co)]; % bounding box
bb = bb + [-0.7 -0.7 0.7 0.5];

if nargin == 5
    wakes = varargin{1};
    k1 = varargin{2};
    k2 = varargin{3};

    N = wakes.m(1);
    % Distinguish jet and outer flow eval by making one wake left-running
    wakes.xo(1:N) = wakes.xo(1:N) + wakes.dx(1:N);
    wakes.yo(1:N) = wakes.yo(1:N) + wakes.dy(1:N);
    wakes.dx(1:N) = -wakes.dx(1:N);
    wakes.dy(1:N) = -wakes.dy(1:N);
    wakes.theta(1:N) = wakes.theta(1:N)+pi - ceil(wakes.theta(1:N)/2/pi)*2*pi;
%    wakes.theta(1:N) = atan2(wakes.dy(1:N),wakes.dx(1:N));

    %--------------------------------------- Mesh jet volume
    k3 = find(wakes.co(  1:N  ,1)>=bb(3),1);
    k4 = find(wakes.co(N+1:2*N,1)>=bb(3),1) + N;
    node = [flipud(wakes.co(N+1:k4,:)); ...
            foils.co(foils.m(1)+(1:k2-1),:); ...
            foils.co(k1+1:foils.m(1),:); ...
            wakes.co(1:k3,:)];

    edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

    [vert{1},~,tria{1},~] = refine2(node,edge,[],opts);

    [U,V] = influence(vert{1},foils,1);
    u{1} = U*foils.gamma + 1;
    v{1} = V*foils.gamma;
    [U,V] = influence(vert{1},wakes,-1); % -1 for jet interior
    u{1} = u{1} + U*wakes.gamma;
    v{1} = v{1} + V*wakes.gamma;

    %------------------------------------- Mesh outer volume
    node = [flipud(wakes.co(1:k3,:)); ...
            foils.co(1:k1+1,:); ...
            foils.co(foils.m(1)+(k2-1:foils.m(2)),:); ...
            wakes.co(N+1:k4,:); ...
            wakes.co(k4,1) bb(4); ...
            bb(1) bb(4);bb(1) bb(2); ...
            wakes.co(k3,1) bb(2)];

    edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

    % if numel(foils.m) > 2
    %     j = cumsum(foils.m(3:end));
    %     edge = [edge;(1:j(end)).'+[0 1]+size(node,1)];
    %     edge(size(node,1)+j,2) = edge(size(node,1)+j,2) - foils.m(3:end).';
    %     node = [node;foils.co(foils.m(1)+foils.m(2)+1:end,:)];
    % end

    [vert{2},~,tria{2},~] = refine2(node,edge,[],opts);

    [U,V] = influence(vert{2},foils,1);
    u{2} = U*foils.gamma + 1;
    v{2} = V*foils.gamma;
    [U,V] = influence(vert{2},wakes,1);
    u{2} = u{2} + U*wakes.gamma;
    v{2} = v{2} + V*wakes.gamma;

    %---------------------------------- Merge the two meshes
    TRI = [tria{1};tria{2}+size(vert{1},1)];
    VTX = [vert{1};vert{2}];
    data.u = [u{1};u{2}]; data.v = [v{1};v{2}];
else
    j = cumsum(foils.m);
    node = [foils.co;4 2;-2 2;-2 -2;4 -2];
    edge = (1:j(end)+4).' + [0 1]; edge(j,2) = edge(j,2) - foils.m.';
    edge(end,2) = j(end) + 1;
    [VTX,~,TRI,~] = refine2(node,edge,[],opts);
    [U,V] = influence(VTX,foils,1);
    data.u = U*foils.gamma + 1;
    data.v = V*foils.gamma;
end

data.q = sqrt(data.u.^2 + data.v.^2);
data.p = 1 - data.q.^2;
if nargin == 5
    CT = 0.5*((abs(wakes.gamma(N+1)) + 1)^2 - 1);
    data.p(1:numel(u{1})) = data.p(1:numel(u{1})) + 2*CT;
end

% Present %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
cmap = crameri(options.Colormap);
colormap(cmap);
hold on; axis image off;

% Choose the line color to contrast with the colormap
LineColor = [0.8 0.8 0.8] - ...
    0.6*round(cmap(ceil(size(cmap,1)/2),:)*[0.299;0.587;0.114]);

if strcmpi(options.CData,'p') || strcmpi(options.CData,'v')
    pivot = 0;
else
    pivot = 1;
end

if strcmpi(options.Mesh,'on')
    EdgeColor = [0.5 0.5 0.5]; % slightly off-color to streamlines
else
    EdgeColor = 'none';
end

% Create the contour plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patch('Faces',TRI,'Vertices',VTX, ...
    'FaceVertexCData',data.(lower(options.CData)), ...
    'FaceColor','interp','EdgeColor',EdgeColor);
cl = get(gca,'CLim');
set(gca,'CLim',max(abs(cl-pivot))*[-1 1]+pivot);

% Draw streamlines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FlowP = tristream(TRI,VTX(:,1),VTX(:,2),data.u,data.v, ...
    zeros(1,51)+bb(1)+0.01,linspace(bb(2)+0.01,bb(4)-0.01,51));
for i = 1:numel(FlowP)
    plot(FlowP(i).x,FlowP(i).y,'-','Color',LineColor);
end

% Configure window
axis([bb(1)+0.25 bb(3)-0.25 bb(2)+0.4 bb(4)-0.2]);
pb = get(gca,'PlotBoxAspectRatio');
pos = get(fig,'Position'); pos(3) = pos(4)*pb(1);
set(fig,'Position',pos);
set(gca,'Position',[0 0 1 1]);

%print(fig,'-r300','multielement-1','-dpng');

visout.fig = fig;
visout.VTX = VTX;
visout.TRI = TRI;
visout.data = data;