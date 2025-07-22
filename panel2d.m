function [Cp,xc,yc,Cl,Cd,Cm,visout,varargout] = panel2d(surfaces,alphaDeg,varargin)
% PANEL2D  Panel method in two dimensions.
%   PANEL2D(SURFACES,ALPHADEG) runs a standard panel method.
%   PANEL2D(SURFACES,ALPHADEG,CT,XDISK) runs the APPAC aeropropulsive analysis
%   PANEL2D(SURFACES,ALPHADEG,CT,XDISK,WAKEOPTIONS) to pass options to APPAC
%   PANEL2D(___,NAME,VALUE) to pass Name-Value arguments with any above syntax
%
%   See also INFLUENCE, SOLVEWAKE, FLOWVIS.
[oper,CT,xDisk,wakeOptions,options] = parseInput(varargin);

nSurfs = numel(surfaces);

if ((oper == 2) && (nSurfs == 1)) || (nSurfs == 0)
    error('panel2d:invalidInput','Not enough surfaces.');
end

R = [cosd(alphaDeg) -sind(alphaDeg);sind(alphaDeg) cosd(alphaDeg)];

% Determine total number of surface panels for memory allocation %%%%%%%%%%%%%
for i = nSurfs:-1:1 % populating the last entry allocates necessary space
    foils.m(i) = size(surfaces{i},1) - 1;
end
M = sum([foils.m]); % total number of panels

% Build unified struct for all lifting surfaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = M;
for i = nSurfs:-1:1
    k = k - foils.m(i); % shift indices to rows in array for writing
    coords = surfaces{i}*R; % rotate coordinates by AoA
    if det(surfaces{i}([2 end-1],:)-surfaces{i}(1,:)) > 0
        coords = flipud(coords); % method requires CW defined coordinates
    end
    foils.xo(1+k:foils.m(i)+k,:) = coords(1:end-1,1);
    foils.yo(1+k:foils.m(i)+k,:) = coords(1:end-1,2);
    foils.dx(1+k:foils.m(i)+k,:) = diff(coords(:,1));
    foils.dy(1+k:foils.m(i)+k,:) = diff(coords(:,2));
end
foils.theta = atan2(foils.dy,foils.dx);
foils.co = [foils.xo+foils.dx/2 foils.yo+foils.dy/2];

% Create aerodynamic influence coefficient matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(M+nSurfs);
for i = 1:nSurfs % Kutta condition
    A(M+i,k+1) = 1;
    A(M+i,k+foils.m(i)+1) = 1;
    k = k + foils.m(i) + 1;
end
[U,V] = influence(foils.co,foils,1);
A(1:M,:) = -U.*sin(foils.theta) + V.*cos(foils.theta); % normal component
B = U.*cos(foils.theta) + V.*sin(foils.theta); % tangent component
RHS = [sin(foils.theta);zeros(nSurfs,1)];

% Solve for global circulation solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if oper == 1
    foils.gamma = A \ RHS;
    Qtan = B*foils.gamma + cos(foils.theta);
    out{1} = foils;
else
    [wakes,foils.gamma,~,~] = solveWake(foils,inv(A),RHS,CT,wakeOptions);
    [U,V] = influence(foils.co,wakes,1);
    D = U.*cos(foils.theta) + V.*sin(foils.theta);
    Qtan = B*foils.gamma + cos(foils.theta) + D*wakes.gamma;
    out = {foils,wakes};
end

% Calculate coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cparr = 1 - Qtan.^2;

xc = mat2cell(foils.co*R(1,:).', foils.m);
yc = mat2cell(foils.co*R(2,:).', foils.m);

if oper == 2
    % Correct Cp aft of the actuator disk where the total pressure is higher
    k1 = find(xc{1} < xDisk, 1, 'last');
    k2 = find(xc{2} < xDisk, 1, 'first');
    idx = getCpCorrectionIds(foils,wakes,k1,k2);
    Cparr(idx) = Cparr(idx) + 2*CT;
end

Cl = -Cparr.'*foils.dx;
Cd =  Cparr.'*foils.dy;
Cm = (Cparr.*foils.dx).'*foils.co(:,1) + (Cparr.*foils.dy).'*foils.co(:,2);
Cm25 = Cm + 0.25*Cl;

Cp = mat2cell(Cparr, foils.m);

% Data visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
visout = -1;
if strcmpi(options.Plot,'on')
    if oper == 1; visout = flowVis(options,foils); end
    if oper == 2; visout = flowVis(options,foils,wakes,k1,k2); end
end

% Print integrated values at the very end
% fprintf(1,'%+4s: %8.5f\n','Cl',Cl,'Cd',Cd,'Cm',Cm,'Cm25',Cm25);

nout = max(nargout,1) - 7;
for i = 1:nout
    varargout{i} = out{i};
end
end



% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [oper,CT,xDisk,wakeOptions,options] = parseInput(args)
% PARSEINPUT  Handle input arguments from the various input syntaxes.
oper = 1; % operating mode (vanilla:1, appac:2)
CT = 0;   % values of CT and xDisk are inconsequential when oper=1
xDisk = 0;

% Default options
wakeOptions.MaxIterations = 50;
wakeOptions.FunctionTolerance = 1e-6;
wakeOptions.RelaxationFactor = 0.5;
wakeOptions.NumPanels = 100;
wakeOptions.WakeLengthChords = 9;
wakeOptions.GammaDecayChords = 3;
wakeOptions.NodeSpacing = 'cosine';
wakeOptions.Display = 'iter';

options.Colormap = 'berlin';
options.Mesh = 'off';
options.Plot = 'off';
options.CData = 'q';

% Perform checks on input arguments and parse accordingly %%%%%%%%%%%%%%%%%%%%
nArgs = numel(args);

if nArgs == 0
    return
end

if nArgs == 1
    error('panel2d:invalidInputSyntax','Invalid input syntax.');
end

if mod(nArgs,2) == 0
    if ~(isnumeric(args{1}) && isscalar(args{1}) && ...
         isnumeric(args{2}) && isscalar(args{2}))&& ...
        ~ischar   (args{1})
        error('panel2d:incorrectInputClass','Incorrect input class.');
    end
    NameValueStartIdx = isnumeric(args{1})*2 + 1;
else
    if (~isnumeric(args{1}) || ~isscalar(args{1}) || ...
        ~isnumeric(args{2}) || ~isscalar(args{2}) || ...
        ~isstruct (args{3}) )
        error('panel2d:incorrectInputClass','Incorrect input class.');
    end
    NameValueStartIdx = 4;
    wakeOptions = updateOptions(wakeOptions, ...
                    [fieldnames(args{3}) struct2cell(args{3})].');
end
if NameValueStartIdx > 2
    oper = 2;
    CT = args{1};
    xDisk = args{2};
end
options = updateOptions(options,args(NameValueStartIdx:nArgs));
end

function options = updateOptions(options,NameValuePairs)
% UPDATEOPTIONS  Update options using pairs of Name-Value arguments.
names = fieldnames(options);
for i = 1:2:numel(NameValuePairs)-1
    k = strcmpi(names,NameValuePairs{i});
    if ~any(k)
        error('panel2d:invalidOptions', ...
            '%s is not a recognized option.',NameValuePairs{i});
    end
    options.(names{k}) = NameValuePairs{i+1};
end
end

function in = getCpCorrectionIds(foils,wakes,k1,k2)
% GETCPCORRECTIONIDS  Get indices of control points that need Cp correction
%   The method of creating the bounding polygon for the domain of increased
%   total pressure is taken from FLOWVIS.
bb = [min(foils.co),max(foils.co)]; % bounding box
bb = bb + [-0.7 -0.7 0.7 0.5];

N = wakes.m(1);
% Distinguish jet and outer flow eval by making one wake left-running
wakes.xo(1:N) = wakes.xo(1:N) + wakes.dx(1:N);
wakes.yo(1:N) = wakes.yo(1:N) + wakes.dy(1:N);
wakes.dx(1:N) = -wakes.dx(1:N);
wakes.dy(1:N) = -wakes.dy(1:N);
wakes.theta(1:N) = wakes.theta(1:N)+pi - ceil(wakes.theta(1:N)/2/pi)*2*pi;

%--------------------------------------- Mesh jet volume
k3 = find(wakes.co(  1:N  ,1)>=bb(3),1);
k4 = find(wakes.co(N+1:2*N,1)>=bb(3),1) + N;
node = [flipud(wakes.co(N+1:k4,:)); ...
        foils.co(foils.m(1)+(1:k2-1),:); ...
        foils.co(k1+1:foils.m(1),:); ...
        wakes.co(1:k3,:)];

edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

in = inpoly2(foils.co,node,edge);
end