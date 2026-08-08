clear
close all

set(0,'defaultAxesFontName','Times')
set(0,'defaultAxesFontSize',8.5)
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

% Use full paths when possible to be robust
filename = mfilename('fullpath');
filepath = fileparts( filename );

addpath([filepath '/mesh2d']); initmsh();

% Medium-scale high-geometric-complexity aeropropulsive problem %%%%%%%%%%%%%%
% surfaceFiles = {'onr-dep/mainVec10.dat','onr-dep/nacelleVec10.dat'};
% surfaceFiles = {'dengwirda/mainElement.dat','dengwirda/nacelle.dat', ...
%    'dengwirda/foreFlap.dat','dengwirda/aftFlap.dat','onr-dep/krueger.dat'};

% surfaceFiles = {'minimaljetfoil/lower-open.dat','minimaljetfoil/upper-open.dat'};
surfaceFiles = {'minimaljetfoil/lower-closed.dat','minimaljetfoil/upper-closed.dat'};

for i = numel(surfaceFiles):-1:1
    fid = fopen([filepath '/airfoils/' surfaceFiles{i}],'r');
    surfaces{i} = cell2mat(textscan(fid,'%f%f','Delimiter',{'\t',','}));
    fclose(fid);
end
% Use offset to make main element TE the rotation point
xoffset = max(surfaces{1}(:,1));
for i = 1:numel(surfaces)
    surfaces{i}(:,1) = surfaces{i}(:,1) - xoffset;
end

% Solve
opts.NumPanels = 200; % optionally pass options to the wake solver
opts.ConvergenceCriterion = 1;
opts.Display = 'final';
% [Cp,xc,Cl,Cd,Cm] = panel2d(surfaces,0,1,0.67-xoffset,opts,'Plot','on','Colormap','vik');
[Cp,xc,Cl,Cd,Cm] = panel2d(surfaces,5,0.2,1,0.8-xoffset,opts,'Plot','on','CData','p');
drawnow;

% Plot Cp distributions
% Sweep across multiple ground heights
% opts.Display = 'none';
h = [0.25 0.5 0.75 1 Inf];
colors = lines(length(h));
f = figure;
hold on;
ax = gca;
for i = 1:length(h)
    if isfinite(h(i))
        hstr = sprintf('%.2f',h(i));
    else
        hstr = '\infty';
    end

    [Cp,xc] = panel2d(surfaces,5,h(i),1,0.8-xoffset,opts);

    % Plot nothing to get correct legend entries and next color
    color = colors(i,:);
    p(2*i-1) = plot(ax,NaN,NaN,'--','Color',color,'DisplayName',sprintf('$h/c=%s$ (lower)',hstr));
    p(2*i) = plot(ax,NaN,NaN,'-','Color',color,'DisplayName',sprintf('$h/c=%s$ (upper)',hstr));
    for j = 1:numel(Cp)
        [~,le] = min(xc{j});
        plot(ax,xc{j}(1:le)+xoffset,Cp{j}(1:le),'--','Color',color);
        plot(ax,xc{j}(le:end)+xoffset,Cp{j}(le:end),'-','Color',color);
    end
end
xlim(ax,[-0.1 1.4])
ylim(ax,[-5 3])
xlabel(ax,'$x/c$','Interpreter','latex','FontSize',10,'Units','Points');
ylabel(ax,'$C_p$','Interpreter','latex','FontSize',10,'Units','Points');
set(ax,'YDir','reverse','Box','off','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on',...
 'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','Units','inches');
tmp = get(ax,'Position');
set(ax,'Position',[tmp(1) tmp(2) 5 2.5])
tmp2 = get(ax,'OuterPosition');
set(f,'Units','inches','Position',[tmp(1) tmp(2) tmp2(3) tmp2(4)]);

legend(p,'Interpreter','latex', ...
    'Location','eastoutside', ...
    'FontSize',8,'Units','Points');

exportgraphics(f,'aeropropulsiveIGE-1.pdf')