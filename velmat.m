function [U,V] = velmat(co,af,part,h)
%VELMAT Compute velocity influence matrices
if nargin == 3
    h = Inf;
end

[U,V] = influence(co,af,part);

if isfinite(h)
    % Compute influence of mirror image if one exists
    [Uimg,Vimg] = influence([co(:,1) -co(:,2)],af,part);
    U = U + Uimg;
    V = V - Vimg;
end