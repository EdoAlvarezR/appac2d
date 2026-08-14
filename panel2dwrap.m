function out = panel2dwrap(surfaces,alphaDeg,h,varargin)

    [Cp,xc,yc,Cl,Cd,Cm,visout,foils,wakes] = panel2d(surfaces,alphaDeg,h,varargin{:});

    out.Cp = Cp;
    out.xc = xc;
    out.yc = yc;
    out.Cl = Cl;
    out.Cd = Cd;
    out.Cm = Cm;
    out.foils = foils;
    out.wakes = wakes;
    out.visout = visout;
    out.alphaDeg = alphaDeg;
    out.h = h;

end