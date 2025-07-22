function out = panel2dwrap(surfaces,alphaDeg,varargin)

    [Cp,xc,yc,Cl,Cd,Cm,visout,foils,wakes] = panel2d(surfaces,alphaDeg,varargin{:});

    out.Cp = Cp;
    out.xc = xc;
    out.yc = yc;
    out.Cl = Cl;
    out.Cd = Cd;
    out.Cm = Cm;
    out.foils = foils;
    out.wakes = wakes;
    out.visout = visout;

end