function out = panel2dwrap(surfaces,alphaDeg,varargin)

    [Cp,xc,Cl,Cd,Cm,foils,wakes] = panel2d(surfaces,alphaDeg,varargin{:});

    out.Cp = Cp;
    out.xc = xc;
    out.Cl = Cl;
    out.Cd = Cd;
    out.Cm = Cm;
    out.foils = foils;
    out.wakes = wakes;

end