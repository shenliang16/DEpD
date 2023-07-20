function Omega = ComputeOmega(DscDis, alpha_l, alpha_p)
if isfield(DscDis, 'line')
    if sum(DscDis.line<16)>10
        Omega.lines = 1 - exp(-alpha_l*DscDis.line.^3);
%     Omega.lines =0.8 * (Omega.lines - min(Omega.lines)) / (max(Omega.lines) - min(Omega.lines)) + 0.1;
    else
        Omega.lines = 0.5;
    end
else
    Omega.lines = 0.5;
end
if isfield(DscDis, 'point')
    Omega.points = 1 - exp(-alpha_p*DscDis.point);
%     Omega.points =0.8 * (Omega.points - min(Omega.points)) / (max(Omega.points) - min(Omega.points)) + 0.1;
else
    Omega.points = 0.5;
end