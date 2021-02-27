function [signals] = SteadyStateSignals(Sspgr, Sbssfp, nSPGR, nbSSFP, u, p)
%[signals] = SteadyStateSignals(Sspgr, Sbssfp, nSPGR, nbSSFP, u, p)
%   Calculates the steady-state signals of SPGR and bSSFP for several
%   acquisition settings

nSS = nSPGR + nbSSFP;

%p = [T1 T2 M0 P0 B0]
%u = [FA TR TE/RF]
signals = zeros(nSS, 1);

for ii=1:nSS
    if ii<=nSPGR
        aux_u(1) = u(ii);
        aux_u(2) = u(ii + nSS);
        aux_u(3) = u(ii + 2*nSS);
        signals(ii) = Sspgr(aux_u,p);
    else
        aux_u(1) = u(ii);
        aux_u(2) = u(ii + nSS);
        aux_u(3) = u(ii + 2*nSS);
        signals(ii) = Sbssfp(aux_u,p);
    end
end

signals = [real(signals(:)); imag(signals(:))];

end

