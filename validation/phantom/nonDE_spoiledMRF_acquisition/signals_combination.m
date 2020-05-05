function [mag_signals] = signals_combination(N, FA, RF, TR, TE, T1, T2, M0, idx)
%[mag_signals] = signals_combination(N, FA, RF, TR, TE, T1, T2, M0, w)
%   Calculates magnitude signals of GRE sequence, returning only the
%   signals at the indexes idx

[signals,~] = cppEPG_GRE_legacy(N, FA, RF, TR, T1, T2, M0);

signals = exp(-TE./T2).* signals;
signals = signals(logical(idx));

mag_signals = abs(signals(:));

end

