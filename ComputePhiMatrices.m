function [Phir,Phic] = ComputePhiMatrices(omega,pr,pc)
%
% Function to compute the frequency response of a Vector Fitting model [1]
%
% Usage:
%   [Phir,Phic] = ComputePhiMatrices(omega,pr,pc)
%
% Input arguments:
%  - omega: frequency samples, column vector. This is angular frequency (omega = 2*pi*f)
%  - pr: real poles
%  - pc: complex conjugate poles (only one pole per pair)
%
% Output arguments:
%  - Phir: coefficient matrix, part associated to real poles
%  - Phic: coefficient matrix, part associated to complex poles
%
% Licensing condition: 
% you can freely use these codes (the "Software") subject to the conditions
% in the LICENSE file. Note that you must cite the following book chapter in the
% publications and product documentation arising from the use of this Software
% 
%  [1] P. Triverio, "Vector Fitting", in P. Benner, S. Grivet-Talocia, A.
%  Quarteroni, G. Rozza, W. H. A. Schilders, L. M. Silveira (Eds.),
%  "Handbook on Model Order Reduction", De Gruyter (to appear).
% 
% Copyright 2019 Piero Triverio, www.modelics.org


    kbar = length(omega);
    nr = length(pr);
    nc = length(pc);

    % Compute the new Phir and Phic matrices
    Phir = zeros(kbar,nr);
    for ii = 1:nr
        Phir(:,ii) = 1./(1j*omega-pr(ii));
    end
    Phic = zeros(kbar,2*nc);
    for ii = 1:nc
        Phic(:,2*ii-1) = 1./(1j*omega-pc(ii)) + 1./(1j*omega-conj(pc(ii)));
        Phic(:,2*ii) = 1j./(1j*omega-pc(ii)) - 1j./(1j*omega-conj(pc(ii)));       
    end
return