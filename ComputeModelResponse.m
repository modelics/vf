function H = ComputeModelResponse(omega,R0,Rr,Rc,pr,pc)
%
% Compute the frequency response of a Vector Fitting model [1]
%
% Usage:
%   H = ComputeModelResponse(omega,R0,Rr,Rc,pr,pc)
%
% Input arguments:
%  - omega: frequency samples, column vector. This is angular frequency (omega = 2*pi*f)
%  - R0: constant coefficient
%  - Rr: residues of real poles, 3D array. First dimension corresponds to system outputs. Second dimension to system inputs. Third dimension corresponds to the various poles. 
%  - Rc: residues of complex conjugate pole pairs (only one per pair)
%  - pr: real poles, column vector
%  - pc: complex poles, column vector. Only one per pair of conjugate poles
%
% Output arguments:
%  - H: model response samples, 3D array. First dimension corresponds to system outputs. Second dimension to system inputs. Third dimension corresponds to frequency. 
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

qbar = size(R0,1);      % number of outputs
mbar = size(R0,2);      % number of inputs
kbar = length(omega);   % number of frequency points

nr = length(pr);        % number of real poles
nc = length(pc);        % number of complex conjugate pairs

% Preallocate space for H
H = zeros(qbar,mbar,kbar);

% This part should be vectorized for higher efficiency
for ik = 1:kbar
    H(:,:,ik) = R0;
    for ir = 1:nr
        H(:,:,ik) = H(:,:,ik) + Rr(:,:,ir)/(1j*omega(ik) - pr(ir));
    end
    for ic = 1:nc
        H(:,:,ik) = H(:,:,ik) + Rc(:,:,ic)/(1j*omega(ik) - pc(ic));
        H(:,:,ik) = H(:,:,ik) + conj(Rc(:,:,ic))/(1j*omega(ik) - conj(pc(ic)));        
    end
end


return