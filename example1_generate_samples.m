%
% Example 1 in [1] - generation of the samples of a rational transfer
% function
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



clear;
close all;

%% Create a test frequency response
kbar = 100;         % number of frequency samples
omega1 = 0.1;         % min frequency
omegamax = 10;        % max frequency
qbar = 1;           % number of outputs
mbar = 1;           % number of inputs

s = 1j*linspace(omega1,omegamax,kbar);
s = s(:);

nbar = 10;                      % number of poles
ncomplpairs = 4;               % number of complex conjugate pairs
nreal = nbar-ncomplpairs*2;     % number of real poles

% Random poles, common to all entries
% Real poles
poles(1:nreal) = -0.3*omegamax*rand(nreal,1);
% Complex conjugate poles
poles(nreal+1:nreal+ncomplpairs) = -0.2*omegamax*rand(ncomplpairs,1) + j*omegamax*rand(ncomplpairs,1);
poles(nreal+ncomplpairs+1:nbar) = conj(poles(nreal+1:nreal+ncomplpairs));

H = zeros(qbar,mbar,kbar);
for q = 1:qbar
    for m = 1:mbar
        % Residues of real poles
        res(1:nreal) = 2*rand(nreal,1)-1;
        % Residues of complex conjugate poles
        res(nreal+1:nreal+ncomplpairs) = rand(ncomplpairs,1) - 1j*rand(ncomplpairs,1);
        res(nreal+ncomplpairs+1:nbar) = conj(res(nreal+1:nreal+ncomplpairs));

        r0= rand;
        H(q,m,:) = r0*ones(1,1,kbar);
        for n = 1:length(poles)
           H(q,m,:) = H(q,m,:) + reshape(res(n)./(s-poles(n)),1,1,kbar); 
        end
    end
end

save example1.mat
