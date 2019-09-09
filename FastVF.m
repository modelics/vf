function Model = FastVF(omega,H,Order,Options)
%
% Fast Vector Fitting, real-valued implementation from [1]
%
% Usage:
%   Model = FastVF(omega,H,Order)
%   Model = FastVF(omega,H,Order,Options)
%
% Inputs:
%  - omega: frequency samples, column vector. This is angular frequency (omega = 2*pi*f)
%  - H: response samples, 3D array. First dimension corresponds to system outputs. Second dimension to system inputs. Third dimension corresponds to frequency. 
%  - Order: desired order of the VF model
%  
% Options (with defaults):
%    - Options.PolesEstimationThreshold: treshold for the first convergence
%    test to stop the iterative estimation of poles (default: 1e-1)
%    - Options.ModelErrorThreshold: treshold for the second convergence test; when modeling error falls below this threshold,
%       the model is returned (default: 1e-3)
%    - Options.MaxIterations: maximum number of iterations. When reached the 
%      algorithm stops even if target accuracy was not reached
%      (default: 5)
%    - Options.EnforceStability: flip model poles into the stable left half plane to enforce model stability 
%      (default: 1 (enabled) )
%    - Options.Debug: print and plot extra information. Pause at each
%    iteration (default: 0 (disabled) )
%
% Outputs:
%    - Model.pr: real poles of the model, column vector
%    - Model.pc: complex poles of the model, one per conjugate pair, column
%    vector
%    - Model.R0: constant term R0 of the model, 2D array
%    - Model.Rr: residues of real poles, 3D array. First and second dimension
%    correspond to outputs and inputs. Third dimension to poles
%    - Model.Rc: residues of complex poles, 3D array. First and second dimension
%    correspond to outputs and inputs. Third dimension to poles. One
%    residue per pair
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



%% Options defaults

% Convergence threshold (poles estimation)
if ~exist('Options','var') || ~isfield(Options,'PolesEstimationThreshold')
    Options.PolesEstimationThreshold = 1e-1;
end

% Convergence threshold (model-samples error)
if ~exist('Options','var') || ~isfield(Options,'ModelErrorThreshold')
    Options.ModelErrorThreshold = 1e-3;
end

% Maximum number of iterations
if ~exist('Options','var') || ~isfield(Options,'MaxIterations')
    Options.MaxIterations = 5;
end

% Enforce stability
if ~exist('Options','var') || ~isfield(Options,'EnforceStability')
    Options.EnforceStability = 1;
end

% Debug mode
if ~exist('Options','var') || ~isfield(Options,'debug')
    Options.debug = 0;
end


%% Initializations

tic;

kbar = length(omega);   % number of samples
qbar = size(H,1);       % number of outputs
mbar = size(H,2);       % number of inputs
nbar = Order;           % desired model order

omega = omega(:);       % Make sure omega is a column vector

Model = [];

% ---- Initial poles
% Number of real poles
nr = mod(nbar,2); % no real poles if nbar is even, one if order is odd
% Number of complex conjugate pairs
nc = (nbar-mod(nbar,2))/2;                              

alpha = 0.01;

% Real pole
if nr == 1 % nbar is odd
    pr = -alpha*max(omega);
else
    pr = [];    
end

% Complex poles
if nc == 1
    pc = (-alpha+1j)*max(omega)/2;
else
    if min(omega) == 0
        pc = (-alpha+1j)*max(omega)*(1:nc)/nc;
    else
        pc = (-alpha+1j)*(min(omega) + (max(omega)-min(omega))/(nc-1)*(0:nc-1));
    end
    pc = pc(:);
end

%% Poles estimation

% Iterative process
iter = 1;
while iter <= Options.MaxIterations
    fprintf('\n -------------- Iteration n.%d --------------\n',iter);
    
    %% QR decompositions
    % Compute the Phir and Phic matrices
    [Phir,Phic] = ComputePhiMatrices(omega,pr,pc);    
    
    % Compute Phi0 and Phi1
%   Phi0 = [ones(kbar,1),Phir, Phic];
%   Phi1 = [Phir, Phic];
     
    % Preallocate matrix used in QR decompositions
    M = zeros(2*kbar,2*nbar+1);
    % Preallocate matrix and right hand side for least squares problem (54)
    Alsq = zeros(nbar*qbar*mbar,nbar);
    blsq = zeros(nbar*qbar*mbar,1);
    
    % Compute the first columns of M, which do not depend on q and m
    M(1:kbar,1) = ones(kbar,1);
    M(1:kbar,2:nr+1) = real(Phir);
    M(1:kbar,nr+2:nbar+1) = real(Phic);
    M(kbar+1:end,2:nr+1) = imag(Phir);
    M(kbar+1:end,nr+2:nbar+1) = imag(Phic);
    
    irow = 0;
    for q = 1:qbar %loop over outputs
        for m = 1:mbar %loop over inputs
            V_Hqm = squeeze(H(q,m,:));
            D_Hqm = spdiags(V_Hqm,0,kbar,kbar);
            % Compute the second part of M
            M(1:kbar,nbar+2:nbar+nr+1) = -real(D_Hqm*Phir);
            M(1:kbar,nbar+nr+2:end) = -real(D_Hqm*Phic);
            M(kbar+1:end,nbar+2:nbar+nr+1) = -imag(D_Hqm*Phir);
            M(kbar+1:end,nbar+nr+2:end) = -imag(D_Hqm*Phic);
            % QR decomposition (53)
            [Qqm,Rqm] = qr(M,0);
            
            % Assemble the coefficient matrix and right hand side of the
            % least squares problem for poles estimation
            Alsq(irow+1:irow+nbar,:) = Rqm(nbar+2:end,nbar+2:end);
            blsq(irow+1:irow+nbar) = Qqm(1:kbar,nbar+2:end)'*real(V_Hqm) + Qqm(kbar+1:end,nbar+2:end)'*imag(V_Hqm);
        end
    end
    
    % Least squares problem (54)
    cw = Alsq\blsq;

    % Evaluate the weigthing function before we compute the new poles estimate
    % This will be used by the first convergence test
    w = 1 + Phir*cw(1:nr) + Phic*cw(nr+1:end);

    % Plot the magnitude of the weighting function
    if Options.debug
        figure(102);
        plot(omega,abs(w));
        xlabel('Omega');
        ylabel('|w(j\omega|');
        title('Weigthing function magnitude');
        grid on
    end

    % Compute the new poles estimate
    A = zeros(nbar,nbar);
    bw = ones(nbar,1);
    for ii = 1:nr
        A(ii,ii) = pr(ii);
        % bw(ii) = 1;
    end
    for ii = 1:nc
        A(nr+2*ii-1:nr+2*ii,nr+2*ii-1:nr+2*ii) = [real(pc(ii)), imag(pc(ii)); -imag(pc(ii)),real(pc(ii))];
        bw(nr+2*ii-1:nr+2*ii,1) = [2;0];
    end
    pnew = eig(A-bw*cw.');

    % Plot the new poles estimate
    if Options.debug
       figure(101); 
       scatter(real(pnew),imag(pnew),40,'ro');
       xlabel('Re\{p_n\}');
       ylabel('Im\{p_n\}');
       title(sprintf('Poles estimate, iteration %d',iter));
       grid on;
    end
    
    % Extract real poles
    ind_rp = find(abs(imag(pnew))<10*eps*abs(pnew));
    pr = real(pnew(ind_rp));
    nr = length(ind_rp);
    
    % Extract complex conjugate pairs of poles
    % Find only the poles with positive imaginary part
    ind_cp = find(imag(pnew)>=10*eps*abs(pnew));
    pc = pnew(ind_cp);
    nc = length(ind_cp);
    
    %% Stability/causality enforcement
    if Options.EnforceStability
        pr = -abs(pr);
        pc = -abs(real(pc)) + 1j*imag(pc);
    end
        
    %% First convergence test
    
    w_minus_one = 1/sqrt(kbar)*norm(abs(w(:)-1));
    % Do a tentative model fitting if either:
    % - the first convergence test is successful 
    % - Options.debug is enabled
    if w_minus_one <= Options.PolesEstimationThreshold
        fprintf('Convergence test (poles estimation): \t\tpassed (%e)\n',w_minus_one);
    else
        fprintf('Convergence test (poles estimation): \t\tfailed (%e)\n',w_minus_one);
    end
    if  w_minus_one <= Options.PolesEstimationThreshold || Options.debug
        %% Tentative final fitting
        [Phir,Phic] = ComputePhiMatrices(omega,pr,pc);
        
        % Compute the matrix of the least squares problem
        Alsq = zeros(2*kbar,nbar+1);
        Alsq(1:kbar,1) = ones(kbar,1);
        Alsq(1:kbar,2:nr+1) = real(Phir);
        Alsq(1:kbar,nr+2:nbar+1) = real(Phic);
        Alsq(kbar+1:end,2:nr+1) = imag(Phir);
        Alsq(kbar+1:end,nr+2:nbar+1) = imag(Phic); 
        
        % Store model coefficients in the output structure Model, in case it
        % will found accurate enough
        Model.pr = pr;
        Model.pc = pc;
        Model.R0 = zeros(qbar,mbar);
        Model.Rr = zeros(qbar,mbar,nr);
        Model.Rc = zeros(qbar,mbar,nc);
        
        % Model-samples error
        err = 0;
        for q = 1:qbar
            for m = 1:mbar
                % Right hand side
                V_Hqm = squeeze(H(q,m,:));
                blsq = [real(V_Hqm);imag(V_Hqm)];
                c_Hqm = Alsq\blsq;
                Model.R0(q,m) = c_Hqm(1);
                Model.Rr(q,m,:) = c_Hqm(2:nr+1);
                Model.Rc(q,m,:) = c_Hqm(nr+2:2:end) + 1j*c_Hqm(nr+3:2:end);
                
                % Plot the given samples vs the model response for the
                % (1,1) entry of the transfer function (if in debug mode)
                if Options.debug && q == 1 && m == 1
                    % Compute model response
                    Htemp = ComputeModelResponse(omega,Model.R0(q,m),Model.Rr(q,m,:),Model.Rc(q,m,:),Model.pr,Model.pc);
                    figure(103);
                        subplot(211);
                        plot(omega,abs(squeeze(H(q,m,:))),'bx');
                        hold on;
                        plot(omega,abs(squeeze(Htemp(q,m,:))),'r-.'); 
                        grid on;
                        xlabel('Omega');
                        ylabel('Magnitude');
                        legend('Samples H_k','Model');
                        hold off; 

                        subplot(212);
                        plot(omega,180/pi*angle(squeeze(H(q,m,:))),'bx');
                        hold on;
                        plot(omega,180/pi*angle(squeeze(Htemp(q,m,:))),'r-.');
                        grid on;
                        xlabel('Omega');
                        ylabel('Phase [deg]');
                        legend('Samples H_k','Model');
                        hold off;                          
                end
                
                err = err+sum(abs(Alsq*c_Hqm-blsq).^2);
            end
        end
        err = sqrt(err)/sqrt(qbar*mbar*kbar);
        
        if err <= Options.ModelErrorThreshold
            fprintf('Convergence test (model-samples error): \tpassed (%e)\n',err);
            fprintf('Model identification successful\n'); 
            fprintf('Modeling time: %f s\n',toc);
            return
        else
            fprintf('Convergence test (model-samples error): \tfailed (%e)\n',err);
        end
    end
    
    if Options.debug
        pause;
    end
    
    iter = iter + 1;
    
end

fprintf('Warning: could not reach the desired modeling error within the allowed number of iterations\n');
fprintf('Modeling time: %f s\n',toc);

return
