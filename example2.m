%
% Example 2 - fitting of the aortic input impedance of a pediatric patient
% [1]
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

% Patient 1 from 
% M. K. Sharp, G. M. Pantalos, L. Minich, L. Y. Tani, E. C. McGough, and
% J. A. Hawkins. Aortic input impedance in infants and children. Journal of
% Applied Physiology, 88(6):2227-2239, 2000.

% Magnitude and phase
Zmag = [3125.90 448.66 340.70 492.55 450.52 906.79 574.12 456.40 570.80 546.01 434.76];
Zphase = pi/180*[0.00 -25.71 5.64 23.14 33.82 -6.50 27.56 14.03 16.24 34.98 25.55];

% Impedance [dyn*s*cm^-5]
Z = Zmag.*exp(1j*Zphase);

% Heart rate [1/minute]
HR = 152.4;

% Period [s]
T = 60/HR;

omega0 = 2*pi/T;

%% Frequency response
omega = (0:10)'*omega0;
Z = reshape(Z,1,1,length(Z));

figure(93);
plot(omega/2/pi,abs(squeeze(Z)),'bo');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dyn\cdots\cdotcm^{-5}]');
grid on
hold on;

figure(94);
plot(omega/2/pi,180/pi*angle(squeeze(Z)),'bo');
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
grid on
hold on;

%% Fitting
Order = 8;
% Because of the noise and poor sampling in the original samples, we need
% to essentially exclude the first convergence test (never satisfied)
Options.PolesEstimationThreshold = 100; 

Model = FastVF(omega,Z,Order,Options);

%% Compute model response over a finer grid
omega_model = linspace(min(omega), max(omega),100);

Z_model = ComputeModelResponse(omega_model,Model.R0,Model.Rr,Model.Rc,Model.pr,Model.pc);

% Plot the frequency response of the model
figure(93);
plot(omega_model/2/pi,abs(squeeze(Z_model)),'r-.','LineWidth',1.5);
legend('Samples', ['Model (order=',num2str(Order),')']);

figure(94);
plot(omega_model/2/pi,180/pi*angle(squeeze(Z_model)),'r-.','LineWidth',1.5);
legend('Samples', ['Model (order=',num2str(Order),')'],'Location','SouthEast');