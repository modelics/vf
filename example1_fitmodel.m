%
% Example 1 in [1] - fitting of the model
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

load example1.mat

% To use the same random response in the book
%load example1_asinbook.mat

% Optionally one can add some noise
% H = H + 0.001*(rand(size(H))+1j*rand(size(H)));


% Transfer function entry to display
q = 1;
m = 1;

%% plot poles and samples
figure(2);
scatter(real(poles),imag(poles),40,'bx');
xlabel('Re\{p_n\}');
ylabel('Im\{p_n\}');
grid on;
hold on;

figure(3);
plot(s/1j,squeeze(abs(H(q,m,:))),'b');
xlabel('\omega [rad/s]');
ylabel('Magnitude');
grid on;
hold on;

figure(4)
plot(s/1j,180/pi*squeeze(angle(H(q,m,:))),'b');
xlabel('\omega [rad/s]');
ylabel('Phase [deg]');
grid on;
hold on;

%% Vector Fitting
Order = nbar;
Options = [];
% Options.debug = 1;

Model = FastVF(s/1j,H,Order,Options);


%% Analyze model
% plot ploes
p_model = [Model.pr;Model.pc; conj(Model.pc)];
figure(2);
scatter(real(p_model),imag(p_model),40,'ro');
legend('Exact','Model','Location','NorthWest');

% Compute frequency response of the model
Hmodel = ComputeModelResponse(s/1j,Model.R0,Model.Rr,Model.Rc,Model.pr,Model.pc);

figure(3);
plot(s/1j,abs(squeeze(Hmodel(q,m,:))),'r-.','LineWidth',1.5);
legend('Samples', 'Model');

figure(4);
plot(s/1j,180/pi*angle(squeeze(Hmodel(q,m,:))),'r-.','LineWidth',1.5);
legend('Samples', 'Model');



