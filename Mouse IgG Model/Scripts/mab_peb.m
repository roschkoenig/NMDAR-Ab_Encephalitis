%% mab_peb
%==========================================================================
% This routine loads the inverted individual time window DCMs and performs
% parametric empirical Bayesian analyiss

clear all

% Housekeeping
%==========================================================================
D           = mab_housekeeping;
fs          = filesep;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fdcm        = D.Fdcm;
Fanalysis   = D.Fanalysis;

mind        = D.mind;
Tplot       = D.Tplot;      Tlabel      = D.Tlabel;
Gplot       = D.Gplot;      Glabel      = D.Glabel;

load([Fanalysis fs 'DCM_All.mat']);
LFP = spm_eeg_load([Fanalysis fs 'LFP_MEEG']);
clear D

% PEB Analysis
%==========================================================================
% Make second level model space
%--------------------------------------------------------------------------
clear M ptz atb itx
k = 1/10;           % Inverse time constant
H = 1 * 1/0.37;     % Maximum height
i       = 0;
ptz_t   = [];

% PTZ time curves
%--------------------------------------------------------------------------
seg     = 1:fix(length(SLIDE)/4);
tim_ax  = linspace(0, 60, length(seg));

for w = tim_ax
    i = i+1;
    ptz_t(i) = H*k*w * exp(-k*w);
end

ptz(seg) = zeros(1,length(seg));
ptz(seg + seg(end)) = ptz_t;
ptz     = [ptz, ptz];

% Antibody definitions
%--------------------------------------------------------------------------
atb     = zeros(1, 2*seg(end));
atb_2   = ones(1, 2*seg(end));
atb     = [atb, atb_2];

% Interaction modelling 
%--------------------------------------------------------------------------
for i = 1:length(atb), 
    if ptz(i) == 0,
        itx(i) = 0;
    else
        if atb(i) > 0
            itx(i) = ptz(i);
        else
            itx(i) = -ptz(i);
        end
    end
end

clist   = condlist(LFP);
for c = 1:length(clist)
    uscores = find(clist{c} == '_');
    for u = uscores, clist{c}(u) = ' '; end
end

subplot(3,1,1), plot(atb); title('Main effect of Antibody');
subplot(3,1,2), plot(ptz); title('Main effect of PTZ');
subplot(3,1,3), plot(itx); title('Interaction');
set(gca, 'XTick', 1:seg(end):length(SLIDE))
set(gca, 'XTickLabel', clist(1:seg(end):length(SLIDE)));

% Run PEB
%--------------------------------------------------------------------------
FCM         = SLIDE';
X           = [atb; ptz; itx; ones(1,length(atb)); ]';
Xnames      = {'Antibody', 'PTZ', 'Interaction', 'Static'};

M.X         = X;
M.Xnames    = Xnames;
M.Q         = 'all';

% The model space below yields the full model as winning
%==========================================================================
% Time constant parameters mab_spm_fx_cmc
%--------------------------------------------------------------------------
% G(:,1)  ss -> ss (-ve self)  4    MOD
% G(:,2)  sp -> ss (-ve rec )  4    INH
% G(:,3)  ii -> ss (-ve rec )  4    INH
% G(:,4)  ii -> ii (-ve self)  4    MOD
% G(:,5)  ss -> ii (+ve rec )  4    EXC
% G(:,6)  dp -> ii (+ve rec )  2    EXC
% G(:,7)  sp -> sp (-ve self)  4    MOD
% G(:,8)  ss -> sp (+ve rec )  4    EXC
% G(:,9)  ii -> dp (-ve rec )  2    INH
% G(:,10) dp -> dp (-ve self)  1    MOD
% 
% The order in the DCM structure is as follows
% j     = [7 2 3 4 5 6 8 9 10 1];
% i.e.     M I I M E E E I M  M
% new   =  1 2 3 4 5 6 7 8 9 10

fields{1}   = {'T(1)', 'T(2)', 'T(3)', 'T(4)'};     % time constants
fields{2}   = {'G(2)', 'G(3)', 'G(8)'};             % inh connections
fields{3}   = {'G(5)', 'G(6)', 'G(7)'};             % exc connections
fields{4}   = {'G(1)', 'G(4)', 'G(9)', 'G(10)'};    % mod connections
fields{5}   = {fields{1}{:}, fields{2}{:}};         % time and inh
fields{6}   = {fields{1}{:}, fields{3}{:}};         % time and exc
fields{7}   = {fields{1}{:}, fields{4}{:}};         % time and mod
fields{8}   = {fields{1}{:}, fields{2}{:}, fields{4}{:}};   % time, inh, and mod
fields{9}   = {fields{1}{:}, fields{2}{:}, fields{3}{:}};   % time, inh, and exc
fields{10}   = {fields{2}{:}, fields{3}{:}};         % inh and exc
fields{11}   = {fields{2}{:}, fields{4}{:}};         % inh and mod
fields{12}  = {fields{3}{:}, fields{4}{:}};         % exc and mod
fields{13}  = {fields{2}{:}, fields{3}{:}, fields{4}{:}}; % all coupling
fields{14}  = {fields{1}{:}, fields{2}{:}, fields{3}{:}, fields{4}{:}}; % all 

labels  = { 't', 'g_i', 'g_e', 'g_m', ...
            't, g_i', 't, g_e', 't, g_m', 't, g_i, g_m', 't, g_i, g_e', ...
            'g_i, g_e', 'g_i, g_m', 'g_e, g_m', 'g_{all}', 'all'};
% 
% % Run PEB across reduced second level model space
% %==========================================================================
% for f = 1:length(fields)
%     [PEB, RCM]  = spm_dcm_peb(FCM, M, fields{f});
%     P(f).PEB = PEB;
%     P(f).F   = PEB.F;
%     F(f) = PEB.F;
% end
% save([Fdcm fs 'PEB.mat'], 'P');
% 
% % Calculate the Free energy difference between winning and second model
% %--------------------------------------------------------------------------
% [Fstd Fstg] = sort(F, 'descend');
% dF          = Fstd(1) - Fstd(2);
% 
% % Plot Bayesian model comparison over second level model space
% %--------------------------------------------------------------------------
% subplot(2,1,1), bar(F - min(F));          title(['Free Energy, dF = ' num2str(dF)]);
% subplot(2,1,2), bar(spm_softmax(F'));     title('Posterior Probability');

% Bayesian model reduction over winning PEB (full model)
%--------------------------------------------------------------------------

[PEB RCM]   = spm_dcm_peb(FCM, M, fields{end});
BMA         = spm_dcm_peb_bmc(PEB);
FEB.BMA     = BMA;
FEB.RCM     = RCM;
FEB.PEB     = PEB;

save([Fdcm fs 'Full Empirical Bayes.mat'], 'FEB');

%% Plot DCM outputs
%==========================================================================
load([Fdcm fs 'Full Empirical Bayes']);
RCM     = FEB.RCM;

% Extract predicted and observed spectra
%--------------------------------------------------------------------------
for r = 1:length(RCM)
    pre(:,r) = log(abs(RCM{r}.Hc{1}));
    obs(:,r) = log(abs(RCM{r}.xY.y{1}));
end

% Define plotting ranges
%--------------------------------------------------------------------------
cbar    = [-2.8 1];
yrange  = [2 15];

fqaxis  = RCM{1}.M.Hz;
tmaxis  = linspace(-45, 45, length(RCM)/2);

first   = 1:length(RCM)/2;
second  = first + length(RCM)/2;

subplot(2,2,1)
    imagesc(tmaxis, fqaxis, obs(:,first), cbar);
    ylim(yrange);
    set(gca, 'ydir', 'normal');
    title('Control, observed');

subplot(2,2,2)
    imagesc(tmaxis, fqaxis, obs(:,second), cbar);
    ylim(yrange);
    set(gca, 'ydir', 'normal');
	title('Patient, observed');

subplot(2,2,3)
    imagesc(tmaxis, fqaxis, pre(:, first), cbar);
    ylim(yrange);
    set(gca, 'ydir', 'normal');
    title('Control, predicted');

subplot(2,2,4)
    imagesc(tmaxis, fqaxis, pre(:, second), cbar);
    ylim(yrange);
    set(gca, 'ydir', 'normal');
    title('Patient, predicted');

cmap = flip(cbrewer('div', 'Spectral', 100));
colormap(cmap);

%% Plot individual Parameter Values
%--------------------------------------------------------------------------
Np      = length(FEB.BMA.Pnames);
ploti   = [Tplot, Gplot+length(Tplot)];

Ep = reshape(FEB.BMA.Ep, [Np,4]);
Cp = reshape(diag(FEB.BMA.Cp), [Np,4]);

for e = 1:3
    subplot(3,1,e)
    spm_plot_ci(Ep(ploti,e), Cp(ploti,e));
    ylim([-1 2]);
end
set(gca, 'XTick', 1:Np, 'XTickLabel', {Tlabel{:}, Glabel{:}});
