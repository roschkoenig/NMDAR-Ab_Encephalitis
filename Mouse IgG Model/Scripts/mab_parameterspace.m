% mab_parameterspace
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

load([Fdcm fs 'Full Empirical Bayes.mat'])

%% Illustrate parameter space using first principal eigenmodes
%==========================================================================
% Extract parameters from reduced first level models
%--------------------------------------------------------------------------
for r = 1:length(FEB.RCM)
    T(r,:) = FEB.RCM{r}.Ep.T;
    G(r,:) = FEB.RCM{r}.Ep.G;
end

% Do principal component decomposition
%--------------------------------------------------------------------------
T = full(T);
G = full(G);

[Tcf Tsc] = pca(T, 'Algorithm', 'eig');
[Gcf Gsc] = pca(G, 'Algorithm', 'eig');

% Map first principal components separately for T and G parameters
%--------------------------------------------------------------------------
seps = linspace(0, length(FEB.RCM), 5);
cols = cbrewer('qual', 'Paired', 10);
cols = cols([7 8 9 10],:);

for s = 2:length(seps)
    plid    = seps(s-1) + 1:seps(s);
    
    subplot(2,5,[2 3 7 8]);
    % Plot
        scatter(Tsc(plid,1), Gsc(plid,1), [], cols(s-1,:), 'filled'); hold on
    % Labels
        xlabel('Time constant component');
        ylabel('Connection strength component');
    % Settings
        xlim([-3 3]); ylim([-5 5]);
        axis square
end
legend({'Control, pre PTZ', 'Control, post PTZ', 'Antibody, pre PTZ', 'Antibody, post PTZ'});

subplot(2,5,1), 
% Plot    
    bar(Tcf(Tplot,1)); 
% Labels    
    title('Time constants: First component')
% Settings
    xlim([0 length(Tlabel)+1]); 
    set(gca, 'XTick', 1:length(Tlabel), 'XTickLabel', Tlabel);

subplot(2,5,6), 
% Plot  
    bar(Gcf(Gplot,1)); 
% Labels    
    title('Connection strengths: First component');
% Settings
    xlim([0 length(Glabel)+1]); 
    set(gca, 'XTick', 1:length(Glabel), 'XTickLabel', Glabel);

subplot(2,5,[5 5 9 10]);
for s = 2:length(seps)
    plid    = seps(s-1) + 1:seps(s);
    scatter3(T(plid,1), T(plid,2), G(plid,7), [], cols(s-1,:), 'filled'); hold on
    xlabel('T(1)');
    ylabel('T(2)');
    zlabel('G(7)');
end

set(gcf, 'Position', [300 300 1200 500]);

%% Forward modelling
%--------------------------------------------------------------------------
range1 = [-4 3];
range2 = [-5 5];
stps   = 200;
steps1 = linspace(range1(1), range1(2), stps);
steps2 = linspace(range2(1), range2(2), stps);
Ts     = Tcf(:,1) * steps1; 
Gs     = Gcf(:,2) * steps2;

OCM    = FEB.RCM{1};
Nc     = length(FEB.BMA.Xnames);
Np     = length(FEB.BMA.Pnames);

BasePs.T    = OCM.M.pE.T + FEB.BMA.Ep([1:4] + (Np * (Nc - 1)))';
BasePs.G    = OCM.M.pE.G + FEB.BMA.Ep([5:Np] + (Np * (Nc - 1)))';

clear delta_all Hc_all

for d1 = 1:length(steps1)-1
d1
for d2 = 1:length(steps2)
    Ps          = OCM.Ep;
    Ps.T        = BasePs.T + Ts(:,d1)';
    Ps.G        = BasePs.G + Gs(:,d2)';
    Hctemp      = spm_csd_mtf(Ps, OCM.M, OCM.xU);
    Hc_all{d1,d2}   = Hctemp{1};    
    delta_all(d2,d1)    = mean(abs(Hc_all{d1,d2}(1:4)));
end
end

%%

figure
imagesc(steps1, steps2, log(delta_all));
set(gca, 'Ydir', 'normal');
axis square; colorbar
colormap gray
xlim([-3 3]); ylim([-5 5]);

%% Calculate and plot heatmaps
%--------------------------------------------------------------------------
clear dens plid

seps = linspace(0, length(FEB.RCM), 5);
for s = 2:length(seps)
    plid{s-1} = seps(s-1) + 1:seps(s);
end

for p = 1:length(plid)
for s1 = 2:length(steps1)-1   
for s2 = 2:length(steps2)
    
    % find parameters in t-range
    tT  = Tsc(plid{p});
    hTi = find(tT >= steps1(s1-1) );
    lTi = find(tT < steps1(s1));
    Ti  = intersect(hTi, lTi);
    
    % find parametrs in g-range
    tG  = Gsc(plid{p});
    hGi = find(tG >= steps2(s2-1));
    lGi = find(tG < steps2(s2));
    Gi  = intersect(hGi, lGi);
    
    % find overlap and save occurance number
    both = intersect(Ti, Gi);
    dens{p}(s2,s1) = length(both);
    
end
end
end

%% Plotting routine
heatcols = flip(cbrewer('div', 'Spectral', 100));

for d = 1:length(dens)
    
    sm          = fspecial('gaussian', 20, 20);
    sdens{d}     = filter2(sm, dens{d}); 
    for sm = 1:3
        sdens{d}     = filter2(sm, sdens{d});    
    end

    subplot(1,4,d), 
        limz = [min(min([sdens{:}])) max(max([sdens{:}])) ];
        imagesc(steps1, steps2, sdens{d}, limz); hold on
        axis square
        set(gca, 'Ydir', 'normal');
        colormap(heatcols)
        xlim([-3 3]); ylim([-5 5]);
        
end

figure
subplot(1,3,1)
    contour(steps1(1:end-1), steps2, delta_all, 4);
    xlim([-3 3]); ylim([-5 5]);
    set(gca, 'Ydir', 'normal');
    
    title('Delta power contours');
    colormap gray; axis square

subplot(1,3,2)
    contour(steps1(1:end-1), steps2, log(delta_all), 4);
    xlim([-3 3]); ylim([-5 5]);
    set(gca, 'Ydir', 'normal');
    
    title('Delta power contours: Log transformed');
    colormap gray; axis square
    
subplot(1,3,3)
    delts       = reshape(delta_all, [size(delta_all,1)*size(delta_all,2), 1]);
    dlimvalues  = [50 75 95 99];
    
    for l = 1:length(dlimvalues)
        deltlims(l)   = prctile(delts,dlimvalues(l));
    end
    
    contour(steps1(1:end-1), steps2, delta_all, deltlims)
    xlim([-3 3]); ylim([-5 5]);
    set(gca, 'Ydir', 'normal');
    colormap gray; axis square
    
    title('Delta power contours: Centiles');




%% Integrating Testing for additional variance based on human data
%==========================================================================
% Find human data file and organise everything in structures
%--------------------------------------------------------------------------
Base    = fileparts(Fbase);
load([Base fs 'Patients' fs 'Matlab Files' fs 'Patient_PCA.mat']);

% Mouse PCA
%--------------------------------------------------------------------------
M.BasePs    = BasePs;
M.Tsc       = Tsc;      M.Tcf   = Tcf;
M.Gsc       = Gsc;      M.Gcf   = Gcf;

% Human PCA
%--------------------------------------------------------------------------
H.Tsc       = Pt.dTsc;  H.Tcf   = Pt.dTcf;

% Collate (M) condition specific indices
%--------------------------------------------------------------------------
seps = linspace(0, length(FEB.RCM), 5);
clear plid

for s = 2:length(seps)
    plid{s-1}    = seps(s-1) + 1:seps(s);
end

pid{1} = [plid{1} plid{2}];
pid{2} = [plid{3} plid{4}];

conds = [1 2];  % Control vs Antibody
clear Hc mDelta

for c = 1:length(conds)

condid  = conds(c);   
M.Tmd   = median(M.Tsc(pid{condid}));
M.Gmd   = median(M.Gsc(pid{condid}));

% Set up 'median' antibody condition for mice
%--------------------------------------------------------------------------
BPs     = M.BasePs;
BPs.T   = BPs.T + M.Tmd * M.Tcf(:,1)';
BPs.G   = BPs.G + M.Gmd * M.Gcf(:,1)';

stps    = 200;
steps   = linspace(-2.5, 2.5, stps); 
cols    = flip(cbrewer('div', 'Spectral', stps));

for s = 1:length(steps)
    Ps          = OCM.Ep;
    Ps.T        = BPs.T + steps(s) * H.Tcf(:,1)';
    Ps.G        = BPs.G;
    Hctemp      = spm_csd_mtf(Ps, OCM.M, OCM.xU);
    Hc{c}(s,:)     = abs(Hctemp{1});
    
    figure(1)
        subplot(2,1,c)
        if ~mod(s,10)
        plot(log(abs(Hctemp{1})), 'color', cols(s,:)); hold on
        xlim([1 60]);
        ylim([-6 8]);
        end
        
    mDelta{c}(s) = mean(Hc{c}(s,1:4));
end
end

mx = max(max(log([Hc{1} Hc{2}])));
mn = min(min(log([Hc{1} Hc{2}])));
sccols = cbrewer('qual', 'Paired', 10);
sccols = sccols([8 10],:);

for h = 1:length(Hc)
    figure(2)
    subplot(2,2,h)
        imagesc(OCM.M.Hz, steps, log(Hc{h}), [mn mx]);
        set(gca, 'Ydir', 'normal')
        
    subplot(2,2,[3 4])
        scplot = fix(linspace(1,length(steps), 50));
        scatter(steps(scplot), log(mDelta{h}(scplot)), [], sccols(h,:), 'filled'); hold on
end




