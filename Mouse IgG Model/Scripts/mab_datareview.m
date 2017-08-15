%% mab_datareview
%==========================================================================
% This routine takes the individual traces, corrects for large 
% artefacts and plots various outputs as well as performing statistical
% tests on the frequency output in the four different conditions

clear all

% Housekeeping
%--------------------------------------------------------------------------
D           = mab_housekeeping;
fs          = filesep;
Fbase       = D.Fbase;
Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fdcm        = D.Fdcm;
Fanalysis   = D.Fanalysis;
mind        = D.mind;
load(D.LFP);

% Routine specific definitions
%--------------------------------------------------------------------------
slideplot   = 0;  % Plots individual sliding window frew plots per subject
frq         = [0.5 100];

%% Loop through subjects - No 3 and 9 excluded because of bad signal
%==========================================================================
for mm = mind

% Extract and denoise the signal
%--------------------------------------------------------------------------
clear preft postft dwin drft d set
h       = m{mm}.hdr;
p.pre   = m{mm}.pre;
p.post  = m{mm}.post;

% Recursively baseline correct and remove large amplitude artefacts
%--------------------------------------------------------------------------
d = p.pre;
for r = 3:-1:1
    d           = d-mean(d);
    ran         = 200 * 2^r; 
    d(d>ran)    = 0;
    d(d<-ran)   = 0;
    d = d-mean(d);
end
p.pre = d;
clear d

d = p.post;
for r = 3:-1:1
    d           = d-mean(d);
    ran         = 200 * 2^r; 
    d(d>ran)    = 0;
    d(d<-ran)   = 0;
    d = d-mean(d);
end
p.post = d;

% Sliding window to calculate frequency components with change over time
%--------------------------------------------------------------------------
l   = length(m{mm}.pre);
win = 60* h.Fs;
stp = 15 * h.Fs;
i   = 0;

for s = 1:stp:l-win
    
    frqi    = frq * win/h.Fs * 0.5;
    i       = i+1;
    
    d       = p.pre;
    dwin    = d(s:s+win);   
    drft    = abs(fft(dwin));
    preft(i,:) = (drft(frqi(1):frqi(2))); 
    clear d dwin drft 
    
    d       = p.post;
    dwin    = d(s:s+win);   
    drft    = abs(fft(dwin));
    postft(i,:) = (drft(frqi(1):frqi(2)));
    
end
frq_axis = linspace(frq(1), frq(2), size(postft,2));
tim_axis = [1:stp:l-win] / h.Fs;

p.preft      = preft;
p.postft     = postft;
m{mm}.preft  = preft;
m{mm}.postft = postft;
m{mm}.prefix = p.pre;
m{mm}.postfix = p.post;

% Plot individual subject's frequency spectra
%--------------------------------------------------------------------------
if slideplot
figure 
subplot(3,3,1:6),
    imagesc(tim_axis, frq_axis, [(preft) ; (postft)]');
    ylim([1 20]);
    title(m{mm}.name);
    set(gca, 'YDir', 'normal');

subplot(3,3,7:9),
    plot([p.pre, p.post]);
    ylim([-400 300])
    xlim([0 Inf]);
end
end 
save([Fanalysis fs 'LFP_concat.mat'], 'm');

%% Plot summary frequency plots for four conditions
%==========================================================================
% 2 by 2 design: pre and post PTZ, w/ and w/o Ab

close all
figure 
clear Ab Co preCo preAb postCo postAb
ab = 0;
co = 0;
colormap jet

% Loop through subjects to extract data from the structure
%--------------------------------------------------------------------------
for mm = mind 

    if m{mm}.id == 'c'
        co = co + 1;
        preCo(co,:,:)   = m{mm}.preft;
        postCo(co,:,:)  = m{mm}.postft;
    else 
        ab = ab + 1;
        preAb(ab,:,:)   = m{mm}.preft;
        postAb(ab,:,:)  = m{mm}.postft;
    end
end

% Save condition spectra into separate struct array
%--------------------------------------------------------------------------
C(1).all = preCo;       cname{1} = 'Control pre PTZ';
C(2).all = postCo;      cname{2} = 'Control post PTZ';
C(3).all = preAb;       cname{3} = 'Patient pre PTZ';
C(4).all = postAb;      cname{4} = 'Patient post PTZ';

% Caluclate standard errors, and plot means and SEMs
%--------------------------------------------------------------------------
for c = 1:4
    C(c).msub = squeeze(mean(C(c).all,1));
    C(c).msubtim = squeeze(mean(C(c).msub));
    for f = 1:size(C(c).all,3)
        matr = C(c).all(:,:,f);
        dots = matr(:);
        sem  = std(dots) / sqrt(length(dots));
        C(c).up(f) = mean(mean(matr,1)) + 2 * sem;
        C(c).lo(f) = mean(mean(matr,1)) - 2 * sem;
    end
end

Fs          = m{1}.hdr.Fs;
l           = length(m{1}.pre);

for c = 1:length(C)
    subplot(2,2,c)
        imagesc(tim_axis, frq_axis, C(c).msub')
        set(gca, 'YDir', 'normal')
        title(cname{c});
        ylim([1 10])
end

%% Plot summary frequency distributions for four conditions
%==========================================================================
% Pre PTZ Plot
%--------------------------------------------------------------------------
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100, 700 1000]); 
subplot(2,1,1)
	c = [1 2 3 4];
    conds = {'Control pre PTZ', 'Control post PTZ', 'Antibody pre PTZ', 'Antibody post PTZ'};
    plot(frq_axis(1),C(c(1)).msubtim(1),'r'); hold on    % Bug req initial plot
    [mc1] = plotshaded(frq_axis, [C(c(1)).up; C(c(1)).msubtim; C(c(1)).lo], 'b'); 
    [mp1] = plotshaded(frq_axis, [C(c(2)).up; C(c(2)).msubtim; C(c(2)).lo], 'r'); 
    
    % Set plotting parameters
    ylim([0 max([C.msubtim])]);
    xlim([0 10]);
 	box off
        
    % Labels and Fonts
    s1_t = title('Power spectra for control mice');
    legend([mc1, mp1], conds{c(1:2)}); 
    yl1 = ylabel('Power');

% Post PTZ Plot
%--------------------------------------------------------------------------
subplot(2,1,2)
    plot(frq_axis(1),C(c(3)).msubtim(1),'b'); hold on    % Bug req initial plot
    [mc2 sc2] = plotshaded(frq_axis, [C(c(3)).up; C(c(3)).msubtim; C(c(3)).lo], 'b');
    [mp2 sp2] = plotshaded(frq_axis, [C(c(4)).up; C(c(4)).msubtim; C(c(4)).lo], 'r'); 
    
    % Set plotting parameters
    set(gca, 'TickDir', 'out');
    ylim([0 max([C.msubtim])]);    
    xlim([0 10])
    box off
    
    % Labels and Fonts
    s2_t = title('Power spectra for antibody-treated mice');
    legend([mc2, mp2], conds{c(3:4)});
    xl2 = xlabel('Frequency [Hz]');
    yl2 = ylabel('Power');

set([s1_t, s2_t], 'FontSize', 12);
set([s1_t, s2_t, xl2, yl1, yl2], 'FontWeight', 'bold');


% Calculate Statistics
%==========================================================================
figure
cutoff = 4 * ceil(length(C(1).all) / 100); 
i      = 1;

for c = 1:length(C)
for s = 1:size(C(c).all, 1)
for t = 1:size(C(c).all, 2)
    C(c).lopo(i) = mean(C(c).all(s,t,1:cutoff));
    C(c).hipo(i) = mean(C(c).all(s,t,cutoff:end));
    i = i + 1;
end
end
i = 1;
end

% Plot dotplots of log power within two freq bands
%--------------------------------------------------------------------------
lab = {'Co pre', 'Co post', 'Pt pre', 'Pt post'};

for c = 1:length(C), dat{c} = log(C(c).lopo); end
mab_dotplot(dat, lab, 0);
    
% Set plotting parameters
set(gcf, 'Color', 'w')
    
% Labels and Fonts
yl = ylabel('Log of bandpower');
ti = title('Windowed estimates of delta-band power (1-5Hz)');
set([yl, ti], 'FontWeight', 'bold');
set(ti, 'FontSize', 12);

% Calculate ANOVA
%--------------------------------------------------------------------------
% eff1 denotes the effect of antibody, eff2 denotes the effect of PTZ
% this test also assesses any interaction between these effects

clear ab ptz

i = 1;
for k = 1:length(C(1).hipo)
    y(i) = C(1).lopo(k);
    ab{i}      = 'Co';
    ptz{i}     = 'pre';
    i = i+1;
end

for k = 1:length(C(2).hipo)
    y(i) = C(2).lopo(k);
    ab{i}      = 'Co';
    ptz{i}     = 'post';
    i = i+1;
end

for k = 1:length(C(3).hipo)
    y(i) = C(3).lopo(k);
    ab{i}      = 'Ab';
    ptz{i}     = 'pre';
    i = i+1;
end

for k = 1:length(C(4).hipo)
    y(i) = C(4).lopo(k);
    ab{i}      = 'Ab'; 
    ptz{i}     = 'post';
    i = i+1;
end

p = anovan(y, {ab' ptz'}, 'model', 'interaction', 'varnames', {'ab', 'ptz'});

%% Plot example traces
%==========================================================================
% Time series plots to demonstrate examples of obvious changes in the
% frequency composition
clear d

no_sub = 14;
start  = 5;
downs  = 5;

h           = m{no_sub}.hdr;
Fs          = h.Fs / downs;
startsamp   = start * h.Fs;
fsa         = 500;

pre    = downsample(m{no_sub}.prefix(startsamp:startsamp + 10*h.Fs), downs); 
post   = downsample(m{no_sub}.postfix(startsamp:startsamp + 10*h.Fs), downs);
prft   = abs(fft(m{no_sub}.pre(startsamp:startsamp + 100*h.Fs), fsa));     prft = (prft(1:end/2));
poft   = abs(fft(m{no_sub}.post(startsamp:startsamp + 100*h.Fs), fsa));    poft = (poft(1:end/2));
frq_axis    = linspace(1, h.Fs/2, fsa / 2);
tim_axis    = linspace(0, 10, length(post)); 

figure
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 1000 500]);
subplot(2,1,1), 
    plot(tim_axis, pre, 'k'); hold on
    plot(tim_axis, post - 300, 'r')
    
    % Set plotting parameters
    box off;
    set(gca, 'YTick', [-300 0], 'YTickLabel', {'Post PTZ', 'Pre PTZ'});
    set(gca, 'XColor', 'w', 'XTick', []);
    
    % Labels and Fonts
	ti = title('Example of antibody treated mouse');
    set(ti, 'FontSize', 12, 'FontWeight', 'bold');
    
    
no_sub = 7;
start  = 5;
downs  = 5;

h           = m{no_sub}.hdr;
Fs          = h.Fs / downs;
startsamp   = start * h.Fs;
fsa         = 500;

pre    = downsample(m{no_sub}.pre(startsamp:startsamp + 10*h.Fs), downs); 
post   = downsample(m{no_sub}.post(startsamp:startsamp + 10*h.Fs), downs);
prft   = abs(fft(m{no_sub}.pre(startsamp:startsamp + 100*h.Fs), fsa));     prft = (prft(1:end/2));
poft   = abs(fft(m{no_sub}.post(startsamp:startsamp + 100*h.Fs), fsa));    poft = (poft(1:end/2));
frq_axis    = linspace(1, h.Fs/2, fsa / 2);
tim_axis    = linspace(0, 10, length(post)); 
    
subplot(2,1,2), 
    plot(tim_axis, pre, 'k'); hold on
    plot(tim_axis, post - 300, 'r')
    
    % Set plotting parameters
    box off;
    set(gca, 'YTick', [-300 0], 'YTickLabel', {'Post PTZ', 'Pre PTZ'});
    set(gca, 'XColor', 'w', 'XTick', []);
    
    % Labels and Fonts
	ti = title('Example of control mouse');
  	set(ti, 'FontSize', 12, 'FontWeight', 'bold');
