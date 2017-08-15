% Housekeeping
%==========================================================================
clear all
D           = nae_housekeeping;
fs          = filesep;
Fdata       = D.Fdata;
Fdcm        = D.Fdcm;
files       = cellstr(spm_select('List', Fdata, '^N.*\.mat$'));

%%
for f = 1:length(files)

% Set up DCM structure and invert baseline
%==========================================================================
DCM = [];

sub = files{f}(1:end-4);
% Fix directory of canonical forward matrix
%--------------------------------------------------------------------------
DCM.xY.Dfile        = [Fdata fs 'M' files{f}];

% Load MEEG object and extract sampling rate and info
%--------------------------------------------------------------------------
LFP                 = spm_eeg_load(DCM.xY.Dfile);
Fs                  = fsample(LFP);
smpls               = size(LFP,2);
timax               = linspace(0, smpls/Fs, smpls);
clist               = condlist(LFP);

for c = 1:length(clist)
% Set up DCM details
%--------------------------------------------------------------------------
DCM.options.analysis    = 'CSD';   	% cross-spectral density 
DCM.options.model       = 'CMC';    % structure cannonical microcircuit (for now)
DCM.options.spatial    	= 'LFP';    % virtual electrode input   
DCM.options.Tdcm        = [timax(1) timax(end)] * 1000;     % time in ms

DCM.options.Fdcm    = [1 60];     	% frequency range  
DCM.options.D       = 1;         	% frequency bin, 1 = no downsampling
DCM.options.Nmodes  = 8;          	% number of eigenmodes
DCM.options.han     = 0;         	% no hanning 
DCM.options.trials  = c;            % index of ERPs within file

DCM.Sname           = chanlabels(LFP);
DCM.M.Hz            = DCM.options.Fdcm(1):DCM.options.D:DCM.options.Fdcm(2);
DCM.xY.Hz           = DCM.M.Hz;

% Create DCM Struct and specify DCM.options 
%--------------------------------------------------------------------------
DCM.A     	= {1 1 1};
DCM.B    	= {};
DCM.C   	= sparse(length(DCM.A{1}),0);

% Reorganise model parameters in specific structure
%==========================================================================
DCM.M.dipfit.Nm     = DCM.options.Nmodes;
DCM.M.dipfit.model 	= DCM.options.model;
DCM.M.dipfit.type   = DCM.options.spatial;
DCM.M.dipfit.Nc     = size(LFP,1);
DCM.M.dipfit.Ns     = length(DCM.A{1});

% Load empirical priors
%--------------------------------------------------------------------------
[pE,pC]             = nae_spm_cmc_priors(DCM.A,DCM.B,DCM.C);
load([Fdcm fs 'Priors' fs 'Priors.mat']);
pE                  = Priors.pE;
pC                  = Priors.pC;

DCM.M.pE    = pE;
DCM.M.pC    = pC;

DCM.name            = [Fdcm fs 'DCM_' sub '_' clist{c} '.mat'];
DCM                 = nae_spm_dcm_csd(DCM);
end
end

%% Review DCM fits
%==========================================================================
subfiles = cellstr(spm_select('List', Fdata, '^N.*\.mat$'));

for s = 1:length(subfiles)
    
sub         = subfiles{s}(1:end-4);
dcmfiles    = cellstr(spm_select('FPList', Fdcm, ['^DCM_' sub '.*.mat$']));
conds       = [];
cols        = flip(cbrewer('qual','Paired', 6));

for d = 1:length(dcmfiles) 
    conds{d} = dcmfiles{d}(end-4); 
    load(dcmfiles{d});
    subplot(4,2,s)
    plot(log(abs(DCM.Hc{1})), 'color', cols(2*d,:)); hold on
    plot(log(abs(DCM.xY.y{1})), 'color', cols(2*d-1,:)); hold on
end
end
legend({'B pred', 'B obs', 'P pred', 'P obs', 'R pred', 'R obs'});
hold off

%% Plot N004 example traces
%--------------------------------------------------------------------------
subfiles = cellstr(spm_select('List', Fdata, '^N.*\.mat$'));
s = 4;

sub         = subfiles{s}(1:end-4);
MEEG        = spm_eeg_load([Fdata fs subfiles{s}]);

trials      = [5 52 81];
chid        = [1 11 15 9 2 12 15 10 3 5 7 13 4 6 8 14 17 18 19];

for t = 1:length(trials)
    subplot(2, length(trials), t)
    for c = 1:size(MEEG,1)
        plot(squeeze(MEEG(chid(c),:,trials(t)))-c*100);
        hold on
        ylim([-2100 200]);
        yvals(c)   = -c*100;
    end
end

subplot(2,length(trials),1)
chlabs = chanlabels(MEEG);
set(gca, 'YTick', flip(yvals), 'YTickLabel', flip(chlabs(chid)));


MEEG        = spm_eeg_load([Fdata fs 'M' subfiles{s}]);

subplot(2, length(trials), [1:length(trials)]+length(trials))

for t = 1:length(trials)
    plot(squeeze(MEEG(1,:,trials(t))) - t*100); hold on
end

