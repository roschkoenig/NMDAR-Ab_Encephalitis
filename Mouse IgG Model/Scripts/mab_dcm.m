%% mab_dcm
%==========================================================================
% This routine implements and runs DCMs on the windowed spectral single
% channel data 

D           = mab_housekeeping;
Fbase       = D.Fbase;    
Fscripts    = D.Fscripts;  
Fanalysis   = D.Fanalysis; 
Fdcm        = D.Fdcm;
mind        = D.mind;
fs          = filesep;
clear D

%% Model specification
%==========================================================================
rng('default')
clear DCM SLIDE

DCM.xY.Dfile    = [Fanalysis fs 'LFP_MEEG'];
LFP             = spm_eeg_load(DCM.xY.Dfile);
Fs              = fsample(LFP);
smpls           = size(LFP,2);
timax           = linspace(0, smpls/Fs, smpls);
LFP_conds   	= condlist(LFP);

for c = 1:length(LFP_conds)
disp(['Currently at ' num2str(c) ' of ' num2str(length(LFP_conds)) ' windows']);

% Set up DCM details
%--------------------------------------------------------------------------
DCM.options.analysis    = 'CSD';   	% cross-spectral density 
DCM.options.model       = 'CMC';    % structure cannonical microcircuit (for now)
DCM.options.spatial    	= 'LFP';    % virtual electrode input   
DCM.options.Tdcm        = [timax(1) timax(end)] * 1000;     % time in ms

DCM.options.Fdcm    = [1 60];     	% frequency range  
DCM.options.D       = 1;         	% frequency bin, 1 = no downsampling
DCM.options.Nmodes  = 8;          	% cosine reduction components used 
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
DCM.name            = [Fdcm fs 'DCM_' LFP_conds{DCM.options.trials} '.mat'];

SLIDE{c}            = mab_spm_dcm_csd(DCM);
SLIDE{c}.xY.R       = diag(SLIDE{c}.xY.R);
save([Fanalysis fs 'DCM_All'], 'SLIDE');
end
