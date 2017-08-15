% nae sources
%==========================================================================
% This function will source localise the abnormal activity (paroxysmal and
% rhythmic) using an IID approach and extract a virtual-electrode response
% at the cortical source of maximal power

clear all
D       = nae_housekeeping;
fs      = filesep;
Fdata   = D.Fdata;
files   = cellstr(spm_select('List', Fdata, '^N.*\.mat$'));

%% Prepare files
%==========================================================================
for f = 1:length(files)
    
% Set default 10/20 EEG sensors
%--------------------------------------------------------------------------
MEEG        = spm_eeg_load([Fdata fs files{f}]);    
S.task      = 'defaulteegsens';
S.D         = MEEG;
MEEG        = spm_eeg_prep(S);
S           = [];
save(MEEG)

% Compute leadfields for inverse solutions
%--------------------------------------------------------------------------
conds = condlist(MEEG);
clear job  
job{1}.spm.meeg.source.invert.D = {[Fdata fs files{f}]};
job{1}.spm.meeg.source.invert.val = 1;
job{1}.spm.meeg.source.invert.whatconditions.condlabel = conds(2:end);
job{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'IID';
job{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
job{1}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];
job{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
job{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
job{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
job{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
job{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
job{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
job{1}.spm.meeg.source.invert.modality = {'EEG'};
spm_jobman('run', job)
end

%% Extract maximum intensity MNI voxel locations from source inversion
%==========================================================================
for f = 1:length(files)
    load(files{f})
    
    clear mMAP xMAP
    for d = 1:length(D.other.inv{1}.inverse.J)
        mMAP{d} = mean(D.other.inv{1}.inverse.J{d},2);
        xMAP(d) = max(mMAP{d});
    end
    [v i] = max(xMAP);
    [v l] = max(mMAP{i});
    
    L(f).xyz = fix(D.other.inv{1}.forward.mesh.vert(l,:) * 1000);
    L(f).name = files{f}(1:end-4);
end
save([Fdata fs 'MIP_Locations.mat'], 'L');

%% Extract source waveforms
%==========================================================================
clear L
load([Fdata fs 'MIP_Locations.mat']);

for l = 1:length(L)
    MEEG            = spm_eeg_load([Fdata fs L(l).name '.mat']);
    xyz             = L(l).xyz;
    scalefactor     = sqrt(xyz(1)^2 + xyz(2)^2 + xyz(3)^2);    
    
    S.D             = [Fdata fs L(l).name '.mat'];
    S.dipoles.pnt   = L(l).xyz; 
    S.dipoles.ori   = L(l).xyz / scalefactor;
    S.dipoles.label = {'LFP'};
    
    sD              = spm_eeg_dipole_waveforms(S);
    save(sD);
end
