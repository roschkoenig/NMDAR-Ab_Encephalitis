function D = nae_housekeeping

fs          = filesep;

if strcmp(computer, 'PCWIN64'), Fbase = 'C:\Users\rrosch\Dropbox\Research\Friston Lab\1605 NMDAR encephalitis\1703 Translational\Patients';
else                            Fbase = '/Users/roschkoenig/Dropbox/Research/Friston Lab/1605 NMDAR encephalitis/1703 Translational/Patients';
end

Fanalysis   = [Fbase fs 'Matlab Files'];
Fscripts    = [Fbase fs 'Scripts'];
Fdata       = [Fbase fs 'Data'];
Fdcm        = [Fanalysis fs 'DCM'];

addpath(genpath(Fscripts));
spm('defaults', 'eeg');

D.Fbase     = Fbase;
D.Fanalysis = Fanalysis;
D.Fscripts  = Fscripts;
D.Fdcm      = Fdcm;
D.Fdata     = Fdata;