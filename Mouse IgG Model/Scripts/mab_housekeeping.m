function D = mab_housekeeping

if strcmp(computer, 'MACI64')
    Fbase    = '/Users/roschkoenig/Dropbox/Research/Friston Lab/1605 NMDAR encephalitis/1703 Translational/Mouse IgG Model';
else 
    Fbase    = 'C:\Users\rrosch\Dropbox\Research\Friston Lab\1605 NMDAR encephalitis\1703 Translational\Mouse IgG Model';
end

fs          = filesep;
Fscripts    = [Fbase fs 'Scripts'];
Fdata       = [Fbase fs 'Data'];
Fanalysis   = [Fbase fs 'Matlab Files'];
Fdcm        = [Fanalysis fs 'DCM Files'];
LFP         = [Fanalysis fs 'LFP_concat'];

Tplot   = [2 1 3 4];                  
Tlabel  = {'sp', 'ss', 'ii', 'dp'}; 
Gplot   = [2 3 8 7 5 6 1 10 4 9];     
Glabel  = {'sp>ss', 'ii>ss', 'ii>dp', 'ss>sp', 'ss>ii', 'dp>ii', 'sp', 'ss', 'ii', 'dp'}; 


spm('defaults', 'eeg');
addpath(genpath(Fscripts));

D.Fbase     = Fbase;
D.Fscripts  = Fscripts;
D.Fdata     = Fdata;
D.Fanalysis = Fanalysis;
D.Fdcm      = Fdcm;
D.LFP       = LFP;
D.mind      = [1 2 4 5 6 7 8 10 11 12 13 14 15];

D.Tplot     = Tplot;
D.Tlabel    = Tlabel;
D.Gplot     = Gplot;
D.Glabel    = Glabel;

