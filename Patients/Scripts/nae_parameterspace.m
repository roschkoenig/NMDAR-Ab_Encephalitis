% Housekeeping
%==========================================================================
clear all
D           = nae_housekeeping;
fs          = filesep;
Fdata       = D.Fdata;
Fdcm        = D.Fdcm;
Fbase       = D.Fbase;
files       = cellstr(spm_select('List', Fdata, '^N.*\.mat$'));

% Extract Parameters
%--------------------------------------------------------------------------
load([Fdcm fs 'Full Empirical Bayes']);
count = 1;

for p = 1:length(FEB)
for r = 1:length(FEB(p).RCM)
    Ts(count,:) = FEB(p).RCM{r}.Ep.T;
    cond(count) = FEB(p).RCM{r}.xY.code{1};
    count = count + 1;
end
end

% Identify specific conditions
%--------------------------------------------------------------------------
bi = find(cond == 'B');
pi = find(cond == 'P');
ri = find(cond == 'R');

% Assess differences between background and paroxysms / rhythmic abnormalities
%--------------------------------------------------------------------------
dT      = [];
count   = 1;

idx     = [bi size(Ts,1)];

for b = 1:(length(idx)-1)
    thislot     = idx(b):idx(b+1)-1;
    
    if length(thislot) == 2         % isolated paroxysms
        dT(count,:)     = Ts(thislot(2),:) - Ts(thislot(1),:);
        count           = count + 1;
        
    elseif length(thislot) == 3     % rhythmic abnormalities
        dT(count,:)     = Ts(thislot(2),:) - Ts(thislot(1),:);
        count           = count + 1;
        dT(count,:)     = Ts(thislot(3),:) - Ts(thislot(1),:);
        count           = count + 1;
    end
end


% Perform Bartlett's test for equal variance
%--------------------------------------------------------------------------
p = vartestn(dT, 'display', 'off');
nae_dotplot({dT(:,1), dT(:,2), dT(:,3), dT(:,4)}, {'ss', 'sp', 'ii', 'dp'}, 0);

% Decompose into components and save
%--------------------------------------------------------------------------
[dTcf dTsc] = pca(dT);

Pt.dTcf = dTcf;
Pt.dTsc = dTsc;
Pt.dT   = dT;
save([Fbase fs 'Matlab Files' fs 'Patient_PCA.mat'], 'Pt');
