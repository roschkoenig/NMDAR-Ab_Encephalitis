% nae_peb
%==========================================================================
% This routine loads the inverted individual time window DCMs and performs
% parametric empirical Bayesian analyiss

clear all
close all

% Housekeeping
%==========================================================================
D           = nae_housekeeping;
fs          = filesep;
Fbase       = D.Fbase;
Fdata       = D.Fdata;
Fscripts    = D.Fscripts;
Fdcm        = D.Fdcm;

clear D

% PEB Analysis
%==========================================================================
% Make second level model space
%--------------------------------------------------------------------------
subfiles = cellstr(spm_select('List', Fdata, '^N.*\.mat$'));

for s = 1:length(subfiles) 
    
sub         = subfiles{s}(1:end-4);
dcmfiles    = cellstr(spm_select('FPList', Fdcm, ['^DCM_' sub '.*.mat$']));
conds       = [];
cols        = flip(cbrewer('qual','Paired', 6));

for d = 1:length(dcmfiles), conds{d} = dcmfiles{d}(end-4); end

X      = [];
Xnames = {'All', conds{:}};
X(:,1) = ones(length(conds),1);

for c = 1:length(conds)
    X(:,1+c) = zeros(1,length(conds));
    X(c,1+c) = 1; 
end

% Run PEB
%--------------------------------------------------------------------------
FCM = spm_dcm_load(dcmfiles);

M.X         = X;
M.Xnames    = Xnames;
M.Q         = 'all';

% Time constant parameters 
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
% G Parameters: The order in the DCM structure is as follows
% j     = [7 2 3 4 5 6 8 9 10 1];
% i.e.     M I I M E E E I M  M
% new   =  1 2 3 4 5 6 7 8 9 10
%
% T Parameters: The order is as follows
% ss sp ii dp
%  1  2  3  4  
%  E  I  E  I

% Model space by parameter type nae_spm_fx_cmc
%--------------------------------------------------------------------------
clear fields F labels
fields{1}   = {'T(1)', 'T(2)', 'T(3)', 'T(4)'};     % time constants
fields{2}   = {'G(2)', 'G(3)', 'G(8)'};             % inh connections
fields{3}   = {'G(5)', 'G(6)', 'G(7)'};             % exc connections
fields{4}   = {'G(1)', 'G(4)', 'G(9)', 'G(10)'};    % mod connections
 
labels  = { 't', 'g_i', 'g_e', 'g_m' };

saveit = 1;

% Run PEB across reduced second level model space
%==========================================================================

for f = 1:length(fields)
    
    [PEB, RCM]  = spm_dcm_peb(FCM, M, fields{f});
    F(f) = PEB.F;
    
    if saveit try load([Fdcm fs 'PEB']); catch P = []; end; end
    P(s,f).PEB    = PEB;
    P(s,f).fields = fields{f};
    P(s,f).F      = PEB.F;
    if saveit, save([Fdcm fs 'PEB'], 'P'); end

end
end

%% Identify overall winning PEB model
%==========================================================================
if saveit, load([Fdcm fs 'PEB']); end

clear Fs Fall
for p = 1:size(P,1)
for m = 1:size(P,2)
    Fs(p,m) = P(p,m).F; 
end
end
Fall    = sum(Fs);

% Plot free energies and model posteriors spm_dcm_bmc
%--------------------------------------------------------------------------
subplot(3,1,1), bar(Fall - min(Fall));
subplot(3,1,2), plot(Fs' - min(Fs'));
subplot(3,1,3), 
    [alpha, exp_r, xp] = spm_BMS(Fs, 1e6, 1, 0, 1);
    bar(xp); title('RFX Analysis');
    
set(gca, 'XTick', 1:length(Fall), 'XTickLabel', labels); 

[v l] = max(Fall);

%% Bayesian model reduction over winning PEB (Time constants and inhibitory)
%--------------------------------------------------------------------------
load([Fdcm fs 'PEB']);
for p = 1:size(P,1)
    clear Snames Sfixed sep macseps winseps
    Snames = P(p,l).PEB.Snames
    
    for s = 1:length(Snames)
        macseps     = find(P(p,l).PEB.Snames{s} == '/');
        winseps     = find(P(p,l).PEB.Snames{s} == '\');
        if length(winseps) > length(macseps), sep = winseps(end);
        else sep = macseps(end); end
        
        Sfixed{s} = [Fdcm fs Snames{s}(sep + 1:end)];
    end
    
    FCM         = spm_dcm_load(Sfixed);
    M.X         = P(p,l).PEB.M.X
    M.Xnames    = P(p,l).PEB.Xnames;
    M.Q         = 'all';
    [PEB RCM]   = spm_dcm_peb(FCM', M, fields{l});
    try PMA         = spm_dcm_peb_bmc(PEB); catch PMA = []; end
    
    FEB(p).PMA    = PMA;
    FEB(p).RCM    = RCM;
    FEB(p).PEB    = PEB;
    
end
save([Fdcm fs 'Full Empirical Bayes'], 'FEB');
