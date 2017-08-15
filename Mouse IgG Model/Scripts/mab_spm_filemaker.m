%% mab_spm_filemaker
%==========================================================================
% This routine generates a single MEEG file from the segmented and cleaned
% data
D           = mab_housekeeping;
Fbase       = D.Fbase;    
Fscripts    = D.Fscripts;  
Fanalysis   = D.Fanalysis; 
Fdcm        = D.Fdcm;
mind        = D.mind;
fs          = filesep;
load(D.LFP);
clear D
 
% Restructure data into array
%==========================================================================
Fs      = m{1}.hdr.Fs;
count   = 0;

for mm = mind
    count = count + 1;
    pre(count,:)    = m{mm}.prefix;
    post(count,:)   = m{mm}.postfix;
    mid(count)      = m{mm}.id;
end

l   = size(pre,2);
win = 60 * Fs;
stp = 15 * Fs;
    
% Sliding window to calculate frequency components with change over time
%--------------------------------------------------------------------------
i   = 0;

for s = 1:stp:l-win
    
    % Saves sliding windows in window no * mouse * data  
    i = i+1;
    d = pre;
    prewin(i,:,:) = d(:,s:s+win);

    d = post;
    postwin(i,:,:) = d(:, s:s+win);
end

% Put into single Cell Array
%--------------------------------------------------------------------------
co_prewin   = prewin(:, find(mid == 'c'), :);
co_postwin  = postwin(:, find(mid == 'c'), :);
pt_prewin   = prewin(:, find(mid == 'p'), :);
pt_postwin  = postwin(:, find(mid == 'p'), :);

cond    = {'CO_pre', 'CO_post', 'PT_pre', 'PT_post'};
data    = {co_prewin, co_postwin, pt_prewin, pt_postwin};
timax   = linspace(0, (size(co_prewin,3) - 1)/Fs, size(co_prewin, 3));

i = 0;
clear ftdata trialname
for dd = 1:length(data)
for ww = 1:size(data{dd},1)
for mm = 1:size(data{dd},2)
    i = i + 1;
    ftdata.trial{i} = squeeze(data{dd}(ww,mm,:))';
    ftdata.time{i}  = timax;
    trialname{i}    = [cond{dd} '_W' num2str(ww)];
end
end
end

ftdata.label        = {'LFP'};

D   = spm_eeg_ft2spm(ftdata, [Fanalysis fs 'LFP_MEEG']);
D   = conditions(D, 1:size(D,3), trialname); 
save(D);