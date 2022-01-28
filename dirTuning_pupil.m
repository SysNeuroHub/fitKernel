load('C:\Users\dshi0006\Downloads\hugo_oephysdata_ch23.mat', ...
    'ch','dd','ephysdata');

nTrials = length(dd.eye);
fs_eye = median([dd.eye.fs]);


eyeData = ephysdata.eye;

respWin = [30 250]; %[ms]
preWin = [-30 30];


%%TODO
%obtain pupil-corrected PSTH (fast and slow timescales)
%sampling rate is same to the original


%% concatenate across trials as in my Cell Rep paper 2018
eyeData_cat = concatenate_eye(eyeData);

psth_cat = [];
spk_cat = [];
for itr = 1:length(eyeData)
    %make and concatenate PSTH
    tbins = cat(1, eyeData(itr).t, eyeData(itr).t(end)+eyeData(itr).dt)-0.5*eyeData(itr).dt;
    hh = histogram(ephysdata.spikes.spk{itr}, tbins);
    psth = hh.Values';
    psth_cat = cat(1, psth_cat, psth);
    
    spk_cat = cat(1, spk_cat, ephysdata.spikes.spk{itr}+t0);
end


%% detect and interpolate blinks
marginSize = 40; %frames
[eyeData_rmblk, blinks] = removeBlinks(eyeData_cat, marginSize);


parea_avg = mean(eyeData_rmblk.parea);


%% decompose back into trials
eyeData_rmblk_tr = eyeData;
t={eyeData.t};
psth_tr = cell(nTrials, 1);
headidx = 1;
for itr = 1:nTrials
    nFrames = length(t{itr});
    thisx = eyeData_rmblk.x(headidx:headidx+nFrames-1);
    thisy = eyeData_rmblk.y(headidx:headidx+nFrames-1);
    thispwdth = eyeData_rmblk.pwdth(headidx:headidx+nFrames-1);
    thisphght = eyeData_rmblk.phght(headidx:headidx+nFrames-1);
    thispsth = psth_cat(headidx:headidx+nFrames-1);
    
    %eyeData_rmblk_tr(itr).t = t{itr};
    eyeData_rmblk_tr(itr).x = thisx;
    eyeData_rmblk_tr(itr).y = thisy;
    eyeData_rmblk_tr(itr).pwdth = thispwdth;
    eyeData_rmblk_tr(itr).phght = thisphght;
    
    psth_tr{itr} = thispsth;
    headidx = headidx+nFrames;
end


%% direction tuning
dirs = unique(dd.targetloc);
mresp = zeros(length(dirs),3);
mparea = [];
for idir = 1:length(dirs)
    
    
    theseTr = intersect(find(dd.successTrials), find(dd.targetloc==dirs(idir)));
    
    resp_c = [];
    parea_c = [];
    for itr = 1:length(theseTr)
        thisTr = theseTr(itr);
        theseTimes  = intersect(find(eyeData(thisTr).t - dd.tOnset(thisTr) > 1e-3*respWin(1)), ...
            find(eyeData(thisTr).t - dd.tOnset(thisTr) < 1e-3*respWin(2)));
        resp_c(itr) = sum(psth_tr{thisTr}(theseTimes))/diff(respWin)/1e-3; %[spikes/s]

        theseTimes_pre  = intersect(find(eyeData(thisTr).t - dd.tOnset(thisTr) > 1e-3*preWin(1)), ...
            find(eyeData(thisTr).t - dd.tOnset(thisTr) < 1e-3*preWin(2)));

       parea_c(itr) = mean(eyeData_rmblk_tr(thisTr).parea(theseTimes_pre));
    end
    mresp(idir,1) = mean(resp_c);
    mresp(idir,2) = nanmean(resp_c(parea_c<parea_avg));
    mresp(idir,3) = nanmean(resp_c(parea_c>parea_avg));
    mparea(idir) = mean(parea_c);
    length(find(parea_c>parea_avg))
end

plot(dirs, mresp);
legend('all trials','parea<mean','parea>mean');
xlabel('tgt direction [deg]');
%polarplot(dirs/180*pi, mresp);
screen2png('dirTuning_pupil_preStim');
