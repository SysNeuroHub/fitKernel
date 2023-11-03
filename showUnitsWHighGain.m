%% show tOnset resp of selected units who has high gain (tonset resp w cue) / (tonset rep wo cue)

load('/mnt/syncitium/Daisuke/cuesaccade_data/figPSTH_pop20231026hugo/fitPSTH_pop20231026hugo.mat',...
    'gainInfo_pop','id_pop','Rsqadj_pop');

avgResp_pop = [];
for iunit = 1:numel(id_pop)
    avgResp_pop(:,:,:,iunit) = squeeze(gainInfo_pop(iunit).avgTonsetByCue(:,1,:,:));
end
Rsqadj_pop_r = Rsqadj_pop([2 3 4],:)./Rsqadj_pop([1],:);

%% select units
load('pickUnitsByClass.mat',"funcClass");
[~,okunits] = intersect(id_pop, funcClass.id_all);

Rsqadj_r_selected = Rsqadj_pop_r(:,okunits);
avgResp_selected = avgResp_pop(:,:,:,okunits);
id_selected = id_pop(okunits);


%% extract units with high gain modulation
tgtDirBin = 1; 
nTop = 50;
respWin = [0.05 0.35];
winSamps = gainInfo_pop(1).winSamps;
respWinIdx = find(winSamps >= respWin(1) & winSamps <= respWin(2)); 
gain_selected = squeeze(mean(abs(avgResp_selected(:,respWinIdx,2,:)./avgResp_selected(:,respWinIdx,1,:)), 2));

[~,idx_bestgain] = sort(gain_selected(tgtDirBin,:),'descend');

id_bestgain = id_selected(idx_bestgain(1:nTop));


showScatterTriplets(1-Rsqadj_r_selected,{'vision','eye speed','eye position'},[-0.1 0.8],idx_bestgain(1:nTop));
screen2png(['unitsWBestGain_' num2str(param.cardinalDir(tgtDirBin))]);
close;

mResp = squeeze(mean(avgResp_selected(:,:,:,idx_bestgain(1:nTop)),4));
seResp = squeeze(ste(avgResp_selected(:,:,:,idx_bestgain(1:nTop)),4));

figure;
boundedline(winSamps, squeeze(mResp(tgtDirBin,:,2)), squeeze(seResp(tgtDirBin,:,2)),'r','transparency',.5);
boundedline(winSamps, squeeze(mResp(tgtDirBin,:,1)), squeeze(seResp(tgtDirBin,:,1)),'transparency',.5);
legend('w cue','wo cue');xlabel('time from target onset');ylabel('firing rate');
title(['target at ' num2str(param.cardinalDir(tgtDirBin)) 'deg, n=' num2str(nTop)] );
screen2png(['tonsetRespWBestGain_' num2str(param.cardinalDir(tgtDirBin))]);
close


