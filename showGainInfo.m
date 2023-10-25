function f = showGainInfo(gainInfo)
% use output of geGainInfo
f = figure;

for ii = 1:2
    ax(ii)=subplot(4,2,ii);
    imagesc(gainInfo.winSamps, gainInfo.cardinalDir, squeeze(gainInfo.avgTonsetByCue(:,1,:,ii)));
    vline(gainInfo.respWin,gca,'-','w');
    if ii==1
        title('wo cue');
        ylabel('dir[deg]');
    elseif ii==2
        title('w cue');
        linkcaxes(ax);
        mcolorbar(ax(2),.5);
    end
    
    subplot(4,2,ii+2)
    imagesc(gainInfo.winSamps, gainInfo.cardinalDir, gainInfo.gain_dir(:,:,ii));
    hline(gainInfo.prefDir);vline(0);caxis([0 2]);
    if ii==1
        ylabel('dir[deg]');
    elseif ii==2
        mcolorbar(gca,.5);
    end
    
    subplot(4,2,ii+4)
    imagesc(gainInfo.winSamps, gainInfo.cardinalDist, gainInfo.gain_distCue(:,:,ii));
    vline(0);caxis([0 2]);
    if ii==1
        ylabel('dist from cue');
    elseif ii==2
        mcolorbar(gca,.5);
    end
    
    subplot(4,2,ii+6)
    imagesc(gainInfo.winSamps, gainInfo.cardinalDist, gainInfo.gain_distPref(:,:,ii));
    vline(0); caxis([0 2]);
    if ii==1
        ylabel('dist from pref dir');
    elseif ii==2
        mcolorbar(gca,.5);
    end
    xlabel('time from stim onset[s]');
end
