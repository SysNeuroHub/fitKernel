
[saveServer, rootFolder] = getReady();
for aa = 2
    switch aa
        case 1
            animal = 'hugo'; yidx=1:3;
        case 2
            animal =  'ollie'; yidx = 3;
    end
    for yy = yidx
        switch yy
            case 1
                year = '2021';
            case 2
                year = '2022';
            case 3
                year = '2023';
        end
        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);
        nData = numel(loadNames);
        thisdata = 1:nData;


        saveFolder = fullfile(saveServer, year,animal);%17/6/23
        % load(fullfile(saveFolder, 'assembly20241212.mat'),'param',...
        %     'id_pop','expval_tgt_pop','corr_tgt_pop','corr_tgt_rel_pop',...
        %     'latency_r_pop','avgAmp_hm_pop','p_hm_pop','spkOk_th_pop',...
        %     'spkOkTrials_pop','spkOkUCueTrials_pop','mFiringRate_pop',...
        %     'PtonsetResp_pop','errorIDs');

        for idata = thisdata

            datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            thisid = [animal '/' year '/' datech];
            disp(thisid);

            saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];

            saveName = fullfile(saveFolder, [saveSuffix '.mat']);
            if ~exist(saveName,'file'); continue; end;
            load(saveName,'expval_tgt','corr_tgt','latencyStats','avgAmp_hm','p_hm','spkOk_th', ...
                'spkOkTrials','spkOkUCueTrials','mFiringRate','cellclassInfo','kernelInfo','ntargetTrials',...
                'corr_tgt_rel','spkNGRate','CueTrRate','nLatencyTrials_pref_success','prefDir',...
                'ranksumval_hm', 'ranksumz_hm');

            if ~exist('mFiringRate','var') || mFiringRate < 5; continue; end

            thisDate = [months{idata} '_' dates{idata}];

            id_pop{idata} = thisid;

            expval_tgt_pop{idata} = expval_tgt;
            corr_tgt_pop{idata} = corr_tgt;
            corr_tgt_rel_pop{idata} = corr_tgt_rel;

            latency_r_pop{idata} = latencyStats.latency_r;
            latency_p_pop{idata} = latencyStats.latency_p; %15/2/25
            if isfield(latencyStats,'difflatency')
                difflatency_pop{idata} = latencyStats.difflatency; %26/3/25
            else
                difflatency_pop{idata} = nan;
            end
            nLatencyTrials_pref_success_pop{idata} = nLatencyTrials_pref_success;

            avgAmp_hm_pop{idata} = avgAmp_hm;
            p_hm_pop{idata} = p_hm;
            ranksumval_hm_pop{idata} = ranksumval_hm;
            ranksumz_hm_pop{idata} = ranksumz_hm;
            spkOk_th_pop{idata} = spkOk_th;
            spkOkTrials_pop{idata} = spkOkTrials;
            spkOkUCueTrials_pop{idata} = spkOkUCueTrials;
            mFiringRate_pop{idata} = mFiringRate;
            PtonsetResp_pop{idata} = cellclassInfo.PtonsetResp;
            ntargetTrials_pop{idata} = ntargetTrials;
            errorIDs{idata} = 0;

            spkNGRate_pop{idata} = spkNGRate;%(numel(t_tr)-numel(spkOkTrials))/numel(t_tr)*100;
            CueTrRate_pop{idata} = CueTrRate;%(numel(spkOkTrials)-numel(spkOkUCueTrials))/numel(spkOkTrials)*100;
            kernel_pop(:,idata) = kernelInfo.kernel;
            tlags_pop{idata} = kernelInfo.tlags;
            prefDir_pop{idata} = prefDir;

            clear mFiringRate
        end

        save(fullfile(saveFolder, 'assembly20250207.mat'),...
            'id_pop','expval_tgt_pop','corr_tgt_pop','corr_tgt_rel_pop',...
            'latency_r_pop','avgAmp_hm_pop','p_hm_pop','spkOk_th_pop',...
            'spkOkTrials_pop','spkOkUCueTrials_pop','mFiringRate_pop',...
            'PtonsetResp_pop','errorIDs','ntargetTrials_pop','spkNGRate_pop',"CueTrRate_pop",...
            'kernel_pop','tlags_pop','nLatencyTrials_pref_success_pop','prefDir_pop', ...
            "ranksumval_hm_pop","ranksumz_hm_pop",'difflatency_pop');

        clear  id_pop expval_tgt_pop corr_tgt_pop corr_tgt_rel_pop ...
            latency_r_pop avgAmp_hm_pop p_hm_pop spkOk_th_pop ...
            spkOkTrials_pop spkOkUCueTrials_pop mFiringRate_pop...
            PtonsetResp_pop errorIDs ntargetTrials_pop spkNGRate_pop CueTrRate_pop...
            kernel_pop tlags_pop nLatencyTrials_pref_success_pop prefDir_pop...
            ranksumval_hm_pop ranksumz_hm_pop difflatency_pop
    end
end