%legacy analysis. copied from fitPSTH_test     

%% decompose back into individual trials
        %eyeData_rmotl_tr = decompose_eye(eyeData_rmotl_cat, t_tr);
        
        
%         %% upsampling
%         psth_predicted_all = interp1(predictorInfo.t_r, predicted_all, eyeData_rmotl_cat.t, 'linear');
%         psth_predicted = zeros(predictorInfo.nPredictors, length(eyeData_rmotl_cat.t));
%         for ivar = 1:predictorInfo.nPredictors
%             psth_predicted(ivar,:) = interp1(predictorInfo.t_r, predicted(:,ivar), eyeData_rmotl_cat.t, 'linear');
%         end
%         PSTH_f_upsample = interp1(predictorInfo.t_r, PSTH_f, eyeData_rmotl_cat.t, 'linear');
%         % psth_predicted =  predicted_slow+predicted_fast+predicted_fast_x+predicted_fast_y;
%         
%         
%         %% triggered by external events TOBE REMOVED
%         psth_predicted_all_tr = cell(nTrials, 1);
%         psth_predicted_tr = cell(nTrials, predictorInfo.nPredictors);
%         psth_denoised_tr = cell(nTrials,1);
%         psth_tr = cell(nTrials,1);
%         headidx = 1;
%         for itr = 1:nTrials
%             nFrames = length(t_tr{itr});
%             
%             psth_predicted_all_tr{itr} = psth_predicted_all(headidx:headidx+nFrames-1);
%             psth_tr{itr} = PSTH_f_upsample(headidx:headidx+nFrames-1); %11/1/22
%             psth_denoised_tr{itr} = psth_tr{itr} - psth_predicted_all_tr{itr};
%             
%             for ivar = 1:predictorInfo.nPredictors
%                 psth_predicted_tr{itr,ivar} = psth_predicted(ivar, headidx:headidx+nFrames-1);
%             end
%             headidx = headidx+nFrames;
%         end
        
        psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
%         psth_all = cat(2,psth_tr, psth_predicted_all_tr, psth_predicted_tr); %cell(#trials, 2+nPredictors)
        
        %% avg trials
        %         [f, psth_snippet, parea_snippet, dist_snippet, taxis_snippet] ...
        %             = pupilFigureAvgSingle(dd, eyeData_rmotl_tr, psth_all, param.evName, param.figTWin);
        %         [f, psth_snippet, pdiam_snippet, dist_snippet, taxis_snippet] ...
        %             = pupilFigure(dd, eyeData_rmotl_tr, psth_all, param.evName, param.figTWin);
        %         legend(psthNames(2:end),'location','northwest');
        %         screen2png(fullfile(saveFigFolder,['pupilPsth_' param.evName saveSuffix]), f);
        %         close;


 %% individual trials
        theseTimes = intersect(find(taxis_snippet > param.respWin(1)), find(taxis_snippet < param.respWin(2)));
        msnippet = squeeze(mean(psth_snippet(theseTimes, :, :), 1));
        Rmsnippet = corrcoef(msnippet);
        
        fsnippet = figure('position',[0 0 1900 1400]);
        ax4 = [];
        for itype = 1:length(psthNames)
            ax4(itype)=subplot(2,length(psthNames),itype);
            thisData = squeeze(psth_snippet(:,:,itype))';
            imagesc(taxis_snippet, 1:size(psth_snippet,2), thisData);
            %ylim([200 300]);
            title(psthNames{itype});
            if itype==1
                ylabel('trial');
                xlabel(['time from ' param.evName]);
            end
            crange = prctile(thisData(:),[1 99]);
            if diff(crange)==0
                crange = [crange(1) crange(1)+1];
            end
            caxis(crange);
            mcolorbar(ax4(itype), .2);
            
            ax4(itype+length(psthNames))=subplot(2,length(psthNames),itype+length(psthNames));
            plot(msnippet(:,itype), msnippet(:,1),  '.');
            axis square;
            title(['vs ' psthNames{itype} ', R: ' num2str(Rmsnippet(1,itype))]);
            xlabel(psthNames{itype})
        end
        screen2png(fullfile(saveFigFolder,['indtrials_' saveSuffix]));
        close;