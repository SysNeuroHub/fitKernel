
 function [fOnset] = getFOnset(d)
 %[outcome, time] = getFOnset(d)
 %    fOnset; % fixation behaviour onset time (seconds)
 %
% this function should be obviated if get.numTrials in mdbase works properly
 %cf. cuesaccade.get.fOnset

  % fixation time, t1
      [t,trial,~,state] = d.meta.fixation.state;

      % we want the start of the FIXATING state, but only for SUCCESSful trials...
      success = strcmp(state,'SUCCESS');
      fixating = strcmp(state,'FIXATING');

      ix = ismember(trial,trial(success)) & fixating;
      ff = d.meta.cic.firstFrame.time; % first frame, one per trial
      
      fOnset = NaN([numel(d.meta.condIds),1]);
      fOnset(trial(ix)) = t(ix) - ff(trial(ix));
      
      
%       [t,trial,~,state] = d.meta.choice.state;
%       
%       stateNEW = state;
%       stateNEW = strjoin(stateNEW);
%       fixedStates = strrep(stateNEW, 'FIXATING FREEVIEWING FIXATING SUCCESS', 'FIXATING FREEVIEWING FIXATING FAIL');
%       fixedStates = strrep(fixedStates, '  ', ' 0 ');
%       fixedStates = strsplit(fixedStates);
% 
%       % we want the start of the FIXATING state, but only for SUCCESSful trials...
%       success = strcmp(fixedStates,'SUCCESS');
%       fixating = strcmp(fixedStates,'FIXATING');
% 
%       outcome = nan([d.numTrials,1]);
%       time = nan([d.numTrials,1]);
%       %ix_s = strcmpi(state,'SUCCESS');
%       ix_s = ismember(trial,trial(success)) & fixating; 
%      
%       outcome(trial(ix_s)) = true;
%       time(trial(ix_s)) = t(ix_s);
%       
%       %       ix_f = strcmpi(fixedStates,'FAIL') & fixating;
%       ix_f = strcmpi(fixedStates,'FAIL');
%       outcome(trial(ix_f)) = false;
%       time(trial(ix_f)) = t(ix_f);
%       
%       firstFrames = d.meta.cic.firstFrame.time(1)';
%       time = time - firstFrames;
%       
%       tOnset = d.tOnset;
%       time(isnan(tOnset)) = nan;
%       outcome(isnan(tOnset)) = nan;
 end