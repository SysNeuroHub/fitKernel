function [outcome, time] = getCueOnset(d)
 %[outcome, time] = getCueOnset(d)
 % returns true for trials ending in SUCCESS, false for trials ending in FAIL

 %11/7/22 created from getChoice
 
      %[t,trial,~,state] = d.meta.choice.state;
      
      outcome = d.meta.cueTarget.disabled('time',Inf).data;
      outcome = ~outcome;
      %time = d.meta.cueTarget.disabled('time',Inf).time';
      %time = d.meta.sTarget.cued('time',Inf).time';
      time =  d.meta.cueTarget.startTime';
      
      firstFrames = d.meta.cic.firstFrame.time(1)';
      time = time - firstFrames;
      %time(outcome==0) = nan;
 end