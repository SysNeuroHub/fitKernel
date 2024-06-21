function hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents)
% hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents)

if nargin < 3
    validEvents = 1:numel(onsets_cat.tOnset);
end

latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, validEvents, 1);
latency_bhv_cOn = getTgtBhvLatency(onsets_cat, dd, validEvents, 0);

hesitant = (latency_bhv_cOn - latency_bhv_srt > 0.1); %multiple saccades
hesitant = hesitant + ~isnan(latency_bhv_cOn).*isnan(latency_bhv_srt); % choice without saccades (slow eye movement)
