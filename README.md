# fitKernel

"fitKernel" is a library of matlab functions and scripts to compute Generalized Linear Model (GLM) for cuesaccade data.


# Requirement

* marmodata (loading eye data)
* marmolab-stimuli (loading cuesaccade data)
* dsbox (ridge regression and some visualization)
* https://github.com/pillowlab/neuroGLM (GLM fitting)

Environments under linux and windows machines are tested.



# Usage
mainScript.m is the main script to to "everything": obtain linear kernel and event-triggered spiking traces

Functions about concatenating trial-based data
* concatenate_eye: mamordata.eye concatenated across trials
* concatenate_spk: spike times conctenated across trials
* decompose_eye: eye data decomposed back into trials from concatenated
* retrieveTrIdx_r: retrieve temporal indices of every trial from the concatenated data

Main function for kernel estimation
* fitPSTH_cv: estimate linear kernel with K-fold cross validation using neuroGLM
* fitPSTH: estimate linear kernel with ridge regression in dsbox (obsolete)
* fitMultiplicative: (under testing)

Functions about preparing predictor signals
* processEyeData: main function to concatenate eye data, extract times of saccades, blinks and outliers
* detectOverlapEvents: detect overlap between the 2nd and 3rd inputs  (cf. trace2Event)
* detectPupilOnsets: start/end times of pupil dilation and constriction according to Einstein 2017 JNS
* getChoice: times when task outcome is registered (success/failure) of every trial
* getCueOnset: cue onset times and the task outcome (success/failure) of every trial
* (getCueDirMtx): matrix of cue direction of all trials (NOT YET FUNCTIONING)
* getEyeDirMtx: matrix of eye position over time
* getEyeSpdDirMtx: matrix of eye velocity over time
* getSaccMtx: matrix of saccade direction over time
* getTgtDirMtx: matrix of target direction over time
* getPupilDiameter: converts pupil area into diameter and its percentile
* getSaccDir: classifies each saccade events into one of quantized bins of directions
* removeBlinksEDF: eye data after removing times of blinks detected by eyelink
* removePareaOutliers: eye data after removing times of outliers that crossed a specified threshold
* selectSaccades: saccade times after eliminating those at specified exclusion periods
* preparePredictors: matrix of predictors of specified modalities

Functions for analysis/visualization
* showKernel: shows kernel and predicted PSTH trace over time
* showSaccOnsetResp: shows observed/predicted PSTH traces triggered by saccade onset OUTSIDE of the task
* showFixCueOnsetResp: shows observed/predicted PSTH triggered by fixation and cue onsets
* showTonsetResp: shows observed/predicted PSTH triggered by target and saccade onsets of preferred target direction. Also computes cell class (under development)
* showTonsetByCue: show observed/predicted PSTH triggered by target onsets with/without cues
* (dirTuning): preferred direction of spikes at specified time window (TO BE DELETED)
* (parea_spike_ms): computes correlation between spikes and pupil diameter at specified temporal frequencies
* (pupilFigure): figure of single-trial traces of pupil position/area and spikes triggered at specified type of events 
* (showPspec_parea_psth NOT CHECKED)
* (showDirTuning: NOT COMPLETED) : preferred direction of spikes at specified time window 

Utility functions
* filtPSTH: (a)causally filtered spike trace
* alignMtxDir: a matrix which 2nd dimension is rotated based on its preferred direction
* event2Trace: convert event times into a trace of {0, 1} (cf. trace2Event)
* trace2Event: converts a trace of {0,1} into event times (cf. Event2Trace)
* getMonthDateCh: list of datafile names, months, dates and channels saved under rootFolder
* getPSTH: convert spike times into a trace
* fitResponse: fit a 1D circular gaussian to single trial data of a given time

Other Analysis Scripts
* classifyUnits_pop: curates results of showTonsetResp
* fitPSTH_pop: curates results of fitPSTH_test.m from all neurons and compute statistics
* parea_spike_ms_test" curates results of parea_spike_ms from all neurons
* pupilOnsetsResp: spikes(recorded and predicted)  triggered by pupil dilation/constriction starts (TOBE MERGED W fitPSTH_test)
* selectFitPeriod_test: estimate kernel only during fOnset-tOnset (TO BE DELETED)

# Note

I don't test environments under Mac.

# Author

* Daisuke Shimaoka

# License

fitKernel is confidential.
