# fitKernel

"fitKernel" is a library of matlab functions and scripts to compute Generalized Linear Model (GLM) for cuesaccade data.


# Requirement

* marmodata (loading eye data)
* marmolab-stimuli (loading cuesaccade data)
* dsbox (ridge regression and some visualization)

Environments under linux and windows machines are tested.



# Usage

Functions about concatenating trial-based data
* concatenate_eye: mamordata.eye concatenated across trials
* concatenate_spk: spike times conctenated across trials
* decompose_eye: eye data decomposed back into trials from concatenated

Main function for kernel estimation
* fitPSTH: applies ridge regression a spkie trace to obtain linear kernel

Functions about preparing predictor signals
* detectOverlapEvents: detect overlap between the 2nd and 3rd inputs  (cf. trace2Event)
* detectPupilOnsets: start/end times of pupil dilation and constriction according to Einstein 2017 JNS
* getCueDirMtx: matrix of cue direction of all trials (NOT YET FUNCTIONING)
* getDirMtx: matrix of target direction of all trials (OBSOLETE?)
* getEyeDirMtx: matrix of eye position over time
* getEyeSpdDirMtx: matrix of eye velocity over time
* getPupilDiameter: converts pupil area into diameter and its percentile
* getSaccDir: classifies each saccade events into one of quantized bins of directions
* getSaccMtx: matrix of saccade direction over time
* getTgtDirMtx: matrix of target direction over time
* removeBlinksEDF: eye data after removing times of blinks detected by eyelink
* removePareaOutliers: eye data after removing times of outliers that crossed a specified threshold
* selectSaccades: saccade times after eliminating those at specified exclusion periods
* preparePredictors: matrix of predictors of specified modalities

Functions for analysis/visualization
* dirTuning: preferred direction of spikes at specified time window (TO BE DELETED)
* parea_spike_ms: computes correlation between spikes and pupil diameter at specified temporal frequencies
* pupilFigure: figure of single-trial traces of pupil position/area and spikes triggered at specified type of events 

Utility functions
* filtPSTH: (a)causally filtered spike trace
* alignMtxDir: a matrix which 2nd dimension is rotated based on its preferred direction
* event2Trace: convert event times into a trace of {0, 1} (cf. trace2Event)
* trace2Event: converts a trace of {0,1} into event times (cf. Event2Trace)
* getMonthDateCh: list of datafile names, months, dates and channels saved under rootFolder
* getPSTH: convert spike times into a trace

Scripts
* fitPSTH_test: main script to obtain linear kernel and event-triggered spiking traces (TO BE RENAMED)
* fitPSTH_pop: curate results of fitPSTH_test.m from all neurons and compute statistics
* parea_spike_ms_test" curate results of parea_spike_ms from all neurons
* pupilOnsetsResp: spikes(recorded and predicted)  triggered by pupil dilation/constriction starts (TOBE MERGED W fitPSTH_test)
* targetResp: spikes(recorded and predicted) triggered by target onsets (TO BE MERGED W fitPSTH_test)
* targetResp_avgFig: avg spikes across trials triggered by target onsets (TO BE REMOVED)
* selectFitPeriod_test: estimate kernel only during fOnset-tOnset (TO BE DELETED)

# Note

I don't test environments under Mac.

# Author

* Daisuke Shimaoka

# License

fitKernel is confidential.
