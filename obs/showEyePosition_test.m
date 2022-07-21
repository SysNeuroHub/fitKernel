saveFolder = 'E:/tmp/cuesaccade_data';

animal = 'hugo';

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);
datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
disp(datech);

saveSuffix = [animal replace(datech,'/','_')];

thisDate = [months{idata} '_' dates{idata}];

loadName = fullfile(saveFolder, [saveSuffix '.mat']);
load(loadNames{idata},'dd');
load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']),'catEvTimes',...
    'blinks','outliers','eyeData_rmotl_cat');


%% detect SUCCESS and FAIL trials
onset = catEvTimes.tOnset;
invalidEvents = intersect(find(~isnan(onset)), find(dd.successTrials==0));
validEvents = intersect(find(~isnan(onset)), find(dd.successTrials==1));

[outcome, time] = getChoice(dd);
choiceTimes = time - time(2) + catEvTimes.cOnset(2);%hack

theseEvents = invalidEvents;
thisEye = eyeData_rmotl_cat;
for ii = 1:9
    ax(ii) = subplot(3,3,ii);
    trange = [catEvTimes.fOnset(theseEvents(ii)) choiceTimes(theseEvents(ii))];
    tidx = intersect(find(thisEye.t>trange(1)), find(thisEye.t<trange(2)));
    [~,tidxTOn] = min(abs(thisEye.t-onset(theseEvents(ii))));
    [~,tidxFOn] = min(abs(thisEye.t-catEvTimes.fOnset(theseEvents(ii))));
    %[~,tidxCOn] = min(abs(thisEye.t-catEvTimes.cOnset(theseEvents(ii))));
    [~,tidxCOn] = min(abs(thisEye.t - choiceTimes(theseEvents(ii))));
    
    %[~,tidxFixOn] = min(abs(thisEye.t-catEvTimes.fixOnset(theseEvents(ii))));
    plot(thisEye.x(tidx), thisEye.y(tidx));
    hold on
    plot(thisEye.x(tidxTOn), thisEye.y(tidxTOn), 'ro');
    plot(thisEye.x(tidxFOn), thisEye.y(tidxFOn), 'co');
    plot(thisEye.x(tidxCOn), thisEye.y(tidxCOn), 'mo');
    plot(5*cos(pi/180*dd.targetloc(theseEvents(ii))), 5*sin(pi/180*dd.targetloc(theseEvents(ii))),'k*')
    xlim([-6 6]);ylim([-6 6]);
    vline(0);hline(0);
    title(['trial: ' num2str(theseEvents(ii))]);
    axis square
end
xlabel('x [deg]');    ylabel('y [deg]');
legend('eye','tOnset','fOnset','(cOnset)','target location');

