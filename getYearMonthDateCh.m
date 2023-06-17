function [loadNames, years, months, dates, channels] = getYearMonthDateCh(animal, rootFolder)
%FIXME

foldersByYear = dir(rootFolder);%FIX
foldersByYear = foldersByYear(3:end);
foldersByYear = foldersByYear([foldersByYear.isdir]);
foldersByYear = foldersByYear((cellfun(@numel, {foldersByYear.name}) == 4));
YearNames = {foldersByYear.name};

foldersByMonth = dir(rootFolder);%FIX
foldersByMonth = foldersByMonth(3:end);
foldersByMonth = foldersByMonth([foldersByMonth.isdir]);

monthNames = {foldersByMonth.name};
dateNames = cell(1,length(monthNames));
for im = 1:length(monthNames)
    foldersByDate = dir(fullfile(rootFolder, monthNames{im}));
    foldersByDate = foldersByDate(3:end);
    foldersByDate = foldersByDate([foldersByDate.isdir]);
    dateNames{im} = {foldersByDate.name};
end

idata = 1;
loadNames = [];
months = [];
dates = [];
channels = [];
for im = 1:length(monthNames)
    for id = 1:length(dateNames{im})
        date = [monthNames{im} '/' dateNames{im}{id}];
        
        loadFolder = fullfile(rootFolder, date, 'saved_oephysdata');
        
        fileNames = {dir(fullfile(loadFolder, [animal '_oephysdata_*.mat'])).name};
        
        for ich = 1:length(fileNames)
            loadName = fullfile(loadFolder, fileNames{ich});
            
            if exist(loadName, 'file')
                loadNames{idata} = loadName;
                months{idata} = monthNames{im};
                dates{idata} = dateNames{im}{id};
                chName_c =  textscan(fileNames{ich}, [animal '_oephysdata_ch%d.mat']);
                channels{idata} = chName_c{1};
                idata = idata + 1;
            end
        end
    end
end