function loadName = getCuesaccadeName(rootFolder, ID)

aaa = split(ID,'/');
animal = aaa{1};
year = aaa{2};
months = aaa{3};
dates = aaa{4};
channels = aaa{5};

loadName = fullfile(rootFolder, year, 'cuesaccade_data', months, dates, ...
    'saved_oephysdata', [animal '_oephysdata_ch' channels '.mat']);
