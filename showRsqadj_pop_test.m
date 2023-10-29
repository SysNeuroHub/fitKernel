load('Z:\Shared\Daisuke\cuesaccade_data\figPSTH_pop20230717hugo\fitPSTH_pop20230717hugo.mat')


theseIDs = {'hugo/2021/08August/25/27',... %vision
    'hugo/2022/07July/26/19',... %eye speed
   'hugo/2022/08August/05/2'} %integrator OK

[~, selectedIDs] = intersect(id_pop, theseIDs);


plot(Rsqadj_pop(2,:)'./Rsqadj_pop(1,:)',Rsqadj_pop(3,:)'./Rsqadj_pop(1,:)','.');
hold on;
plot(Rsqadj_pop(2,selectedIDs)'./Rsqadj_pop(1,selectedIDs)',...
Rsqadj_pop(3,selectedIDs)'./Rsqadj_pop(1,selectedIDs)','o');
title('Rsq adjusted');
xlabel('omit eye speed/full mdl');
ylabel('omit eye position/full mdl');


