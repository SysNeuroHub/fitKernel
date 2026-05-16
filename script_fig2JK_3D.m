fig3D_avg = showScatterTriplets3D(corr_tgt_avg_pop(2:4,:), param.predictorNames, [-0.4 1], selectedIDs,'linear',animalid_pop);
screen2png('fig2_3D',fig3D_avg);