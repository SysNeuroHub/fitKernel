%2022 09September/16/4

plot(-onsets_cat.tOnset +onsets_cat.fixOnset, 1:215, 'ro');
hold on;
plot(-onsets_cat.tOnset +onsets_cat.fOnset, 1:215, 'bo');

plot(-onsets_cat.tOnset +onsets_cat.cueOnset, 1:215, 'go');

ylabel('trial id');
xlabel('time to tOnset');
legend('fixOnset','fOnset','cueOnset');
xlim([-5 0]);