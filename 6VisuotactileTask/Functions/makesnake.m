function makesnake(FRz_all, tempsms)
win = [-400 1000]; %-40 1000
sortwin = [50 500]; %50 500
iw = [find(tempsms>=win(1),1) find(tempsms>=win(2),1)]; %figure window
ii = [find(tempsms>=sortwin(1),1) find(tempsms>=sortwin(2),1)] ;       %sorting window
meanFRz_allw = nanmean(FRz_all(:,ii(1):ii(2)),2);
[~,sortmi] = sort(meanFRz_allw);
%nans are at the end
crap = sortmi((end-sum(isnan(meanFRz_allw))+1) :end);
sortmi((end-sum(isnan(meanFRz_allw))+1) :end) = []; sortmi = flipud(sortmi);
sortmi = [sortmi; crap];
figure()
imagesc(win,[1 size(FRz_all,1)],FRz_all(sortmi,iw(1):iw(2)))
cb = colorbar;
cb.Label.String = 'Z-scored firing rate';
caxis([-0.5 1])
line([0 0],[0 size(FRz_all,1)+0.5],'Color',[.9 .5 .4],'LineWidth',4)

end