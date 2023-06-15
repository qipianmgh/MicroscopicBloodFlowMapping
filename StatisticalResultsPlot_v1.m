clear
clc
close all

load('Data_StatisticalResults.mat')

fg_100 = figure();
set(fg_100,'units','inches','position',[5,5,2.15,2]);
ax_acv = axes(fg_100, 'Position', [0.14 0.25 0.83 0.711]);
hold(ax_acv, 'on')

mean_allmice_art = mean([data.set0502.art.meanFlowV, data.set0508.art.meanFlowV, data.set0511.art.meanFlowV]);
STD_allmice_art = std([data.set0502.art.meanFlowV, data.set0508.art.meanFlowV, data.set0511.art.meanFlowV]);

pt_art = [data.set0502.art.meanFlowV, data.set0508.art.meanFlowV, data.set0511.art.meanFlowV, mean_allmice_art];
STD_art = [data.set0502.art.flowVSTD, data.set0508.art.flowVSTD, data.set0511.art.flowVSTD, STD_allmice_art];
errorbar(ax_acv, 0.85, pt_art(1), STD_art(1), 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r')
errorbar(ax_acv, 0.95, pt_art(2), STD_art(2), 'd', 'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g')
errorbar(ax_acv, 1.05, pt_art(3), STD_art(3), '^', 'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b')
errorbar(ax_acv, 1.15, pt_art(4), STD_art(4), 's', 'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')

mean_allmice_cap = mean([data.set0502.cap.meanFlowV, data.set0508.cap.meanFlowV, data.set0511.cap.meanFlowV]);
STD_allmice_cap = std([data.set0502.cap.meanFlowV, data.set0508.cap.meanFlowV, data.set0511.cap.meanFlowV]);

pt_cap = [data.set0502.cap.meanFlowV, data.set0508.cap.meanFlowV, data.set0511.cap.meanFlowV, mean_allmice_cap];
STD_cap = [data.set0502.cap.flowVSTD, data.set0508.cap.flowVSTD, data.set0511.cap.flowVSTD, STD_allmice_cap];
errorbar(ax_acv, 1.85, pt_cap(1), STD_cap(1), 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r')
errorbar(ax_acv, 1.95, pt_cap(2), STD_cap(2), 'd', 'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g')
errorbar(ax_acv, 2.05, pt_cap(3), STD_cap(3), '^', 'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b')
errorbar(ax_acv, 2.15, pt_cap(4), STD_cap(4), 's', 'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')

mean_allmice_vein = mean([data.set0502.vein.meanFlowV, data.set0508.vein.meanFlowV, data.set0511.vein.meanFlowV]);
STD_allmice_vein = std([data.set0502.vein.meanFlowV, data.set0508.vein.meanFlowV, data.set0511.vein.meanFlowV]);

pt_vein = [data.set0502.vein.meanFlowV, data.set0508.vein.meanFlowV, data.set0511.vein.meanFlowV, mean_allmice_vein];
STD_vein = [data.set0502.vein.flowVSTD, data.set0508.vein.flowVSTD, data.set0511.vein.flowVSTD, STD_allmice_vein];
errorbar(ax_acv, 2.85, pt_vein(1), STD_vein(1), 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r')
errorbar(ax_acv, 2.95, pt_vein(2), STD_vein(2), 'd', 'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g')
errorbar(ax_acv, 3.05, pt_vein(3), STD_vein(3), '^', 'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b')
errorbar(ax_acv, 3.15, pt_vein(4), STD_vein(4), 's', 'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')

h_legend = legend(ax_acv, 'Mouse #1', 'Mouse #2', 'Mouse #3','Average across mice');
set(h_legend,'position',[0.40 0.58 0.46 0.25],'Fontname','Calibri','FontSize',8,'Box','off');
xticks(ax_acv, [1,2,3]);
xticklabels(ax_acv, {'Arterioles', 'Capillaries', 'Venules'})
set(ax_acv,'Fontname','Calibri','FontSize',9);
h_ylabel = ylabel(ax_acv, 'Center-line flow velocity (mm/s)');
set(h_ylabel,'Fontname','Calibri','FontSize',9,'position',[0.55 5 -1.0000]);
xlim(ax_acv, [0.75, 3.25]);
ylim(ax_acv, [-1, 12]);
print(fg_100,'-dtiffn','-r600','Fig7_a.tiff'); 




%%%
%%- branching order (flow velocity) boxplot
pt_cap_art = [];
pt_cap_vein = [];
ptSTD_cap_art = [];
ptSTD_cap_vein = [];
O = [];


for ordr = 1:4
    %- from arteriol side
    flowV_art = [mean(abs(data.set0502.table(data.set0502.table(:,4) == 2 & data.set0502.table(:,5) == ordr, 2))),...
        mean(abs(data.set0508.table(data.set0508.table(:,4) == 2 & data.set0508.table(:,5) == ordr, 2))),...
        mean(abs(data.set0511.table(data.set0511.table(:,4) == 2 & data.set0511.table(:,5) == ordr, 2)))];
    flowVSTD_art = [std(abs(data.set0502.table(data.set0502.table(:,4) == 2 & data.set0502.table(:,5) == ordr, 2))),...
    std(abs(data.set0508.table(data.set0508.table(:,4) == 2 & data.set0508.table(:,5) == ordr, 2))),...
    std(abs(data.set0511.table(data.set0511.table(:,4) == 2 & data.set0511.table(:,5) == ordr, 2)))];
   
    flowV_art_mean_allmice = mean(flowV_art);
    flowV_art_STD_allmice = std(flowV_art);

    
    flowV_vein = [mean(abs(data.set0502.table(data.set0502.table(:,4) == 2 & data.set0502.table(:,6) == ordr, 2))),...
        mean(abs(data.set0508.table(data.set0508.table(:,4) == 2 & data.set0508.table(:,6) == ordr, 2))),...
        mean(abs(data.set0511.table(data.set0511.table(:,4) == 2 & data.set0511.table(:,6) == ordr, 2)))];
    flowVSTD_vein = [std(abs(data.set0502.table(data.set0502.table(:,4) == 2 & data.set0502.table(:,6) == ordr, 2))),...
        std(abs(data.set0508.table(data.set0508.table(:,4) == 2 & data.set0508.table(:,6) == ordr, 2))),...
        std(abs(data.set0511.table(data.set0511.table(:,4) == 2 & data.set0511.table(:,6) == ordr, 2)))];
    
    flowV_vein_mean_allmice = mean(flowV_vein);
    flowV_vein_STD_allmice = std(flowV_vein);    
    
    pt_cap_art = [pt_cap_art flowV_art flowV_art_mean_allmice];
    ptSTD_cap_art = [ptSTD_cap_art flowVSTD_art flowV_art_STD_allmice];
    pt_cap_vein = [pt_cap_vein flowV_vein flowV_vein_mean_allmice];
    ptSTD_cap_vein = [ptSTD_cap_vein flowVSTD_vein flowV_vein_STD_allmice];
    

    O = [O, ordr, ordr, ordr, ordr];  

end


fg_200 = figure();
set(fg_200,'units','inches','position',[5,5,2.15,2]);
ax_bror = axes(fg_200, 'Position', [0.14 0.25 0.83 0.711]);
hold(ax_bror, 'on');
for ordr = 1:4
    errorbar(ax_bror, ordr-0.15, pt_cap_art(4*(ordr-1)+1), ptSTD_cap_art(4*(ordr-1)+1), 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r')
    errorbar(ax_bror, ordr-0.05, pt_cap_art(4*(ordr-1)+2), ptSTD_cap_art(4*(ordr-1)+2), 'd', 'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g')
    errorbar(ax_bror, ordr+0.05, pt_cap_art(4*(ordr-1)+3), ptSTD_cap_art(4*(ordr-1)+3), '^', 'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b')
    errorbar(ax_bror, ordr+0.15, pt_cap_art(4*(ordr-1)+4), ptSTD_cap_art(4*(ordr-1)+4), 's', 'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
end

xticks(ax_bror, [1,2,3,4])
xticklabels(ax_bror, {'1', '2', '3', '4'})
set(ax_bror,'Fontname','Calibri','FontSize',9);
h_ylabel = ylabel(ax_bror, 'Capillary flow velocity (mm/s)');
set(h_ylabel,'Fontname','Calibri','FontSize',9,'position',[0.2 3.5 -1.00]);
h_xlabel = xlabel({'Branching order from precapillary';'arterioles'});
set(h_xlabel,'Fontname','Calibri','FontSize',9,'position',[2.5 -2.2 -1.00]);
xlim(ax_bror, [0.5, 4.5]);
ylim(ax_bror, [-1, 8]);
print(fg_200,'-dtiffn','-r600','Fig7_b.tiff'); 



fg_300 = figure();
set(fg_300,'units','inches','position',[5,5,2.15,2]);
ax_bror2 = axes(fg_300, 'Position', [0.14 0.25 0.83 0.711]);
hold(ax_bror2, 'on');
for ordr = 1:4
    errorbar(ax_bror2, ordr-0.15, pt_cap_vein(4*(ordr-1)+1), ptSTD_cap_vein(4*(ordr-1)+1), 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r')
    errorbar(ax_bror2, ordr-0.05, pt_cap_vein(4*(ordr-1)+2), ptSTD_cap_vein(4*(ordr-1)+2), 'd', 'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g')
    errorbar(ax_bror2, ordr+0.05, pt_cap_vein(4*(ordr-1)+3), ptSTD_cap_vein(4*(ordr-1)+3), '^', 'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b')
    errorbar(ax_bror2, ordr+0.15, pt_cap_vein(4*(ordr-1)+4), ptSTD_cap_vein(4*(ordr-1)+4), 's', 'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
end
xticks(ax_bror2, [1,2,3,4]);
xticklabels(ax_bror2, {'1', '2', '3', '4'});
set(ax_bror2,'Fontname','Calibri','FontSize',9);
h_ylabel = ylabel(ax_bror2, 'Capillary flow velocity (mm/s)');
set(h_ylabel,'Fontname','Calibri','FontSize',9,'position',[0.18 1.5 -1.0000]);
h_xlabel = xlabel({'Branching order from postcapillary';'venules'});
set(h_xlabel,'Fontname','Calibri','FontSize',9,'position',[2.5 -0.38 -1.0000]);
xlim(ax_bror2, [0.5, 4.5])
ylim(ax_bror2, [0,3])
print(fg_300,'-dtiffn','-r600','Fig7_c.tiff'); 



