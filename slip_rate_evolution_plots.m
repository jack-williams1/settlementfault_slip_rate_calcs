%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displacement-time plots for Settlement Fault %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SEE SEPERATE PLOT FOR RUNNING DISPLACEMENT-TIME SIMULATIONS

load slip_rate_results

%plot (a) single slip rate measurement example, and (b) 50 succesful paths
%use to recreate Fig 12 in manuscript
figure(1);

tiledlayout(2,1,'TileSpacing','Compact'); nexttile

indx=101;%select step sequence for example

ax_xlim=([0 134]); ax_ylim=([0 5]); fntsize=10.5;

p1=stairs(stored_age{indx}/1000,stored_sevd{indx},'Color',[0.7 0.9 0.7],'LineWidth',0.7); hold on

for ii=2:height(time)

    plot([(time(ii,1)-time(ii,2))/1000,(time(ii,1)+time(ii,2))/1000],[displacement(ii,1),displacement(ii,1)],'k-','LineWidth',1.2); hold on
    plot([time(ii,1)/1000,time(ii,1)/1000],[(displacement(ii,1)-displacement(ii,2)),(displacement(ii,1)+displacement(ii,2))],'k-','LineWidth',1.2); hold on
end

upper_age=max(stored_age{ii}); %recalc oldest event age for sequecne
%plot slip rate paths
p2=plot([upper_age/1000,stored_age{indx}(2)/1000],[max(stored_sevd{indx}),stored_sevd{indx}(2)],'r--','LineWidth',1.2);%plot long term slip rate path
p3=plot([max(tmp_stored_age{2,indx})/1000,stored_age{indx}(2)/1000],[max(tmp_stored_sevd{2,indx}),stored_sevd{indx}(2)],'b--','LineWidth',1.2);%plot short term slip rate path
xlim([ax_xlim(1) ax_xlim(2)]),ylim([ax_ylim(1) ax_ylim(2)]); pbaspect([2 1 1]);
xlabel('Time (ka)'); ylabel('Vertical Displacement (m)');

legend([p1 p2 p3],{'Simulated step sequence','late Quaternary uplift rate','short-term uplift rate'},'location','northeast','FontSize',fntsize-1.5)
set(gca,'FontSize',fntsize);

t1=title('(a)','fontsize',fntsize,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.1 h1(2) h1(3)]);

nexttile

%Plot 50x subsample displacement-time paths using stair function
for ii=1:50
    stairs(stored_age{ii}/1000,stored_sevd{ii},'Color',[0.7 0.9 0.7],'LineWidth',0.5); hold on
end

%Plot constraints
for ii=2:height(time)

    plot([(time(ii,1)-time(ii,2))/1000,(time(ii,1)+time(ii,2))/1000],[displacement(ii,1),displacement(ii,1)],'k-','LineWidth',1.5); hold on
    plot([time(ii,1)/1000,time(ii,1)/1000],[(displacement(ii,1)-displacement(ii,2)),(displacement(ii,1)+displacement(ii,2))],'k-','LineWidth',1.5); hold on
end

xlim([ax_xlim(1) ax_xlim(2)]),ylim([ax_ylim(1) ax_ylim(2)]); pbaspect([2 1 1])
xlabel('Time (ka)'); ylabel('Vertical Displacement (m)');

t1=title('(b)','fontsize',11,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.1 h1(2) h1(3)]);

set(gcf,'Position',[340 499 563 600]);

%Fit and plot Kernel distributions to uplift rate distributions
%used to recreate Fig S7 in manuscript
[ks_stsr,xi_stsr]=ksdensity(st_slip_rate); [ks_lqsr,xi_lqsr,]=ksdensity(lq_slip_rate);

%plot slip rate variability
figure(2);

p_stsr=plot(xi_stsr,ks_stsr,'b-','LineWidth',0.9); hold on %currently only plottig short term estimate
%uncomment to also plot slip rate variability of late Quaternary slip rate estimate
%p_lqsr=plot(xi_lqsr,ks_lqsr,'r-','LineWidth',0.9); 
%legend([p_lqsr,p_stsr],{['late Quaternary' char(10) 'uplift rate'],['short-term' char(10) 'uplift rate']},'FontSize',fntsize-2);
xlim([0 0.8]); xlabel('uplift rate (mm/yr)'); ylabel('f'); axis square
set(gca,'FontSize',fntsize);

%% Analyse mean and std uplift and net slip rate

median_lq_slip_rate=[median(lq_slip_rate),iqr(lq_slip_rate)]; %late quaternary uplift rate
median_st_slip_rate=[median(st_slip_rate),iqr(st_slip_rate)]; %short term uplift rate

median_net_lq_slip_rate=[median(lq_net_slip_rate),iqr(lq_net_slip_rate)]; %late quaternary uplift rate
median_net_st_slip_rate=[median(st_net_slip_rate),iqr(st_net_slip_rate)]; %short term uplift rate

%K-S test for whether the distributions are from the same parent distribution
[h,p]=kstest2(st_slip_rate,lq_slip_rate); %ks test for whether the slip rates are statistically distinct

%find total number of events between Owaka trench constraint and height of Catlins lake terrace 
lq_events=length(find(event_num(:,3)>0));

%derive slip rate variability (Cowie et al 2012) of step sequences
intersr=zeros(500,1); time_window=10000; time_bins=[0:time_window:110000]; 

for ii=1:500 %for each set of earthquake time-displacement paths
    
    tmp_sr=zeros(length(time_bins)-1,1);

    for kk=1:length(time_bins)-1
        
        if kk==1
            tmp_disp1=0;
        else
            tmp_indx1=find(stored_age{ii}>time_bins(kk),1,'first');%indx for displacement at start of window
            tmp_disp1=stored_sevd{ii}(tmp_indx1-1);%total displacement at start of time window
        end

            tmp_indx2=find(stored_age{ii}>time_bins(kk)+time_window,1,'first');%indx for displacement at end of window
            tmp_disp2=stored_sevd{ii}(tmp_indx2-1);%total displacement at end of time window
    
            tmp_sr(kk)=(tmp_disp2-tmp_disp1)/time_window;%time bin slip rate
    end
    
    mean_sr=stored_sevd{ii}(tmp_indx2-1)/stored_age{ii}(tmp_indx2);

    %Slip rate variability for 10000 year return window
    % for that time-displacement path (Cowie et al 2012)
    intersr(ii)=std(tmp_sr)/mean_sr;

end

intersr_stats=[mean(intersr),std(intersr)];

%% Plot all paths using heatmap [not used in manuscript]

figure (3);

%Plot heatmap for all paths (following 'makeHeatmap.m' from Hatem et al (2021))

X=[cell2mat(stored_age)'./1000,cell2mat(stored_sevd)'];
Xp=[ax_xlim(1):1:ax_xlim(2)]; Yp=[ax_ylim(1):0.25:ax_ylim(2)];%Set bin sizes for histogram
N = hist3(X,'Ctrs',{Xp, Yp}); %Create bivariate histogram

s1 = surf(Xp, Yp, N','edgecolor','none'); view([0 90]); caxis([0 7])

run orbl_cm_interp.m %specialty color map
colormap(orbl_interp)
xlim([ax_xlim(1) ax_xlim(2)]),ylim([ax_ylim(1) ax_ylim(2)]); pbaspect([2 1 1]);
xlabel('Time (ka)'); ylabel('Vertical Displacement (m)');
colorbar; h = colorbar; ylabel(h, 'Count in bin');

set(gcf,'Position',[340 499 563 298]);

%% How many sequences have two events younger than 1300 cal. yrs BP (consistent with Akatore) 

for kk=1:500
    tmp(kk)=length(find(stored_age{kk}<1300));
end

length(find(tmp>2))

