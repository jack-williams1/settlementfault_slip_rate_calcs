%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displacement-time paths for Settlement Fault %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SEE SEPERATE SCRIPT FOR MAKING FIGURE PLOTS %%%

%Displacement-time paths created that pass through three known constraints
%on Settlement Fault activity. Displacement paths created for a single
%event displacement that follows a normal distribution

%time constraints from Holocene estuarine terrace, Owaka trench, and last interglacial terrace
%Ages with respect to 1950 CE (i.e. cal. yeas BP)
time=vertcat([0,0],[3654,190],[19700,1900],[125000,7000]);
%displacement constraints (vertical uplift) on the Settlement Fault for each time constraint
displacement=vertcat([0,0],[1.9,0.4],[2.5,0.2],[2,2]);

%assign a truncated normal distriburtion to represent the single event displacement
pd=makedist('normal',0.8,0.5); sevd_dist=truncate(pd,0.4,2.0);

%create various empty variables and set target number of valid paths
target_simu=500; stored_sevd={}; tmp_stored_sevd=cell(2,1); inc_sevd=0;
event_num=zeros(target_simu,height(displacement)-1);  
attempt_count=zeros(height(displacement)-1,1);
lq_slip_rate=zeros(target_simu,1); lq_net_slip_rate=zeros(target_simu,1);
st_slip_rate=zeros(target_simu,1); st_net_slip_rate=zeros(target_simu,1);

%Uniform distribution for Settlement Fault dip
f_dip=makedist('Uniform','Upper',60,'Lower',30);

%% Create displacement histories


for hh=2:height(displacement) %Run loop for each displacement constraint
    tmp_count1=0;
    count=1; tmp_count=0;
    
    disp_target_dist=makedist('Uniform','Upper',displacement(hh,1)+displacement(hh,2),...
        'Lower',displacement(hh,1)-displacement(hh,2));
    
    %Create displacement paths for that marker until sufficient number of fitted paths are created
    while count<=target_simu 
    
        tmp_count=tmp_count+1;
        
        kk=1; 

        if hh==2 %start new displacement path 
            inc_sevd=0; %if first event, displacement is 0
        else
            inc_sevd=max(tmp_stored_sevd{hh-2,count});%if later event, start with displacement from last event
        end   
        
        %randomly sample displacement target constraint
        disp_target=sort(random(disp_target_dist,[2,1]),'ascend');
        
        %if disp_target is already included by previous displacement
        if disp_target(2)<=inc_sevd
            tmp_stored_sevd{hh-1,count}=inc_sevd;
            tmp_count1=tmp_count1+1;
            count=count+1;%succesful path, exit first while loop
        
        else    
            %Successively add displacement until displacement exceeds dispacement marker max constraint
            while inc_sevd(kk)<disp_target(2)
            
                sevd=random(sevd_dist);%sample a random single event vertical displacement (sevd)
                inc_sevd(kk+1)=sevd+inc_sevd(kk);%add to incremental sevd in simulation
        
                %if inc_sevd within displacement range, store path
                if inc_sevd(kk+1)>disp_target(1) & inc_sevd(kk+1)<disp_target(2)
                    tmp_stored_sevd{hh-1,count}=inc_sevd;%store path
                    event_num(count,hh-1)=kk; %count number of events between displacements within each path
                
                    %force path to be independent of the one previously created
                    if count>1 & inc_sevd(2)==tmp_stored_sevd{hh-1,count-1}(end)
                        tmp_stored_sevd{hh-1,count-1}=tmp_stored_sevd{hh-1,count};
                        event_num(count-1,hh-1)=kk; 
                    else
                        count=count+1;%succesful path, exit first while loop
                    end
                end  
            kk=kk+1; 
            end
        end
 
    end %target number of paths developed, move onto next displacement constraint
    
    attempt_count(hh-1)=tmp_count;
    
end


% Create time histories

for hh=2:height(displacement) 
    
    %Time histories must lie between two constraints, which are randomly
    %sampled from dating constraints
    if hh==2 
        %force most recent displacement to be pre-European history (1840 or 110
        %cal. years BP)
        time_low=110; time_up=makedist('normal',time(hh,1),time(hh,2)/2);
        time_dist=makedist('Uniform','Lower',time_low,'Upper',random(time_up));
    elseif hh==3
        time_low=makedist('normal',time(hh-1,1),time(hh-1,2)/2); time_up=makedist('normal',time(hh,1),time(hh,2));
        time_dist=makedist('Uniform','Lower',random(time_low),'Upper',random(time_up));
    elseif hh==4
        time_low=makedist('normal',time(hh-1,1),time(hh-1,2)); time_up=makedist('Uniform','Lower',time(hh,1)-time(hh,2),'Upper',time(hh,1)+time(hh,2));
        time_dist=makedist('Uniform','Lower',random(time_low),'Upper',random(time_up));
    end

    for ii=1:target_simu

        tmp_len=length(tmp_stored_sevd{hh-1,ii});
        tmp_stored_age{hh-1,ii}=zeros(tmp_len,1);

        if hh==2
            %if tmp_len=2, 1 event between 0-3.5 ka
            %if tmp_len>2, >1 event between 0-3.5 ka
            tmp_stored_age{hh-1,ii}=zeros(tmp_len,1);
            tmp_stored_age{hh-1,ii}(1)=time_low;
          
             %if >1 event between 0-3.5 Ka, simulate some events timings in this period
             if tmp_len>2     
                tmp_stored_age{hh-1,ii}(2:end-1)=random(time_dist,[tmp_len-2,1]);
             end

             %enforce last event in this displacement step to be ~3.5 Ka
             tmp_stored_age{hh-1,ii}(end)=random(time_up); 
             tmp_stored_age{hh-1,ii}=sort(tmp_stored_age{hh-1,ii},'ascend');
            
        elseif hh==3
            %if at least 1 event between 3.5-20 Ka, simulate some event timings
             if tmp_len>1     
                tmp_stored_age{hh-1,ii}=random(time_dist,[tmp_len-1,1]);
             else %if no events, maintain date from previous state
                tmp_stored_age{hh-1,ii}=tmp_stored_age{hh-2,ii}(end); 
             end
             
              
        elseif hh==4

             if tmp_len==1 %no event occurred in time displacement history
                tmp_stored_age{hh-1,ii}=random(time_up);
             else
                tmp_stored_age{hh-1,ii}(1:end-1)=sort(random(time_dist,[tmp_len-1,1]),'ascend');
                tmp_stored_age{hh-1,ii}(end)=random(time_up); %force time displacement history to end ~125 ka
             end
        end

        tmp_stored_age{hh-1,ii}=tmp_stored_age{hh-1,ii}';
    end %end ii loop
end %end hh loop


%concatenate event displacements and times from each path
for ii=1:target_simu 
    stored_sevd{ii}=[unique(horzcat(tmp_stored_sevd{:,ii})),max(horzcat(tmp_stored_sevd{:,ii}))];
    stored_age{ii}=[horzcat(tmp_stored_age{:,ii})];
    
    total_disp=max(stored_sevd{ii})-stored_sevd{ii}(2);%Cumulative displacement minus MRE
    upper_age=max(stored_age{ii}); 
    total_age = upper_age-stored_age{ii}(2); %time period record by time-displacement history
    lq_slip_rate(ii)=(total_disp*1000)/total_age;%long term slip rate for each sequence

    st_disp=max(tmp_stored_sevd{2,ii})-stored_sevd{ii}(2);%Cumulative displacement since Owaka trench minus MRE
    st_age = max(tmp_stored_age{2,ii})-stored_age{ii}(2);
    st_slip_rate(ii)=(st_disp*1000)/st_age; %short term slip rate for each sequence
    
    flt_dip=random(f_dip);
    lq_net_slip_rate(ii)=lq_slip_rate(ii)/sin(deg2rad(flt_dip)); 
    st_net_slip_rate(ii)=st_slip_rate(ii)/sin(deg2rad(flt_dip)); 
end


%% Save analysis

save('slip_rate_results','time','displacement','stored_age','stored_sevd',...
    'tmp_stored_age','tmp_stored_sevd','st_slip_rate','lq_slip_rate',...
    'lq_net_slip_rate','st_net_slip_rate','event_num');

    
