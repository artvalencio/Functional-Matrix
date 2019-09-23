function [EEG,pulse,marker]=preprocessfunc(filename,chanlocs)
%PREPROCESSFUNC EEG pre-processing function for analysis with NeuroMat-HC experiment data
%Requires EEGLAB
%----------------------------------------
%Inputs:
%- filename: the *.edf EEG data file
%- chanlocs: the *.sfp file of the electrode positions
%----------------------------------------
%Outputs:
%- EEG:     the output pre-processed EEG
%- pulse:   the data of the pulse
%- marker:  the information of the instant when the visual stimulus is
%           applied
%---------------------------------------
%Usage example
%[EEG,pulse,marker]=preprocess('FB0127J9.edf','cap10-20.sfp')
%-----------------------------------------------------------------------------------
%  (CC-BY-4.0) Arthur Valencio
%  Institute of Computing, State University of Campinas
%  Research, Innovation and Dissemination Center for Neuromathematics (RIDC NeuroMat)
%  FAPESP fellowship #2018/0900-8. RIDC NeuroMat also supported by FAPESP #2013/07699-0

eeglab redraw;
%import EEG data
EEG=pop_biosig(filename,'channels',1:19);
pulse=pop_biosig(filename,'channels',27);
%get marker of visual stimulus
marker=pop_biosig(filename,'channels',36);
a=zeros(1,length(marker.data));
a(marker.data<-2e4)=1;
dur=duration(a);
a(dur<25*500)=2;
a(dur>40*500)=2;
marker.data=a;
clear a;
%trim data to the actual experiment
k=0;
time=[];
for i=1:length(marker.data)
    if and(i==1,marker.data(i)==2)
        k=k+1;
        time(k,1)=marker.times(i)/1000;
        idx=find(marker.data(i+1:end)~=marker.data(i),1,'first');
        if ~isempty(idx)
            time(k,2)=marker.times(idx)/1000;
        else
            time(k,2)=marker.times(end)/1000;
        end
    elseif and(and(i>1,marker.data(i)==2),marker.data(i)~=marker.data(i-1))
        k=k+1;
        time(k,1)=marker.times(i)/1000;
        idx=find(marker.data(i+1:end)~=marker.data(i),1,'first');
        if ~isempty(idx)
            time(k,2)=marker.times(idx)/1000;
        else
            time(k,2)=marker.times(end)/1000;
        end
    end
end
disp(time)
if ~isempty(time)
        EEG=pop_select(EEG,'notime',time);
        marker=pop_select(marker,'notime',time);
        pulse=pop_select(pulse,'notime',time);
end
%import electrode positions
EEG=pop_chanedit(EEG,'load',chanlocs);
%filter data
EEG=pop_eegfiltnew(EEG,1,[]);%low-pass: 1Hz
EEG=pop_eegfiltnew(EEG,[],50);%high-pass:50Hz
%re-reference
EEG=pop_reref(EEG,[]);



%extract epochs



[ALLEEG EEG]=eeg_store(EEG,EEG);
eeglab redraw;

end

function dur=duration(x)
%gets the duration of a visual marker pulse
for i=1:length(x)
    if i==1
         temp=x(i+1:length(x));
            idx=find(temp~=x(i),1,'first');
            if ~isempty(idx)
                dur(i)=idx;
            else
                dur(i)=length(temp);
            end 
    elseif i==length(x)
        if x(i)~=x(i-1)
            dur(i)=1;
        else
            dur(i)=dur(i-1);
        end
    else
        if x(i)==x(i-1)
            dur(i)=dur(i-1);
        else
            temp=x(i+1:length(x));
            idx=find(temp~=x(i),1,'first');
            if ~isempty(idx)
                dur(i)=idx;
            else
                dur(i)=length(temp);
            end 
        end
    end
end
end