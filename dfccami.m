function [cami_xy,cami_yx,mutual_info,diridx,te_xy,te_yx] = ...
    dfccami(data,shift,windowlen,step,ini_part,L,units)
%DFCCAMI Retrieves the Dynamic Functional Connectivity matrix from EEG data
%using information-theoretical measures
%-----------------------------------------------------------------------------------
%Inputs
%- data:        the EEG data, given as a matrix, with each column being a
%               EEG channel
%- shift:       the maximum shift allowed for analysis considering 
%               time-delay between locations. Given in time-series points
%               units. Select shift=0 for disregard time-delay effects
%- windowlen:   the interval (window length) on which the infromation measure 
%               is calculated
%- step:        step (jump) given between calculation of correlation in one window
%               and on the other. Select step=1 for better precision,
%               step>1 for faster computation
%- ini_part:    initial partition to encode the data.
%               e.g.: ini_part=-1:0.5:1 generates 6 symbols, i.e. '0' for values
%               below -1, '1' for those in (-1,-0.5),'2' for (-0.5,0),'3'
%               for(0,0.5),'4' for (0.5,1) and '5' for values greater than
%               1.
%- L:           length of symbolic sequence to be considered. 
%               It is fixed that L_past(x)=L_past(y)=L_fut(y).
%               L=1, sequence is the initial partition. 
%               L=2, sequence is string with two points of the past leading
%               to two points of the future (e.g. '00' in x and '11' in
%               y->'01' in future of y), and so on
%-units:        'bits' or 'nats'. In case of typos, 'bits' are assumed as standard.      
%
%NOTE: TO USE CORRELATION MEASURES INSTEAD, USE DFC
%-----------------------------------------------------------------------------------
%Outputs
% Dynamic Functional Connectivity matrices of
%- cami_xy:     Causal Mutual Information in the direction X->Y
%- cami_yx:     Causal Mutual Information in the direction Y->X
%- mutual_info: Mutual Information shared between X and Y
%- diridx:      Directionality Index (direction of net flow of causal information),
%               positive for X->Y, negative for Y->X
%- te_xy:       Transfer Entropy in the direction X->Y
%- te_yx:       Transfer Etropy in the direction Y->X
%All the outputs are given in the form of structs where:
%- .value:      maximal absolute value of the information-theoretical measure
%               considering time-delay
%- .lag:        the time-delay where this maximal absolute value is found
%-----------------------------------------------------------------------------------
%Usage example
%[cami_xy,cami_yx,mutual_info,diridx,te_xy,te_yx] = dfccami(eeg,20,1000,10,-60:10:60,1,'bits')
%-----------------------------------------------------------------------------------
%  (CC-BY-4.0) Arthur Valencio
%  Institute of Computing, State University of Campinas
%  Research, Innovation and Dissemination Center for Neuromathematics (RIDC NeuroMat)
%  FAPESP fellowship #2018/0900-8. RIDC NeuroMat also supported by FAPESP #2013/07699-0
    
    %initializing
    n_chans=length(data(1,:));
    
    %find appropriate Taken's theorem tau of the system (zero autocorrelation)
    tau=findtau(data,n_chans);
    
    disp('Calculating for electrode pair x,y:');
    
    %calculating values and lags
    for i=1:n_chans
        for j=1:n_chans
            if i~=j
                disp(strcat(num2str(i),',',num2str(j)))
                selecttau=max(tau(i),tau(j));
                [cami_xy{i,j},cami_yx{i,j},mutual_info{i,j},diridx{i,j},te_xy{i,j},te_yx{i,j}]=...
                    slidingcorr(data(:,i),data(:,j),shift,windowlen,step,selecttau,ini_part,L,units);
            end
        end
    end
end

function [cami_xy,cami_yx,mutual_info,diridx,te_xy,te_yx] = ...
         slidingcorr(x,y,shift,windowlen,step,tau,ini_part,L,units)
    tslen=length(x);
    %calculates values:
    for shift2=-shift:shift
        if shift2>=0
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                [cami_xyval(k,shift2+shift+1),cami_yxval(k,shift2+shift+1),...
                    mutual_infoval(k,shift2+shift+1),diridxval(k,shift2+shift+1),...
                    te_xyval(k,shift2+shift+1),te_yxval(k,shift2+shift+1)]...
                   = serialmultithreadcami(x(i:i+windowlen),y(i+shift2:i+shift2+windowlen),...
                   L,L,ini_part,ini_part,tau,units);          
            end
        else
            k=0;
            for i=shift+1:step:tslen-windowlen-shift
                k=k+1;
                [cami_xyval(k,shift2+shift+1),cami_yxval(k,shift2+shift+1),...
                    mutual_infoval(k,shift2+shift+1),diridxval(k,shift2+shift+1),...
                    te_xyval(k,shift2+shift+1),te_yxval(k,shift2+shift+1)]...
                   = serialmultithreadcami(x(i+shift2:i+shift2+windowlen),y(i:i+windowlen),...
                   L,L,ini_part,ini_part,tau,units);  
            end
        end
    end
    %calculates peak of measures and respective lags:
    cami_xy=calcpeak(cami_xyval,shift);
    cami_yx=calcpeak(cami_yxval,shift);    
    mutual_info=calcpeak(mutual_infoval,shift);
    diridx=calcpeak(diridxval,shift);
    te_xy=calcpeak(te_xyval,shift);
    te_yx=calcpeak(te_yxval,shift);   
end

function out=calcpeak(x,shift)
    absx=abs(x);
    xlen=length(x(:,1));
    xval(1:xlen)=NaN;
    lag(1:xlen)=NaN;
    for i=1:xlen
        if ~isempty(find(absx(i,:),1)) %if measure is not totally zero
            %try
            [~,lag(i)]=getpeak(absx(i,:),-shift:shift);
            xval(i)=x(i,lag(i)+shift+1);
            %catch
            %    lag(i)=NaN;
            %    xval(i)=NaN;
            %end
        else %else measure is zero and lag undefined
            lag(i)=NaN;
            xval(i)=0;
        end
    end
    out.value=xval;
    out.lag=lag;
end

function [peak,lag]=getpeak(x,idx)
%find peak of absolute value of correlation within -shift:shift interval
    if idx==0
        peak=x;
        lag=0;
    else
        [peak,index]=max(x);
        lag=idx(index);
    end
end

function tau=findtau(x,n_chans)
    for chan=1:n_chans
        for i=1:1000
            r=corrcoef(x(1:end-i,chan),x(i+1:end,chan));
            corrval(i,chan)=r(1,2);
        end
        tau(chan)=find(corrval(:,chan)<=0,1);
    end
end