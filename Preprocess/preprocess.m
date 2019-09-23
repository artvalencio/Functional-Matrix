filename=input('EEG filename: ','s');
[EEG,pulse,marker]=preprocessfunc(strcat(filename,'.edf'),'cap10-20.sfp');
%run ica
EEG=pop_runica(EEG,'icatype','binica');
%select ica components
EEG=pop_selectcomps(EEG,1:19);
[ALLEEG EEG]=eeg_store(ALLEEG,EEG);
eeglab redraw;