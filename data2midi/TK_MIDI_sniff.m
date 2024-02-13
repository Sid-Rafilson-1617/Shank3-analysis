clear
% userd=getuserdir;%function I downloaded on matlab central
jj= (['\\F-moving-data\shnk3 (a)\Sniff']);%folder where your data are
cd(jj)
Fs=800;%sampling frequency
smoo=25;%smoothing in ms
smoo=round(1e3*smoo/Fs);%convert to samples
fn='NiDAQ_sniff.dat';%sniff analog recording
sniffpile=dir(fn);%finds files that match a condition, e.g., filename



sniff_locs=[];%defining an array to hold the inhalation times

sniff_file =sniffpile.name;
sniff_file = fopen(sniff_file,'r');
sniff=fread(sniff_file,'float64');
sniff=sniff(1:4:length(sniff));%analog sniff trace
sniff=smooth(sniff,smoo,'sgolay');% smoothing the signal with the sgolay method


figure %making a figure to display the test
hold on
tax=(1:length(sniff))./Fs./60;%defines an x axis in minutes

subplot(2,1,1)
plot(tax,sniff,'k-','LineWidth',1)

ylabel('Temp (AU)')

% yticks([])
% xlim([0 30])
ylim([prctile(sniff,1) prctile(sniff,99)])%sets the y-axis range 
title('Temp signal')
set(gcf,'Position',[400 400 1000 400])%controls position and size of the figure

if ~isempty(sniff)%make sure there is signal there

    %The temperature signal has slow drift, and that makes it harder to
    %peak-detect

    %break the signal into smaller chunks within which to normalize
    
    scanner=round(1:(5*Fs):length(sniff));% dividing the time series into 5 second windows
    i_locs=[];%;zeros(length(in_locs),1);
    e_locs=[];%zeros(length(in_locs),1);


    for scan=2:length(scanner)

        znindow=(scanner(scan-1)):scanner(scan);%positions within that window

        zniff=zscore(sniff(znindow));%z-scoring within that window



        %MinPeakProminence is important but I figured it out empirically,
        %this can differ for different signals, so definitely play with
        %that parameter for your particular signal
        [sniff_pks,in_locs] = findpeaks(zniff,'MinPeakDistance',smoo,'MinPeakProminence',0.5); %find inhalation points

        [ex_pks,ex_locs] = findpeaks(-zniff,'MinPeakDistance',smoo,'MinPeakProminence',0.5); %find exhalation points
        % %doesn't work great for the exhalations!
        if mod(scan,100)==2
            %this code lets you check if the peakfinding is ok
            figure
            hold on
            plot(zniff)
            plot(in_locs,zniff(in_locs),'ro')
            plot(e_locs,zniff(e_locs),'go')
            close
        end
        %in_locs are the inhalation times for that window
        for peax= 1:length(in_locs)
            this=in_locs(peax);

            xes=ex_locs(ex_locs>this);
            if ~isempty(xes)
                %                 if (xes(1)-this)<(1000)  && (xes(1)-this)>(50)

                %scanner(scan-1) is the start time of the window
                %this is the inhalation time within that window
                %add them to get the absolute time in the recording.
                i_locs=[i_locs; this+scanner(scan-1)];


%                 e_locs=[e_locs; xes+scanner(scan-1)];
            end
        end

    end





    subplot(2,1,2)
    hold on
    fsniff=(Fs./(diff(i_locs))); %get the sniff frequencies in per seconds

    llocs=i_locs(2:end)./Fs./60;%converts the inhalation times to minutes

    scatter(llocs,(fsniff),8, 'k','filled' )%plots inhalation times and instantaneous sniff frequencies 



    %     bwindist=imresize(windist,[size(windist,1) 100],'nearest');
    %     hold on
    %     imagesc(flipud(rot90(bwindist)))
    %             colormap(ax(1),colorcet('R2'))
    set(gcf,'Position',[400 400 1000 400])
    xlim([1 llocs(end)])
%     ylim([0 16])

    xlabel('Time (minutes)')
    title('Sniff dot plot')
    subplot(2,1,1)
    xlim([1 llocs(end)])

    %             sv=split(sv,'raw_sniffs');
    %             sv=sv{2};
%     print(gcf, '-dsvg',sv);

end
%%
scale=[2 3 2 2 3 2 3 2 2 3 2 3 2 2 3 2 3 2 2 3 2 3 2 2 3 2 3 2 2 3 2 3 2 2 3 2 3 2 2 3 ];
%inelegant way to make a pentatonic 
scale=cumsum(scale);

sniffsec=i_locs./Fs;%converts inhalation times to seconds

%midi events have 6 parameters, track, ____ , pitch, velocity, note on
%time, note off time
midsn=ones(length(i_locs),6);%this array holds these 6 parameters for each event
fsniff=ceil(Fs./(diff(i_locs)));%make the instantaneous sniff frequencies integers

nsniff=fsniff.*0;%an array for the pitch sequence
%in this loop, convert the frequency into a note that's in the scale.
for no=1:length(fsniff)
    nsniff(no)=scale(fsniff(no));


end

% fsniff(fsniff>12)=12;
%         fsniff(fsniff<3)=0;
%         fsniff(fsniff>=3 & fsniff<6)=5;
%         fsniff(fsniff>5 & fsniff<10)=7;
%         fsniff(fsniff>7)=12;


midsn(:,3)=[nsniff; nsniff(end)];%whole pitch sequence
midsn(:,4)=127;%sets all the velocities to max

midsn(:,5)=sniffsec-bounds(1)./Fs;%note on times
midsn(:,6)=midsn(:,5)+[diff(sniffsec)-0.025; 0.1];%note off times


midi_new = matrix2midi_bpm(midsn,5e2,960);%this code formats a matrix to convert to midi
midfile=[fn '_sniff_midi.mid'];
% s2mfile_dir = [userd '/Dropbox/MATLAB/FM_Ephys/spike2midi/'];

% cd(s2mfile_dir)
writemidi(midi_new, midfile);
