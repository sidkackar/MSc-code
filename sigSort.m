function [sig] = sigSort(obj, spikeTimes, clu, eventTimes, timeWinPre, timeWinPost)
        
        % sigSort takes the data block file (obj) as the minimum required input and outputs the IDs of neurons/clusters which are responsive to stimuli. It compares the 
        % neural activity after the stimulus to the baseline neuronal activity before the stimulus, in the time window of 1 sec for each as default. I have carried out 
        % a paired t-test between the neuronal activity on both sides to test for significance.
        
        % It is pointed out that methods based on the paired t test require defining a time window, and so, shuffle tests are preferred for this task. But, given that I
        % am taking a broad default time window of 1 sec both before and after the stimulus, is it enough to overcome this deficiency?
        
        % sig = output variable, containing IDs of neurons/clusters which are responsive to specified stimuli
        % obj = the processed experimental data, can be specified by date/day/session/experimental conditions
        % spikeTimes = all time points at which neuronal activity is measured
        % clu = contains IDs/amplitudes (not sure as I wrote this 2 years ago) of clusters
        % eventTimes = the time points of the event (eg. specific stimuli - auditory/visual/both etc.) the response to which we want to measure
        % timeWinPre/timeWinPost = the specified time window before (Pre) / after (Post) the specified event, each is later resolved into an upper and lower limit
        
        % Removing trials where the test subject made no response
        
        if strcmp(obj.blocks{1,1}.expDef, 'multiSpaceWorld')
            fltblk = prc.filtStruct(obj.blocks{1,1}, obj.blocks{1,1}.responseMade~=0); 
        else
            fltblk = obj.blocks;
        end
                
        % Assigning default values
                
        if ~exist('spikeTimes', 'var'); spikeTimes = fltblk.ephSpikeTimes; end
                
        if ~exist('clu', 'var'); clu = fltblk.ephSpikeTemplates; end
                
        if ~exist('eventTimes', 'var'); eventTimes = fltblk.stimPeriodStart; end
        
        % timePreLow/timePostLow is the lower limit of time window and timePreUp/timePostUp is the upper limit, set as 0-1 sec by default
                
        if ~exist('timeWinPre', 'var'); timePreLow = 0; timePreUp = 1;
        elseif length(timeWinPre) ~= 1; timePreLow = timeWinPre(1,1); timePreUp = timeWinPre(end);
        else
        timePreLow = 0; timePreUp = timeWinPre;
        end

        if ~exist('timeWinPost', 'var'); timePostLow = 0; timePostUp = 1;
        elseif length(timeWinPost) ~= 1; timePostLow = timeWinPost(1,1); timePostUp = timeWinPost(end);
        else
        timePostLow = 0; timePostUp = timeWinPost;
        end

        eventTimes = eventTimes(~isnan(eventTimes));
        
        %Selecting unique clusters
        
        cluList = unique(clu);
        sig = zeros(length(cluList),1);
        
        %Testing each cluster for significance of response
        
        for i = 1:length(cluList)
                
            % sigSpikes contains all the spike times of the cluster under consideration    
                
            sigSpikes = spikeTimes(clu == cluList(i));
            eventIndex = 1:length(eventTimes);

            % allSpikes contains all spike times in the specified time window, and is later split into preSpikes (before event) and postSpikes (after event)

            allSpikes = arrayfun(@(x) sigSpikes(sigSpikes >= (eventTimes(x)-timePreUp) & sigSpikes <= (eventTimes(x)+timePostUp))-eventTimes(x), eventIndex, 'uni',0);
            preSpikes = cellfun(@(x) x(x <= 0 - timePreLow), allSpikes, 'uni',0);
            postSpikes = cellfun(@(x) x(x >= 0 - timePostLow), allSpikes, 'uni',0);

            % preSpikeIdx and postSpikeIdx are a sort of measure of no. of spikes per unit time before and after the event, which are compared through a paired t-test to 
            % determine statistical significance of response
            
            preSpikeIdx = cellfun(@length, preSpikes)/(timePreUp-timePreLow);
            postSpikeIdx = cellfun(@length, postSpikes)/(timePostUp-timePostLow);
            z = ttest(postSpikeIdx, preSpikeIdx, 'alpha', 0.05/sqrt(length(cluList)))
            if z == 1; sig(i) = cluList(i); end    
        end
        
        %Storing and sorting significant clusters
        
        sig = sort(sig(sig~=0));
        end
