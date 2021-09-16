%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameData(WhichSimonBlock)
%Pull the data from the Game file for plotting activity during the Simon game

   GameDataDir = '*DataDirectory*\Data\GameData';
        GameDataFiles = what(GameDataDir);
        GameDataFiles = GameDataFiles.mat;
        UseTheseBlocks = [2 3 4 5 7 9 10 11 12 13];    %first and last blocks are rest, blocks 6 and 8 were corrupted  
        GameDataFiles = GameDataFiles(UseTheseBlocks);
  
   GameData_filename = GameDataFiles{WhichSimonBlock};

   GData = load(GameData_filename);
   
   GameSequences = GData.sData.params.game.targetsPerTrial;
   GameSequences = GameSequences(:,1:4);
   GameTrialTimes = GData.sData.results.trialTime;
   
   GameTrialOutcomes = GData.sData.results.trialSuccess;
   
   GameClock = datetime(datevec(GData.sData.gameClock));
   GameClock = datenum(GameClock);
   
   GamePosition = GData.sData.cursorPos;
  
   GameNewTrial = GData.sData.currentGameState;
   dGameNewTrial = diff(GameNewTrial);
   GameNewTrial(1) = [];
   NewTrialStarts = find((GameNewTrial == 4) & (dGameNewTrial == 1));
   GameNewTrial = NewTrialStarts;
   
   
    Xlim1 = datetime(datevec(GameClock(1)));
   Xlim2 = datetime(datevec(GameClock(end)));
   
 %
   clear PutativeGameNewTrialTimes PutativeGameNewTrialEndTimes GameNewTrialTimes GameNewTrialEndTimes

   %Identify putative trial times here; then use the movement in the
   %GamePosition vectors to more precisely define the trial start and stop
   %times as some period around the movement:
   
   %Identify the Max Amplitude in the GamePosition 2xN vector to standin as
   %a marker of movement:
   
   MaxGamePosition = sum(abs(GamePosition'));
   MaxGamePosition(MaxGamePosition < 1e-6) = 0;
   
   for i = 1:length(GameNewTrial)
   
       PutativeGameNewTrialTimes(i) = GameClock(GameNewTrial(i)-50);
      
       if (GameNewTrial(i)+200) <= length(GameClock)
            PutativeGameNewTrialEndTimes(i) = GameClock(GameNewTrial(i)+200);
       else 
           PutativeGameNewTrialEndTimes(i) = GameClock(end);
       end
   end
   %%% Now use those boundaries to define real bocks based on time before
   %%% (100 time steps) and after (50 time steps) movement of the cursor
   %%% (which is locked during target display)
     for i = 1:length(GameNewTrial)
   
        startIndex = find(GameClock == PutativeGameNewTrialTimes(i));
        endIndex = find(GameClock == PutativeGameNewTrialEndTimes(i));
        
        ConsiderPositionSegment = MaxGamePosition(startIndex:endIndex);
        
        Movement = find(ConsiderPositionSegment>0);
        MovementStart = Movement(1);
        MovementEnd = Movement(end);
        
       GameNewTrialTimes(i) = GameClock(startIndex + MovementStart - 100); 
       
       if length(GameClock) >= (startIndex+MovementEnd+50)
            GameNewTrialEndTimes(i) = GameClock(startIndex + MovementEnd + 50);
       else
            GameNewTrialEndTimes(i) = GameClock(end);
       end
   end
   

end

