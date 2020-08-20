function PlotDispersion(session,ManipName,Init,Final,Pas,NbFrame,Norm,P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the pair dispersion computed with the PairDisp.m function 
% / ! \ Before doing this algorithm you have to compute the pair dispersion
% with PairDisp.m. You can use the compile version of this function locate
% in /4d-ptv/Post-Processing/Pair_dispersion_EBernard/
%
%% Input:
% session       : the path where the directory of the experiment is 
% ManipName     : the name of your experiment
% Init          : the first initial separation you want to plot 
% Final         : the last initial separation you want to plot 
% Pas           : the difference between two initial separation 
% NbFrale       : the minimum size you used to compute the Dinitial
% Norm          : a boolean (true of false) you want to compensate the
% curve by the time at a certain power in order to see some regime 
% P             : the power of the compensation 
%% Output 
% A beautifull graph with all your pair dispersion in loglog scale 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Norm','var')
    Norm=false;
end

folderint = sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/PairDispFinal/',session.input_path,ManipName,NbFrame);



if Norm==true
    for i=1:((Final-Init)/Pas+1)
        filename=sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/PairDispFinal/pair_disp_%d.mat',session.input_path,ManipName,NbFrame,Init+(i-1)*Pas);
        load(filename);
        hold on;
        BinSize=Result.BinSize;
        MD=Result.MeanD;
        Time=Result.Time;
        Tstar=Result.Tstar;
        TimeNormTs=Time./Tstar;
        TimeNorm=Time.^P;
        MdNorm=MD./TimeNorm;
        loglog(TimeNormTs,MdNorm,'o','DisplayName',sprintf('Initial separation %d -%.2f mm',Init+Pas*(i-1),Init+Pas*(i-1)+BinSize));   
    end
    
%% Graphic parameters 
xlabel('t/t*');

ylabel(sprintf('<D_0^2>/t^{%d}',P));
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('location','best');
grid on;

else 
    for i=1:((Final-Init)/Pas+1)
        filename=sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/PairDispFinal/pair_disp_%d.mat',session.input_path,ManipName,NbFrame,Init+(i-1)*Pas);
        load(filename);
        hold on;
        BinSize=Result.BinSize;
        MD=Result.MeanD;
        Time=Result.Time;
        Tstar=Result.Tstar;
        TimeNormTs=Time./Tstar;
        loglog(TimeNormTs,MD,'o','DisplayName',sprintf('Initial separation %d - %.2f mm',Init+Pas*(i-1),Init+Pas*(i-1)+BinSize));   
    end

%% Graphic parameters 
xlabel('t/t*');

ylabel('<D_0^2>');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('location','best');
grid on;
end

end
