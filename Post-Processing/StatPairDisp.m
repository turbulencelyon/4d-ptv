function [Dist,Index]=StatPairDisp(track,MinFrame,MaxFrame,range,seuil1,seuil2)
Dist=[];
Index=[];
for i=MinFrame:range:(MaxFrame-range)
    [I,D]=SelectInit(i,track,i+range-1,seuil1,seuil2);
    Dist=vertcat(Dist,D);
    Index=vertcat(Index,I);
end
end