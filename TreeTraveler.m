%% This function calculate center charge position and total charge at a given node/bin

function [BinTable] = TreeTraveler(bin,BinTable,BinNum,Mytree,DataChar)

    Children = find(Mytree.BinParents == bin);
    ChildNum = size(Children,2);

    if ChildNum == 0 % is a leave
        ptIndex = find(Mytree.PointBins == bin);
        
        if isempty(ptIndex) ~= 1 % not an empty leave
            BinTable(bin, 1) = Mytree.Points(ptIndex, 1);
            BinTable(bin, 2) = Mytree.Points(ptIndex, 2);
            BinTable(bin, 3) = Mytree.Points(ptIndex, 3);
            BinTable(bin, 4) = DataChar(ptIndex);
        end
    else         
    
        Pos = [0,0,0]; % initial value of position
        Char = 0;      % initial value of charge
        
        for n=1:1:8 % we know all node has 0 or 8 children since it's an octree
          
          ChildrenR = find(Mytree.BinParents == bin);
          chIndex = ChildrenR(n);
          BinTable = TreeTraveler(chIndex,BinTable,BinNum,Mytree,DataChar);
          
          Char = Char + BinTable(chIndex, 4);
          Pos = Pos + BinTable(chIndex, 1:3).*BinTable(chIndex, 4);
          
        end
        
        if Char ~=0;
            BinTable(bin, 1) = Pos(1)/Char;
            BinTable(bin, 2) = Pos(2)/Char;
            BinTable(bin, 3) = Pos(3)/Char;
            BinTable(bin, 4) = Char;
        end
        
    end
end
    
  
    
 
