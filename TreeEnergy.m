%This function will perform the energy traversal for a given point, pt1 
%travels throughout the tree from root to the leaves
%until s/d<theta is satisfied

function [E, Count] = TreeEnergy(point,bin,E,Count,DataChar,BinTable,Mytree)

pt1 = point;
PtChar1 = DataChar(pt1);
PtPos1 = Mytree.Points(pt1,:);
children = find(Mytree.BinParents==bin);
PtIndex = find(Mytree.PointBins == bin);
if isempty(children)==1 
    
    if isempty(PtIndex)~=1
       pt2  = PtIndex;
        
       if pt2 ~= pt1
% current bin is a leaf and the point in the bin is not the point we chose
        
        pt2 = PtIndex(1);
        PtChar2 = DataChar(pt2);
        PtPos2 = Mytree.Points(pt2,:);

        R12 = norm(PtPos1-PtPos2);
        E = E + (PtChar1*PtChar2./(2*R12));
        Count = Count + 1;
       end
    end
else
    PtChar2 = BinTable(bin,4);
    PtPos2 = BinTable(bin, 1:3);
    
    % Now we need to calculate s and d;
    BDY = Mytree.BinBoundaries(bin, :);
    BDYpt1 = [BDY(1),BDY(2),BDY(3)]; % 1,2,3 in the BinBoundery corrispond to the upperleft point coodinates
    BDYpt2 = [BDY(4),BDY(5),BDY(6)]; % 4,5,6 in the BinBoundery corrispond to the bottomright point coodinates
    
    s = norm(BDYpt1-BDYpt2); %maxium dimension is the diagonal dimension
    d = norm(PtPos1-PtPos2);
    theta = 0.01;
    
    if (s/d) < theta
        E = E + (PtChar1*PtChar2/(2*d));
        Count = Count + 1;
    else
        for i=1:1:8  % again, 8 children since it's an octree
            children = find(Mytree.BinParents==bin);
            child = children(i);
            [E, Count] = TreeEnergy(point,child,E,Count,DataChar,BinTable,Mytree);
        end 
    end
end


end