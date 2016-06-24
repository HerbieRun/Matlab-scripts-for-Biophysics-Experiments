%%Fist part, read the data

pts = 1000;
DataPos = csvread('points.csv');
DataChar = csvread('charges.csv');

%% Build an OctTree (Maxium 1 point allowed in each bin/node)

Mytree = OcTree(DataPos,'binCapacity',1);
BinNum = Mytree.BinCount; %how many bins in the tree

%% After a tree is built, we need to store the position of center charge and total charge to each bin

BinTable = zeros(BinNum,4); 

%BinTable is our final output in this part, it's a BinNum*4 dimension Table
%BinTable stores the x,y,z coordinate of each charge and the total charge
%in each bin number. 


for i=1:1:BinNum
    BinTable = TreeTraveler(i,BinTable,BinNum,Mytree,DataChar); %TreeTraveler function was separately built
end
 
%% Now we can start calculate the energy by traversal the tree

E = 0; %initial energy
Count = 0; % Count will keep counting how many times the energy calculation was done

for i=1:1:pts
        bin = 1; % always search from the root;
        [E,Count] = TreeEnergy(i,bin,E,Count,DataChar,BinTable,Mytree);
end

%% Output Section
E
Count