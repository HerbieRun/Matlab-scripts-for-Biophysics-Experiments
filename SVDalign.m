% find rigid boday transformation from M1 to M2
% R, optimal rotation matrix
% t, transformation

%% Construct the x,y,z position sets of M1 and M2

M1 = pdbread ('model1.pdb');
M2 = pdbread ('model2.pdb');

M1set = zeros (56:3); 
M2set = zeros (56:3); 
%Based on sequence alignment, we only compare 56 residues


for i=1:56
    M1set (i,1) = M1.Model.Atom(i+4).X; 
    %starting from the 5th residue in Model1
    M2set (i,1) = M2.Model.Atom(i).X;
    %starting from the 1st residue in Model2 
    i=i+1;
end


for j=1:56
    M1set (j,2) = M1.Model.Atom(j+4).Y; 
    M2set (j,2) = M2.Model.Atom(j).Y; 
    j=j+1;
end

for k=1:56
    M1set (k,3) = M1.Model.Atom(k+4).Z; 
    M2set (k,3) = M2.Model.Atom(k).Z; 
    k=k+1;
end

M1mean = mean(M1set);
M2mean = mean(M2set);

for i=1:56
     M1zm(i,:) = M1set(i,:) - M1mean;
     M2zm(i,:) = M2set(i,:) - M2mean;
end

C = M1zm'*M2zm;
[U,S,V] = svd(C);
D = diag(ones(1,3));
D(3,3) = det(U*V');
R = U*D*V';
t = M1set(1,:)'-R*M2set(1,:)';

%% Now R solved, apply the rotation to M2set

clear i;
NewM2 = zeros (56,3);
for i=1:56
    NewM2 (i,:) = R*M2set(i,:)'+t;
    i=i+1;
end

%% Creat a new pdb file called M3
M3 = M2;

clear i j k;

for i=1:56
    M3.Model.Atom(i).X = NewM2 (i,1); 
    i = i+1;
end


for j=1:56
    M3.Model.Atom(j).Y = NewM2 (j,2); 
    j = j+1;
end

for k=1:56
    M3.Model.Atom(k).Z = NewM2 (k,3); 
    k = k+1;
end

pdbwrite('newmodel2.pdb',M3);

%% RMSD 

clear i;
t = zeros(1,56);
for i=1:56
    t(i)= (norm(M1set(i,:)'-NewM2(i,:)'))^2;
    i = i+1;
end

rmsd = sqrt(sum(t)/56)
