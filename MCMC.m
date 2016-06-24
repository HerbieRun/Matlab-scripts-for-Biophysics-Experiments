%% Read the atom positions and extract all the coordinates
M = pdbread('methane.pdb');

% By prereading the pdb file, we know the 2nd atom was carbon 
% and it's at th origin of the coordinates.
C = [0,0,0];

% The 1,3,4,5 atom are Hydrogen
H1x = M.Model.HeterogenAtom(1).X;
H1y = M.Model.HeterogenAtom(1).Y;
H1z = M.Model.HeterogenAtom(1).Z;
H1 = [H1x,H1y,H1z];

H2x = M.Model.HeterogenAtom(3).X;
H2y = M.Model.HeterogenAtom(3).Y;
H2z = M.Model.HeterogenAtom(3).Z;
H2 = [H2x,H2y,H2z];

H3x = M.Model.HeterogenAtom(4).X;
H3y = M.Model.HeterogenAtom(4).Y;
H3z = M.Model.HeterogenAtom(4).Z;
H3 = [H3x,H3y,H3z];

H4x = M.Model.HeterogenAtom(5).X;
H4y = M.Model.HeterogenAtom(5).Y;
H4z = M.Model.HeterogenAtom(5).Z;
H4 = [H4x,H4y,H4z];

    % now, calculate energy of initial configuration
    kr = 367; ro = 1.08;
    BondSum = (norm(H1-C)-ro)^2 + (norm(H2-C)-ro)^2 + (norm(H3-C)-ro)^2 + (norm(H4-C)-ro)^2;
    BondSum = 0.5*kr*BondSum;
    
    Angle12 = acos(((H1-C)*(H2-C)')/(norm(H1-C)*norm(H2-C)));
    Angle12 = Angle12/pi*180;
    Angle13 = acos(((H1-C)*(H3-C)')/(norm(H1-C)*norm(H3-C)));
    Angle13 = Angle13/pi*180;
    Angle14 = acos(((H1-C)*(H4-C)')/(norm(H1-C)*norm(H4-C)));
    Angle14 = Angle14/pi*180;
    Angle23 = acos(((H2-C)*(H3-C)')/(norm(H2-C)*norm(H3-C)));
    Angle23 = Angle23/pi*180;
    Angle24 = acos(((H2-C)*(H4-C)')/(norm(H2-C)*norm(H4-C)));
    Angle24 = Angle24/pi*180;
    Angle34 = acos(((H3-C)*(H4-C)')/(norm(H3-C)*norm(H4-C)));
    Angle34 = Angle34/pi*180;
    
    ktheta = 35; theta0 = 109.5;
    AngleSum = (Angle12-theta0)^2 + (Angle13-theta0)^2 + (Angle14-theta0)^2 + (Angle23-theta0)^2 + (Angle24-theta0)^2 + (Angle34-theta0)^2;
    AngleSum = 0.5*ktheta*AngleSum;
    
    Eo = BondSum + AngleSum;
    Einitial = Eo; % Record Eo before any simulation

%% Simulation and Energy calculation
n = 10000000; % n, repeat times
Etrack = zeros(1,n); % Etrack will keep track of energy

for i=1:1:n
    
    % now, generate a random model 
    % Note that when generating a new model, the carbon position remains
    % constant at the origian
    % x,y,z are within (-1.5A to +1.5A) this gives a reasonable bond length
    
    H1n = -1.5+3*rand(1,3);
    H2n = -1.5+3*rand(1,3);
    H3n = -1.5+3*rand(1,3);
    H4n = -1.5+3*rand(1,3);
   
    
    % now, calculate energy of new configuration
    kr = 367; ro = 1.08;
    BondSum = (norm(H1n-C)-ro)^2 + (norm(H2n-C)-ro)^2 + (norm(H3n-C)-ro)^2 + (norm(H4n-C)-ro)^2;
    BondSum = 0.5*kr*BondSum;
    
    Angle12 = acos(((H1n-C)*(H2n-C)')/(norm(H1n-C)*norm(H2n-C)));
    Angle12 = Angle12/pi*180;
    Angle13 = acos(((H1n-C)*(H3n-C)')/(norm(H1n-C)*norm(H3n-C)));
    Angle13 = Angle13/pi*180;
    Angle14 = acos(((H1n-C)*(H4n-C)')/(norm(H1n-C)*norm(H4n-C)));
    Angle14 = Angle14/pi*180;
    Angle23 = acos(((H2n-C)*(H3n-C)')/(norm(H2n-C)*norm(H3n-C)));
    Angle23 = Angle23/pi*180;
    Angle24 = acos(((H2n-C)*(H4n-C)')/(norm(H2n-C)*norm(H4n-C)));
    Angle24 = Angle24/pi*180;
    Angle34 = acos(((H3n-C)*(H4n-C)')/(norm(H3n-C)*norm(H4n-C)));
    Angle34 = Angle34/pi*180;
    
    ktheta = 35; theta0 = 109.5;
    AngleSum = (Angle12-theta0)^2 + (Angle13-theta0)^2 + (Angle14-theta0)^2 + (Angle23-theta0)^2 + (Angle24-theta0)^2 + (Angle34-theta0)^2;
    AngleSum = 0.5*ktheta*AngleSum;
    
    Enew = BondSum + AngleSum; % Record the newly calculated energy 
    
    kT = 1;
    alpha = exp(-(Enew-Eo)/kT); % alpha, acceptance ratio
    
    if alpha >1,
        H1 = H1n; H2 = H2n; H3 = H3n; H4 = H4n; 
        Eo = Enew;
        
    else rnumber = rand;
        if alpha >= rnumber
            H1 = H1n; H2 = H2n; H3 = H3n; H4 = H4n; 
            Eo = Enew;
        else 
            % Eo = Eo;
        end
    end
    
    Etrack(i) = Eo; % Record the newly calculated energy each time
            
       
end

%% Final part, output, gives a plot and all the bond lengths and bond angles

plot(Etrack)

bond1 = norm(H1)
bond2 = norm(H2)
bond3 = norm(H3)
bond4 = norm(H4)

angle1R = acos((H1*H2')/(norm(H1)*norm(H2)));
angle1 = angle1R/pi*180

angle2R = acos((H1*H3')/(norm(H1)*norm(H3)));
angle2 = angle2R/pi*180

angle3R = acos((H1*H4')/(norm(H1)*norm(H4)));
angle3 = angle3R/pi*180

angle4R = acos((H2*H3')/(norm(H2)*norm(H3)));
angle4 = angle4R/pi*180

angle5R = acos((H2*H4')/(norm(H2)*norm(H4)));
angle5 = angle5R/pi*180

angle6R = acos((H3*H4')/(norm(H3)*norm(H4)));
angle6 = angle5R/pi*180





