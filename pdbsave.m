Mnew = M;
Mnew.Model.HeterogenAtom(1).X = H1(1);
Mnew.Model.HeterogenAtom(1).Y = H1(2);
Mnew.Model.HeterogenAtom(1).Z = H1(3);

Mnew.Model.HeterogenAtom(3).X = H2(1);
Mnew.Model.HeterogenAtom(3).Y = H2(2);
Mnew.Model.HeterogenAtom(3).Z = H2(3);

Mnew.Model.HeterogenAtom(4).X = H3(1);
Mnew.Model.HeterogenAtom(4).Y = H3(2);
Mnew.Model.HeterogenAtom(4).Z = H3(3);

Mnew.Model.HeterogenAtom(5).X = H4(1);
Mnew.Model.HeterogenAtom(5).Y = H4(2);
Mnew.Model.HeterogenAtom(5).Z = H4(3);

pdbwrite('newMethane.pdb',Mnew);