N = 100;

C = cauchy_rnd(-1,1,2,N);
A = normrnd(0,1,2,N);

Test_Knn = [1, 2, 3, 2, 4, 3, 4; 1, 2, 3, 5, 3, 6, 8]

A = [Test_Knn(1,1), Test_Knn(2,1)]
B = [Test_Knn(1,2), Test_Knn(2,2)]
C = [Test_Knn(1,3), Test_Knn(2,3)]
D = [Test_Knn(1,4), Test_Knn(2,4)]
E = [Test_Knn(1,5), Test_Knn(2,5)]
F = [Test_Knn(1,6), Test_Knn(2,6)]
G = [Test_Knn(1,7), Test_Knn(2,7)]

function [e_dist] = distance(A, B)
  e_dist = sqrt((B(1,2) - A(1, 2))**2 + (A(1,1) - B(1, 1))**2);
endfunction

function y = mink(A,k)
  A = sort(A);
  y = A(1:k)
endfunction  

A_D = struct ("Dist_2_B", {{distance(B, A), 'B'}},"Dist_2_C", {{distance(C, A), 'C'}}, "Dist_2_D", {{distance(D, A), 'D'}}, "Dist_2_E", {{distance(E, A), 'E'}}, "Dist_2_F", {{distance(F, A), 'F'}}, "Dist_2_G", {{distance(G, A), 'G'}});
B_D = struct ("Dist_2_A", {{distance(A, B), 'A'}},"Dist_2_C", {{distance(C, B), 'C'}}, "Dist_2_D", {{distance(D, B), 'D'}}, "Dist_2_E", {{distance(E, B), 'E'}}, "Dist_2_F", {{distance(F, B), 'F'}}, "Dist_2_G", {{distance(G, B), 'G'}});
C_D = struct ("Dist_2_A", {{distance(A, C), 'A'}},"Dist_2_B", {{distance(C, B), 'B'}}, "Dist_2_D", {{distance(D, C), 'D'}}, "Dist_2_E", {{distance(E, C), 'E'}}, "Dist_2_F", {{distance(F, C), 'F'}}, "Dist_2_G", {{distance(G, C), 'G'}});
D_D = struct ("Dist_2_A", {{distance(A, D), 'A'}},"Dist_2_B", {{distance(B, D), 'B'}}, "Dist_2_C", {{distance(C, D), 'C'}}, "Dist_2_E", {{distance(E, D), 'E'}}, "Dist_2_F", {{distance(F, D), 'F'}}, "Dist_2_G", {{distance(G, D), 'G'}});
E_D = struct ("Dist_2_A", {{distance(A, E), 'A'}},"Dist_2_B", {{distance(B, E), 'B'}}, "Dist_2_C", {{distance(C, E), 'C'}}, "Dist_2_D", {{distance(D, E), 'D'}}, "Dist_2_F", {{distance(F, E), 'F'}}, "Dist_2_G", {{distance(G, E), 'G'}});
F_D = struct ("Dist_2_A", {{distance(A, F), 'A'}},"Dist_2_B", {{distance(B, F), 'B'}}, "Dist_2_C", {{distance(C, F), 'C'}}, "Dist_2_D", {{distance(D, F), 'D'}}, "Dist_2_E", {{distance(E, F), 'E'}}, "Dist_2_G", {{distance(G, F), 'G'}});
G_D = struct ("Dist_2_A", {{distance(A, G), 'A'}},"Dist_2_B", {{distance(B, G), 'B'}}, "Dist_2_C", {{distance(C, G), 'C'}}, "Dist_2_D", {{distance(D, G), 'D'}}, "Dist_2_E", {{distance(E, G), 'E'}}, "Dist_2_G", {{distance(F, G), 'F'}});


ANums = fieldnames(A_D);
A_C = struct2cell(A_D);
B_C = struct2cell(B_D);
C_C = struct2cell(C_D);
D_C = struct2cell(D_D);
E_C = struct2cell(E_D);
F_C = struct2cell(F_D);
G_C = struct2cell(G_D);

AEE = (A_C(:,1));
AE = cell2mat(AEE);
AE = AE(:)
BEE = (B_C(:,1));
BE = cell2mat(BEE);
BE = BE(:)
CEE = (C_C(:,1));
CE = cell2mat(CEE);
CE = CE(:)
DEE = (D_C(:,1));
DE = cell2mat(DEE);
DE = DE(:)
EEE = (E_C(:,1));
EE = cell2mat(EEE);
EE = EE(:)
FEE = (F_C(:,1));
FE = cell2mat(FEE);
FE = FE(:)
GEE = (G_C(:,1));
GE = cell2mat(GEE);
GE = GE(:)
                      
                      
%Min_A = mink(A_C(1, 2, :)(:), 3)
%Min_B = mink(B_C(1, 2, :)(:), 3)
%Min_C = mink(C_C(1, 2, :)(:), 3)
%Min_D = mink(D_C(1, 2, :)(:), 3)
%Min_E = mink(E_C(1, 2, :)(:), 3)
%Min_F = mink(F_C(1, 2, :)(:), 3)
%Min_G = mink(G_C(1, 2, :)(:), 3)



%figure(); hold on;
%plot(Test_Knn(1,:), Test_Knn(2,:), 'bo')
%plot(A(1,:), A(2,:), 'ro')
%plot(C(1,:), C(2,:), 'bx')