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
             
             
distance_matrix = [distance(B, A), distance(C, A), distance(D, A), distance(E, A), distance(F, A), distance(G, A);
                   distance(A, B), distance(C, B), distance(D, B), distance(E, B), distance(F, B), distance(G, B);
                   distance(A, C), distance(B, C), distance(D, C), distance(E, C), distance(F, C), distance(G, C);
                   distance(A, D), distance(B, D), distance(C, D), distance(E, D), distance(F, D), distance(G, D);
                   distance(A, E), distance(B, E), distance(C, E), distance(D, E), distance(F, E), distance(G, E);
                   distance(A, F), distance(B, F), distance(C, F), distance(D, F), distance(E, F), distance(G, F);
                   distance(A, G), distance(B, G), distance(C, G), distance(D, G), distance(E, G), distance(F, G)];
                   

Min_A = mink(distance_matrix,3,2)
figure(); hold on;
plot(Test_Knn(1,:), Test_Knn(2,:), 'bo')
%plot(A(1,:), A(2,:), 'ro')
%plot(C(1,:), C(2,:), 'bx')