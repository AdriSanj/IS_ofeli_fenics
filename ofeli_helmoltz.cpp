//AUTHORS: Xose Garcia Quelle, Monica Loureiro Varela, Adrian Sanjurjo Garcia.

#include "OFELI.h"
#include "math.h"
#include "Electromagnetics.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
cout << "Example of using Ofeli with a Helmholtz problem" << endl;
real_t L = Pi; // Length of the bar.
size_t N = 10; // Number of elements.
double h = L/double(N);

// Building the mesh.
Mesh malla(L, N);

//Here we define the number of nodes.
size_t numero_nodos;      
numero_nodos = malla.getNbNodes();

// Making the vector at the RHS.
Vect<double> b(malla);
for (size_t i = 0 ; i < numero_nodos + 1 ; i++ ) {
    b(i) = 1;
      };
b *= h;

// Here we define the matrix in the LHS and its components, based in the FeNics example.
TrMatrix<double> K(numero_nodos);

Vect<double> u(malla);
double k = 1.5;   // Wave number.
double Ke_11, Ke_12, Me_11, Me_12;

// Elements of the matrix K.
Ke_11 = 1/h;
Ke_12 = -1/h;

// lements of the matrix M.
Me_11 = h/3;
Me_12 = h/6;

// We define here the elements of the matrix K.
double A1, A2;
A1 = Ke_11 - k*k*Me_11;
A2 = Ke_12 - k*k*Me_12;

// Boundary conditions.
K(1,1) = 1; K(1,2) = 0; b(1) = 0;
K(numero_nodos,numero_nodos) = 1; K(numero_nodos-1,numero_nodos) = 0; b(numero_nodos) = 0;

// Here we assign the elements in K.
for (size_t i = 2; i < numero_nodos; i++) {
  K(i,i  ) =  2*A1;
  K(i,i+1) = A2;
  K(i,i-1) = A2;
}

// Uncomment if you want to see the K matrix. 
//cout << K << endl;

K.Solve(b);

cout << "Solution of the problem: " << endl;
cout << "\nSolution:\n" << b;

return 0;
}
