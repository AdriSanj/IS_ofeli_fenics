#include "OFELI.h"
#include "math.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
cout << "Ejemplo de uso de Ofeli Helmholtz. " << endl;

// Construimos la malla en 1D.
real_t L = Pi; // L = Longitud de la barra.
size_t N = 10; // N = Numero de ELEMENTOS // Nodos = 11
double h = L/double(N);
Mesh malla(L, N);

double k = 1.5; // Numero de onda.

// Vector b y termino fuente.
Vect<double> b(malla);
size_t numero_nodos;
numero_nodos = malla.getNbNodes();
for (size_t i = 0 ; i < numero_nodos + 1 ; i++ ) {
    b(i) = 1;
      };
b *= h;
//cout << b << endl;

TrMatrix<double> K(numero_nodos); // Matriz de rigidez tridiagonal.

Vect<double> u(malla); // Vector u (solucion).

double Ke_11, Ke_12, Me_11, Me_12;

// Realizamos las contribuciones y ensamblaje de las matrices de rigidez:
  // Elementos para K.
Ke_11 = 1/h;
Ke_12 = -1/h;
  // Elementos para M.
Me_11 = h/3;
Me_12 = h/6;
  // Sumas para elementos:
double A1, A2;
A1 = Ke_11 - k*k*Me_11;
A2 = Ke_12 - k*k*Me_12;



// Condiciones de contorno.
K(1,1) = 1; K(1,2) = 0; b(1) = 0;
K(numero_nodos,numero_nodos) = 1; K(numero_nodos-1,numero_nodos) = 0; b(numero_nodos) = 0;

for (size_t i = 2; i < numero_nodos; i++) {
  K(i,i  ) =  2*A1;
  K(i,i+1) = A2;
  K(i,i-1) = A2;
}

//cout << K << endl;

// Resolvemos y printeamos la solucion.
K.Solve(b);
cout << "\nSolution:\n" << b;
return 0;
}
