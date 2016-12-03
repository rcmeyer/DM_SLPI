#include "iterativ.h"
#include <iostream>
#include<complex>
#include "Dense"
#include "Sparse"
using namespace Eigen;
using namespace std;

int main()
{
  cout << "Choisissez la valeur de k, le nombre d'itérations" << endl;
  int k;
  cin >> k;
  cout << "Choisissez la valeur de n, la taille de la matrice sera N=n²" << endl;
  int n;
  cin >> n;
  int N(pow(n,2));
  cout << "Choisissez la valeur de eps, la tolérance" << endl;
  double eps;
  cin >> eps;

  cout << "Quelle matrice A souhaitez-vous ?" << endl;
  cout << "1/ Matrice alphaI+BtB" << endl;
  cout << "2/ Matrice Laplacien" << endl;
  cout << "3/ BCSSTR16 Problème d'un barrage" <<endl;
  cout << "4/ NOS6 Problème de Poisson"<< endl;
  cout << "5/ NOS7 Equation diffusion en 3D"<< endl;
  cout << "6/ BCSSTR14 Structure d'un toit" << endl;
  int choix;
  cin >> choix;
  if ((choix !=1) && (choix!=2)&& (choix!=3)
  && (choix!=4) && (choix!=5) && (choix!=6) )
  {
    cout << "Le programme va s'arreter" << endl;
  return 0;
}
  iterativ iterativ(N,k,eps,choix);
  VectorXd x;
  x.resize(N);
/*
  cout << "Choisissez la méthode de résolution" << endl;
  cout << "1) Gradient à Pas Optimal " <<endl;
  cout << "2) Résidu minimum" <<endl;
  cout << "3) Gradient conjugué " << endl;
  int choix_method;
  cin >> choix_method;

 switch(choix_method)
  {
    case 1:
      x=iterativ.GPO();
      cout << "Avec le GPO x est"<< x << endl;
    break;

    case 2:
      x=iterativ.ResMin();
      cout << "Avec le Résidu Minimum x est"<< x << endl;
    break;

    case 3:
      x=iterativ.GC();
      cout << "Avec le GC x est"<< x << endl;
    break;

    default:
      cout  << "Ce choix n'est pas valable, le programme va s'arreter"<< endl;
      exit(0);
  }*/

  x=iterativ.GPO();
  cout << "Avec le GPO x est"<< x << endl;
  x=iterativ.ResMin();
  cout << "Avec le Résidu Minimum x est"<< x << endl;
  x=iterativ.GC();
  cout << "Avec le GC x est"<< x << endl;

 x=iterativ.GCEigen();
  cout << "x eigen est"<< x << endl;


    return 0;
}
