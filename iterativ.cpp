#include "iterativ.h"
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;
using namespace Eigen;


// Constructeur de la matrice A, de B et x0
iterativ::iterativ(int N, int k, double eps, int choix)
// k est le nombre d'itaration
:_N(N), _k(k), _eps(eps), _choix(choix)
{
  _b.resize(_N);
  for (int i = 0; i < _N; i++)
  {
    _b(i)=1.;
  }
  _x0.resize(_N);
  srand(time(0));
  _x0.setRandom(_N);
  system (("mkdir -p ./" + to_string(_N)).c_str());
}

//Fonction construction matrice
MatrixXd iterativ::Matrice()
{
  _A.setZero(_N,_N);
  if (_choix==1)
  {
    _A.setIdentity();
    _A=0.1*_A; /// Ici alpha=0.1
    MatrixXd B;
    B.resize(_N,_N);
    srand(time(0));
    B.setRandom(_N,_N);
    for (int i = 0; i < _N; i++) // a optimiser !!!!
    {
      for (int j = 0; j < _N; j++)
      {
        B(i,j)=abs(B(i,j));
      }
    }
    _A=_A+B.transpose()*B;
  }
  else if (_choix==2) //Laplacien
  {
    int n(sqrt(_N));
    MatrixXd B, I;
    B.resize(n,n);
    I.resize(n,n);
    I.setIdentity();
    I=-1*I;
    SparseMatrix<double> BSparse(n,n);
    vector<Triplet<double> > triplets;
  for (int i=0; i<n; ++i)
  {
    triplets.push_back({i,i,4.});
    if (i > 0)
    triplets.push_back({i,i-1,-1.});
    if (i < n-1)
    triplets.push_back({i,i+1,-1.});
  }

  BSparse.setFromTriplets(triplets.begin(), triplets.end());
  B=BSparse;

  for( int i = 0; i < n; i++)
   {
    _A.block(i*n,i*n,n,n)=B; // démarre au coeff i*n, i*n et de taille n*n

    if (i<n-1)
    {
      _A.block(i*n,(i+1)*n,n,n)=I;
      _A.block((i+1)*n,i*n,n,n)=I;
    }
  }
  _A=_A/(pow(1./(_N+1),2));
  //le pas h=1/(N+1)
  }
  else if (_choix==3)
  {
    SparseMatrix<double> ASparse;
    VectorXd vecteur;//vecteur triplet ????
    ofstream bcsstk16;

    for (int i = 0; i < 4884; i++)
    {
       getline(bcsstk16, vecteur[i]);
       for(int j=0;j<vecteur.size();j++)
            tab[nLigne][i]=ligne[i];
    }
  }
  else //// On pourra l'enlever plus tard
  {cout << "Le choix de matrice n'est pas valable" << endl;}

  return _A;
}

//Gradient à pas optimal
VectorXd iterativ::GPO()
{
  //_A.resize(_N,_N); // on a peut etre pas besoin de cette ligne
  _A=Matrice(); //On construit la matrice
  _r.resize(_N);
  _r=_b-_A*_x0;
  ofstream mon_flux;
  mon_flux.open(to_string(_N)+ "/GPO.txt",ios::out);
  _z.resize(_N);
  _x.resize(_N);
  _x=_x0;
  double alpha(0);
  int i(0);

  while( (i<_k+1)&&(_r.norm()>_eps))
  {
    if(mon_flux)
    {
      mon_flux << i << "  "<< _r.norm() <<endl;
    }
    _z=_A*_r;
    alpha=(_r.dot(_r))/(_z.dot(_r)) ;// Produit scalaire
    _x=_x+ alpha*_r;
    _r=_r-alpha*_z;
    i++;
  }
  if (i==_k+1)
  {
    cout << "Tolérance non atteinte : "<< _r.norm() << endl;
  }
  return _x;
}

//Résidu Minimum
VectorXd iterativ::ResMin()
{
  //_A.resize(_N,_N);
  _A=Matrice();
  _r.resize(_N);
  _z.resize(_N);
  _x.resize(_N);
  double alpha(0);
  int i(0);
  _r=_b-_A*_x0;
  ofstream mon_flux;
  mon_flux.open(to_string(_N)+ "/ResMin.txt",ios::out);
  _x=_x0;

  while( (i<_k+1)&&(_r.norm()>_eps))
  {
    if(mon_flux)
    {
      mon_flux << i << "  "<< _r.norm() <<endl;
    }
    _z=_A*_r;
    alpha=(_r.dot(_z))/(_z.dot(_z)) ;// Produit scalaire
    _x=_x+ alpha*_r;
    _r=_r-alpha*_z;
    i++;
  }
  if (i>_k)
  {
    cout << "Tolérance non atteinte : "<< _r.norm() << endl;
  }
  return _x;
}

//Gradient Conjugué
VectorXd iterativ::GC()
{
  //_A.resize(_N,_N);
  _A=Matrice();
  _r.resize(_N);
  _p.resize(_N);
  _z.resize(_N);
  _x.resize(_N);
  double alpha(0);
  double beta(0);
  double gamma(0);
  int i(0);

  _r=_b-_A*_x0;
  ofstream mon_flux;
  mon_flux.open(to_string(_N)+ "/GC.txt",ios::out);
  _x=_x0;
  _p=_r;
  beta=_r.norm();
  while( (i<_k+1)&&(beta>_eps))
  {
    if(mon_flux)
    {
      mon_flux << i << "  "<< _r.norm() <<endl;
    }
    _z=_A*_p;
    alpha=(_r.dot(_r))/(_z.dot(_p));// Produit scalaire
    _x=_x+alpha*_p;
    VectorXd r(_r); //r=rk et _r=rk+1
    _r=_r-alpha*_z;
    gamma=(_r.dot(_r))/(r.dot(r));
    _p=_r+gamma*_p;
    beta=r.norm();
    i++;
  }
  if (i>_k)
  {
    cout << "Tolérance non atteinte : "<< _r.norm() << endl;
  }
    return _x;
}

//Gradient conjugué Eigen
VectorXd iterativ::GCEigen()
{
  _x.resize(_N); //_A.resize(_N,_N);
   _A=Matrice();
  SparseMatrix<double> ASparse(_A.rows(),_A.rows());

  vector<Triplet<double>> triplets;
  for (int i=0; i<ASparse.rows(); ++i)
  {
    for (int j=0; j<ASparse.rows(); j++)
    {
      triplets.push_back({i,j,_A(i,j)});
    }
	}
	ASparse.setFromTriplets(triplets.begin(), triplets.end());

	triplets.clear();

  ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
  cg.compute(ASparse);
  _x = cg.solve(_b);
/*
  cout << "matrice A GCEigen" << endl ;
  cout << ASparse << endl;
  cout << "vecteur b GCEigen" << endl ;
  cout << _b << endl;
*/
  return _x;

}

//Résidu minimum Eigen
/*VectorXd iterativ::ResMinEigen()
{
  _x.resize(_N); //_A.resize(_N,_N);
   _A=Matrice();
  SparseMatrix<double> ASparse(_A.rows(),_A.rows());

  vector<Triplet<double>> triplets;
  for (int i=0; i<ASparse.rows(); ++i)
  {
    for (int j=0; j<ASparse.rows(); j++)
    {
      triplets.push_back({i,j,_A(i,j)});
    }
	}
	ASparse.setFromTriplets(triplets.begin(), triplets.end());

	triplets.clear();

  GMRES<SparseMatrix<double>, Lower|Upper > solver;
  solver.compute(ASparse);
  _x = solver.solve(_b);

  return _x;

} //*/
