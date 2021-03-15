/*
Fichier header avec l'ensemble des fonctions définies pour le projet POOCNS du voyageur de commerce
Par Ayman Makki, Jean Hiouani et Baya Khiter
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <vector>

using namespace std;

//____________________________________
// FONCTIONS POUR LE CALCUL DES DISTANCES
//_______________________________________

// Conversion des angles en radian
double deg_to_rad(double deg)
{
  double rad = deg*(M_PI/180);
  return(rad);
}

// Calcul de la ditance entre 2 villes en km à partir de leurs coordonnées GPS
double km_dist(double lat1, double lon1, double lat2, double lon2)
{
  double R = 6371;
  double dlat = deg_to_rad(lat2-lat1);
  double dlon = deg_to_rad(lon2-lon1);
  double a = sin(dlat/2)*sin(dlat/2) + cos(deg_to_rad(lat1))*cos(deg_to_rad(lat2))*sin(dlon/2)*sin(dlon/2);
  double c = 2*atan(sqrt(a)/sqrt(1-a));
  double d = R*c;
  return(d);
}


// Calcul de la distance total pour un chemin particulier
double calcul_energie_km(vector<vector<double>> coord, vector<int> trajet)
{
  double energie =0;
  int n = trajet.size()-1;
  for (int i=0; i<n; i++)
  {
    // On calcule à chaque fois la distance à l'aide de la fonction précédente en indiquant les coordonnées de chaque ville
    energie += km_dist(coord[1][trajet[i]],coord[0][trajet[i]],coord[1][trajet[i+1]],coord[0][trajet[i+1]]);
  }
  // Rajout de la distance entre la dernière ville et la ville de départ
  energie += km_dist(coord[1][trajet[n]],coord[0][trajet[n]],coord[1][trajet[0]],coord[0][trajet[0]]);
  return(energie);
}

// Calcul de la distance lorsqu'on se sert de points générés aléatoirement
double calcul_energie(vector<vector<double>> coord, vector<int> trajet)
{
  double energie =0;
  int n = trajet.size()-1;
  for (int i=0; i<n; i++)
  {
    energie += pow((coord[0][trajet[i]] - coord[0][trajet[i+1]]),2) + pow((coord[1][trajet[i]] - coord[1][trajet[i+1]]),2);
  }
  energie += pow((coord[0][trajet[n]] - coord[0][trajet[0]]),2) + pow((coord[1][trajet[n]] - coord[1][trajet[0]]),2);
  return(pow(energie,0.5));
}

//____________________________________
// FONCTIONS POUR l'ALGORITHME DU RECUIT SIMULE
//_______________________________________


// Fonction qui permutent l'ensemble des points entre l'indice i et l'indice j sur un chemin
vector<int> permutation(vector<int> trajet, int i, int j)
{
  int imin = min(i,j);
  int imax = max(i,j)+1;
  // Définition d'un sous vecteur avec les indices à inverser
  vector<int> sub(&trajet[imin],&trajet[imax]);
  // Inversion du vecteur en utilisant les pointeurs
  reverse(sub.begin(),sub.end());
  // Réintégration du vecteur inverser dans le chemin
  copy(sub.begin(),sub.end(), trajet.begin()+imin);
  return trajet;
}

// Fonction de décision implémentant l'algorithme de Metropolis-Hasting
bool metropolis(double e_old, double e_new, double T)
{
  // Calcul du différentiel de distance
  double var_e = e_new-e_old;
  if(var_e < 0) // Si jamais le nouveau trajet est meilleur
  {
    return(true); // on garde la permutation
  }
  else // Autrement, on la garde avec une certaine probabilité
  {
    if ((double) rand() /RAND_MAX > exp(-(var_e)/T))
    {
      return(false); // on rejette la permutation
    }
    else
    {
      return(true); // on garde la permutation
    }
  }
}

// Fonction de décision pour le schéma de température décrit dans l'énoncé
bool actu_temp(int n, int k, double C)
{
  if ( n<= exp((C*(k+1))) )
  {
    return(false);
  }
  else
  {
    return(true);
  }
}

//______________________________________________
// FONCTIONS UTILITAIRES POUR MANIPULER LES DONNEES
//_______________________________________________

// Permet d'afficher un vecteur
void print(vector<int> p)
{
  for (int i=0; i<p.size(); i++)
  cout << p[i] << " " ;
  cout << endl;
}

// Permet d'écrire les coordonnées dans un fichier texte
// Les facteurs xscale et yscale permettent de mettre à l'échelle pour la représentation sur une carte
void write(vector<vector<double>> coord, vector<int> trajet, string filename, double xscale=180, double yscale=90)
{
  ofstream file(filename);
  int n = trajet.size();
  for (int i=0; i<n ; i++)
  {
    file << coord[0][trajet[i]]*xscale/180 << " " << coord[1][trajet[i]]*yscale/90<< endl;
  }
  // On rajoute le dernier trajet entre la dernière ville et le point de départ
  file << coord[0][trajet[0]]*xscale/180 << " " << coord[1][trajet[0]]*yscale/90<< endl;
  file.close();
}

// Permet d'écrire les données des vecteurs de température et d'énergie
void write_vec(vector<double> vec, string filename)
{
  ofstream file(filename);
  int n = vec.size();
  for (int i=0; i<n ; i++)
  {
    file << i << " " << vec[i]<< endl;
  }
  file.close();
}

//____________________________________
// FONCTIONS POUR TRACER LES GRAPHIQUES
//_______________________________________

// Permet de représenter un chemin sur la carte et de sauvegarder le résultat en png
void plot_map(string data, string output)
{
  ofstream file("plot.gnuplot");
  file << "set terminal png size 1600,1200" << endl;
  file << "set style line 1 lc rgb \"#ffffff\" lt 1 lw 2 pt 7 ps 1" << endl;
  file << "set pointintervalbox 3" << endl ;
  file << "set output \"" << output << "\"" << endl;
  file << "plot \"projection.jpg\" binary filetype=jpg center=(0,0) with rgbimage , \"" << data << "\" with linespoints ls 1" << endl ;
  file.close();
  system("gnuplot plot.gnuplot");
}

// Permet de représenter un chemin dans le cas de points tirés aléatoirement
void plot_random(string data, string output, double xmax, double ymax)
{
  ofstream file("plot.gnuplot");
  file << "set terminal png size 1600,1200" << endl;
  file << "set style line 1 lc rgb \"#000000\" lt 1 lw 2 pt 7 ps 1" << endl;
  file << "set pointintervalbox 3" << endl ;
  file << "set output \"" << output << "\"" << endl;
  file << "set xrange[-0.5:" << xmax+0.5 << "]" << endl;
  file << "set yrange[-0.5:" << ymax+0.5 << "]" << endl;
  file << "plot \"" << data << "\" with linespoints ls 1" << endl ;
  file.close();
  system("gnuplot plot.gnuplot");
}

// Permet de représenter les vecteurs d'évolution de la température et d'énergie
void plot_vec(string data, string output, bool logx=false, bool logy=false)
{
  ofstream file("plot.gnuplot");
  file << "set terminal png" << endl;
  if (logx || logy)
  {
    file << "set logscale ";
    if (logx)
      file << "x";
    if (logy)
      file << "y";
    file << endl;
  }
  file << "set output \"" << output << "\"" << endl;
  file << "plot \"" << data << "\" with lines" << endl;
  file.close();
  system("gnuplot plot.gnuplot");
}

// Permet de lire le fichier csv avec les coordonnées des villes et de les rajouter dans une matrice (2,238)
vector<vector<double>> open_csv(string filename)
{
  vector<vector<double>> data(2);
  string x,y;
  ifstream file(filename);
  if (!file.is_open())
    cout << "Fichier non trouvé" << endl ;

  while (file.good())
  {
    file >> x >> y; // lecture d'une ligne
    // rajout des coordonnées converties en 'double'
    data[1].push_back((stod(x)));
    data[0].push_back(stod(y));
  }
  return(data);
}

// Permet d'initialiser un vecteur aléatoire de coordonnées pour obtenir un échantillon
// de villes générées aléatoirement
vector<vector<double>> random_init(int n, int xmax, int ymax)
{
  vector<vector<double>> data(2, vector<double>(n));
  for (int i=0; i<n; i++)
  {
    data[0][i] = rand()%xmax;
    data[1][i] = rand()%ymax;
  }
  return(data);
}

bool check_perms(vector<int> perm_j, vector<int> perm_i)
{
  for (int i = 0; i< perm_j.size(); i++)
  {
    if (perm_i[i]==perm_j[i])
      return(false);
  }
  return(true);
}


// Fonction répérant l'algorithme principal afin de pouvoir tester différents paramètres
void toto(double T0, double Tmin, double C)
{
  srand(time(NULL));

  vector<vector<double>> coord = open_csv("dta.csv");

  int N = coord[0].size();
  int n=1;
  double e_new, e_min, e_old;
  double T = T0;
  int i,j;
  int k=1;
  vector<int> trajet(N);
  vector<double> temp;
  vector<double> energie;
  vector<double> low_e;
  vector<int> best(N);

  for (int i=0; i<N; i++)
  {
    trajet[i] = i;
  }

  e_old = calcul_energie_km(coord, trajet);
  e_min = e_old;

  temp.push_back(T);
  energie.push_back(e_old);
  low_e.push_back(e_min);

  cout << T0 << " " << Tmin << " " << C << endl;

  auto start = chrono::steady_clock::now();

  while (T>Tmin)
  {
    i = rand()%N;
    j = rand()%N;
    if (i==j)
      continue;

    trajet = permutation(trajet, i, j);
    e_new = calcul_energie_km(coord, trajet);

    if (metropolis(e_old, e_new, T))
    {
      e_old = e_new;
    }
    else
    {
      trajet = permutation(trajet, i,j);

    }

    // Schéma de température de l'énconcé
    if (actu_temp(n,k,C))
    {
      k = k+1;
      T = T0/k;
    }

    n=n+1;
    if (e_old < e_min)
    {
      e_min = e_old;
      best = trajet;
    }

    temp.push_back(T);
    energie.push_back(e_old);
    low_e.push_back(e_min);
  }
  auto end = chrono::steady_clock::now();
  double total = double(chrono::duration_cast<chrono::seconds> (end-start).count());
  cout << "total time : " << total <<  endl ;
  cout << calcul_energie_km(coord, best) << endl;
  cout << n << endl ;

}
