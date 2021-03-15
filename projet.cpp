using namespace std;
#include "projet.h"
/*
PROGRAMME PRINCIPAL DU RECUIT SIMULE
Il y a deux possibilités :
- Lancer l'algorithme sur un ensemble de ville générées aléatoirement
- Lancer l'algorithme sur les 238 villes de l'énoncé
Pour passer de l'un à l'autre, il suffit de commenter/décommenter les parties indiquées le long de l'algorithme (3 enddroits différents)
// PARTIE ALEATOIRE pour le premier choix
// PARTIE ENONCE pour le deuxième
*/
int main()
{

    srand(time(NULL));
    // Définit un pointeur de fonction pour pouvoir choisir la fonction
    // de calcul de distances en fonction du type de ville
    double (*calcul)(vector<vector<double>> coord, vector<int> trajet);

    // PARTIE ENONCE
    double xscale = 960, yscale = 480;            // Facteurs de mise en forme pour la représentation sur la carte
    vector<vector<double>> coord = open_csv("dta.csv");   // Lecture des coordonnées des villes
  int N = coord[0].size();                      // On récupère le nombre de villes
    calcul = calcul_energie_km ;                // On prend la fonction correspondant au calcul avec longitude et lattitude
    // PARTIE ALEATOIRE
    /*
    double xscale = 180, yscale = 90;
    double xmax = 100, ymax = 100;            // Limites du domaine pour générer les villes aléatoires
    int N = 200                               // Nombre de villes à générer
    vector<vector<double>> coord = random_init(N, xmax, ymax);  // On génère les villes aléatoirement
    calcul = calcul_energie;                // On prend la fonction de distance euclidienne
    //FIN PARTIE ALEATOIRE
    */

    // On définit les variables nécessaires et paramètres du modèle
    double T0 = 10, C=0.01, q=0.6;
    double Tmin = 1, T= T0 ;
    double e_new, e_min, e_old;
    int k=1, n=1,i,j;

    // On définit les vecteurs qui nous servent à stocker nos résultats
    vector<int> trajet(N), best(N);
    vector<double> temp, energie, low_e;

    // On initialise le chemin avec les villes dans l'ordre canonique
    for (int i=0; i<N; i++)
      trajet[i] = i;

    // Calcul de la distance iniitale
    e_old = calcul(coord, trajet);
    e_min = e_old;

    // On remplit les vecteurs de température et d'éngergie
    temp.push_back(T);
    energie.push_back(e_old);
    low_e.push_back(e_min);

    cout << "Début de l'algorithme" << endl ;
    cout << "Distanc de départ : " << e_old << endl ;

    // On écrit les données dans un fichier pour pouvoir effectuer le graphique
    write(coord, trajet, "avant.txt", xscale, yscale);
    // PARTIE ENONCE
    plot_map("avant.txt", "avant.png");
    //  PARTIE ALEATOIRE
    //plot_random("avant.txt", "avant.png", xmax, ymax);
    // FIN PARTIE ALEATOIRE

    // On début le chrono pour connaître le temps d'exécution
    auto start = chrono::steady_clock::now();

    while (T0 > 5*Tmin) // Condition d'arrêt || A modifier en T > Tmin si on utilise le schéma de température par palier !!
    {
      // On tire deux indices aléatoires distincts
      i = rand()%N;
      j = rand()%N;
      if (i==j)
        continue;
      cout << "Itération " << n << "    ||    " ;

      // On permure le chemin et on calcule la nouvelle distance
      trajet = permutation(trajet, i, j);
      e_new = calcul(coord, trajet);

      // On se sert de la fonction de décision de Metropolis Hasting pour savoir si on conserve ou non la permutation
      if (metropolis(e_old, e_new, T))
      {
        e_old = e_new;   // Si True, alors on la conserve et on enregistre la nouvelle énergie
      }
      else
      {
        trajet = permutation(trajet, i,j);    // Sinon on repermute le chemin pour retourner au chemin précédent

      }
      // Système de température par palier (énoncé)
      /*
      if (actu_temp(n,k,C))
      {
        k = k+1;
        T = T0/k;
      }*/

      // Système de tempérautre avec remontée (plus efficace)
      if (T<Tmin)
      {
        // Remontée quand on atteint la température minimale
        T0*=q;
        T=T0;
        k=1;
      }
      else
      {
        // Décroissance exponentielle autrement
        T = T0*exp(-k/1e5);
        k=k+1;
      }

      n=n+1;

      // On sauvegarde toujours le meilleur résultat pour pouvoir obtenir le meilleur chemin à la fin
      if (e_old < e_min)
      {
        e_min = e_old;
        best = trajet;
      }

      // Affichage des résultats
      cout << "Score : " << e_old << "    || \t Temp : " << T << "    ||    " ;
      cout << "Meilleur score : " << e_min << endl ;

      // Sauvegarde des données de l'itéation
      temp.push_back(T);
      energie.push_back(e_old);
      low_e.push_back(e_min);
    }
    // Fin du chrono
    auto end = chrono::steady_clock::now();
    double total = double(chrono::duration_cast<chrono::seconds> (end-start).count());


    cout << "Temps total : " << total <<  endl ;
    cout << "Distance optimale : " << calcul(coord, best) << endl;

    // On écrit le meilleur chemin pour pouvoir le représenter visuellement
    write(coord, best, "apres.txt", xscale, yscale);
    // PARTIE ENONCE
    plot_map("apres.txt", "apres.png");
    // PARTIE ALEATOIRE
    //plot_random("apres.txt", "apres.png", xmax, ymax);
    // FIN PARTIE ALEATOIRE

    // On écrit les vecteurs de température et d'énergie pour pouvoir les tracer
    write_vec(temp, "temp.txt");
    write_vec(energie, "energie.txt");
    write_vec(low_e, "low_e.txt");
    plot_vec("temp.txt", "temp.png", true, false);
    plot_vec("energie.txt", "energie.png", true, true);
    plot_vec("low_e.txt", "low_e.png", true, true);

    // On suprime tous les fchiers temporaires non nécessaires
    system("rm *.txt");
    system("rm *.gnuplot");


}
