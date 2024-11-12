
Les feuilles : 1,2,3,4   
Les noeuds internes : 5,6   
La racine  : 7  

---------------------------- Parametres -----------------------------------------------

#define P     30    // Nombre de panneau, doit etre >= 7
#define alpha 2    // Tau d'EPs generé par panneau

#define N	7	   // Nombre de cellule dans le réseau
#define K	4	   // Nombre de Cox (Erlang dans ce cas) phase par cellule 

int Gama[N] = {5,5,5,5,5,5,5}; // Le taux de perte dans chaque file


int Lambda[N] = {1,2,3,4,0,0,0};  // Le taux d'arrivée des DP par dans chaque file
int Mu = 50 // Le taux de service identique dans chaque file, aussi identtique dans chaque phase (cas d'erlang) 

#pi = le nombre de panneaux affectés a la cellule i


 ---------------------------- Les 10 configurations selectionnées ---------------------------------


	  #p1       #p2       #p3		#p4		  #p5		#p6		  #p7       P          E[N]       E[T] 
	  
1 	2.00000   3.00000   4.00000   5.00000   4.00000   6.00000   6.00000   30.00000   15.93892   1.59389   <--- Configuration optimale

2  	3.00000   5.00000   4.00000   3.00000   4.00000   5.00000   6.00000   30.00000   19.00862   1.90086   

3   2.00000   5.00000   5.00000   5.00000   2.00000   5.00000   6.00000   30.00000   21.00173   2.10017   

4   6.00000   2.00000   2.00000   6.00000   2.00000   6.00000   6.00000   30.00000   24.12136   2.41214   

5	6.00000   6.00000   2.00000   3.00000   2.00000   5.00000   6.00000   30.00000   26.85180   2.68518   

6	2.00000   3.00000   5.00000   5.00000   5.00000   4.00000   6.00000   30.00000   39.39571   3.93957   

7 	3.00000   4.00000   5.00000   3.00000   5.00000   4.00000   6.00000   30.00000   41.00554   4.10055   

8   3.00000   4.00000   2.00000   6.00000   5.00000   4.00000   6.00000   30.00000   43.05624   4.30562   

9   5.00000   5.00000   2.00000   3.00000   5.00000   4.00000   6.00000   30.00000   45.02914   4.50291   

10  6.00000   6.00000   2.00000   4.00000   2.00000   4.00000   6.00000   30.00000   47.64287   4.76429   



