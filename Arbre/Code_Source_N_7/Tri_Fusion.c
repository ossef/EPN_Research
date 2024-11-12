#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


void fusionner(int T[10],int deb,int mil,int fin)
 {
        int *table1;
        int deb2=mil+1;
        int compt1=deb;
        int compt2=deb2;
        int i;
        
        //Tableau intermediaire, pour ne pas perdre les données de la partie1 (car on va "remplacer" les elements selon le bon placement)
        table1=malloc((mil-deb+1)*sizeof(int));

        //on recopie les éléments de la partie1 du tableau
        for(i=deb;i<=mil;i++)
            {
            table1[i-deb]=T[i];
            }
                        
        for(i=deb;i<=fin;i++)
            {        
            if (compt1==deb2) 	// L'algo s'arréte quand les element de table1 ont tous été placés, donc compt1 arrive au debut de la 2éme partie
                {
                break; 			//tous les éléments ont donc été classés
                }
            else if (compt2==(fin+1)) // les elements de la 2éme partie ont tous été placés et qu'il reste des elements de la partie1
                {
                T[i]=table1[compt1-deb]; //on ajoute les éléments restants de la partie1
                compt1++;                
                }
            else if (table1[compt1-deb]<T[compt2]) //Si l'element dans table1 est inferieur a l'element dans T
                {
                T[i]=table1[compt1-deb]; // Alors on prend l'element de table1
                compt1++;				   // On incremente le compteur de la partie1
                }
            else                     	     //Si l'element dans T est inferieur a l'element dans table1
                {
                T[i]=T[compt2];       // Alors on prend l'element de T
                compt2++;				  // On incremente le compteur le la partie2
                }  
            }
        free(table1);
 }
        

void tri_fusion_bis(int T[10],int deb,int fin)
        {										 // si (deb == fin) il n y'a rien a faire car il y'a qu'un seul element dans le tableau :)
        if (deb!=fin)						     // Sinon ...
            {
				int milieu=(fin+deb)/2;
				tri_fusion_bis(T,deb,milieu);   //decoupage a gauche
				tri_fusion_bis(T,milieu+1,fin);  // decoupage a droite
				fusionner(T,deb,milieu,fin);     //fusion des deux parties
            }
        }

void tri_fusion(int T[10], int n)
     {
												// Si le tableau contient ZERO elements, "fusionner" le vide avec lui méme ? hmmm ...
     if (n>0)								    // Sinon  
            {
				tri_fusion_bis(T,0,n-1);     // Pour declencher le processus de tri
            }
     }
     
int Alea_Int( int a, int b) // Generer un entier aleatoire dans [a,b]
{
    return (int)(a + (b - a) * ((double) rand() / RAND_MAX));
}
     

int main() {
	
	int i;
	srandom(getpid());
	
	int T[10];
	
	for (i=0; i<10 ; i++)
	{
		T[i] = Alea_Int(1,100);
	}
	
	printf("Avant tri : [");
	for (i=0; i<10 ; i++)
	{
		printf("%5d", T[i]);
	}
	printf("] \n");
	
	tri_fusion(T,10);
	
	printf("Apres tri : [");
	for (i=0; i<10 ; i++)
	{
		printf("%5d", T[i]);
	}
	printf("] \n");
	
		
	
	exit(0);
}
