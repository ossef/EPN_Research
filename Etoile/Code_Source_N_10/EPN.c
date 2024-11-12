#include <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include  <math.h>

#define DEBUG 2

#define EPSILON 1e-3  //condition d'arret , cas du Rho iteratif
#define NITER 50	  //Stabilité pour 50 cas , cas du Rho  iteratif

#define P     40    // Nombre de panneau, doit être >= 7
#define alpha 2    // Tau d'EPs generé par panneau

#define N	10	   // Nombre de cellule dans le réseau
#define K	4	   // Nombre de Cox phase par cellule 

#define Infty 50
#define T_Infty Infty*Infty*Infty*Infty

//int Alpha[N] = {10,10,10,10,10,10,10};
int Gama[N] = {4,4,4,4,4,4,4,4,4,2}; 

double Rewards[T_Infty][13]; //stocker les rewards afin de les trier 

int Lambda[N] = {1,1,1,2,2,2,3,3,3,0};  //{4,5,6};
int Mu[N][K] = {{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{40,40,40,40},
				{100,100,100,100}};
				
double ProbaCox[N][K] =   {{1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0},
						  {1, 1, 1, 0}};
			
double ProbaRout[N][N] =  {  { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},  
							 { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},   
							 { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},   
							 { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},  
							 { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},  
							 { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},  
							 { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},      
						     { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},  
						     { 0,    0,    0,    0,    0,    0,    0,    0,    0,    1},  
						     { 0.01,   0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0} };  


typedef struct etat{
	int Y[K]; // l'etat (Y1,Y2,Y3,Y4)
   } Etat;

typedef struct tab_pi{
	Etat tab [N][T_Infty];
   } Tab_Pi;

Tab_Pi PI;


						 
/* - Par ordre de calul- */

double Beta[N][K];	  // 1
double RhoTheta[N];   // 2 - en inversant la matrice
double Theta[N];      // 3 - Utilise Beta(i,n)
double Rho[N]; 		  // 4 - Utilise RhoTheta et Theta

/* - Pour le caclul des pertes - */

double Loss1[N];
double Loss2[N];
double Loss[N];

/* - Pour l'heuristique 2 - */

int DInf[N];
int DSup[N];


double **inverse;
int Cst;

void Lecture_Matrice_Inverse() // Lecture de (I-P)^-1 où P est la matrice re routage
{
	FILE *f = fopen("Matrice_Routage.data","r");
	if(f == NULL)
	{
		printf("Problème : Le ficher 'Matrice_Routage.data' n'existe pas dans ce dossier ! \n");
		exit(0);
	}
	
	int Nelement;
	fscanf(f,"%d",&Nelement);
	if(Nelement != N)
	{
		printf("Problème : Le ficher 'Matrice_Routage.data' contient une matrice de taille %d, alors que ça doit être de taille %d ! \n",Nelement,N);
		exit(0);
	}
	
	int i,j;
	inverse =  (double **)malloc(N*sizeof(double));
		if( inverse == NULL) {printf("Pas assez de memoire de stockage !"); exit(0);}
	for(i=0 ; i<N ; i++)
	{
		inverse[i] = (double *)malloc(N*sizeof(double));
			if( inverse[i] == NULL) {printf("Pas assez de memoire de stockage !"); exit(0);}
	}
	for(i=0 ; i<N ; i++)
	{
		for(j=0 ; j<N ; j++)
		{
			fscanf(f,"%lf",&inverse[i][j]);
		}
	}
	
	if(DEBUG == 1)
	{
		printf("\nThe matrix (I - P)^-1 is : \n"); 
		   for (i = 0;i < N; i++)
			{
			 for (j = 0;j < N; j++)
			   {
				 printf("\t%f", inverse[i][j]);
			   }
				printf("\n");
			}
	}
	
}

void Init_Tabs2()
{
	int  j, i1, i2, i3, i4, c;
	c = 0;	
		for (i1=0; i1<Infty ; i1++) 
		{
			for (i2=0; i2<Infty ; i2++) 
			{
				for (i3=0; i3<Infty ; i3++) 
				{
					for (i4=0; i4<Infty ; i4++) 
					{
						for (j=0; j<N ; j++) 
						{
							PI.tab[j][c].Y[0] = i1;
							PI.tab[j][c].Y[1] = i2;
							PI.tab[j][c].Y[2] = i3;
							PI.tab[j][c].Y[3] = i4;
						}
						c++;
					}
				}
			}
		}
}

void Init_Tabs1()
{
	int i, n;
	for ( i=0; i<N ; i++) 
	{
		Rho[i] = 0;
		Theta[i] = 0;
		RhoTheta[i] = 0;
		Loss1[i] = 0;
		Loss2[i] = 0;
		Loss[i] = 0;
	}
	
	for ( i=0; i<N ; i++) 
	{
		for ( n=0; n<K ; n++) 
			{
				Beta[i][n] = 0;
			}
	}
	
}

long int factorielle(long int n)
{
	long int s = 1;
	int i;
	
	if (n >= 2)
	{
		for (  i=1; i<=n ; i++) 
		{
				s *= i;
		}
	}
	return s;
}


void Calcul_Beta(int Alpha[N])
{
	int i,n;
	for ( i=0; i<N ; i++) 
	{
		Beta[i][0] = (Alpha[i]*1.0)/(Gama[i] + Mu[i][0]); 
	}
	
	for ( i=0; i<N ; i++) 
	{
		for ( n=1; n<K ; n++)
		{
			Beta[i][n] = (Mu[i][n-1]* ProbaCox[i][n-1]* Beta[i][n-1]*1.0)/Mu[i][n];
			
		}
	}
}


void Calcul_Rho()
{
	int i,j,n;
	
	for (i = 0;i < N; i++)
			{
				 for (j = 0;j < N; j++)
				 {
					 RhoTheta[i] += Lambda[j]*inverse[j][i]*1.0;
				 }
			}
	
	if(DEBUG == 1)
	{
		printf("RHOTHETA :");
		for (i = 0;i < N; i++)
			printf("%f ", RhoTheta[i]);
		printf("\n");
	}
	
	
	for ( i=0; i<N ; i++) 
	{
		for ( n=0; n<K ; n++)
		{
			Theta[i] += (1-ProbaCox[i][n])*Mu[i][n]*Beta[i][n]*1.0;
		}
	}
	
	if(DEBUG == 1)
	{
		printf("THETA :");
		for (i = 0;i < N; i++)
			printf("%f ", Theta[i]);
		printf("\n");
	}
	
	for ( i=0; i<N ; i++) 
	{
		Rho[i] = RhoTheta[i]/Theta[i];
	}
	
	if(DEBUG == 1)
	{
		printf("Rho :");
		for (i = 0;i < N; i++)
			printf("%.10lf ", Rho[i]);
		printf("\n");
	}
	
}

void Affiche_Beta()
{
 
 if( DEBUG == 1 )
 {
	 int i,n;
	 printf("Affichage Beta(i,n) \n");
	 for ( i=0; i<N ; i++) 
		{
			printf("%d: ",i+1);
			for ( n=0; n<K ; n++) 
				{
					printf("%f  ",Beta[i][n]);
				}
			printf("\n"); 
		}
 }
}

void Affiche_Rho()
{
	if(DEBUG == 1)
	{
		int i;
		printf("THETA :");
		for (i = 0;i < N; i++)
			printf("%f ", Theta[i]);
		printf("\n");
		
		printf("Rho :");
		for (i = 0;i < N; i++)
			printf("%.10lf ", Rho[i]);
		printf("\n");
	}
}

int Test_Stabilite()
{
	int i,n;
	double s;
	int bool1 = 0;
	int bool2 = 0;
	
	for ( i=0; i<N ; i++) 
	{
		s = 0;
		for ( n=0; n<K ; n++)
		{
			s += Beta[i][n];
		}
		if( s >= 1 )
		{
			if(DEBUG == 1) {printf("Instabilité au niveau de Beta(%d,*) ! \n",i+1);}
			bool1 = 1;
		}
	}
	if(bool1 == 0)
		if(DEBUG == 1) {printf("Stabilité des Beta(i,*) OK !\n");}
	for ( i=0; i<N ; i++) 
	{
		if( Rho[i] >= 1 )
		{
			if(DEBUG == 1) {printf("Instabilité au niveau de Rho(%d) ! \n",i+1);}
			bool2 = 1;
		}
	}
	if(bool2 == 0)
		if(DEBUG == 1) {printf("Stabilité des Rho(i) OK !\n");}
	if(bool1 == 1 || bool2 == 1)
		return 1; //Instable
	return 0;     //Stable		
}

double Calcul_N_Number() /* -----  Lemma 4.2 -------- E[X] -------*/
{
	int i,j;
	double s2 , s3;
	
	
	s3 = 0; // Totale waiting time
	for ( i=0; i<N ; i++) 
	{
		s2 = 0; // Somme intermediare
		for ( j=0; j<N ; j++) 
		{
			s2 += Theta[j]*Rho[j]*ProbaRout[j][i]; 
		}
		
		s3 += ( ( Lambda[i] + s2 ) / ( Theta[i] - (Lambda[i] + s2) ) );
	}
	
	if(DEBUG == 1)
	{
		printf("Mean number of DPs in the EPN : %.10lf \n",s3);
	}
	return s3;
}

double Calcul_T_WaitingTime() /* ----  lemma 4.4 ----- E[T] --------*/
{
	int i;
	double s1, s2;
	
	s2 = Calcul_N_Number();

	s1 = 0; // Les arrivées externes
	for ( i=0; i<N ; i++) 
	{
		s1 += Lambda[i]*1.0; 
	}
		
	s2 = s2/s1 ; 
	if(DEBUG == 1)
	{
		printf("Mean Waiting time of a DP in the EPN : %.10lf \n",s2);
	}
	return s2;
}


int norme(int i, int c)
{
	int s = 0;
	for (int n=0; n<K ; n++) 
		{
			s += PI.tab[i][c].Y[n];
		}
	return s;
}


double Calcul_Proba_Etat(int i, int c ) // Pi (Y1,Y2,Y3)
{
	int n;
	double s1 , p1;
	
	s1 = 0;
	for (n=0; n<K ; n++) 
		{
			s1 += Beta[i][n];
		}
	
	p1 = (1-s1)*1.0;
	
	for (n=0; n<K ; n++) 
		{
			p1 *= ( pow(Beta[i][n],PI.tab[i][c].Y[n]) / (factorielle(PI.tab[i][c].Y[n])*1.0) );
		}
	
	return p1;
}

double Calcule_reward1(int i, int c)
{
	double s1 = 0;
	for (int n=0; n<K ; n++) 
		{
			s1 += (1 - ProbaCox[i][n])*Mu[i][n]*PI.tab[i][c].Y[n]*1.0;
		}
	return s1;	
}

double Calcule_reward2(int i, int c)
{
	return 	PI.tab[i][c].Y[0]*Gama[i]*1.0;
}

void Calcul_Loss1_Loss2_All()
{
	int i,c,x;
	double s = 0;
	
	for (i=0; i<N ; i++)  
	{		
		for (c=0; c<T_Infty ; c++)
		{
			if( norme(i,c) > 0 )
			{
				Loss1[i] += ( Calcul_Proba_Etat(i,c) * Calcule_reward1(i,c) ) ;
				Loss2[i] += ( Calcul_Proba_Etat(i,c) * Calcule_reward2(i,c) ) ;
			}	
		}
		Loss1[i] *= (1 - Rho[i]);
		
		for (x=0; x<Infty ; x++)
		{
			s += (1 - Rho[i]) * pow(Rho[i],x) * 1.0 ;
		}
		
		Loss2[i] *= s;
	}
	
	if(DEBUG == 1)
	{
		printf("Loss1 :");
		for (i = 0;i < N; i++)
			printf("%.10lf ", Loss1[i]);
		printf("\n");

		printf("Loss2 :");
		for (i = 0;i < N; i++)
			printf("%.10lf ", Loss2[i]);
		printf("\n");
	}
}

void Calcul_Loss_All()
{
	Calcul_Loss1_Loss2_All();
	
	int i;
	for (i=0; i<N ; i++)  
	{
		Loss[i] = Loss1[i] + Loss2[i];
	}
	
	if(DEBUG == 1)
	{
		printf("Loss :");
		for (i = 0;i < N; i++)
			printf("%.10lf ", Loss[i]);
		printf("\n");
	}
}

double Calcul_EP_Lost() /* ------- Lemma 4.1 ------  Loss^t*/
{
	Calcul_Loss_All();

	double s = 0;
	for (int i=0; i<N ; i++)  
	{
		s += Loss[i];
	}
	
	if(DEBUG == 1)
	{
		printf("Tau de perte des EPs dans l'EPN : %.10lf \n",s);
	}
	
	return s;
}
	

void Print_Fichier(FILE*f, int combine)
{

int j;

	//for (i=0; i<combine ; i++)
	//{
		for (j=0; j<13 ; j++)
		{
			fprintf(f,"%.5lf   ",Rewards[0][j]); // Rewards[0][j]) pour afficher que le meilleur cas
		}
		fprintf(f,"\n");
	//}
	
}


void Tri_Selection(int n)
{
	int i,j;
	double tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12;
	
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{	
			if ( Rewards[i][12] > Rewards[j][12] ) 
			{
				tmp0 = Rewards[i][0];
				Rewards[i][0] = Rewards[j][0];
				Rewards[j][0] = tmp0;
				
				tmp1 = Rewards[i][1];
				Rewards[i][1] = Rewards[j][1];
				Rewards[j][1] = tmp1;
				
				tmp2 = Rewards[i][2];
				Rewards[i][2] = Rewards[j][2];
				Rewards[j][2] = tmp2;
				
				tmp3 = Rewards[i][3];
				Rewards[i][3] = Rewards[j][3];
				Rewards[j][3] = tmp3;
				
				tmp4 = Rewards[i][4];
				Rewards[i][4] = Rewards[j][4];
				Rewards[j][4] = tmp4;
				
				tmp5 = Rewards[i][5];
				Rewards[i][5] = Rewards[j][5];
				Rewards[j][5] = tmp5;
				
				tmp6 = Rewards[i][6];
				Rewards[i][6] = Rewards[j][6];
				Rewards[j][6] = tmp6;
				
				tmp7 = Rewards[i][7];
				Rewards[i][7] = Rewards[j][7];
				Rewards[j][7] = tmp7;
				
				tmp8 = Rewards[i][8];
				Rewards[i][8] = Rewards[j][8];
				Rewards[j][8] = tmp8;
				
				tmp9 = Rewards[i][9];
				Rewards[i][9] = Rewards[j][9];
				Rewards[j][9] = tmp9;
				
				tmp10 = Rewards[i][10];
				Rewards[i][10] = Rewards[j][10];
				Rewards[j][10] = tmp10;
				
				tmp11 = Rewards[i][11];
				Rewards[i][11] = Rewards[j][11];
				Rewards[j][11] = tmp11;
				
				tmp12 = Rewards[i][12];
				Rewards[i][12] = Rewards[j][12];
				Rewards[j][12] = tmp12;
			}
		}
	}
}


void Calcul_Rho_Anaytique(int Alpha[N])
{
	int i,n;
	double s1,s2;

	Calcul_Beta(Alpha);
	
	s1 = 0;
	s2 = 0;
	
	for ( i=0; i<N ; i++)
	{
		for ( n=0; n<K ; n++)
		{
			Theta[i] += (1-ProbaCox[i][n])*Mu[i][n]*Beta[i][n]*1.0;
		}
		s1 += Lambda[i];
		s2 += ProbaRout[N-1][i];
	}
	
	Rho[N-1] = s1/(Theta[N-1]*( 1 - s2)*1.0);
	
	for ( i=0; i<N-1 ; i++)
	{
		Rho[i] = (Lambda[i] + (ProbaRout[N-1][i]*s1) /((1 - s2)*1.0))/(Theta[i]*1.0);
	}
	
	/*printf("Rho Analytique :");
		for (i = 0;i < N; i++)
			printf("%.10lf ", Rho[i]);
		printf("\n"); */
}

void Brute_Force()
{
	
int i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,combine;
int Alpha[N];

printf("\n-------------------- Brute Force Go ---------------------\n");

FILE *f = fopen("Rewards_Brute_Force.res","w");
fprintf(f,"\nTopologie d'etoile (LORA) avec N = %d cellules et K = %d Cox phases  \n",N,K);
fprintf(f,"Les feuilles : 1,2,3,4,5,6,7,8,9   \n");
fprintf(f,"Le noeud centrale : 10   \n");
fprintf(f,"\n ---------------------------- Algorithme Brute Force ---------------------------------\n");


	combine = 0;
	
	fprintf(f,"\nEntrée: La puissance des panneaux    : Alpha = %d ",alpha);
	fprintf(f,"\nEntrée: Nombre de panneau à affecter : P = %d \n",P);
	for (int p=P; p<=P ; p++)			//for (int p=N; p<=P ; p++)
    {
		for (i1=1; i1<=p ; i1++)
			{
				for (i2=1; i2<=p ; i2++)
				{
					for (i3=1; i3<=p ; i3++)
					{
						for (i4=1; i4<=p ; i4++)
						{
							for (i5=1; i5<=p ; i5++)
							{
								for (i6=1; i6<=p ; i6++)
								{
									for (i7=1; i7<=p ; i7++)
									{
										for (i8=1; i8<=p ; i8++)
										{
											for (i9=1; i9<=p ; i9++)
											{
												for (i10=1; i10<=p ; i10++)
												{
													if( i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 == p )
													{
														Alpha[0] = i1*alpha; Alpha[1] = i2*alpha;Alpha[2] = i3*alpha;Alpha[3] = i4*alpha;
														Alpha[4] = i5*alpha; Alpha[5] = i6*alpha;Alpha[6] = i7*alpha;Alpha[7] = i8*alpha;
														Alpha[8] = i9*alpha; Alpha[9] = i10*alpha;
														
														Init_Tabs1();
														Calcul_Rho_Anaytique(Alpha);
														
														if( Test_Stabilite() == 0 )
														{
															Rewards[combine][0] = i1; Rewards[combine][1] = i2; Rewards[combine][2] = i3; Rewards[combine][3] = i4;
															Rewards[combine][4] = i5; Rewards[combine][5] = i6; Rewards[combine][6] = i7; Rewards[combine][7] = i8;
															Rewards[combine][8] = i9; Rewards[combine][9] = i10;
															Rewards[combine][10] = i1+i2+i3+i4+i5+i6+i7+i8+i9+i10;
															Rewards[combine][11] = Calcul_N_Number(); 
															Rewards[combine][12] = Calcul_T_WaitingTime();
															
															combine++;
															
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
	
			
	 //fprintf(f,"==> Nombre de cas stable est : %d \n \n",combine);
	 printf("==> Nombre de cas stable est : %d \n",combine);

	 
	 if(combine != 0 )
	 {	
		 Tri_Selection(combine);
	 }
	 
	 Print_Fichier(f,combine);
	 
	 /*if(combine != 0 )
	 {	
		fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d, la ditribution optimal des panneaux est ci-dessus (première ligne)",alpha,P);
		printf("==> Sortie: Pour la valeur de alpha = %d et P = %d, la ditribution optimal est dans le fichier 'Rewards_Force_Brute.res' ! \n",alpha,P);
	 }
	 else
	 {
		fprintf(f,"\nSortie: Aucun cas stable pour alpha = %d et P = %d , Veuillez diminuer P ou alpha !",alpha,P);
		printf("==> Sortie: Aucun cas stable pour alpha = %d et P = %d , Veuillez diminuer P ou alpha ! \n",alpha,P);
	 } */
   }
	 fclose(f);
	
}

int Tour_Affectation(int Alpha[N])
{
	int i;
	int t = 0;
	int max = Alpha[6];
	
	for (i=0; i<6 ; i++)
	{
		if( Alpha[i] < max)
			t++;
	}
	
	return t;
}

void Heuristique1()
{
	
int i1,i2,i3,i4,i5,i6,i7,combine,combinep,p;
int Alpha[N];

printf("\n-------------------- Heuristique1 Go ---------------------\n");
FILE *f = fopen("Rewards_Heuristique1.res","w");
fprintf(f,"\nTopologie d'Arbre avec N = %d cellules et K = %d Cox phases  \n",N,K);
fprintf(f,"Les feuilles : 1,2,3,4   \n");
fprintf(f,"Les noeuds internes : 5,6   \n");
fprintf(f,"La racine  : 7   \n");
fprintf(f,"\n ---------------------------- Algorithme Heuristique1 ---------------------------------\n");


	p = N;
	combine = 2e9;
	combinep = 1e9;
	
	
	fprintf(f,"\nEntrée: La puissance des panneaux: Alpha = %d \n",alpha);

	while (  (((combine != combinep && (combine>0 || combinep>0) ) || (combine == 0 && combinep == 0))) && p <= P )
	{	
	combinep = combine;
	combine = 0;

		for (i1=1; i1<=p ; i1++)
			{
				for (i2=1; i2<=p ; i2++)
				{
					for (i3=1; i3<=p ; i3++)
					{
						for (i4=1; i4<=p ; i4++)
						{
							for (i5=1; i5<=p ; i5++)
							{
								for (i6=1; i6<=p ; i6++)
								{
									for (i7=1; i7<=p ; i7++)
									{
										if( i1 + i2 + i3 + i4 + i5 + i6 + i7 <=p )
										{
											Alpha[0] = i1*alpha; Alpha[1] = i2*alpha;Alpha[2] = i3*alpha;Alpha[3] = i4*alpha;
											Alpha[4] = i5*alpha; Alpha[5] = i6*alpha;Alpha[6] = i7*alpha;
											
											
											Init_Tabs1();
											Calcul_Beta(Alpha);
											Affiche_Beta();

											Calcul_Rho();
											Affiche_Rho();
											
											if( Test_Stabilite() == 0 )
											{
												Rewards[combine][0] = i1; Rewards[combine][1] = i2; Rewards[combine][2] = i3; Rewards[combine][3] = i4;
												Rewards[combine][4] = i5; Rewards[combine][5] = i6; Rewards[combine][6] = i7; 
												Rewards[combine][7] = i1+i2+i3+i4+i5+i6+i7;
												Rewards[combine][8] = Calcul_N_Number(); 
												Rewards[combine][9] = Calcul_T_WaitingTime();
												
												combine++;
												
											}
										}
									}
								}
							}
						}
					}
				}
			}	
			
	 fprintf(f,"\n------------ Nbre de panneau pouvant être utiliser = '%d'  ------------------ \n\n",p);
	 Tri_Selection(combine);
	 fprintf(f,"==> Nombre de combinaisons stable: combine = %d \n",combine);
	 
	 if( combine != 0 )
	 {
		    fprintf(f,"==> Meilleur configuration intiale (ci-dessous) : \n");
		    fprintf(f,"==> Correspond à un nombre de panneau = %.2f  \n \n",Rewards[0][7]);
			Print_Fichier(f,combine);
			fprintf(f,"\n");
			printf("Iteration P = %d , combinaisons = %d\n",p,combine);
			break;
	 }
	 printf("Iteration P = %d , combinaisons = %d\n",p,combine);
	 p++;
	} //fin while
	if( combine == 0) // aucune solution trouvée inferieur à P
	{
		fprintf(f,"\nSortie: Aucun cas stable pour alpha = %d et P = %d , Veuillez diminuer P ou alpha !",alpha,P);
		printf("==> Sortie: Aucun cas stable pour alpha = %d et P = %d , Veuillez diminuer P ou alpha ! \n",alpha,P);
	}
	else  //----------------- Debut de construction des meilleurs solutions -------------------
	{
	int arret,i,k;
    arret = 0;
	int tour;
	
	while( arret == 0 )
	{
	//printf("TEST1\n");
	     for (k=0; k<7 ; k++)
				 {
					 Alpha[k] = Rewards[0][k]*alpha;
				 }	 
		tour = Tour_Affectation(Alpha);
		
		if( tour == 0)
		{
				//printf("TEST2\n");
				 		
			     Alpha[6] += alpha;
				 
				 Init_Tabs1();
				 Calcul_Beta(Alpha);
				 Affiche_Beta();

				 Calcul_Rho();
				 Affiche_Rho();
				 
				 //Pavant = Papres;
				 if( Test_Stabilite() == 0 )
				 {

					Rewards[0][6] += 1; 
					Rewards[0][7] = Rewards[0][1]+Rewards[0][1]+Rewards[0][2]+Rewards[0][3]+Rewards[0][4]+Rewards[0][5]+Rewards[0][6];
					Rewards[0][8] = Calcul_N_Number(); 
					Rewards[0][9] = Calcul_T_WaitingTime();
					//Papres = Rewards[0][7];
					
					fprintf(f,"==> Nombre de panneau à utiliser = %.2f \n",Rewards[0][7]); 
					fprintf(f,"==> Meilleur configuration (ci-dessous): \n\n");
					Print_Fichier(f,combine);
					fprintf(f,"\n");
					
					printf("Iteration P = %.4f, OK \n",Rewards[0][7]);
				 }
				 else
				 {
					printf("Iteration P = %.4f, KO ! \n",Rewards[0][7]);

				 } 
				 if(Rewards[0][7] == P)
				 	{arret = 1;break;}
				 //if( Rewards[0][9] <= ARRET )
				 //	arret = 1;
		}
		else
		{	
			//printf("TEST3\n");
		 for (i=5; i>=0 ; i--)
		 {
			 //if(Rewards[0][i] < Rewards[0][6])
			 if(Alpha[i] < Alpha[6])
			 {

				 Alpha[i] += alpha;
				 
				 Init_Tabs1();
				 Calcul_Beta(Alpha);
				 Affiche_Beta();

				 Calcul_Rho();
				 Affiche_Rho();
				 
				 //Pavant = Papres;
				 if( Test_Stabilite() == 0 )
				 {
					//printf("TEST4\n");
					Rewards[0][i] += 1; 
					Rewards[0][7] = Rewards[0][0]+Rewards[0][1]+Rewards[0][2]+Rewards[0][3]+Rewards[0][4]+Rewards[0][5]+Rewards[0][6];
					Rewards[0][8] = Calcul_N_Number(); 
					Rewards[0][9] = Calcul_T_WaitingTime();
					//Papres = Rewards[0][7];
					
					fprintf(f,"==> Nombre de panneau effectif à utiliser = %.2f \n",Rewards[0][7]); 
					fprintf(f,"==> Meilleur configuration (ci-dessous): \n\n");
					Print_Fichier(f,combine);
					fprintf(f,"\n");

					printf("Iteration P = %.4f, OK \n",Rewards[0][7]);
				 }
				 else
				 {
					printf("Iteration P = %.4f, KO ! \n",Rewards[0][7]);
				 }
				 if(Rewards[0][7] == P)
				 {arret = 1; break;}
				 //if( Rewards[0][9] <= ARRET )
				 //	arret = 1;
			 }
		}
	   }//fin else
	} // fin while
		fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d, la ditribution optimal des panneaux est ci-dessus (première ligne)",alpha,P);
		printf("==> Sortie: Pour la valeur de alpha = %d et P = %d, la ditribution optimal est dans le fichier 'Rewards_Heuristique1.res' ! \n",alpha,P);
	}
fclose(f);
	
}

void Test_configuration()
{
	int i1, i2, i3 , i4, i5, i6, i7;
	int Alpha[N];
	
	i1 = 2; i2 = 3; i3 = 3; i4 = 5; i5 = 3; i6 = 6; i7 = 8;
	Alpha[0] = i1*alpha; Alpha[1] = i2*alpha;Alpha[2] = i3*alpha;Alpha[3] = i4*alpha;
	Alpha[4] = i5*alpha; Alpha[5] = i6*alpha;Alpha[6] = i7*alpha;
											
	//if( i1+i2+i3+i4+i5+i6+i7 == P )
	//{										
		Init_Tabs1();
		Calcul_Beta(Alpha);
		Affiche_Beta();

		Calcul_Rho();
		//Affiche_Rho();
												
		if( Test_Stabilite() == 0 )
		{
			printf("Nombre moyen de DPs  = %f \n", Calcul_N_Number());
			printf("Temps d'attente moyen d'un DP  = %f \n", Calcul_T_WaitingTime());		
		}
	//}
	//else
	//{
		//printf("Attention le nombre de panneau doit être P = %d\n",P);
	//}
}

int Condition_Arret(double tab1[N], double tab2[N])
{	
	int i;
	double s = 0;
	for (i=0; i<N ; i++)
	{
		s += fabs(tab1[i] - tab2[i]);
	}
	if( s < EPSILON )
	{
		Cst++;
	}
	else
	{
		Cst = 0;
	}
	
	return Cst == NITER ? 1:0;
}
void Calcul_Rho_Iterative(int Alpha[N])
{
	int i,j,n;
	double s1, RhoAVT[N];
	Cst = 0;

	Calcul_Beta(Alpha);
	
	for ( i=0; i<N ; i++)
	{
		RhoAVT[i] = 0;
		Rho[i] = 1; // initialisation à n'imorte quelle valeur differente de 0
		for ( n=0; n<K ; n++)
		{
			Theta[i] += (1-ProbaCox[i][n])*Mu[i][n]*Beta[i][n]*1.0;
		}
	}
	
	
	int c = 0;
	while( Condition_Arret(RhoAVT,Rho) == 0 )
	{
		for ( i=0; i<N ; i++)
		{
			RhoAVT[i] = Rho[i];
		}

		for ( i=0; i<N ; i++)
		{
			s1 = 0;
			for ( n=0; n<K ; n++)
			{
				for ( j=0; j<N ; j++)
				{
					s1 += (1-ProbaCox[j][n])*Mu[j][n]*Beta[j][n]*Rho[j]*ProbaRout[j][i]*1.0;
				}
			}
			
			Rho[i] = (Lambda[i] + s1)/(Theta[i]*1.0);
		}
		
		c++;
	}
	/*printf("Aprés %d iterations \n Rho Itertative :",c);
		for (i = 0;i < N; i++)
			printf("%.10lf ", Rho[i]);
		printf("\n");*/
		//printf("E[N] = %.5lf  and E[T] = %.5lf \n ",Calcul_N_Number(),Calcul_T_WaitingTime());		

}



void Calcul_Bornes()
{
	int i;
	double s,s1,s2;
	
	for ( i=0; i<N ; i++)
	{
		s1 += Lambda[i];
		s2 += ProbaRout[N-1][i];
	}
	s2 -= ProbaRout[N-1][N-1];
	
	for (i=0; i<N-1 ; i++)
			{
				s = ((Gama[i] + Mu[i][0])/(alpha*Mu[i][0]*1.0)) * (Lambda[i] + ((ProbaRout[N-1][i]*s1) /((1 - s2)*1.0)) );
				
				//printf("%.20lf \n", s);
				DInf[i] = ceil(s);
				if(DInf[i] == s)	
					DInf[i]++;
			}
			
			s = ((Gama[N-1] + Mu[N-1][0])/(alpha*Mu[N-1][0]*1.0)) * (s1 /((1 - s2)*1.0)) ;
				
				//printf("%.20lf \n", s);
				DInf[N-1] = ceil(s);
				if(DInf[N-1] == s)	
					DInf[N-1]++;
	
	int n,r;
			
	for (i=0; i<N ; i++)
			{
				s = (Gama[i] + Mu[i][0])/(alpha*1.0);
				s1 = 1;
				for (n=1; n<K ; n++)
				{
					s2 = 1;
					for (r=0; r<n ; r++)
					{
						s2 *= ProbaCox[i][r]*1.0;
					}
					s2 *= Mu[i][0]/(Mu[i][n]*1.0);
					
					s1 += s2;
				}
				s *= (1.0/s1);
				DSup[i] = floorl(s);
				if(DSup[i] == s)	
					DSup[i]--;
				
			}
			
	if(DEBUG == 2)
	{
		printf("Borne Inf: [");
		for (i=0; i<N ; i++)
			{
				printf("%d - ",DInf[i]);
			}
		printf("]\n");
		printf("Borne Sup: [");
		for (i=0; i<N ; i++)
			{
				printf("%d - ",DSup[i]);
			}
		printf("]\n");
	}
				
}

void Test_Bornes(FILE *f)
{
 int i;
 for (i=0; i<N ; i++)
 {
	 if(DSup[i] < DInf[i])
	 {
		 break;
	 }
 }
 if( i == N )
 {
	 printf("Verification des bornes, Ok \n");
 }
 else
 {
	 printf("Problème Bornes: Borne Sup inferieur à la Borne Inf, cellule %d  KO \n",i+1);
	 exit(0);
 }
 
 int s1, s2;
 s1 = 0;
 s2 = 0;
 for (i=0; i<N ; i++)
 {
	 s1 += DInf[i];
	 s2 += DSup[i];
 } 
 
 if( P < s1 )
 {
	 printf("Problème Panneaux: Il faut au moins %d panneaux pour ces paramètres, KO \n",s1);
	 exit(0);
 }
 else
 {
	 printf("Verification des panneaux, Ok \n");
 }
 if( P >= s2 )
 {
	int Alpha[N];
	for (i=0; i<N ; i++)
	{
		Alpha[i] =  DSup[i]*alpha; 
	}
	Init_Tabs1();
	Calcul_Rho_Anaytique(Alpha);


	
	if( Test_Stabilite() == 0 )
	{
		s1 = 0;
		for (i=0; i<N ; i++)
					{fprintf(f,"%5d    ",DSup[i]); s1 += DSup[i];}
					fprintf(f,"%5d    ",s1);
					
			double ss1, ss2, ss3;
			ss1 = 0; ss2 = 0;
			//ss3 = Calcul_EP_Lost();
			
			for (i=0; i<N ; i++)
			{
					ss1 += Loss1[i]; 
					ss2 += Loss2[i];
			}
			
			fprintf(f,"%.5lf    ",Calcul_N_Number());
			fprintf(f,"%.5lf    ",Calcul_T_WaitingTime());	
			fprintf(f,"%.5lf    ",ss1);
			fprintf(f,"%.5lf    \n",ss2);
			//fprintf(f,"%.5lf  \n",ss3);	

					
		 fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d: %d panneaux ne seront pas utilisé (voir ci-dessus)",alpha,P,P-s1);
		 printf("\nSortie: Pour la valeur de alpha = %d et P = %d: %d panneaux ne seront pas utilisé 'Rewards_Heuristique2.res'\n",alpha,P,P-s1);
		 fclose(f);
		 exit(0);
	 }
	 else
	 { printf("EUH problème, la borne SUP n'est pas autorisé !? \n");fclose(f);exit(0); }
  }
 
}

int Max_Rho(int tmp[N], FILE*f) // Je retourne l'indice de la cellule la plus chargée qui peut recevoir un panneau
{
	int i;
	int bool[N];
	
	
	for ( i=0; i<N ; i++)
	{
		if(tmp[i] + 1 <= DSup[i]) // OK
			bool[i] = 1;
		else                      // KO, Max atteint
			bool[i] = -1;
	}
	
	int test = 0;
	int s = 0;
	for ( i=0; i<N ; i++)
	{
		s += DSup[i];
		if( bool[i] == -1)
			test++;
	}
	if( test == N ){ printf("La borne sup est atteinte ! On ne peux plus rajouter de panneaux solaire ! \n ");
					fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d: %d panneaux ne seront pas utilisé (voir ci-dessus)",alpha,P,P-s);
					printf("\nSortie: Pour la valeur de alpha = %d et P = %d: %d panneaux ne seront pas utilisé 'Rewards_Heuristique2.res'\n",alpha,P,P-s);
					exit(0);
					}
	double max;
	int indice;
	
	for ( i=0; i<N ; i++)
	{
		if( bool[i] == 1 )
		{
			 max = Rho[i];
			 indice = i;	
			 break;
		}	
	}
	
	for ( i=0; i<N ; i++)
	{
		if(max < Rho[i] && bool[i] == 1 ) 
		{
			max = Rho[i];
			indice = i;
		}
	}
	return indice;
}


void Heuristique2(){
	
	printf("\n-------------------- Heuristique2 Go ---------------------\n\n");
	FILE *f = fopen("Rewards_Heuristique2.res","w");
	fprintf(f,"\nTopologie d'etoile (LORA) avec N = %d cellules et K = %d Cox phases  \n",N,K);
	fprintf(f,"Les feuilles : 1,2,3,4,5,6,7,8,9   \n");
	fprintf(f,"Le noeud centrale : 10   \n");
	fprintf(f,"\n ---------------------------- Algorithme Heuristique2 ---------------------------------\n\n");
	fprintf(f,"Entrée: Puissance d'un panneau alpha = %d et  le nombre de panneaux  P = %d panneaux \n\n",alpha,P);
	printf("\nEntrée: Puissance d'un panneau alpha = %d et  le nombre de panneaux  P = %d panneaux \n\n",alpha,P);

	Calcul_Bornes();
    Test_Bornes(f);
    
    int i;
    int tmp[N], Alpha[N], s;
    double s1, s2, s3;
    
    s = 0;
    for (i=0; i<N ; i++)
	{
		tmp[i] = DInf[i];
		s += tmp[i];
	}
	printf("Iteration P = %d \n",s);
   
    // --------- On commence par prendre la configuration de la borne inf ---------------------
    for (i=0; i<N ; i++)
	{
		Alpha[i] =  tmp[i]*alpha;
	}
																				
	Init_Tabs1();
	Calcul_Rho_Anaytique(Alpha);

												
	if( Test_Stabilite() == 0 )
		{
			s1 = 0; s2 = 0;
			s3 = Calcul_EP_Lost();
			
			for (i=0; i<N ; i++)
			{
					s1 += Loss1[i]; 
					s2 += Loss2[i];
			}
			
			for (i=0; i<N ; i++)
					fprintf(f,"%5d    ",tmp[i]);
			fprintf(f,"%5d    ",s);
			fprintf(f,"%.5lf    ",Calcul_N_Number());
			fprintf(f,"%.5lf    ",Calcul_T_WaitingTime());	
			fprintf(f,"%.5lf    ",s1);
			fprintf(f,"%.5lf    ",s2);
			fprintf(f,"%.5lf  \n",s3);	
		}
	else{ printf("EUH problème, la borne Inf n'est pas autorisé !? \n");fclose(f);exit(0); }
	if( s == P)
	 { 	
		fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée des panneaux est ci-dessus",alpha,P);
		printf("==> Sortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée est dans le fichier 'Rewards_Heuristique2.res' ! \n",alpha,P);
		fclose(f);
		exit(0);
	 }
	 
	 // --------- On debute la construction a partir de la borne inf  ---------------------
	 int max;
	 while( s < P )
	 {
		max = Max_Rho(tmp,f);
		tmp[max]++;
		Alpha[max] += alpha;
		s++;
		printf("Iteration P = %d \n",s);

		
		Init_Tabs1();
		Calcul_Rho_Anaytique(Alpha);

												
		if( Test_Stabilite() == 0 )
		{
			s1 = 0; s2 = 0;
			s3 = Calcul_EP_Lost();
			
			for (i=0; i<N ; i++)
			{
					s1 += Loss1[i]; 
					s2 += Loss2[i];
			}
			
			for (i=0; i<N ; i++)
					fprintf(f,"%5d    ",tmp[i]);
			fprintf(f,"%5d    ",s);
			fprintf(f,"%.5lf    ",Calcul_N_Number());
			fprintf(f,"%.5lf    ",Calcul_T_WaitingTime());	
			fprintf(f,"%.5lf    ",s1);
			fprintf(f,"%.5lf    ",s2);
			fprintf(f,"%.5lf  \n",s3);	
		}
		else {  printf("EUH problème, On avais deja verifié que ça allais être stable  !? \n"); fclose(f);exit(0); }
	}
	
	fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée des panneaux est ci-dessus ",alpha,P);
	printf("==> Sortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée est dans le fichier 'Rewards_Heuristique2.res' ! \n",alpha,P);
	
	
fclose(f);	
}


double Function(int phi,int i) //  fonction calcule f(phi_i)
{
	int j;
	double s1,s2,s3,s4;
	
	s1 = 0;
	s2 = 0;
	for (j=0; j<N ; j++)
	{
		s1 += Lambda[j];
		s2 += ProbaRout[N-1][j];
	}
	

	if( i == N-1 )
	{
		s3 = (Gama[N-1]+Mu[N-1][0])*( s1/((1-s2)*1.0) )*1.0;
		s4 = s3*( 1.0/((alpha*phi*Mu[N-1][0]) - s3 ))*( 1.0/s1 )*1.0;
	}
	else
	{
		s3 = (Gama[i]+Mu[i][0])*( (ProbaRout[N-1][i]*s1)/((1-s2)*1.0) + Lambda[i] )*1.0;
		s4 = s3*( 1.0/((alpha*phi*Mu[i][0]) - s3 ))*( 1.0/s1 )*1.0;
		
	}
	
	
	if( s4 < 0 ) {printf("OULAH, probleme, F(PHI(%d)) = %f  est negatif !! \n",i,s4); exit(0);}
	return s4;
}

double Objective_Function(int tmp[N])
	{
		int i;
		double s = 0;
		
		for (i=0; i<N ; i++)
		{
			s += Function(tmp[i],i);
		}
		
		return s;
	}


void Gradient(){
	
	printf("\n-------------------- Gradient Go ---------------------\n\n");
	FILE *f = fopen("Rewards_Gradient.res","w");
	fprintf(f,"\nTopologie d'etoile (LORA) avec N = %d cellules et K = %d Cox phases  \n",N,K);
	fprintf(f,"Les feuilles : 1,2,3,4,5,6,7,8,9   \n");
	fprintf(f,"Le noeud centrale : 10   \n");
	fprintf(f,"\n ---------------------------- Algorithme Gradient ---------------------------------\n\n");
	fprintf(f,"Entrée: Puissance d'un panneau alpha = %d et  le nombre de panneaux  P = %d panneaux \n\n",alpha,P);
	printf("\nEntrée: Puissance d'un panneau alpha = %d et  le nombre de panneaux  P = %d panneaux \n\n",alpha,P);

	Calcul_Bornes();
    Test_Bornes(f);
    
    int i,s;
    int tmp[N],Alpha[N];
    double s1, s2, s3;
    
    
     // --------- On commence par prendre la configuration de la borne inf ---------------------
    s = 0;
    for (i=0; i<N ; i++)
	{
		tmp[i] = DInf[i];
		s += tmp[i];
		Alpha[i] = tmp[i]*alpha;
	}
	printf("Iteration P = %d \n",s);
	
	Init_Tabs1();
	Calcul_Rho_Anaytique(Alpha);
	
	if(Test_Stabilite() == 1) {printf("Hmmmm un problème est survenue, cas instable \n");fclose(f);exit(0);}
   
			s1 = 0; s2 = 0;
			s3 = Calcul_EP_Lost();
			
			for (i=0; i<N ; i++)
			{
					s1 += Loss1[i]; 
					s2 += Loss2[i];
			}
			
			for (i=0; i<N ; i++)
					fprintf(f,"%5d    ",tmp[i]);
			fprintf(f,"%5d    ",s);
			fprintf(f,"%.5lf    ",Calcul_N_Number());
			fprintf(f,"%.5lf    ",Calcul_T_WaitingTime());	
			fprintf(f,"%.5lf    ",s1);
			fprintf(f,"%.5lf    ",s2);
			fprintf(f,"%.5lf  \n",s3);	
			
	if( s == P)
	 { 	
		fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée des panneaux est ci-dessus",alpha,P);
		printf("==> Sortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée est dans le fichier 'Rewards_Gradient.res' ! \n",alpha,P);
		fclose(f);
		exit(0);
	 }
	 
	 // --------- On debute la construction a partir de la borne inf  ---------------------
	 int indice;
	 double some, functions1[N], functions2[N],min;
	 int bool[N];
	 
		
	 while( s < P )
	 {
		for (i=0; i<N ; i++)
		{
			bool[i] = 0;
		}
		some = Objective_Function(tmp);

		
		for (i=0; i<N ; i++)
		{
			if(tmp[i] + 1 <= DSup[i])
			{
				functions2[i] = some - Function(tmp[i],i) + Function(tmp[i]+1,i);
				bool[i] = 1;
			}
			else
			{
				bool[i] = -1;
			}
		}
		indice = -1;
		for (i=0; i<N ; i++)
		{
			if(bool[i] == 1)
			{
				indice = i;
				min = functions2[i];
				break;
			}
		}
		
		for (i=0; i<N ; i++)
		{
			if(bool[i] == 1 && functions2[i] < min )
			{
				indice = i;
				min = functions2[i];
			}
		}
		if(indice == -1) {printf("Toute les cellules on atteints leur borne sup \n"),fclose(f);exit(0);}
		
		tmp[indice]++;
		//functions1[indice] =  Function(tmp[indice],indice);
		for (i=0; i<N ; i++)
		{
			Alpha[i] = tmp[i]*alpha;
		}
		
		Init_Tabs1();
		Calcul_Rho_Anaytique(Alpha);

		
		s++;
		printf("Iteration P = %d \n",s);
		
		if(Test_Stabilite() == 1) {printf("Hmmmm un problème est survenue, cas instable \n");fclose(f);exit(0);}

			s1 = 0; s2 = 0;
			s3 = Calcul_EP_Lost();
			
			for (i=0; i<N ; i++)
			{
					s1 += Loss1[i]; 
					s2 += Loss2[i];
			}
			
			for (i=0; i<N ; i++)
					fprintf(f,"%5d    ",tmp[i]);
			fprintf(f,"%5d    ",s);
			fprintf(f,"%.5lf    ",Calcul_N_Number());
			fprintf(f,"%.5lf    ",Calcul_T_WaitingTime());	
			fprintf(f,"%.5lf    ",s1);
			fprintf(f,"%.5lf    ",s2);
			fprintf(f,"%.5lf  \n",s3);	
	} //fin While
	
	fprintf(f,"\nSortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée des panneaux est ci-dessus ",alpha,P);
	printf("==> Sortie: Pour la valeur de alpha = %d et P = %d, la ditribution proposée est dans le fichier 'Rewards_Gradient2.res' ! \n",alpha,P);
	
	
fclose(f);	
}






int main(int argc,char *argv[])
{

if( P < N )
{
	printf("Erreur ! le nombre de panneaux doit être au moins egale au nombre de cellules: P >= N donc P >= %d \n",N);
	exit(0);
}
Init_Tabs2();
//Lecture_Matrice_Inverse();


//Calcul_Rho_Iterative(Alpha);
//Calcul_Rho_Anaytique(Alpha);

//Brute_Force();
//Heuristique1();
//Heuristique2();
Gradient();
//Test_configuration();


free(inverse);

exit(0);
}
