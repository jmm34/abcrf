
#include <Rmath.h>
#include <R.h>


void findweightsAmelioree(double *ONv, int *inbag, double *NNv, double *WEv, int *nobs, int *nnew, int *ntree, double *thres, int *counti, int *normalise){
  
  int k,i,j;
  double absol;
  int meancount;
	
/* Une boucle pour les arbres */
  for ( k = 1; k <= (*ntree); ++k ){
	  /* Une boucle pour parcourir les prédictions du test set */
    for ( i=1; i<=(*nnew); ++i ){
      meancount = 0;
		/* Une boucle pour parcourir les prédictions du training set */
      for ( j=1; j<=(*nobs); ++j ){
		absol = ONv[ j+(k-1)* (*nobs)-1 ] - NNv[ i+(k-1)* (*nnew)-1 ];
		/* On regarde si on a égalité des prédictions (donc tombent dans la même feuille) et si l'obs a bien servi à la construction de l'arbre */
		if( ( absol <= (*thres) ) && ( absol >= -(*thres) ) && ( inbag[ j+(k-1)* (*nobs)-1 ] != 0 ) ){
			counti[ j-1 ] = 1;
			meancount = meancount + 1;
		}else{
			counti[ j-1 ] = 0;
		}
      }
      if( meancount >= 1 ){
		if( ( *normalise ) >= 1){
			for ( j=1 ; j <= (*nobs); ++j ){
				/* On fait la somme des indicatrice sur le nombre d'éléments dans la feuille */
				WEv[ j+(i-1)* (*nobs)-1 ] = WEv[ j+(i-1)* (*nobs)-1 ] + counti[ j-1 ]/ (double) (meancount)  ;
			}
		}else{
			for ( j=1; j<=(*nobs); ++j ){
				/* Ou bien on incrémente de zéro */
				WEv[ j+(i-1)* (*nobs)-1 ] = WEv[ j+(i-1)* (*nobs)-1 ] + counti[ j-1 ] ;
			}
		}

      }
      

    }
  
  }    
} 
