
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef NEW_H_
#define NEW_H_

#ifndef MAIN_H_
#include <infoli.h>
#endif /* MAIN_H_ */

/*** FUNCTION PROTOTYPES ***/
//void ComputeNetwork_new(bool,bool ,cellState* , mod_prec*, int , int, mod_prec*,int,mod_prec* );
//void ComputeNetwork_old(bool,bool ,cellState* , mod_prec*, int , int, mod_prec*,int,mod_prec* );
void ComputeOneCellTimeMux0_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux1_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux2_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux3_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux4_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux5_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux6_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux7_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);
void ComputeOneCellTimeMux8_newJH(int,  mod_prec*,mod_prec* , int, int, mod_prec*);

//void ComputeOneCell( int,int, mod_prec, mod_prec* );
void ComputeOneCell0_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell1_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell2_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell3_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell4_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell5_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell6_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);
void ComputeOneCell7_newJH( int,int, mod_prec ,mod_prec*, int, mod_prec*);


Dend CompDend_newJH(struct Dend, mod_prec , mod_prec,mod_prec*, int , mod_prec*,int);
mod_prec DendHCurr_newJH(struct channelParams, mod_prec );
mod_prec DendCaCurr_newJH(struct channelParams, mod_prec );
mod_prec DendKCurr_newJH(struct channelParams );
mod_prec DendCal_newJH(struct channelParams );
dendCurrVoltPrms DendCurrVolt_newJH(struct dendCurrVoltPrms );
mod_prec IcNeighbors_newJH(mod_prec*, mod_prec, int, mod_prec*,int);

Soma CompSoma_newJH(struct Soma, mod_prec ,mod_prec );
channelParams SomaCalcium_newJH(struct channelParams );
channelParams SomaSodium_newJH(struct channelParams );
channelParams SomaPotassium_newJH(struct channelParams );
mod_prec SomaPotassiumX_newJH(struct channelParams );
mod_prec SomaCurrVolt_newJH(struct somaCurrVoltPrms );

Axon CompAxon_newJH(struct Axon , mod_prec );
channelParams AxonSodium_newJH(channelParams );
mod_prec AxonPotassium_newJH(channelParams );
mod_prec AxonCurrVolt_newJH(axonCurrVoltPrms );

#endif /* NEW_H_ */
