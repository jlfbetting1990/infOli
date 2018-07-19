
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef NEW_H_
#define NEW_H_

#ifndef MAIN_H_
#include <infoli.h>
#endif /* MAIN_H_ */

/*** CONSTANTS ***/

#define DIV_FOUR 0.25f
#define DIV_FIVE 0.2f
#define DIV_THIRTEEN 1/13.9f
#define FOUR 1/4.2f
#define EIGHT 1/8.5f
#define THIRTY 1/30.00f
#define SEVEN 1/7.3f
#define FIVE_FIVE 1/5.5f
#define FIVE_EIGHT -1/5.8f
#define THIRTY_THREE 1/33.0f
#define TEN 0.1f
#define TWELVE -1/12.0f
#define NINE_HUNDRED 1/900.0f
#define HUNDRED -0.01

/*** FUNCTION PROTOTYPES ***/

//void ComputeOneCell( int,int, mod_prec, mod_prec* );
void ComputeOneCell0_newJH(int, mod_prec ,mod_prec*, int, mod_prec*);


Dend CompDend_newJH(struct Dend, mod_prec , mod_prec,mod_prec*, int , mod_prec*);
mod_prec DendHCurr_newJH(struct channelParams, mod_prec );
mod_prec DendCaCurr_newJH(struct channelParams, mod_prec );
mod_prec DendKCurr_newJH(struct channelParams );
mod_prec DendCal_newJH(struct channelParams );
dendCurrVoltPrms DendCurrVolt_newJH(struct dendCurrVoltPrms );
mod_prec IcNeighbors_newJH(mod_prec*, mod_prec, int, mod_prec*);

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
