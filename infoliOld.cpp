#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
//#include <sys/time.h>
#include <time.h>
#include "infoli.h"


//cellState local_state[IO_NETWORK_DIM1*IO_NETWORK_DIM2][TIME_MUX_FACTOR];
cellState local_state0[MAX_TIME_MUX];
cellState local_state1[MAX_TIME_MUX];
cellState local_state2[MAX_TIME_MUX];
cellState local_state3[MAX_TIME_MUX];
cellState local_state4[MAX_TIME_MUX];


cellState new_state0[MAX_TIME_MUX];
cellState new_state1[MAX_TIME_MUX];
cellState new_state2[MAX_TIME_MUX];
cellState new_state3[MAX_TIME_MUX];
cellState new_state4[MAX_TIME_MUX];

mod_prec IAppin0[MAX_TIME_MUX];
mod_prec IAppin1[MAX_TIME_MUX];
mod_prec IAppin2[MAX_TIME_MUX];
mod_prec IAppin3[MAX_TIME_MUX];
mod_prec IAppin4[MAX_TIME_MUX];
mod_prec Connectivity_Matrix0[CONN_MATRIX_MAX/5];
mod_prec Connectivity_Matrix1[CONN_MATRIX_MAX/5];
mod_prec Connectivity_Matrix2[CONN_MATRIX_MAX/5];
mod_prec Connectivity_Matrix3[CONN_MATRIX_MAX/5];
mod_prec Connectivity_Matrix4[CONN_MATRIX_MAX/5];


void ComputeNetwork_old(bool ini,bool new_matrix,cellState IniArray[IO_NETWORK_SIZE], mod_prec iAppin[IO_NETWORK_SIZE] , int N_Size, int Mux_Factor,mod_prec Connectivity_Matrix[CONN_MATRIX_SIZE], int Conn_Matrix_Size, mod_prec cellOut[IO_NETWORK_SIZE]){

	int i,j,k;
//	for(j=0; j<IO_NETWORK_SIZE; ++j) {
//		printf("cell %d: Ca2Plus: %.9f\n",j,IniArray[j].dend.Ca2Plus);
//	}
//	return;

	//returnState AxonOut;
	mod_prec neighVdend0[MAX_N_SIZE];
	mod_prec neighVdend1[MAX_N_SIZE];
	mod_prec neighVdend2[MAX_N_SIZE];
	mod_prec neighVdend3[MAX_N_SIZE];
	mod_prec neighVdend4[MAX_N_SIZE];

	//Initialize cell clusters
	if(ini){

		ComputeNetwork_label0:	for(j=0;j<Mux_Factor;j++){

			local_state0[j]= IniArray[j];
			local_state1[j]= IniArray[j+Mux_Factor];
			local_state2[j]= IniArray[j+(Mux_Factor*2)];
			local_state3[j]= IniArray[j+(Mux_Factor*3)];
			local_state4[j]= IniArray[j+(Mux_Factor*4)];
		}

	}


	if(new_matrix){
		int max_set = (IO_NETWORK_SIZE*IO_NETWORK_SIZE)/HW_CELLS; 
		int its = max_set*IO_NETWORK_SIZE;
		ComputeNetwork_connectivity: for (j=0;j<IO_NETWORK_SIZE*Mux_Factor;++j){
			Connectivity_Matrix0[j] = Connectivity_Matrix[((j*IO_NETWORK_SIZE) % (IO_NETWORK_SIZE*IO_NETWORK_SIZE-1))];
			Connectivity_Matrix1[j] = Connectivity_Matrix[((j*IO_NETWORK_SIZE) % (IO_NETWORK_SIZE*IO_NETWORK_SIZE-1))+Mux_Factor];
			Connectivity_Matrix2[j] = Connectivity_Matrix[((j*IO_NETWORK_SIZE) % (IO_NETWORK_SIZE*IO_NETWORK_SIZE-1))+2*Mux_Factor];
			Connectivity_Matrix3[j] = Connectivity_Matrix[((j*IO_NETWORK_SIZE) % (IO_NETWORK_SIZE*IO_NETWORK_SIZE-1))+3*Mux_Factor];
			Connectivity_Matrix4[j] = Connectivity_Matrix[((j*IO_NETWORK_SIZE) % (IO_NETWORK_SIZE*IO_NETWORK_SIZE-1))+4*Mux_Factor];
		}
//		printf("Matrix 0\n");
//		for (j=0;j<IO_NETWORK_SIZE*Mux_Factor;++j) {
//			printf("%.2f ",Connectivity_Matrix0[j]);
//		}
//		printf("\n");
//		
//		printf("Matrix 1\n");
//		for (j=0;j<IO_NETWORK_SIZE*Mux_Factor;++j) {
//			printf("%.2f ",Connectivity_Matrix1[j]);
//		}
//		printf("\n");
//
//		printf("Matrix 2\n");
//		for (j=0;j<IO_NETWORK_SIZE*Mux_Factor;++j) {
//			printf("%.2f ",Connectivity_Matrix2[j]);
//		}
//		printf("\n");
//
//		printf("Matrix 3\n");
//		for (j=0;j<IO_NETWORK_SIZE*Mux_Factor;++j) {
//			printf("%.2f ",Connectivity_Matrix3[j]);
//		}
//		printf("\n");
//
//		printf("Matrix 4\n");
//		for (j=0;j<IO_NETWORK_SIZE*Mux_Factor;++j) {
//			printf("%.2f ",Connectivity_Matrix4[j]);
//		}
//		printf("\n");
//		return;
	}
	//Save invoked inputs on Block RAM and pick up the neighboring inputs for each cluster.
	ComputeNetwork_label1:	for(j=0;j<Mux_Factor;j++){


		IAppin0[j]=iAppin[j+(Mux_Factor*0)];
		IAppin1[j]=iAppin[j+(Mux_Factor*1)];
		IAppin2[j]=iAppin[j+(Mux_Factor*2)];
		IAppin3[j]=iAppin[j+(Mux_Factor*3)];
		IAppin4[j]=iAppin[j+(Mux_Factor*4)];
	}

	//Neighboring inputs preparation
	k=0;
	ComputeNetwork_label2:	for(i=0;i<Mux_Factor;i++){

		neighVdend0[i] = local_state0[i].dend.V_dend ;
		neighVdend0[i+Mux_Factor] = local_state1[i].dend.V_dend ;
		neighVdend0[i+2*Mux_Factor] = local_state2[i].dend.V_dend ;
		neighVdend0[i+3*Mux_Factor] = local_state3[i].dend.V_dend ;
		neighVdend0[i+4*Mux_Factor] = local_state4[i].dend.V_dend ;
	    	
		k=k+HW_CELLS;
	}
	k=0;
	ComputeNetwork_label3:	for(i=0;i<Mux_Factor;i++){

		neighVdend1[i] = local_state0[i].dend.V_dend ;
		neighVdend1[i+Mux_Factor] = local_state1[i].dend.V_dend ;
		neighVdend1[i+2*Mux_Factor] = local_state2[i].dend.V_dend ;
		neighVdend1[i+3*Mux_Factor] = local_state3[i].dend.V_dend ;
		neighVdend1[i+4*Mux_Factor] = local_state4[i].dend.V_dend ;
		
		k=k+HW_CELLS;
	}
	k=0;
	ComputeNetwork_label4:	for(i=0;i<Mux_Factor;i++){

		neighVdend2[i] = local_state0[i].dend.V_dend ;
		neighVdend2[i+Mux_Factor] = local_state1[i].dend.V_dend ;
		neighVdend2[i+2*Mux_Factor] = local_state2[i].dend.V_dend ;
		neighVdend2[i+3*Mux_Factor] = local_state3[i].dend.V_dend ;
		neighVdend2[i+4*Mux_Factor] = local_state4[i].dend.V_dend ;

		k=k+HW_CELLS;
	}
	k=0;
	ComputeNetwork_label5:	for(i=0;i<Mux_Factor;i++){

		neighVdend3[i] = local_state0[i].dend.V_dend ;
		neighVdend3[i+Mux_Factor] = local_state1[i].dend.V_dend ;
		neighVdend3[i+2*Mux_Factor] = local_state2[i].dend.V_dend ;
		neighVdend3[i+3*Mux_Factor] = local_state3[i].dend.V_dend ;
		neighVdend3[i+4*Mux_Factor] = local_state4[i].dend.V_dend ;

		k=k+HW_CELLS;
	}

	k=0;
	ComputeNetwork_label6:	for(i=0;i<Mux_Factor;i++){

		neighVdend4[i] = local_state0[i].dend.V_dend ;
		neighVdend4[i+Mux_Factor] = local_state1[i].dend.V_dend ;
		neighVdend4[i+2*Mux_Factor] = local_state2[i].dend.V_dend ;
		neighVdend4[i+3*Mux_Factor] = local_state3[i].dend.V_dend ;
		neighVdend4[i+4*Mux_Factor] = local_state4[i].dend.V_dend ;

		k=k+HW_CELLS;
	}

	ComputeOneCellTimeMux0(0,IAppin0,neighVdend0,N_Size,Mux_Factor,Connectivity_Matrix0);
	ComputeOneCellTimeMux1(1,IAppin1,neighVdend1,N_Size,Mux_Factor,Connectivity_Matrix1);
	ComputeOneCellTimeMux2(2,IAppin2,neighVdend2,N_Size,Mux_Factor,Connectivity_Matrix2);
	ComputeOneCellTimeMux3(3,IAppin3,neighVdend3,N_Size,Mux_Factor,Connectivity_Matrix3);
	ComputeOneCellTimeMux4(4,IAppin4,neighVdend4,N_Size,Mux_Factor,Connectivity_Matrix4);


	ComputeNetwork_label8:	for(j=0;j<Mux_Factor;j++){
		local_state0[j]=new_state0[j];
		local_state1[j]=new_state1[j];
		local_state2[j]=new_state2[j];
		local_state3[j]=new_state3[j];
		local_state4[j]=new_state4[j];
	}
	//Create system output(all cell axon currents)
	ComputeNetwork_label9:	for(j=0;j<Mux_Factor;j++){

		cellOut[j+(TIME_MUX_FACTOR*0)] = local_state0[j].axon.V_axon;
		cellOut[j+(TIME_MUX_FACTOR*1)] = local_state1[j].axon.V_axon;
		cellOut[j+(TIME_MUX_FACTOR*2)] = local_state2[j].axon.V_axon;
		cellOut[j+(TIME_MUX_FACTOR*3)] = local_state3[j].axon.V_axon;
		cellOut[j+(TIME_MUX_FACTOR*4)] = local_state4[j].axon.V_axon;
	}



	//return cellOut;

}

void ComputeOneCellTimeMux0(int cluster,  mod_prec iAppin[MAX_TIME_MUX], mod_prec neighVdend[MAX_N_SIZE], int N_Size, int Mux_Factor,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){
	int j,i;
	
	mod_prec * cm_p = Connectivity_Matrix;

	// call cell execution
	ComputeOneCellTimeMux0_label10:for(j=0;j<Mux_Factor;j++){


//		printf("\nCalculating cell %d#%d:\n",cluster, j);
		ComputeOneCell0(cluster,j,iAppin[j],neighVdend,N_Size,cm_p);
		cm_p += IO_NETWORK_SIZE;
		}



}

void ComputeOneCellTimeMux1(int cluster,  mod_prec iAppin[MAX_TIME_MUX], mod_prec neighVdend[MAX_N_SIZE], int N_Size, int Mux_Factor,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){
	int j,i;

	mod_prec * cm_p = Connectivity_Matrix;


	// call cell execution
	ComputeOneCellTimeMux1_label11:for(j=0;j<Mux_Factor;j++){

//		printf("\nCalculating cell %d#%d:\n",cluster, j);
		ComputeOneCell1(cluster,j,iAppin[j],neighVdend,N_Size,cm_p);
		cm_p += IO_NETWORK_SIZE;

		}



}


void ComputeOneCellTimeMux2(int cluster,  mod_prec iAppin[MAX_TIME_MUX], mod_prec neighVdend[MAX_N_SIZE], int N_Size, int Mux_Factor,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){
	int j,i;

	mod_prec * cm_p = Connectivity_Matrix;


	// call cell execution
	ComputeOneCellTimeMux2_label12:for(j=0;j<Mux_Factor;j++){
//		printf("\nCalculating cell %d#%d:\n",cluster, j);
		ComputeOneCell2(cluster,j,iAppin[j],neighVdend,N_Size,cm_p);
		cm_p += IO_NETWORK_SIZE;

		}



}

void ComputeOneCellTimeMux3(int cluster,  mod_prec iAppin[MAX_TIME_MUX], mod_prec neighVdend[MAX_N_SIZE], int N_Size, int Mux_Factor,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){
	int j,i;


	mod_prec * cm_p = Connectivity_Matrix;


	// call cell execution
	ComputeOneCellTimeMux3_label13:for(j=0;j<Mux_Factor;j++){



//		printf("\nCalculating cell %d#%d:\n",cluster, j);
		ComputeOneCell3(cluster,j,iAppin[j],neighVdend,N_Size,cm_p);
		cm_p += IO_NETWORK_SIZE;

		}



}

void ComputeOneCellTimeMux4(int cluster,  mod_prec iAppin[MAX_TIME_MUX], mod_prec neighVdend[MAX_N_SIZE], int N_Size, int Mux_Factor,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){
	int j,i;


	mod_prec * cm_p = Connectivity_Matrix;


	// call cell execution
	ComputeOneCellTimeMux4_label14:for(j=0;j<Mux_Factor;j++){



//		printf("\nCalculating cell %d#%d:\n",cluster, j);
		ComputeOneCell4(cluster,j,iAppin[j],neighVdend,N_Size,cm_p);
		cm_p += IO_NETWORK_SIZE;

		}



}

//Top Inferior Olive Cell compute function including the 3 computational compartments.
void ComputeOneCell0(int i,int j,mod_prec iAppin, mod_prec neighVdend[MAX_N_SIZE], int N_Size,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){

	cellState prevCellState;
	prevCellState = local_state0[j];

    //The three compartments can be computed concurrently but only across a single sim step
	new_state0[j].dend = CompDend(prevCellState.dend, prevCellState.soma.V_soma , iAppin, neighVdend,N_Size,Connectivity_Matrix,j);
	new_state0[j].soma = CompSoma(prevCellState.soma,prevCellState.dend.V_dend, prevCellState.axon.V_axon);
	new_state0[j].axon = CompAxon(prevCellState.axon ,prevCellState.soma.V_soma);



}

void ComputeOneCell1(int i,int j,mod_prec iAppin, mod_prec neighVdend[MAX_N_SIZE], int N_Size,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){

	cellState prevCellState;
	prevCellState = local_state1[j];

    //The three compartments can be computed concurrently but only across a single sim step
	new_state1[j].dend = CompDend(prevCellState.dend, prevCellState.soma.V_soma , iAppin, neighVdend,N_Size,Connectivity_Matrix,j);
	new_state1[j].soma = CompSoma(prevCellState.soma,prevCellState.dend.V_dend, prevCellState.axon.V_axon);
	new_state1[j].axon = CompAxon(prevCellState.axon ,prevCellState.soma.V_soma);



}

void ComputeOneCell2(int i,int j,mod_prec iAppin, mod_prec neighVdend[MAX_N_SIZE], int N_Size,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){

	cellState prevCellState;
	prevCellState = local_state2[j];

    //The three compartments can be computed concurrently but only across a single sim step
	new_state2[j].dend = CompDend(prevCellState.dend, prevCellState.soma.V_soma , iAppin, neighVdend,N_Size,Connectivity_Matrix,j);
	new_state2[j].soma = CompSoma(prevCellState.soma,prevCellState.dend.V_dend, prevCellState.axon.V_axon);
	new_state2[j].axon = CompAxon(prevCellState.axon ,prevCellState.soma.V_soma);



}

void ComputeOneCell3(int i,int j,mod_prec iAppin, mod_prec neighVdend[MAX_N_SIZE], int N_Size,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){

	cellState prevCellState;
	prevCellState = local_state3[j];

    //The three compartments can be computed concurrently but only across a single sim step
	new_state3[j].dend = CompDend(prevCellState.dend, prevCellState.soma.V_soma , iAppin, neighVdend,N_Size,Connectivity_Matrix,j);
	new_state3[j].soma = CompSoma(prevCellState.soma,prevCellState.dend.V_dend, prevCellState.axon.V_axon);
	new_state3[j].axon = CompAxon(prevCellState.axon ,prevCellState.soma.V_soma);



}
void ComputeOneCell4(int i,int j,mod_prec iAppin, mod_prec neighVdend[MAX_N_SIZE], int N_Size,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5]){

	cellState prevCellState;
	prevCellState = local_state4[j];

    //The three compartments can be computed concurrently but only across a single sim step
	new_state4[j].dend = CompDend(prevCellState.dend, prevCellState.soma.V_soma , iAppin, neighVdend,N_Size,Connectivity_Matrix,j);
	new_state4[j].soma = CompSoma(prevCellState.soma,prevCellState.dend.V_dend, prevCellState.axon.V_axon);
	new_state4[j].axon = CompAxon(prevCellState.axon ,prevCellState.soma.V_soma);



}

Dend CompDend(Dend prevDend, mod_prec prevSoma , mod_prec iAppIn,mod_prec neighVdend[MAX_N_SIZE], int N_Size,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5],int j){

	struct Dend d_output;
    struct channelParams chPrms;
    struct dendCurrVoltPrms chComps, newchComps;



    //H current calculation
     chPrms.prevComp1 = prevDend.Hcurrent_q;
     d_output.Hcurrent_q = DendHCurr(chPrms, prevDend.V_dend);
	 //printf("Hcurrent_q: %.9f\n",d_output.Hcurrent_q);
     //Ca current calculation
     chPrms.prevComp1 = prevDend.Calcium_r;
     d_output.Calcium_r  = DendCaCurr(chPrms, prevDend.V_dend);
	 //printf("Calcium_r: %.9f\n",d_output.Calcium_r);
     //K plus current calculation
     chPrms.prevComp1 = prevDend.Potassium_s;
     chPrms.prevComp2 = prevDend.Ca2Plus;
     d_output.Potassium_s = DendKCurr(chPrms );
	 //printf("Potassium_s: %.9f\n",d_output.Potassium_s);
     //K plus current calculation
     chPrms.prevComp1 = prevDend.Ca2Plus;
     chPrms.prevComp2 =prevDend.I_CaH;
     d_output.Ca2Plus  = DendCal(chPrms);
	 //printf("Ca2Plus: %.9f\n",d_output.Ca2Plus);

     //Neighbors Ic accumulation
     chComps.iC = IcNeighbors(neighVdend, prevDend.V_dend,N_Size,Connectivity_Matrix,j);
     //Compartment output calculation
     chComps.iApp = iAppIn;
     chComps.vDend = prevDend.V_dend;
     chComps.vSoma = prevSoma;
     chComps.q = d_output.Hcurrent_q;
     chComps.r = d_output.Calcium_r;
     chComps.s = d_output.Potassium_s;
     newchComps = DendCurrVolt(chComps);
     d_output.I_CaH = newchComps.newI_CaH;
	 //printf("I_CaH: %.9f\n",d_output.I_CaH);
     d_output.V_dend = newchComps.newVDend;
	 //printf("V_dend: %.9f\n",d_output.V_dend);

    return d_output;
}

mod_prec DendHCurr(struct channelParams chPrms , mod_prec prevV_dend){

    mod_prec q_inf, tau_q, dq_dt, q_local;

    //opt for division by constant
    mod_prec const div_four = 0.25;


    //Get inputs
   // mod_prec prevV_dend = chPrms.v;
    mod_prec prevHcurrent_q = chPrms.prevComp1;

    // Update dendritic H current component
    q_inf = 1 /(1 + expf((prevV_dend + 80) * div_four));
    tau_q = 1 /(expf(-0.086 * prevV_dend - 14.6) + expf(0.070 * prevV_dend - 1.87));
    dq_dt = (q_inf - prevHcurrent_q) / tau_q;
    q_local = DELTA * dq_dt + prevHcurrent_q;
    //Put result
    //chPrms.newComp1 = q_local;

    return q_local;
}
mod_prec DendCaCurr(struct channelParams chPrms, mod_prec prevV_dend){

    mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt, r_local;

    //opt for division by constant
    mod_prec const div_five = 0.2;
    mod_prec const div_thirteen = 1/13.9;
    //Get inputs
    //mod_prec prevV_dend = chPrms.v;
    mod_prec prevCalcium_r = chPrms.prevComp1;

    // Update dendritic high-threshold Ca current component
    alpha_r = 1.7 / (1 + expf( -(prevV_dend - 5) *div_thirteen));
    beta_r = 0.02 * (prevV_dend + 8.5) / (expf((prevV_dend + 8.5) * div_five) - 1);
    r_inf = alpha_r / (alpha_r + beta_r);
    tau_r = 5 / (alpha_r + beta_r);
    dr_dt = (r_inf - prevCalcium_r) / tau_r;
    r_local = DELTA * dr_dt + prevCalcium_r;
    //Put result
    //chPrms.newComp1 = r_local;

    return r_local;
}
mod_prec DendKCurr(struct channelParams chPrms){

    mod_prec  alpha_s, beta_s, s_inf, tau_s, ds_dt, s_local;

    //Get inputs
    mod_prec prevPotassium_s = chPrms.prevComp1;
    mod_prec prevCa2Plus = chPrms.prevComp2;

    // Update dendritic Ca-dependent K current component
    alpha_s = min((0.00002*prevCa2Plus), 0.01);
    beta_s = 0.015;
    s_inf = alpha_s / (alpha_s + beta_s);
    tau_s = 1 / (alpha_s + beta_s);
    ds_dt = (s_inf - prevPotassium_s) / tau_s;
    s_local = DELTA * ds_dt + prevPotassium_s;
    //Put result
    //chPrms.newComp1 = s_local;

    return s_local;
}
//Consider merging DendCal into DendKCurr since DendCal's output doesn't go to DendCurrVolt but to DendKCurr
mod_prec DendCal(struct channelParams chPrms){

    mod_prec  dCa_dt, Ca2Plus_local;

    //Get inputs
    mod_prec prevCa2Plus = chPrms.prevComp1;
	//printf("prevCa2Plus: %.9f\n",prevCa2Plus);
    mod_prec prevI_CaH = chPrms.prevComp2;
	//printf("prevI_CaH: %.9f\n",prevI_CaH);

    // update Calcium concentration
    dCa_dt = -3 * prevI_CaH - 0.075 * prevCa2Plus;
    Ca2Plus_local = DELTA * dCa_dt + prevCa2Plus;
    //Put result
    //chPrms.newComp1 = Ca2Plus_local;//This state value is read in DendKCurr

    return Ca2Plus_local;
}

dendCurrVoltPrms DendCurrVolt(struct dendCurrVoltPrms chComps){

    //Local variables
    mod_prec I_sd, I_CaH, I_K_Ca, I_ld, I_h, dVd_dt;

    //Get inputs
    mod_prec I_c = chComps.iC;
//	printf("I_c: %.9f\n",I_c);
    mod_prec I_app = chComps.iApp;
//	printf("I_app: %.9f\n",I_app);
    mod_prec prevV_dend = chComps.vDend;
//	printf("prevV_dend: %.9f\n",prevV_dend);
    mod_prec prevV_soma = chComps.vSoma;
//	printf("prevV_soma: %.9f\n",prevV_soma);
    mod_prec q = chComps.q;
//	printf("q: %.9f\n",q);
    mod_prec r = chComps.r;
//	printf("r: %.9f\n",r);
    mod_prec s = chComps.s;
//	printf("s: %.9f\n",s);
    // DENDRITIC CURRENTS

    // Soma-dendrite interaction current I_sd
    I_sd   = (G_INT / (1 - P1)) * (prevV_dend - prevV_soma);
    // Inward high-threshold Ca current I_CaH
    I_CaH  =  G_CAH * r * r * (prevV_dend - V_CA);
    // Outward Ca-dependent K current I_K_Ca
    I_K_Ca =  G_K_CA * s * (prevV_dend - V_K);
    // Leakage current I_ld
    I_ld   =  G_LD * (prevV_dend - V_L);
    // Inward anomalous rectifier I_h
    I_h    =  G_H * q * (prevV_dend - V_H);

    dVd_dt = (-(I_CaH   + I_sd  + I_ld + I_K_Ca + I_c + I_h) + I_app) / C_M;

    //Put result (update V_dend)
    chComps.newVDend = DELTA * dVd_dt + prevV_dend;
    chComps.newI_CaH = I_CaH;//This is a state value read in DendCal
    return chComps;
}
mod_prec IcNeighbors(mod_prec neighVdend[MAX_N_SIZE], mod_prec prevV_dend, int N_Size ,mod_prec Connectivity_Matrix[CONN_MATRIX_MAX/5],int j){

    int i,Bit_Indicator,pos, Integer_Indicator, Array_Fragment;
    mod_prec f, V, I_c, V_acc, F_acc;
    //printf("Ic[0]= %f\n", neighVdend[0]);

    //opt division by constant
    mod_prec const hundred = -1/100.0;

    I_c = 0;
    Array_Fragment = MAX_N_SIZE*j;

	//printf("Slice of connectivity matrix:\n");
IcNeighbors_label0:	for(i=0;i<N_Size;i++){
				//		printf("%.2f ",Connectivity_Matrix[i]);
						V = prevV_dend - neighVdend[i];
						//printf("V%d: %.9f\n",i,V);
						//printf("%f - %f\n",prevV_dend, neighVdend[i]);
						//f = 0.8 * exp(-1*pow(V,2)/100) + 0.2;    // SCHWEIGHOFER 2004 VERSION
        				//I_c = I_c + (Connectivity_Matrix[i] * f * V);
        				f= V*expf(V*V*hundred);
        				F_acc += f*Connectivity_Matrix[i] ;
        				V_acc += V*Connectivity_Matrix[i] ;
			//		printf("USING: %.2f\n",Connectivity_Matrix[i]);
					}
	//				printf("\n");

	I_c = (0.8 * F_acc + 0.2 * V_acc);

    return I_c;
}


Soma CompSoma(Soma prevSoma , mod_prec prevDend, mod_prec prevAxon){

	struct Soma s_output;
    struct channelParams chPrms, chPrmsNew;
    struct somaCurrVoltPrms chComps;

    // update somatic components
    // SCHWEIGHOFER:


    chPrms.v = prevSoma.V_soma;
    chPrms.prevComp1 = prevSoma.Calcium_k;
    chPrms.prevComp2 = prevSoma.Calcium_l;

    //Ca current calculation
    chPrmsNew = SomaCalcium(chPrms);

    s_output.Calcium_k = chPrmsNew.newComp1;
    s_output.Calcium_l = chPrmsNew.newComp2;
//	printf("Calcium_k: %.9f\n",s_output.Calcium_k);
//	printf("Calcium_l: %.9f\n",s_output.Calcium_l);

    chPrms.v = prevSoma.V_soma;
    chPrms.prevComp1 = prevSoma.Sodium_m;
    chPrms.prevComp2 = prevSoma.Sodium_h;

    //Sodium current calculation
    chPrmsNew = SomaSodium(chPrms);

    s_output.Sodium_m = chPrmsNew.newComp1;
    s_output.Sodium_h = chPrmsNew.newComp2;
//	printf("Sodium_m: %.9f\n",s_output.Sodium_m);
//	printf("Sodium_h: %.9f\n",s_output.Sodium_h);


    chPrms.v = prevSoma.V_soma;
    chPrms.prevComp1 = prevSoma.Potassium_n;
    chPrms.prevComp2 = prevSoma.Potassium_p;

    //Potassium current calculation
    chPrmsNew = SomaPotassium(chPrms);

    s_output.Potassium_n = chPrmsNew.newComp1;
    s_output.Potassium_p = chPrmsNew.newComp2;

//	printf("Potassium_n: %.9f\n",s_output.Potassium_n);
//	printf("Potassium_p: %.9f\n",s_output.Potassium_p);

    chPrms.v = prevSoma.V_soma;
    chPrms.prevComp1 = prevSoma.Potassium_x_s;

    //Potassium X current calculation
    chPrmsNew.newComp1 = SomaPotassiumX(chPrms);

    s_output.Potassium_x_s = chPrmsNew.newComp1;
//	printf("Potassium_x_s: %.9f\n",s_output.Potassium_x_s);
    //Compartment output calculation
    chComps.g_CaL = prevSoma.g_CaL;
    s_output.g_CaL = prevSoma.g_CaL;
    chComps.vDend = prevDend;
    chComps.vSoma = prevSoma.V_soma;
    chComps.vAxon = prevAxon;
    chComps.k = s_output.Calcium_k;
    chComps.l = s_output.Calcium_l;
    chComps.m = s_output.Sodium_m;
    chComps.h = s_output.Sodium_h;
    chComps.n = s_output.Potassium_n;
    chComps.x_s = s_output.Potassium_x_s;
    chComps.newVSoma = SomaCurrVolt(chComps);

    s_output.g_CaL = prevSoma.g_CaL;
    s_output.V_soma = chComps.newVSoma;
	//printf("V_soma: %.9f\n",s_output.V_soma);
    return s_output;
}
channelParams SomaCalcium(struct channelParams chPrms){

    mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, k_local, l_local;

    //Get inputs
    mod_prec prevV_soma = chPrms.v;
    mod_prec prevCalcium_k = chPrms.prevComp1;
    mod_prec prevCalcium_l = chPrms.prevComp2;


    //opt for division with constant
    mod_prec const four = 1/4.2;
    mod_prec const eight = 1/8.5; //no change in area by synthesizer on V7
    mod_prec const thirty = 1/30.0;  //no change in area by synthesizer on V7
    mod_prec const seven = 1/7.3; //no change in area by synthesizer on V7


    k_inf = (1 / (1 + expf(-1 * (prevV_soma + 61)  * four)));
    l_inf = (1 / (1 + expf((     prevV_soma + 85.5) * eight)));
    tau_k = 1;
    tau_l = ((20 * expf((prevV_soma + 160) * thirty) / (1 + expf((prevV_soma + 84) * seven ))) +35);
    dk_dt = (k_inf - prevCalcium_k) / tau_k;
    dl_dt = (l_inf - prevCalcium_l) / tau_l;
    k_local = DELTA * dk_dt + prevCalcium_k;
    l_local = DELTA * dl_dt + prevCalcium_l;
    //Put result
    chPrms.newComp1= k_local;
    chPrms.newComp2= l_local;

    return chPrms;
}
channelParams SomaSodium(struct channelParams chPrms){

    mod_prec m_inf, h_inf, tau_h, dh_dt, m_local, h_local;

    //Get inputs
    mod_prec prevV_soma = chPrms.v;
    //mod_prec prevSodium_m = *chPrms->prevComp1;
    mod_prec prevSodium_h = chPrms.prevComp2;

    //opt for division by constant
    mod_prec const five_five = 1/5.5; //no change in area by synthesizer on V7
    mod_prec const five_eight = -1/5.8; //no change in area by synthesizer on V7
    mod_prec const thirty_three = 1/33.0;  //no change in area by synthesizer on V7

    // RAT THALAMOCORTICAL SODIUM:


    m_inf   = 1 / (1 + (expf((-30 - prevV_soma)* five_five)));
    h_inf   = 1 / (1 + (expf((-70 - prevV_soma)* five_eight)));
    tau_h   =       3 * expf((-40 - prevV_soma)* thirty_three);
    dh_dt   = (h_inf - prevSodium_h)/tau_h;
    m_local       = m_inf;
    h_local       = prevSodium_h + DELTA * dh_dt;
    //Put result
    chPrms.newComp1 = m_local;
    chPrms.newComp2 = h_local;

    return chPrms;
}
channelParams SomaPotassium(struct channelParams chPrms){

    mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt, n_local, p_local;

    //Get inputs
    mod_prec prevV_soma = chPrms.v;
    mod_prec prevPotassium_n = chPrms.prevComp1;
    mod_prec prevPotassium_p = chPrms.prevComp2;

    //opt for division by constant
    mod_prec const ten = 0.1; //precision error if 1/10 is used
    mod_prec const twelve = -1/12.0; //no change in area by synthesizer on V7
    mod_prec const nine_hundred = 1/900.0;  //no change in area by synthesizer on V7


    // NEOCORTICAL
    n_inf = 1 / (1 + expf( ( -3 - prevV_soma) * ten));
    p_inf = 1 / (1 + expf( (-51 - prevV_soma) * twelve));
    tau_n =   5 + (  47 * expf( -(-50 - prevV_soma) * nine_hundred));
    tau_p = tau_n;
    dn_dt = (n_inf - prevPotassium_n) / tau_n;
    dp_dt = (p_inf - prevPotassium_p) / tau_p;
    n_local = DELTA * dn_dt + prevPotassium_n;
    p_local = DELTA * dp_dt + prevPotassium_p;
    //Put result
    chPrms.newComp1 = n_local;
    chPrms.newComp2 = p_local;

    return chPrms;
}
mod_prec SomaPotassiumX(struct channelParams chPrms){

    mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, x_s_local;

    //Get inputs
    mod_prec prevV_soma = chPrms.v;
    mod_prec prevPotassium_x_s = chPrms.prevComp1;

    //opt for division by constant
     mod_prec const ten = 0.1; //no change in area by synthesizer on V7

    // Voltage-dependent (fast) potassium
    alpha_x_s = 0.13 * (prevV_soma + 25) / (1 - expf(-(prevV_soma + 25) * ten));
    beta_x_s  = 1.69 * expf(-0.0125 * (prevV_soma + 35));
    x_inf_s   = alpha_x_s / (alpha_x_s + beta_x_s);
    tau_x_s   =         1 / (alpha_x_s + beta_x_s);
    dx_dt_s   = (x_inf_s - prevPotassium_x_s) / tau_x_s;
    x_s_local       = 0.05 * dx_dt_s + prevPotassium_x_s;
    //Put result
    chPrms.newComp1 = x_s_local;

    return chPrms.newComp1;
}
mod_prec SomaCurrVolt(struct somaCurrVoltPrms chComps){

    //Local variables
    mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

    //Get inputs
    mod_prec g_CaL = chComps.g_CaL;
    mod_prec prevV_dend = chComps.vDend;
    mod_prec prevV_soma = chComps.vSoma;
    mod_prec prevV_axon = chComps.vAxon;
    mod_prec k = chComps.k;
    mod_prec l = chComps.l;
    mod_prec m = chComps.m;
    mod_prec h = chComps.h;
    mod_prec n = chComps.n;
    mod_prec x_s = chComps.x_s;

    // SOMATIC CURRENTS

    // Dendrite-soma interaction current I_ds
    I_ds  = (G_INT / P1) * (prevV_soma - prevV_dend);
    // Inward low-threshold Ca current I_CaL
    I_CaL = g_CaL * k * k * k * l * (prevV_soma - V_CA); //k^3
    // Inward Na current I_Na_s
    I_Na_s  = G_NA_S * m * m * m * h * (prevV_soma - V_NA);
    // Leakage current I_ls
    I_ls  = G_LS * (prevV_soma - V_L);
    // Outward delayed potassium current I_Kdr
    I_Kdr_s = G_KDR_S * n * n * n * n * (prevV_soma - V_K); // SCHWEIGHOFER
    // I_K_s
    I_K_s   = G_K_S * x_s*x_s*x_s*x_s * (prevV_soma - V_K);
    // Axon-soma interaction current I_as
    I_as    = (G_INT / (1 - P2)) * (prevV_soma - prevV_axon);

    dVs_dt = (-(I_CaL   + I_ds  + I_as + I_Na_s + I_ls   + I_Kdr_s + I_K_s)) / C_M;
    chComps.newVSoma = DELTA * dVs_dt + prevV_soma;

    return chComps.newVSoma;
}
Axon CompAxon(Axon prevAxon, mod_prec prevSoma){

	struct Axon a_output;
    struct channelParams chPrms, chPrmsNew;
    struct axonCurrVoltPrms chComps;

    // update somatic components
    // SCHWEIGHOFER:

    chPrms.v = prevAxon.V_axon;
    chPrms.prevComp1 =prevAxon.Sodium_h_a;

    //Sodium current calculation
    chPrmsNew = AxonSodium(chPrms);

    a_output.Sodium_h_a =  chPrmsNew.newComp1;
    a_output.Sodium_m_a =  chPrmsNew.newComp2;
	//printf("Sodium_h_a: %.9f\n",a_output.Sodium_h_a);
	//printf("Sodium_m_a: %.9f\n",a_output.Sodium_m_a);

    chPrms.v = prevAxon.V_axon;
    chPrms.prevComp1 = prevAxon.Potassium_x_a;

    //Potassium current calculation
    chPrmsNew.newComp1 = AxonPotassium(chPrms);

    a_output.Potassium_x_a =  chPrmsNew.newComp1;
	//printf("Potassium_x_a: %.9f\n",a_output.Potassium_x_a);

    chComps.vSoma = prevSoma;
    chComps.vAxon = prevAxon.V_axon;
    //Compartment output calculation
    chComps.m_a = a_output.Sodium_m_a;
    chComps.h_a = a_output.Sodium_h_a;
    chComps.x_a = a_output.Potassium_x_a;
    chComps.newVAxon = AxonCurrVolt(chComps);


    a_output.V_axon =  chComps.newVAxon;
	//printf("V_axon: %.9f\n",a_output.V_axon);


    return a_output;
}

channelParams AxonSodium(struct channelParams chPrms){

    mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a, m_a_local, h_a_local;

    //Get inputs
    mod_prec prevV_axon = chPrms.v;
    mod_prec prevSodium_h_a = chPrms.prevComp1;

    //opt for division by constant
        mod_prec const five_five = 1/5.5; //no change in area by synthesizer on V7
        mod_prec const five_eight = -1/5.8; //no change in area by synthesizer on V7
        mod_prec const thirty_three = 1/33.0;  //no change in area by synthesizer on V7


    // Update axonal Na components
    // NOTE: current has shortened inactivation to account for high
    // firing frequencies in axon hillock


    m_inf_a   = 1 / (1 + (expf((-30 - prevV_axon) * five_five)));
    h_inf_a   = 1 / (1 + (expf((-60 - prevV_axon) * five_eight)));
    tau_h_a   =     1.5 * expf((-40 - prevV_axon) * thirty_three);
    dh_dt_a   = (h_inf_a - prevSodium_h_a)/tau_h_a;
    m_a_local = m_inf_a;
    h_a_local = prevSodium_h_a + DELTA * dh_dt_a;
    //Put result
    chPrms.newComp1 = h_a_local;
    chPrms.newComp2 = m_a_local;

    return chPrms;
}
mod_prec AxonPotassium(struct channelParams chPrms){

    mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, x_a_local;

    //Get inputs
    mod_prec prevV_axon = chPrms.v;
    mod_prec prevPotassium_x_a = chPrms.prevComp1;

    //opt for division by constant
    mod_prec const ten = 0.1; //no change in area by synthesizer on V7

    // D'ANGELO 2001 -- Voltage-dependent potassium
    alpha_x_a = 0.13 * (prevV_axon + 25) / (1 - expf(-(prevV_axon + 25) * ten));
    beta_x_a  = 1.69 * expf(-0.0125 * (prevV_axon + 35));
    x_inf_a   = alpha_x_a / (alpha_x_a + beta_x_a);
    tau_x_a   =         1 / (alpha_x_a + beta_x_a);
    dx_dt_a   = (x_inf_a - prevPotassium_x_a) / tau_x_a;
    x_a_local = 0.05 * dx_dt_a + prevPotassium_x_a;
    //Put result
    chPrms.newComp1 = x_a_local;

    return chPrms.newComp1;
}
mod_prec AxonCurrVolt(struct axonCurrVoltPrms chComps){

    //Local variable
    mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

    //Get inputs
    mod_prec prevV_soma = chComps.vSoma;
    mod_prec prevV_axon = chComps.vAxon;
    mod_prec m_a = chComps.m_a;
    mod_prec h_a = chComps.h_a;
    mod_prec x_a = chComps.x_a;

    // AXONAL CURRENTS
    // Sodium
    I_Na_a  = G_NA_A  * m_a * m_a * m_a * h_a * (prevV_axon - V_NA);
    // Leak

    I_la    = G_LA * (prevV_axon - V_L);
    // Soma-axon interaction current I_sa
    I_sa    = (G_INT / P2) * (prevV_axon - prevV_soma);
    // Potassium (transient)
    I_K_a   = G_K_A * x_a *x_a *x_a *x_a * (prevV_axon - V_K);
    dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;
    chComps.newVAxon = DELTA * dVa_dt + prevV_axon;

    return chComps.newVAxon;
}



inline mod_prec min(mod_prec a, mod_prec b){
    return (a < b) ? a : b;
}

bool Check_bit(int var, int pos){
	bool connection;

	connection = ((var) & (1<<(pos)));

	return connection;
}
