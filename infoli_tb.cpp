#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include "infoli.h"

int main(int argc, char *argv[]){

    char *inFileName;
    char *conFileName;
    char *outFileName = "InferiorOlive_Output.txt";
    FILE *pInFile;
    FILE *pOutFile;
    FILE *conInFile;
    char *iAppBuf;
    const int iAppBufSize =  IAPP_MAX_CHARS*IO_NETWORK_SIZE+1;
    mod_prec iAppArray[IO_NETWORK_SIZE];
    int i, j, k, p, q, n;
    bool ini , new_matrix;
    int simSteps = 0;
    int simTime = 0;
    int inputFromFile = 0;
	bool connectivityMatrixInput = false;
    int initSteps;
    returnState cellOut_1, cellOut_2;
	float result_1, result_2, percDifference;
	float cumulDiff = 0;
	int cumulCounter = 0;
    cellState IniArray[IO_NETWORK_SIZE];
    cellCompParams cellCompParamsPtr;
    int seedvar;
	float temp_cond_value;
	mod_prec cond_value;
    char temp[100];//warning: this buffer may overflow
    mod_prec iApp;
    mod_prec Connectivity_Matrix[CONN_MATRIX_SIZE];
    //timestamp_t t0, t1, secs;
    //double secs;


    printf("Inferior Olive Model (%d cell network)\n", IO_NETWORK_SIZE);

    //Process command line arguments
	switch(argc) {
		case 1 :
        	inputFromFile = 0;
        	printf("Warning: No input file has been specified. A one-pulse input will be used.\n");
			break;
 		case 3 : 
			connectivityMatrixInput = true;
			conFileName = argv[2];
	    	conInFile = fopen(conFileName,"r");
			if (conInFile==NULL) {
				printf("Error: Could not open %s\n", conFileName);
				exit(EXIT_FAILURE);
			}
			printf("Using file %s as connectivity matrix\n", conFileName);
			// no break
		case 2 :
        	inputFromFile = 1;
        	inFileName = argv[1];//comment out for a hardcoded name
        	pInFile = fopen(inFileName,"r");
        	if(pInFile==NULL){
            	printf("Error: Couldn't open %s\n", inFileName);
            	exit(EXIT_FAILURE);
        	}
			printf("Using file %s as input file\n", inFileName);
			break;
		default :
        	printf("Error: Too many arguments.\nUsage: /InferiorOlive.x <Iapp_input_file> <Connection_matrix_file>, ./InferiorOlive.x <Iapp_input_file>, or ./InferiorOlive.x\n");
        	exit(EXIT_FAILURE);
    }

    //Open output file
    pOutFile = fopen(outFileName,"w");
    if(pOutFile==NULL){
        printf("Error: Couldn't create %s\n", outFileName);
        exit(EXIT_FAILURE);
    }
    sprintf(temp, "#simSteps Time(ms) Input(Iapp) Output(V_axon)");
    fputs(temp, pOutFile);
	for (j=0; j<IO_NETWORK_SIZE; ++j) { // analogy with original code
		sprintf(temp, "%d ", j);
		fputs(temp,pOutFile);
	}
	fputs("\n",pOutFile);
    //Malloc for iAppBuffer holding iApp arrays, one 2D array (a single line in the file though) at the time
    printf("Malloc'ing memory...\n");
    printf("iAppBuf: %dB\n", iAppBufSize);
    iAppBuf = (char *)malloc(iAppBufSize);
    //printf("%p\n",iAppBuf);
    if(iAppBuf==NULL){
        printf("Error: Couldn't malloc for iAppBuf\n");
        exit(EXIT_FAILURE);
    }

    for(j=0;j<IO_NETWORK_SIZE;j++){


    	IniArray[j] = InitState();

        }

    //Initialize g_CaL
    seedvar = 1;
    for(j=0;j<IO_NETWORK_SIZE;j++){
            srand(seedvar++);   // use this for debugging, now there is difference
            IniArray[j].soma.g_CaL = 0.68;
    }

    //initialize connection Matrix
	int cnt = 0;
    for (j=0;j<CONN_MATRIX_SIZE; j++){
		if (connectivityMatrixInput) {
			fscanf(conInFile, "%f", &temp_cond_value);
			Connectivity_Matrix[j] = (mod_prec)temp_cond_value;
		} else {
    		Connectivity_Matrix[j] = CONDUCTANCE;
		}
		printf("%.2f ",Connectivity_Matrix[j]);
		++cnt;
		if (cnt >= IO_NETWORK_SIZE) {
			printf("\n");
			cnt = 0;
		}
	}
	j = 2;
    simTime = SIMTIME; // in miliseconds
    if(inputFromFile){
	printf("Reading from file...");
        simSteps = 0;
        //Read full lines until end of file. Every iteration (line) is one simulation step.
        while(ReadFileLine(iAppBuf, iAppBufSize, pInFile, iAppArray)){
            if (simSteps==0){
              	ini=1;
               	new_matrix = "true";
            } else {
               	ini=0;
               	new_matrix = "false";
            }
            //Compute one sim step for all cells
            ComputeNetwork_old(ini,new_matrix, IniArray, iAppArray,IO_NETWORK_SIZE,TIME_MUX_FACTOR,Connectivity_Matrix,CONN_MATRIX_SIZE,cellOut_1.axonOut);
            ComputeNetwork_new(ini,new_matrix, IniArray, iAppArray,IO_NETWORK_SIZE,TIME_MUX_FACTOR,Connectivity_Matrix,CONN_MATRIX_SIZE,cellOut_2.axonOut);
            for(j=0;j<IO_NETWORK_SIZE;j++){
					result_1 = cellOut_1.axonOut[j];
					result_2 = cellOut_1.axonOut[j];
					percDifference = abs((result_2-result_1)/result_1)*100;
					sprintf(temp,"%.16f <==> %.16f (%.16f%%)\n", cellOut_1.axonOut[j], cellOut_2.axonOut[j], percDifference);
                    fputs(temp, pOutFile);
					cumulDiff += percDifference;
					++cumulCounter;
            }
            simSteps++;
        }
    }else{
	printf("No input file. Using one-pulse input...\n");
        simSteps = ceil(simTime/DELTA);

        for(i=0;i<simSteps;i++){
            //Compute one sim step for all cells
            //printf("simSteps: %d\n", i);
            if(i>20000-1 && i<20500-1){ iApp = 6;} // start @ 1 because skipping initial values
            else{ iApp = 0;}


           for(j=0;j<IO_NETWORK_SIZE;j++){

                	 iAppArray[j] = iApp;

           }
           sprintf(temp, "%d %.2f %.1f ", i+1, i*0.05,  iAppArray[0]); // start @ 1 because skipping initial values
           fputs(temp, pOutFile);

                    n = 0;

                    //Compute Network...
                    if (i==0){
                    	ini=1;
                    	new_matrix = "true";
                    }
                    else{
                    	ini=0;
                    	new_matrix = "false";
                    }
                    ComputeNetwork_old(ini,new_matrix, IniArray, iAppArray,IO_NETWORK_SIZE,TIME_MUX_FACTOR,Connectivity_Matrix,CONN_MATRIX_SIZE,cellOut_1.axonOut);

                   for(j=0;j<IO_NETWORK_SIZE;j++){
                	   sprintf(temp, "%d: %.8f ",j,cellOut_1.axonOut[j]);
                	   fputs(temp, pOutFile);
                   }

                    //sprintf(temp, "%.8f  %.8f  %.8f %.8f  %.8f  %.8f %.8f  %.8f  %.8f", cellCompParamsPtr.newCellState.soma.V_soma , cellCompParamsPtr.newCellState.soma.Calcium_k, cellCompParamsPtr.newCellState.soma.Calcium_l, cellCompParamsPtr.newCellState.soma.Potassium_n ,cellCompParamsPtr.newCellState.soma.Potassium_p, cellCompParamsPtr.newCellState.soma.Potassium_x_s, cellCompParamsPtr.newCellState.soma.Sodium_h, cellCompParamsPtr.newCellState.soma.Sodium_m, cellCompParamsPtr.newCellState.soma.g_CaL);
                    //fputs(temp, pOutFile);
                    //sprintf(temp, "\n%.8f  %.8f  %.8f ", cellCompParamsPtr.newCellState.axon.Potassium_x_a , cellCompParamsPtr.newCellState.axon.Sodium_h_a,cellCompParamsPtr.newCellState.axon.Sodium_m_a);
                    //fputs(temp, pOutFile);
                    //cellCompParamsPtr[j][k].prevCellState = cellCompParamsPtr[j][k].newCellState;

                   sprintf(temp, "\n");
               	   fputs(temp, pOutFile);
        }
    }

    printf("%d ms of brain time in %d simulation steps\n", simTime, simSteps);
    free(iAppBuf);
    fclose (pOutFile);
    if(inputFromFile){ fclose (pInFile);}
	printf("Average difference between results: %.16f%%\n",cumulDiff/((float)cumulCounter));
    return 0;
}

int ReadFileLine(char *iAppBuf, int iAppBufSize, FILE *pInFile, mod_prec *iAppArray){
    //FIXME: make this function more robust
    char *strNumber;
    int i = 0;
    //Get one line
    if(fgets(iAppBuf, iAppBufSize, pInFile)){
        //Convert the ASCII string of one element to a double precision floating point value
		//printf("p: %p::%s\n",&iAppBuf,iAppBuf);
        strNumber = strtok(iAppBuf," ");
        i = 0;
        //printf("Line:\n");
        while ((strNumber != NULL) && (i<IO_NETWORK_SIZE)){
            iAppArray[i] = atof(strNumber);//atof() should change if using integers or fixed point
            //printf ("(%s) %0.2f ", strNumber, iAppArray[i]);
            strNumber = strtok(NULL, " ");
            i++;
        }
        if(i<IO_NETWORK_SIZE){
            //BUG: if only one element is missing but the line ends in a space, the error is not detected
            printf("Error: Input line doesn't have enough elements, only %d\n", i);
            exit(EXIT_FAILURE);
        }
        return 1;//success
    }else{
        if(!feof(pInFile)){
        printf("Error: Reading from input file didn't finish successfully\n");
        exit(EXIT_FAILURE);
        }
        return 0;//end of file
    }
}

cellState InitState(){
    //int j, k;
    cellState initState;
    //Initial dendritic parameters
    initState.dend.V_dend = -60;
    initState.dend.Calcium_r = 0.0112788;// High-threshold calcium
    initState.dend.Potassium_s = 0.0049291;// Calcium-dependent potassium
    initState.dend.Hcurrent_q = 0.0337836;// H current
    initState.dend.Ca2Plus = 3.7152;// Calcium concentration
    initState.dend.I_CaH   = 0.5;// High-threshold calcium current
    //Initial somatic parameters
    initState.soma.g_CaL = 0.68; //default arbitrary value but it should be randomized per cell
    initState.soma.V_soma = -60;
    initState.soma.Sodium_m = 1.0127807;// Sodium (artificial)
    initState.soma.Sodium_h = 0.3596066;
    initState.soma.Potassium_n = 0.2369847;// Potassium (delayed rectifier)
    initState.soma.Potassium_p = 0.2369847;
    initState.soma.Potassium_x_s = 0.1;// Potassium (voltage-dependent)
    initState.soma.Calcium_k = 0.7423159;// Low-threshold calcium
    initState.soma.Calcium_l = 0.0321349;
    // Initial axonal parameters
    initState.axon.V_axon = -60;
    //sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
    initState.axon.Sodium_m_a = 0.003596066;// Sodium (thalamocortical)
    initState.axon.Sodium_h_a = 0.9;
    initState.axon.Potassium_x_a = 0.2369847;// Potassium (transient)

    //Copy init sate to all cell states
  //  for(j=0;j<IO_NETWORK_DIM1;j++){
  //     for(k=0;k<IO_NETWORK_DIM2;k++){
  //      	cellCompParamsPtr[j][k].prevCellState = initState;
  //      }
  //  }

    return (initState);
}
