#include "phaseCodeCircuit.h"
#include "PerfectMatching.h"

//--------------------------------------------------------------
//---------------------- START OF main( )  ----------------------
//--------------------------------------------------------------
int main ( int argc, char *argv[] )
{
    
    time_t rSeed;
    rSeed = time (NULL);
    unsigned long labelNumber=0;
    
    int numQ=196;
    long int numCycle=1000;
    double errorRateGate=0.014;
    double errorRateMeasure=0.014;
    double errorRateDetection=0.014;
    double weight=4.5;
    double weightTime=0.85;
    double thresholdr=3;
    double weightZ=1.1;
    int frequency=2;
    double biasratio=0;
    if (argc<3){ //this message will print if there is not at least one command line parameter
        printf("\nINFORMATION ON COMMAND LINE PARAMETERS:\nA usage example is './thisProg -r=5 -i=200 -em=0.1 -eg=0.1 -ep=0.1 -t=200'\n");
        exit(1);
    }
    
    int i;
    for (i=1; i<argc;  i++){
        if (argv[i][0]=='-'){
            int j=1;
            while (argv[i][j]=='-') j++;
            int failedToParse=1;
            char thisFlag=argv[i][j]; j++;
            if (argv[i][j]=='='){
                while (j>=0){ 		//this loop erases the '-b=' in a string like '-b=1.3' so that only '   1.3' remains, so it can be converted to a number cleanly
                    argv[i][j]=' ';
                    j--;
                }
                if (thisFlag=='i'){
                    numCycle=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='u'){
                    numQ=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='g'){
                    errorRateGate=strtof(argv[i], NULL);
                    failedToParse=0;
                }
//                if (thisFlag=='d'){
//                    errorRateDetection=strtof(argv[i], NULL);
//                    failedToParse=0;
//                }
                if (thisFlag=='r'){
                    labelNumber=atoi(argv[i])*1000;
                    failedToParse=0;
                }
                if (thisFlag=='f'){
                    frequency=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='w'){
                    weight=strtof(argv[i], NULL);
                    failedToParse=0;
                }
                if (thisFlag=='z'){
                    weightZ=strtof(argv[i], NULL);
                    failedToParse=0;
                }
                if (thisFlag=='v'){
                    weightTime=strtof(argv[i], NULL);
                    failedToParse=0;
                }
                if (thisFlag=='b'){
                    biasratio=strtof(argv[i], NULL);
                    failedToParse=0;
                }
                if (thisFlag=='e'){
                    thresholdr=strtof(argv[i], NULL);
                    failedToParse=0;
                }
            }
            if (failedToParse==1){ printf("\nCannot interpret command line argument %s. (Note use '=' as in '-s=100'). ABORTING.\n",argv[i]); exit(1); }
            
        }else{ // an argument not starting with - detected
            printf("A command line parameter %s seen; cannot interpret this. ABORTING.\n",argv[i]); exit(1);
        }
    }
    double ra=(labelNumber+rSeed)%10000;
    srand(ra);
    errorRateMeasure=errorRateGate;
    errorRateDetection=errorRateGate/3;
    double flipRatio=1;
    
    //main simulation
    int n= sqrt(numQ);
    int numRun=2*n;
    double threshold=thresholdr;
    int index; int jndex; int p; int q; int t;
    clock_t clock(void) ;
    clock_t startTime,endTime;
    
    //space allocation
    //Store global information
    int *globalStore;
    globalStore=(int *)calloc(numQ*4,sizeof(int));
    if (globalStore == NULL){
        printf("space failed to be allocated！1\n"); exit(1);}
    int *copyX;
    copyX=(int *)calloc(numQ,sizeof(int));
    if (copyX == NULL){
        printf("space failed to be allocated！12\n"); exit(1);}
    int *copyZ;
    copyZ=(int *)calloc(numQ,sizeof(int));
    if (copyZ == NULL){
        printf("space failed to be allocated！3\n"); exit(1);}
    
    int rightX1=0; int rightZ1=0; int right1=0;
    int rightX2=0; int rightZ2=0; int right2=0;
    
    //-----------------------------------------------------
    long int cycle;
    for (cycle=0;cycle<numCycle;cycle++){
        
        int labelTableX[8000]={0};//store labels of syndrome changes
        int labelTableZ[8000]={0};
        int labelX=0; //the number of syndrome changes
        int labelZ=0;
        int errorListZ[8000]={0};
        
        int *matrixX;//store all possible errors
        matrixX=(int *)calloc(numQ*(numRun+1),sizeof(int));
        if (matrixX == NULL){
            printf("space failed to be allocated！4\n"); exit(1);}
        
        for (p=0;p<4*numQ;p++)
            globalStore[p]=0;
        for (p=0;p<numQ;p++){
            copyZ[p]=0;
            copyX[p]=0;
        }
       
        //begin stabiliser check
        int count=0; int countx=0; int countz=0;
        while (count<numRun){
            
//            while((count+1)%frequency!=0 && count<numRun){
                q=0; p=0;
                while (p<(n/2)){
                    for (q=0;q<n/2;q++){
                        t=n+1+2*n*p+2*q;
                        checkNewTwoGate(t,0,globalStore,copyX,numQ,&labelX,labelTableX,countx,1,matrixX,errorRateGate,errorRateMeasure,errorRateDetection,biasratio,flipRatio);
                        matrixX[countx*numQ+t-1]=2;
                    }
                    p++;
                }
                count++;
                countx++;
            
//            }
//
//            if (count<numRun){
//                
//                q=0; p=0;
//                while (p<(n/2)){
//                    for (q=0;q<n/2;q++){
//                        t=2*n*p+2*q;
//                        checkNewTwoGate(t,1,globalStore,copyZ,numQ,&labelZ,labelTableZ,countz,1,matrixX,errorRateGate,errorRateMeasure,errorRateDetection,biasratio,flipRatio);
//                    }
//                    p++;
//                }
//                count++;
//                countz++;
//            }
        }
        
        
        //the last cycle
        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2;q++){
                t=n+1+2*n*p+2*q;
                checkNewTwoGate(t,0,globalStore,copyX,numQ,&labelX,labelTableX,countx,0,matrixX,errorRateGate,errorRateMeasure,0,biasratio,flipRatio);
                matrixX[countx*numQ+t-1]=2;
            }
            p++;
        }
//        q=0; p=0;
//        while (p<(n/2)){
//            for (q=0;q<n/2;q++){
//                t=2*n*p+2*q;
//                checkNewTwoGate(t,1,globalStore,copyZ,numQ,&labelZ,labelTableZ,countz,0,matrixX,errorRateGate,errorRateMeasure,errorRateDetection,biasratio,flipRatio);
//            }
//            p++;
//        }
    
        if (labelX>0){
            
            double *tableX;
            tableX=(double *)calloc(labelX*labelX*3,sizeof(double));
            if (tableX == NULL){
                printf("space failed to be allocated！5\n"); exit(1);}
            
            int *path;
            path=(int *)calloc(labelX*(labelX-1)*15,sizeof(int));
            if (path == NULL){
                printf("space failed to be allocated！6\n"); exit(1);}
            for (index=0;index<labelX*(labelX-1)*15;index++)
                path[index]=-1;
            
            int numNodesX=0;

            startTime += clock();
            int numEdgesX=creatTableNewFourAdapted(labelTableX, labelX, tableX, numQ,matrixX,weight,weightTime,countx+1,path,threshold);
            endTime += clock();
            
            numNodesX=labelX;
            int workEdgeX=numEdgesX*2;

            //the lines below initialise the system for perfect matching
            struct PerfectMatching::Options options;
            PerfectMatching *pmX;
            pmX= new PerfectMatching(numNodesX, workEdgeX);
            
            //the lines below switch off the information that is otherwise printed out by the perfect matching system
            options.verbose = false;
            pmX->options = options;
            //now we just load the graph information into the system
            for (index=0; index<numEdgesX; index++){
                pmX->AddEdge(tableX[3*index], tableX[3*index+1], tableX[3*index+2]);
            }

            pmX->Solve(); //this instruction tell the perfect matching system to find the solution
            for (index=0; index<numNodesX; index++){
                int jndex = pmX->GetMatch(index); //this line just gets the partner for node i from the system
                if (jndex>index){
//                    printf("di=%d,jn=%d\n",index,jndex);
                    correctionXnewTwo(index,jndex,globalStore,path,n,labelX);
//                    correctionXnewTwoSimple(index,jndex,n, labelTableX,globalStore,numQ);

                }
            }

            for (p=0;p<4*numQ;p++)
                globalStore[p]=globalStore[p]%2;
            
            
            rightZ1=sortLogicalZTwo1(globalStore,numQ);
            rightZ2=sortLogicalZTwo2(globalStore,numQ);
            
            free(tableX);
            free(path);
            delete pmX;
        }
        else{
            rightZ1=0;
            rightZ2=0;
        }
        
        if (labelZ>0){
            
            double *tableZ;
            tableZ=(double *)calloc(labelZ*labelZ*3,sizeof(double));
            if (tableZ == NULL){
                printf("space failed to be allocated！7\n"); exit(1);}
            
            int numNodesZ=0;
            
            int numEdgesZ=creatTableSimpleTwo(labelTableZ, labelZ, tableZ, numQ,1,weightZ);
            numNodesZ=labelZ;
            int workEdgeZ=numEdgesZ*2;
            
            //----------------------------------------------------------
            
            struct PerfectMatching::Options options;
            PerfectMatching *pmZ;
            pmZ = new PerfectMatching(numNodesZ, workEdgeZ);
            //the lines below switch off the information that is otherwise printed out by the perfect matching system
            options.verbose = false;
            pmZ->options = options;
            //now we just load the graph information into the system
            for (index=0; index<numEdgesZ; index++){
                pmZ->AddEdge(tableZ[3*index], tableZ[3*index+1], tableZ[3*index+2]);
            }
            pmZ->Solve(); //this instruction tell the perfect matching system to find the solution
            for (index=0; index<numNodesZ; index++){
                int jndex = pmZ->GetMatch(index); //this line just gets the partner for node i from the system
                if (jndex>index)
                    correctionZnewTwo(index,jndex,n,labelTableZ,globalStore,numQ);
            }
            
            for (p=0;p<4*numQ;p+=2)
                globalStore[p]=globalStore[p]%2;
            
            rightX1=sortLogicalXTwo1(globalStore,numQ);
            rightX2=sortLogicalXTwo2(globalStore,numQ);
            
            delete pmZ;
            free(tableZ);
        }
        else{
            rightX1=0;
            rightX2=0;
        }
    
        if ((rightX1+rightZ1)==0)
            right1++;
        if ((rightX2+rightZ2)==0)
            right2++;
        
//        if (rightX1==1)
//            right1++;
//        if (rightX1==1)
//            right1++;
        
//        if (rightX1==0 && rightZ1==1)
//            right1++;
//        if (rightX2==0 && rightZ2==1)
//            right1++;
//        if (rightX1==1 && rightZ1==0)
//            right2++;
//        if (rightX1==1 && rightZ1==0)
//            right2++;
        
        free(matrixX);
        
//        for (p=0;p<n;p++){
//            for (q=0;q<n;q++)
//                printf("%d ",(globalStore[4*(p*n+q)]+globalStore[4*(p*n+q)+2])%2);
//            printf("\n");
//        }
//        printf("\n");
    }
    
    double ratio=(double)(right1+right2)/(double)numCycle/2;
    printf("%f\n",ratio);

    double duration=(double)(endTime-startTime)/CLOCKS_PER_SEC;
    printf("%f\n",duration/numCycle);
    free(globalStore);
    free(copyX);
    free(copyZ);
    
    return 0;
}
