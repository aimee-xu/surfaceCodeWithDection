#include "simpleCodeNewv2.h"
#include "PerfectMatching.h"

//--------------------------------------------------------------
//---------------------- START OF main( )  ----------------------
//--------------------------------------------------------------
int main ( int argc, char *argv[] )
{
    time_t rSeed;
    rSeed = time (NULL);
    unsigned long labelNumber=0;
    
    int numQ=64;
    int numCycle=1000;
    
    double weight=10;
    double biasratio=0.1;
    double thresholdr=3;
    double errorRateGate=0.14;
    // begin paste -------------------------------
    if (argc<3){ //this message will print if there is not at least one command line parameter
        printf("\nINFORMATION ON COMMAND LINE PARAMETERS:\nA usage example is './thisProg -r=5 -i=200 -em=0.1 -eg=0.1 -ep=0.1 -t=200'\n");
        printf("  -i is the numcycle\n");
        printf("  -c is the number of iterations\n");
        printf("  -a is errorrate ancillar\n");
        printf("  -d is error rate data qubit\n");
        printf("  -u is number of qubit\n");
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
                if (thisFlag=='r'){
                    labelNumber=atoi(argv[i])*1000;
                    failedToParse=0;
                }
                if (thisFlag=='w'){
                    weight=strtof(argv[i], NULL);
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

   //main simulation
    
    int n= sqrt(numQ);
    int index; int jndex; int p; int q; int t;
    
    //space allocation
    int *globalStore;
    globalStore=(int *)calloc(numQ*2,sizeof(int));
    if (globalStore == NULL){
        printf("space failed to be allocated！2\n"); exit(1);}
    
    int rightX1=0; int rightZ1=0; int right1=0;
    int rightX2=0; int rightZ2=0; int right2=0;
    
    //-----------------------------------------------------
    long int cycle;
    for (cycle=0;cycle<numCycle;cycle++){
        
        int labelTableX[4000]={0};//store labels of syndrome changes
        int labelTableZ[4000]={0};
        int labelX=0; //the number of syndrome changes
        int labelZ=0;
        int errorListZ[4000]={0};
        
        int *matrixX;//store all possible errors
        matrixX=(int *)calloc(numQ,sizeof(int));
        if (matrixX == NULL){
            printf("space failed to be allocated！2\n"); exit(1);}
        
        for (p=0;p<2*numQ;p++){
            globalStore[p]=0;
        }
        
        //put errors to the data qubit
        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2;q++){
                int t1=2*n*p+2*q;
                int t2=n+1+2*n*p+2*q;
                
                int tem=0;
                int tempx=0; int tempz=0; int temp2x=0; int temp2z=0;
                for (tem=0;tem<2;tem++){
                    if (((double)rand()/(double)RAND_MAX)<errorRateGate*biasratio)
                        tempx++;
                    if (((double)rand()/(double)RAND_MAX)<errorRateGate*(1-biasratio))
                        tempz++;
                    if (((double)rand()/(double)RAND_MAX)<errorRateGate*biasratio)
                        temp2x++;
                    if (((double)rand()/(double)RAND_MAX)<errorRateGate*(1-biasratio))
                        temp2z++;
                }
                
                if (tempx==1)
                    globalStore[2*t1]=1;
                if (temp2x==1)
                    globalStore[2*t2]=1;
                if (tempz==2)
                    globalStore[2*t1+1]=1;
                if (temp2z==2)
                    globalStore[2*t2+1]=1;
                if (tempz==1){
                    matrixX[t1]=1;
                    if (((double)rand()/(double)RAND_MAX)<0.5)
                        globalStore[2*t1+1]=1;
                }
                if (temp2z==1){
                    matrixX[t2]=1;
                    if (((double)rand()/(double)RAND_MAX)<0.5)
                        globalStore[2*t2+1]=1;
                }
            }
            p++;
        }
        
        //begin stabiliser check
        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2;q++){
                int t1=2*n*p+2*q;
                int t2=n+1+2*n*p+2*q;
                checkNewTwo(t1,1,globalStore,numQ,&labelZ,labelTableZ,matrixX,0);
                checkNewTwo(t2,0,globalStore,numQ,&labelX,labelTableX,matrixX,1);
            }
            p++;
        }
        
        if (labelX>0){
            
            double *tableX;
            tableX=(double *)calloc(labelX*labelX*3,sizeof(double));
            if (tableX == NULL){
                printf("space failed to be allocated！2\n"); exit(1);}
            
            int *path;
            path=(int *)calloc(labelX*(labelX-1)*15,sizeof(int));
            if (path == NULL){
                printf("space failed to be allocated！6\n"); exit(1);}
            for (index=0;index<labelX*(labelX-1)*15;index++)
                path[index]=-1;
            
            int numNodesX=0;
    
            int numEdgesX=creatTableNewFourAdapted(labelTableX, labelX, tableX, numQ,matrixX,weight,thresholdr,path);
            
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
                if (jndex>index)
                    correctionXnewTwo(index,jndex,globalStore,path,n,label);
            }
            
            for (p=0;p<2*numQ;p++)
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
            
            int *tableZ;
            tableZ=(int *)calloc(labelZ*labelZ*3,sizeof(int));
            if (tableZ == NULL){
                printf("space failed to be allocated！2\n"); exit(1);}
            
            int numNodesZ=0;
            int numEdgesZ=creatTableSimple(labelTableZ, labelZ, tableZ, numQ,1);
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
                    correctionZnew(index,jndex,n,labelTableZ,globalStore,numQ);
            }
            
            for (p=0;p<2*numQ;p++)
                globalStore[p]=globalStore[p]%2;
                            sortError(globalStore,numQ);
                            for (p=0;p<n;p++){
                                for (q=0;q<n;q++)
                                    printf("%d ",globalStore[2*(p*n+q)]);
                                printf("\n");
                            }
                            printf("\n");
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
        
        free(matrixX);
    }
    
    double ratio=(double)(right1+right2)/(double)numCycle/2;
    printf("%f\n",ratio);
    free(globalStore);
    
    return 0;
}
