#ifndef _SIMPLECODENEW_H_
#define _SIMPLECODENEW_H_

# include "allFunctions.h"
# include "Dijv2.h"

void measureNoiseTwoGate(int q, int *Q, double gateError, double biasRatio)
{
    double ran1; double ran2;
    ran1 =((double)rand()/(double)RAND_MAX);
    ran2 =((double)rand()/(double)RAND_MAX);
    double ran3 =((double)rand()/(double)RAND_MAX);
    
    if (ran1<gateError){
        double t1;
        if (biasRatio>0)
            t1=biasRatio/(1+biasRatio);
        else
            t1=0;
        
        if (ran2<t1){
            if (ran3<=oneoverthree)
                *(Q+2*q)+=1;
            else if (ran3>oneoverthree && ran3<= 2*oneoverthree){
                *(Q+2*q)+=1;*(Q+2*q+1)+=1;
            }
            else if (ran3>2*oneoverthree)
                *(Q+2*q+1)+=1;
        }
        else if (ran2>=t1)
        {*(Q+2*q+1)+=1;}
        
        *(Q+2*q+1)=*(Q+2*q+1)%2; *(Q+2*q)=*(Q+2*q)%2;
    }
}

void hadamardNoiseTwoGate(int q, int *Q, double gateError, double biasRatio)
{
    sortHadamardError(q,Q);
    measureNoiseTwoGate(q, Q, gateError, biasRatio);
}

void cnotTwoGate(int q1, int q2, int *Q, double gateError, double biasRatio)
{
    
    sortCnotError(q1,q2,Q);
    double ran;
    ran =((double)rand()/(double)RAND_MAX);
    if (ran < gateError)
    {
        double ran2;
        ran2 =((double)rand()/(double)RAND_MAX);
        double ran3;
        ran3 =((double)rand()/(double)RAND_MAX);
        double t1;
        if (biasRatio>0)
            t1=biasRatio/(1+biasRatio);
        else
            t1=0;
        
        if (ran2<t1){
            
            int myRnd=1 + rand()%15;
            int firstError=myRnd/4;
            int secndError=myRnd%4;
            
            if (firstError==1) {*(Q+2*q1)+=1;}
            if (firstError==2) {*(Q+2*q1)+=1; *(Q+2*q1+1)+=1;}
            if (firstError==3) {*(Q+2*q1+1)+=1;}
            if (secndError==1) {*(Q+2*q2)+=1;}
            if (secndError==2) {*(Q+2*q2)+=1; *(Q+2*q2+1)+=1;}
            if (secndError==3) {*(Q+2*q2+1)+=1;}
            
        }
        if (ran2>=t1){
            
            if (ran3<=oneoverthree)
                *(Q+2*q1+1)+=1;
            else if (ran3>oneoverthree && ran3<= 2*oneoverthree){
                *(Q+2*q1+1)+=1;*(Q+2*q2+1)+=1;
            }
            else if (ran3>2*oneoverthree)
                *(Q+2*q2+1)+=1;
        }
        
        *(Q+2*q1)=*(Q+2*q1)%2; *(Q+2*q1+1)=*(Q+2*q1+1)%2; *(Q+2*q2)=*(Q+2*q2)%2; *(Q+2*q2+1)=*(Q+2*q2+1)%2;
        
    }
}

void errorDetection(int q, int numQ, int *globalStore, int *matrix, int count, double errorRateDetection, double errorRateGate, double biasRatio, double flipRatio)
{
    int patternD[6]={0};
    measureNoiseTwoGate(0,patternD,errorRateDetection,biasRatio);
    cnotTwoGate (0,1, patternD, errorRateDetection,biasRatio);
    cnotTwoGate (0,2, patternD, errorRateDetection,biasRatio);
    measureNoiseTwoGate(0,patternD, errorRateDetection,biasRatio);
    
    if ((globalStore[4*q+1]+globalStore[4*q+3]+patternD[1])%2==1){
        double ran;
        ran =((double)rand()/(double)RAND_MAX);
        if (ran<flipRatio){
            *(globalStore+4*q+1)+=1;
            *(globalStore+4*q+1)%=2;
        }
        else{
            *(globalStore+4*q+3)+=1;
            *(globalStore+4*q+3)%=2;
        }
        
        matrix[numQ*count+q]=1;
    }
    
    int index;
    for (index=0;index<4;index++)
        *(globalStore+4*q+index)+=patternD[2+index];
    
}

void checkNewTwoGate(int q, int type, int *globalStore, int *copy,int numQ, int *label, int *labelTable, int count, int addError, int *matrix, double errorRateGate, double errorRateMeasure, double errorRateDetection, double biasRatio,double flipRatio)
{
   
    int n = sqrt(numQ);
    int i;
    int actual=0;
    int ql[5]={0}; int index; int jndex;
    
    if (type==0){
        
        ql[3] = q;
        ql[4] = q-1;
        ql[0] = q-n-1;
        ql[1] = q+n-1;
        if (ql[1]>=numQ) ql[1]=ql[1]-numQ;
        if (q%(2*n)==(n+1))
            ql[2]=q+n-2;
        else
            ql[2]=q-2;
        
        //error detection
        
        for (index=0;index<4;index++)
            errorDetection(ql[index], numQ, globalStore, matrix, count,errorRateDetection,errorRateGate,biasRatio,flipRatio);
        
        //initialise ancillar
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
        //perfect check
        int sum=0;
        for (index=0;index<4;index++){
            if (globalStore[4*ql[index]+1]%2==1)
//                printf("ql=%d,ql4=%d\n",ql[index],ql[4]);
            sum+=globalStore[4*ql[index]+1];
        }
        
        if (sum%2==1)
            *(globalStore+4*ql[4]+1)=1;
        else
            *(globalStore+4*ql[4]+1)=0;
//        printf("gl=%d\n",globalStore[4*ql[4]+1]);
        //add error
        if (addError==1){
            
            int patternX[10]={0};
            measureNoiseTwoGate(4,patternX, errorRateMeasure,biasRatio);
//            hadamardNoiseTwoGate(4,patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,0, patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,1, patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,2, patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,3, patternX, errorRateGate,biasRatio);
//            hadamardNoiseTwoGate(4,patternX, errorRateGate,biasRatio);
            measureNoiseTwoGate(4,patternX, errorRateMeasure,biasRatio);

            for (index=0;index<4;index++){
                for (jndex=0;jndex<2;jndex++)
                    *(globalStore+4*ql[index]+jndex)+=patternX[2*index+jndex];
            }
            
            if (patternX[9]==1)
                *(globalStore+4*ql[4]+1)+=1;
        }
        
        actual = *(globalStore+4*ql[4]+1)%2;
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
    }
    if (type==1){
        
        ql[2] = q;
        ql[4] = q+1;
        ql[0] = q-n+1;
        if (ql[0]<0) ql[0]+=numQ;
        ql[1] = ql[2]+n+1;
        if (ql[2]%(2*n)==(n-2))
            ql[3]=ql[2]-n+2;
        else
            ql[3]=ql[2]+2;
        
        //initialise ancillar
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
        //perfect check
        int sum=0;
        for (index=0;index<4;index++)
            sum+=globalStore[4*ql[index]]+globalStore[4*ql[index]+2];
        if (sum%2==1)
            *(globalStore+4*ql[4])=1;
        else
            *(globalStore+4*ql[4])=0;
        
        //add error
        if (addError==1){
        
            int patternZ[20]={0};
            measureNoiseTwoGate(8,patternZ,errorRateMeasure,biasRatio);
            measureNoiseTwoGate(9,patternZ,errorRateMeasure,biasRatio);
            cnotTwoGate (0,8, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (1,9, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (2,8, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (3,9, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (4,8, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (5,9, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (6,8, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (7,9, patternZ, errorRateGate,biasRatio);
            measureNoiseTwoGate(8,patternZ,errorRateMeasure,biasRatio);
            measureNoiseTwoGate(9,patternZ,errorRateMeasure,biasRatio);
            
            for (index=0;index<4;index++){
                for (jndex=0;jndex<4;jndex++)
                    *(globalStore+4*ql[index]+jndex)+=patternZ[4*index+jndex];
            }
            
            if ((patternZ[16]+patternZ[18])==1)
                *(globalStore+4*ql[4])+=1;
        }
        actual = *(globalStore+4*ql[4])%2;
        
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
    }
    
    if (*(copy+ql[4])!=actual){
        
//        printf("%d\n",ql[4]);
        *(labelTable+*label*2)=count;
        *(labelTable+*label*2+1)=ql[4];
        *label+=1;
        if (*label>4000) {printf("label beyond array size\n"); exit(1);}
    }
    *(copy+ql[4])=actual;

}

void checkNewTwoGateSimple(int q, int type, int *globalStore, int *copy,int numQ, int *label, int *labelTable, int count, int addError, double errorRateGate, double errorRateMeasure, double biasRatio)
{
    
    int n = sqrt(numQ);
    int i;
    int actual=0;
    int ql[5]={0}; int index; int jndex;
    
    if (type==0){
        
        ql[3] = q;
        ql[4] = q-1;
        ql[0] = q-n-1;
        ql[1] = q+n-1;
        if (ql[1]>=numQ) ql[1]=ql[1]-numQ;
        if (q%(2*n)==(n+1))
            ql[2]=q+n-2;
        else
            ql[2]=q-2;
        
        //initialise ancillar
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
        //perfect check
        int sum=0;
        for (index=0;index<4;index++)
            sum+=globalStore[4*ql[index]+1];
        
        if (sum%2==1)
            *(globalStore+4*ql[4]+1)=1;
        else
            *(globalStore+4*ql[4]+1)=0;
        
        //add error
        if (addError==1){
            
            int patternX[10]={0};
            
            measureNoiseTwoGate(4,patternX, errorRateMeasure,biasRatio);
            hadamardNoiseTwoGate(4,patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,0, patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,1, patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,2, patternX, errorRateGate,biasRatio);
            cnotTwoGate (4,3, patternX, errorRateGate,biasRatio);
            hadamardNoiseTwoGate(4,patternX, errorRateGate,biasRatio);
            measureNoiseTwoGate(4,patternX, errorRateMeasure,biasRatio);

            for (index=0;index<4;index++){
                for (jndex=0;jndex<2;jndex++)
                    *(globalStore+4*ql[index]+jndex)+=patternX[2*index+jndex];
            }
      
            if (patternX[8]==1)
                *(globalStore+4*ql[4]+1)+=1;
        }
        
        actual = *(globalStore+4*ql[4]+1)%2;
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
    }
    if (type==1){
        
        ql[2] = q;
        ql[4] = q+1;
        ql[0] = q-n+1;
        if (ql[0]<0) ql[0]+=numQ;
        ql[1] = ql[2]+n+1;
        if (ql[2]%(2*n)==(n-2))
            ql[3]=ql[2]-n+2;
        else
            ql[3]=ql[2]+2;
        
        //initialise ancillar
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
        //perfect check
        int sum=0;
        for (index=0;index<4;index++)
            sum+=globalStore[4*ql[index]];
        if (sum%2==1)
            *(globalStore+4*ql[4])=1;
        else
            *(globalStore+4*ql[4])=0;
        
        //add error
        if (addError==1){
            
            int patternZ[10]={0};
            measureNoiseTwoGate(4,patternZ,errorRateMeasure,biasRatio);
            cnotTwoGate (0,4, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (1,4, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (2,4, patternZ, errorRateGate,biasRatio);
            cnotTwoGate (3,4, patternZ, errorRateGate,biasRatio);
            measureNoiseTwoGate(4,patternZ,errorRateMeasure,biasRatio);
            
            for (index=0;index<4;index++){
                for (jndex=0;jndex<2;jndex++)
                    *(globalStore+4*ql[index]+jndex)+=patternZ[2*index+jndex];
            }
            
            if (patternZ[8]==1)
                *(globalStore+4*ql[4])+=1;
            }
        
        actual = *(globalStore+4*ql[4])%2;
        
        for (index=0;index<4;index++)
            *(globalStore+4*ql[4]+index)=0;
    }
    
    if (*(copy+ql[4])!=actual){
        
        //        printf("%d\n",ql[4]);
        *(labelTable+*label*2)=count;
        *(labelTable+*label*2+1)=ql[4];
        *label+=1;
        if (*label>4000) {printf("label beyond array size\n"); exit(1);}
    }
    *(copy+ql[4])=actual;
    
}

int creatTableSimpleTwo (int *labelTable, int label, double *table, int numQ, int type,double weight)
{
    int n=sqrt(numQ);
    int i=0; int j=0;
    int q1x; int q1y; int q2x; int q2y; int qa1; int qa2;
    
    int distanceY; int distanceX;
    int numError=0;
    int distanceTime=0;
    
    for (i=0;i<2*label;i+=2){
        for(j=i+2;j<2*label;j+=2){
            
            distanceTime=abs(labelTable[i]-labelTable[j]);
            qa1=labelTable[i+1];qa2=labelTable[j+1];
            
            if (type==0){
                q1x=(qa1-n)/(n*2); q2x=(qa2-n)/(n*2);
                q1y=qa1%n/2; q2y=qa2%n/2;
            }
            if (type==1){
                q1x=(qa1-1)/(n*2); q2x=(qa2-1)/(n*2);
                q1y=(qa1-1)%n/2; q2y=(qa2-1)%n/2;
            }
            
            if ((q2x-q1x)>=0){
                if ((q2x-q1x)> (n/4))
                    distanceX = q1x+n/2-q2x;
                else
                    distanceX = q2x-q1x;
            }
            else{
                if ((q1x-q2x)> (n/4))
                    distanceX = q2x+n/2-q1x;
                else
                    distanceX = q1x-q2x;
            }
            
            if ((q2y-q1y)>=0){
                if ((q2y-q1y)> (n/4))
                    distanceY = q1y+n/2-q2y;
                else
                    distanceY = q2y-q1y;
            }
            else{
                if ((q1y-q2y)> (n/4))
                    distanceY = q2y+n/2-q1y;
                else
                    distanceY= q1y-q2y;
            }

            double distance=distanceX+distanceY+weight*distanceTime;
            *(table+3*numError+2)=distance;
            *(table+3*numError)=i/2;
            *(table+3*numError+1)=j/2;
            numError++;
            //            printf("%d,%d,%f\n",i,j,distance);
        }
    }
    return numError;
}

int creatTableNewFourAdapted (int *labelTable, int label, double *table, int numQ, int *matrixModel, double weight, double weightTime, int timeD, int *path, double threshold)
{
    
    int i=0; int j=0; int k=0;
    int n=sqrt(numQ);
    int size = n;
    int Maxn = numQ;
    int index; int jndex;
    int numError=0;
    
//    int temp;
//    for (temp=0;temp<=2;temp++){
//    
//        for (index=0;index<n;index++){
//            for (jndex=0;jndex<n;jndex++)
//                printf("%d ",matrixModel[temp*numQ+index*n+jndex]);
//            printf("\n");
//        }
//        printf("\n");
//    }
//    printf("\n");
//    printf("\n");
    
    
    Gragh gragh = Creategragh(size, matrixModel, weight, timeD, weightTime);
//    printf("3\n");
    
    double *dist;
    dist= (double *)calloc(Maxn*timeD/4,sizeof(double));
    
    int *prev;
    prev = (int *)calloc(Maxn*timeD/4, sizeof(int));
    
    for (index=0;index<2*label;index+=2){
        int temp=0;
        for (temp=0;temp<Maxn*timeD/4;temp++){
            dist[temp]=0;
            prev[temp]=0;
        }
        int name=labelTable[index]*Maxn/4+(labelTable[index+1]/n-1)*n/4+labelTable[index+1]%n/2;
//        printf("%d,name=%d\n",labelTable[index+1],name);
        dijkstra(gragh, name, prev, dist, Maxn*timeD/4,weight,threshold);
//        printf("finish-2\n");
        
//        for (temp=0;temp<Maxn*timeD/4;temp++){
//            printf("%d %d   ",temp,prev[temp]);
//        }
//        printf("\n");
        
        for (jndex=index+2;jndex<2*label;jndex+=2){
            
            int templabel=labelTable[jndex+1];
//            printf("%d,%d\n",templabel,labelTable[jndex]);
            
            double distance=dist[labelTable[jndex]*Maxn/4+(templabel/n-1)*n/4+templabel%n/2];
            
//            if (distance<100)
//                printf("%f\n",distance);
            
            *(path+numError*30)=labelTable[index+1];//the first two store information about which pair
            *(path+numError*30+1)=labelTable[jndex+1];
//            printf("%d %d\n",labelTable[index+1],labelTable[jndex+1]);
            int previous=prev[labelTable[jndex]*Maxn/4+(labelTable[jndex+1]/n-1)*n/4+labelTable[jndex+1]%n/2];
//            *(path+numError*20+19)=distance;
            
            int temp2=0;
            while(previous!=name){
                
                int temp3=((previous%(n/2))*2+n*((previous/(n/2))*2+1))%(n*n);//change back the label
                if (temp3!=labelTable[index+1] && temp3!=labelTable[jndex+1] && temp3!=*(path+numError*30+1+temp2)){
                    *(path+numError*30+2+temp2)=temp3;
                    temp2++;
//                    printf("%d ",temp3);
                }
                previous=prev[previous];
                
                if (temp2>28){
                    printf("exceed!");
                    exit(1);
                }
            }
//            printf("\n");
            *(table+3*numError+2)=distance;
            *(table+3*numError)=index/2;
            *(table+3*numError+1)=jndex/2;
            numError++;
        }
        
    }
//    printf("4\n");
    for (i = 0; i<Maxn*timeD/4; i++)
        free(gragh->G[i]);
    free(gragh);
    free(prev);
    free(dist);
    
    return numError;
}

void correctionZnewTwo(int i, int j, int n, int *labelTableZ, int *globalStore, int numQ){
    
    int p; int qa1=labelTableZ[2*i+1]; int qa2=labelTableZ[2*j+1];
    if (qa1<0 || qa2<0){
        printf("ancilla wrong\n");
        exit(1);
    }
    int irow=(qa1-1)/(n*2);
    int icoll=(qa1-1)%n/2;
    int jrow=(qa2-1)/(n*2);
    int jcoll=(qa2-1)%n/2;
    
    if (irow!=jrow || icoll!=jcoll){
        
        
        int rowL = findLarger(irow,jrow);
        
        if (rowL!=jrow){
            int temprow=irow;
            int tempcoll=icoll;
            irow=jrow; icoll=jcoll;
            jrow=temprow; jcoll=tempcoll;
        }
        //printf("irow=%d,icoll=%d,jrow=%d,jcoll=%d\n",irow,icoll,jrow,jcoll);
        if ((jrow-irow)> n/4){
            for (p=0;p<(irow+n/2-jrow);p++){
                int axis=irow*2*n+1+2*icoll-n*(2*p+1);
                if (axis<0) axis+=numQ;
                globalStore[4*axis]++;
            }
        }
        else{
            for (p=0;p<(jrow-irow);p++)
                globalStore[4*(irow*2*n+1+2*icoll+n*(2*p+1))]++;
        }
        
        if (icoll>jcoll){
            if ((icoll-jcoll)> n/4){
                for (p=0;p<(jcoll+n/2-icoll);p++){
                    int axis=jrow*2*n+1+2*jcoll-2*p-1;
                    if (axis<(jrow*2*n)) axis+=n;
                    globalStore[4*axis]++;
                }
            }
            else{
                for (p=0;p<(icoll-jcoll);p++)
                    globalStore[4*(jrow*2*n+1+2*jcoll+2*p+1)]++;
            }
        }
        if (icoll<jcoll){
            if ((jcoll-icoll)> n/4){
                for (p=0;p<(icoll+n/2-jcoll);p++){
                    int axis=jrow*2*n+1+2*jcoll+2*p+1;
                    if (axis>(jrow*2*n+n-2)) axis-=n;
                    globalStore[4*axis]++;
                }
            }
            else{
                for (p=0;p<(jcoll-icoll);p++)
                    globalStore[4*(jrow*2*n+1+2*jcoll-2*p-1)]++;
            }
        }
        
    }
    
}

void correctionXnewTwoSimple(int i, int j, int n, int *labelTableX, int *globalStore, int numQ){
    
    int p; int qa1=labelTableX[2*i+1]; int qa2=labelTableX[2*j+1];
//    printf("qa1=%d,qa2=%d\n",qa1,qa2);
    if (qa1<0 || qa2<0){
        printf("ancilla wrong\n");
        exit(1);
    }
    int irow=(qa1-n)/(n*2);
    int icoll=qa1%n/2;
    int jrow=(qa2-n)/(n*2);
    int jcoll=qa2%n/2;
    
    if (irow!=jrow || icoll!=jcoll){
        
        
        int rowL = findLarger(irow,jrow);
        
        if (rowL!=jrow){
            int temprow=irow;
            int tempcoll=icoll;
            irow=jrow; icoll=jcoll;
            jrow=temprow; jcoll=tempcoll;
        }
        
        if ((jrow-irow)> n/4){
            for (p=0;p<(irow+n/2-jrow);p++){
                int axis=irow*2*n+n+2*icoll-n*(2*p+1);
                if (axis<0) axis+=numQ;
                globalStore[4*axis+1]++;
//                printf("%d ",axis);
            }
        }
        else{
            for (p=0;p<(jrow-irow);p++){
                globalStore[4*(irow*2*n+n+2*icoll+n*(2*p+1))+1]++;
//                printf("%d ",irow*2*n+n+2*icoll+n*(2*p+1));
            }
        }
        
        if (icoll>jcoll){
            if ((icoll-jcoll)> n/4){
                for (p=0;p<(jcoll+n/2-icoll);p++){
                    int axis=jrow*2*n+n+2*jcoll-2*p-1;
                    if (axis<(jrow*2*n+n))
                        axis+=n;
                    globalStore[4*axis+1]++;
//                    printf("%d ",axis);
                }
            }
            else{
                for (p=0;p<(icoll-jcoll);p++){
                    globalStore[4*(jrow*2*n+n+2*jcoll+2*p+1)+1]++;
//                    printf("%d ",jrow*2*n+n+2*jcoll+2*p+1);
                }
            }
        }
        if (icoll<jcoll){
            if ((jcoll-icoll)> n/4){
                for (p=0;p<(icoll+n/2-jcoll);p++){
                    int axis=jrow*2*n+n+2*jcoll+2*p+1;
                    if (axis>(jrow*2*n+2*n-1)) axis-=n;
                    globalStore[4*axis+1]++;
//                    printf("%d ",axis);
                }
            }
            else{
                for (p=0;p<(jcoll-icoll);p++){
                    globalStore[4*(jrow*2*n+n+2*jcoll-2*p-1)+1]++;
//                    printf("%d ",jrow*2*n+n+2*jcoll-2*p-1);
                }
            }
        }
        
    }
//    printf("\n");

}

void correctionXnewTwo(int i, int j, int *globalStore, int *path, int n, int label)
{
    
    int numError=(2*label-i-1)*i/2+(j-i-1);
//    printf("%d\n",path[30*numError+19]);
    
//    printf("i=%d,j=%d\n",path[30*numError],path[30*numError+1]);
    if (path[30*numError]!=path[30*numError+1]){
        
        if (path[30*numError+2]<0){//to connect the only two ancilla
            
            int avg=(path[30*numError+1]+path[30*numError])/2;
            if (abs(avg-path[30*numError+1])==n || abs(avg-path[30*numError+1])==1){
                *(globalStore+4*avg+1)+=1;
//                printf("%d\n",avg);
            }
            
            else if (abs(path[30*numError+1]-path[30*numError])<n){
                *(globalStore+4*(avg+n/2)+1)+=1;
//                printf("%d\n",avg+n/2);
            }
            
            else if (abs(path[30*numError+1]-path[30*numError])>n){
                *(globalStore+4*(avg-n/2*n)+1)+=1;
//                printf("%d\n",avg-n/2*n);
            }
            
        }
        else {
        
            int temp=2;
            while (path[30*numError+temp]>0){//to connect the last ancilla to the second first
//                printf("%d  ",path[30*numError+temp]);
                int avg=(path[30*numError+temp-1]+path[30*numError+temp])/2;
                if (abs(avg-path[30*numError+temp])==n || abs(avg-path[30*numError+temp])==1){
                    *(globalStore+4*avg+1)+=1;
//                    printf("%d ",avg);
                }
                
                else if (abs(path[30*numError+temp]-path[30*numError+temp-1])<n){
                    *(globalStore+4*(avg+n/2)+1)+=1;
//                    printf("%d ",avg+n/2);
                }
                
                else if (abs(path[30*numError+temp]-path[30*numError+temp-1])>n){
                    *(globalStore+4*(avg-n/2*n)+1)+=1;
//                    printf("%d ",avg-n/2*n);
                }
                
                
                temp++;
        
                if (temp>=20){
                    printf("path wrong\n");
                    exit(1);
                }
            }
            //to connect the second to the first
            int avg=(path[30*numError+temp-1]+path[30*numError])/2;
            if (abs(avg-path[30*numError])==n || abs(avg-path[30*numError])==1){
                *(globalStore+4*avg+1)+=1;
//                printf("%d\n",avg);
            }
            
            else if (abs(path[30*numError]-path[30*numError+temp-1])<n){
                *(globalStore+4*(avg+n/2)+1)+=1;
//                printf("%d\n",avg+n/2);
            }
            
            else if (abs(path[30*numError]-path[30*numError+temp-1])>n){
                *(globalStore+4*(avg-n/2*n)+1)+=1;
//                printf("%d\n",avg-n/2*n);
            }
            
        }
    }
    
}

int sortLogicalXTwo2(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp2=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp2+=globalStore[4*(n*(n-1)+p+1)]+globalStore[4*(n*(n-1)+p+1)+2];
    
    flag=temp2%2;
    
    return flag;
    
}

int sortLogicalXTwo1(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp1=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp1+=globalStore[4*(n*p)]+globalStore[4*(n*p)+2];
    
    flag=temp1%2;
    
    return flag;
    
}

int sortLogicalZTwo1(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp1=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp1+=globalStore[4*p+1];
    
    flag=temp1%2;
    return flag;
    
}

int sortLogicalZTwo2(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp2=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp2+=globalStore[4*(2*n-1+n*p)+1];
    
    flag=temp2%2;
    
    return flag;
}
#endif
