#ifndef _SIMPLECODENEW_H_
#define _SIMPLECODENEW_H_

# include "allFunctions.h"
# include "Dijv2.h"

void checkNew(int q, int type, int *globalStore, int numQ, int *label, int *labelTable, int *matrix, int *matrixNext)
{
    int qa;
    int n = sqrt(numQ);
    int i;
    int actual=0;

    if (type==0){
        qa = q-1;
        measurePerfectX(globalStore,q,numQ);
        actual = *(globalStore+2*qa+1);
    }
    if (type==1){
        qa = q+1;
        measurePerfectZ(globalStore,q,numQ);
        actual = *(globalStore+2*qa);
    }
    if (actual==1){
        *(labelTable+*label)=qa;
        *label+=1;
    }
    matrix[qa]=1;
    matrixNext[qa]=0;
}


void checkNewTwo(int q, int type, int *globalStore, int numQ, int *label, int *labelTable, int *matrix, int flag)
{
    int qa;
    int n = sqrt(numQ);
    int i;
    int actual=0;
    
    if (type==0){
        qa = q-1;
        measurePerfectX(globalStore,q,numQ);
        actual = *(globalStore+2*qa+1);
    }
    if (type==1){
        qa = q+1;
        measurePerfectZ(globalStore,q,numQ);
        actual = *(globalStore+2*qa);
    }
    if (actual==1){
        *(labelTable+*label)=qa;
        *label+=1;
        
    }
    if (*label>4000) {printf("label beyond array size\n"); exit(1);}
    if (flag==1)
        matrix[qa]=2;
    if (flag==0)
        matrix[qa]=0;
}

int creatTableSimple (int *labelTable, int label, int *table, int numQ, int type)
{
    int n=sqrt(numQ);
    int i=0; int j=0;
    int q1x; int q1y; int q2x; int q2y; int qa1; int qa2;
    
    int distanceY; int distanceX;
    int numError=0;
    
    for (i=0;i<label;i++){
        for(j=i+1;j<label;j++){
            
            qa1=labelTable[i];qa2=labelTable[j];
            
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

            double distance=distanceX+distanceY;
            *(table+3*numError+2)=distance;
            *(table+3*numError)=i;
            *(table+3*numError+1)=j;
            numError++;
            //            printf("%d,%d,%f\n",i,j,distance);
        }
    }
    return numError;
}

int creatTableNewFourAdapted (int *labelTable, int label, double *table, int numQ, int *matrixModel, double weight, double threshold, int *path)
{
    
    int n = sqrt(numQ);
    int size = n;
    int Maxn = numQ;
    int index; int jndex;
    int numError=0;
    
    Gragh gragh = Creategragh(size, matrixModel, weight, 1, 1000);
    
    double *dist;
    dist= (double *)calloc(Maxn/4,sizeof(double));
    
    int *prev;
    prev = (int *)calloc(Maxn/4, sizeof(int));
    
    for (index=0;index<label;index++){

        int temp=0;
        for (temp=0;temp<Maxn/4;temp++){
            dist[temp]=0;
            prev[temp]=0;
        }
        int name=(labelTable[index]/n-1)*n/4+labelTable[index]%n/2;
        
        dijkstra(gragh, name, prev, dist, Maxn/4,weight,threshold);
        
        for (jndex=index+1;jndex<label;jndex++){
            
            int templabel=labelTable[jndex];
            
            double distance=dist[(templabel/n-1)*n/4+templabel%n/2];
            
            *(path+numError*30)=labelTable[index];//the first two store information about which pair
            *(path+numError*30+1)=labelTable[jndex];
            
            int previous=prev[(labelTable[jndex]/n-1)*n/4+labelTable[jndex]%n/2];
            
            int temp2=0;
            while(previous!=name){
                
                int temp3=((previous%(n/2))*2+n*((previous/(n/2))*2+1))%(n*n);//change back the label
                if (temp3!=labelTable[index] && temp3!=labelTable[jndex] && temp3!=*(path+numError*30+1+temp2)){
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
            
            *(table+3*numError+2)=distance;
            *(table+3*numError)=index;
            *(table+3*numError+1)=jndex;
            numError++;
        }
        
    }
    
    return numError;
}

void correctionXnew(int i, int j, int n, int *labelTableX, int *globalStore, int numQ){
    
    int p; int qa1=labelTableX[i]; int qa2=labelTableX[j];
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
                globalStore[2*axis+1]++;
            }
        }
        else{
            for (p=0;p<(jrow-irow);p++)
                globalStore[2*(irow*2*n+n+2*icoll+n*(2*p+1))+1]++;
        }
        
        if (icoll>jcoll){
            if ((icoll-jcoll)> n/4){
                for (p=0;p<(jcoll+n/2-icoll);p++){
                    int axis=jrow*2*n+n+2*jcoll-2*p-1;
                    if (axis<(jrow*2*n+n)) axis+=n;
                    globalStore[2*axis+1]++;
                }
            }
            else{
                for (p=0;p<(icoll-jcoll);p++)
                    globalStore[2*(jrow*2*n+n+2*jcoll+2*p+1)+1]++;
            }
        }
        if (icoll<jcoll){
            if ((jcoll-icoll)> n/4){
                for (p=0;p<(icoll+n/2-jcoll);p++){
                    int axis=jrow*2*n+n+2*jcoll+2*p+1;
                    if (axis>(jrow*2*n+2*n-1)) axis-=n;
                    globalStore[2*axis+1]++;
                }
            }
            else{
                for (p=0;p<(jcoll-icoll);p++)
                    globalStore[2*(jrow*2*n+n+2*jcoll-2*p-1)+1]++;
            }
        }
        
    }
}

void correctionXnewTwo(int i, int j, int *globalStore, int *path, int n, int label)
{
    
    int numError=(2*label-i-1)*i/2+(j-i-1);
   
    if (path[30*numError]!=path[30*numError+1]){
        
        if (path[30*numError+2]<0){//to connect the only two ancilla
            
            int avg=(path[30*numError+1]+path[30*numError])/2;
            if (abs(avg-path[30*numError+1])==n || abs(avg-path[30*numError+1])==1){
                *(globalStore+2*avg+1)+=1;
                //                printf("%d\n",avg);
            }
            
            else if (abs(path[30*numError+1]-path[30*numError])<n){
                *(globalStore+2*(avg+n/2)+1)+=1;
                //                printf("%d\n",avg+n/2);
            }
            
        }
        else {
            
            int temp=2;
            while (path[30*numError+temp]>0){//to connect the last ancilla to the second first
                //                printf("%d  ",path[30*numError+temp]);
                int avg=(path[30*numError+temp-1]+path[30*numError+temp])/2;
                if (abs(avg-path[30*numError+temp])==n || abs(avg-path[30*numError+temp])==1){
                    *(globalStore+2*avg+1)+=1;
                    //                    printf("%d ",avg);
                }
                
                else if (abs(path[30*numError+temp]-path[30*numError+temp-1])<n){
                    *(globalStore+2*(avg+n/2)+1)+=1;
                    //                    printf("%d ",avg+n/2);
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
                *(globalStore+2*avg+1)+=1;
                //                printf("%d\n",avg);
            }
            
            else if (abs(path[30*numError]-path[30*numError+temp-1])<n){
                *(globalStore+2*(avg+n/2)+1)+=1;
                //                printf("%d\n",avg+n/2);
            }
            
            else if (abs(path[30*numError]-path[30*numError+temp-1])>n){
                *(globalStore+2*(avg-n/2*n)+1)+=1;
                //                printf("%d\n",avg-n/2*n);
            }
            
        }
    }
    
}

void correctionZnew(int i, int j, int n, int *labelTableZ, int *globalStore, int numQ){
    
    int p; int qa1=labelTableZ[i]; int qa2=labelTableZ[j];
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
                globalStore[2*axis]++;
            }
        }
        else{
            for (p=0;p<(jrow-irow);p++)
                globalStore[2*(irow*2*n+1+2*icoll+n*(2*p+1))]++;
        }
        
        if (icoll>jcoll){
            if ((icoll-jcoll)> n/4){
                for (p=0;p<(jcoll+n/2-icoll);p++){
                    int axis=jrow*2*n+1+2*jcoll-2*p-1;
                    if (axis<(jrow*2*n)) axis+=n;
                    globalStore[2*axis]++;
                }
            }
            else{
                for (p=0;p<(icoll-jcoll);p++)
                    globalStore[2*(jrow*2*n+1+2*jcoll+2*p+1)]++;
            }
        }
        if (icoll<jcoll){
            if ((jcoll-icoll)> n/4){
                for (p=0;p<(icoll+n/2-jcoll);p++){
                    int axis=jrow*2*n+1+2*jcoll+2*p+1;
                    if (axis>(jrow*2*n+n-2)) axis-=n;
                    globalStore[2*axis]++;
                }
            }
            else{
                for (p=0;p<(jcoll-icoll);p++)
                    globalStore[2*(jrow*2*n+1+2*jcoll-2*p-1)]++;
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
        temp2+=globalStore[2*(n*(n-1)+p+1)];
    
    flag=temp2%2;
    
    return flag;
    
}

int sortLogicalXTwo1(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp1=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp1+=globalStore[2*(n*p)];
    
    flag=temp1%2;
    
    return flag;
    
}

int sortLogicalZTwo1(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp1=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp1+=globalStore[2*p+1];
    
    flag=temp1%2;
    return flag;
    
}

int sortLogicalZTwo2(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp2=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp2+=globalStore[2*(2*n-1+n*p)+1];
    
    flag=temp2%2;
    
    return flag;
}

#endif
