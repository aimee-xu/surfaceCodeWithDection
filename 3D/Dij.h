#ifndef _DIJ_H_
#define _DIJ_H_

typedef struct node * Gragh;

struct node
{
    int Nv;
    int Ne;
    double **G;
    
};

int getConnection(int ***matrix, double **G, int size, double w, int Tp, double tw)
{
    
    int count = 0; int m = size; int n = size;
    int i;
    int j;
    int t;
    int Maxn = m * n;
    
    for (t = 0; t < Tp; t++)
    {
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                
                int itmp = (i + 1) % m;
                int jtmp = (j + 1) % n;
                
                int ttmp = (t + 1) % Tp;
                if (matrix[t][i][j] == 2)
                {
                    G[t*Maxn + i * n + j][ttmp*Maxn + i * n + j] = tw;
                    G[ttmp*Maxn + i * n + j][t*Maxn + i * n + j] = tw;
                    
                }
                
                
                if ((matrix[t][i][j] == 1)||( matrix[t][i][j] == 2))
                {
                    
                    if ((matrix[t][itmp][j] == 1) || (matrix[t][itmp][j] == 2)) {
                        G[t*Maxn+i*n + j][t*Maxn + itmp*n + j] = 1;
                        G[t*Maxn + itmp*n + j][t*Maxn + i*n + j] = 1;
                        count++;
                    }
                    else {
                        G[t*Maxn + i*n + j][t*Maxn + itmp*n + j] = w;
                        G[t*Maxn + itmp*n + j][t*Maxn + i*n + j] = w;
                    }
                    
                    if ((matrix[t][i][jtmp] == 1) ||(matrix[t][i][jtmp] == 2)) {
                        G[t*Maxn + i*n + j][t*Maxn + i*n + jtmp] = 1;
                        G[t*Maxn + i*n + jtmp][t*Maxn + i*n + j] = 1;
                        count++;
                    }
                    else {
                        G[t*Maxn + i*n + j][t*Maxn + i*n + jtmp] = w;
                        G[t*Maxn + i*n + jtmp][t*Maxn + i*n + j] = w;
                    }
                    
                }
                else {
                    G[t*Maxn + i*n + j][t*Maxn + itmp*n + j] = w;
                    G[t*Maxn + itmp*n + j][t*Maxn + i*n + j] = w;
                    G[t*Maxn + i*n + j][t*Maxn + i*n + jtmp] = w;
                    G[t*Maxn + i*n + jtmp][t*Maxn + i*n + j] = w;
                }
                
                
            }
        }
        
    }
    
    
    for (i = 0; i <Maxn*Tp; i++)
    {
        for (j = 0; j < Maxn*Tp; j++) {
            if (G[i][j] == 0)
                G[i][j] = 9999;
            if (i == j)
                G[i][j] = 0;
        }
    }
    
    return count;
}

Gragh Creategragh(int size, int *matrixModel, double w, int Tp, double tw)
{
    Gragh T = (Gragh)malloc(sizeof(struct node));
    if (T == NULL) {
        printf("space failed to be allocated！1\n"); exit(1);
    }
    
    int m = size, n = size;
    int Maxn = size * size;
    int ***matrix;
    int i; int j = 0;
    int t = 0;
    matrix = (int ***)malloc(Tp * sizeof(int **));
    for (t = 0; t < Tp; t++) {
        matrix[t] = (int **)malloc(m * sizeof(int *));
        for (i = 0; i < m; i++) {
            matrix[t][i] = (int *)malloc(m * sizeof(int));
        }
    }
    
    T->G = (double **)calloc(Maxn*Tp, sizeof(double *));
    for (i = 0; i < Maxn*Tp; i++) {
        T->G[i] = (double *)calloc(Maxn*Tp, sizeof(double));
        if (T->G[i] == NULL) {
            printf("space failed to be allocated！2\n"); exit(1);
        }
    }
    
    for (t = 0; t < Tp; ++t)
    {
        for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j)
                matrix[t][i][j] = matrixModel[t*Maxn + i * m + j];
        }
    }
    
    T->Nv = m * n*Tp;
    T->Ne = getConnection(matrix, T->G, size,w,Tp,tw);
    
    for (t = 0; t < Tp; t++){
        for (i = 0; i < m; i++)
            free(matrix[t][i]);
    }
    return T;
}


void dijkstra(Gragh g, int vs, int*prev, double*dist, int Maxn)
{
    int i, j, k;
    double min;
    double tmp;
    int *flag;      // flag[i]=1表示"顶点vs"到"顶点i"的最短路径已成功获取。
    flag = (int *)malloc(Maxn * sizeof(int));
    
    int INF = 9999;
    // 初始化
    for (i = 0; i < g->Nv; i++)
    {
        flag[i] = 0;              // 顶点i的最短路径还没获取到。
        prev[i] = 0;              // 顶点i的前驱顶点为0。
        dist[i] = g->G[vs][i];// 顶点i的最短路径为"顶点vs"到"顶点i"的权。
    }
    
    // 对"顶点vs"自身进行初始化
    flag[vs] = 1;
    dist[vs] = 0;
    
    // 遍历G.vexnum-1次；每次找出一个顶点的最短路径。
    for (i = 0; i < g->Nv; i++)
    {
        // 寻找当前最小的路径；
        // 即，在未获取最短路径的顶点中，找到离vs最近的顶点(k)。
        min = INF;
        for (j = 0; j < g->Nv; j++)
        {
            if (flag[j] == 0 && dist[j]<min)
            {
                min = dist[j];
                k = j;
            }
        }
        // 标记"顶点k"为已经获取到最短路径
        flag[k] = 1;
        
        // 修正当前最短路径和前驱顶点
        // 即，当已经"顶点k的最短路径"之后，更新"未获取最短路径的顶点的最短路径和前驱顶点"。
        for (j = 0; j < g->Nv; j++)
        {
            tmp = (g->G[k][j] == INF ? INF : (min + g->G[k][j])); // 防止溢出
            if (flag[j] == 0 && (tmp  < dist[j]))
            {
                dist[j] = tmp;
                prev[j] = k;
            }
        }
    }
    
    free(flag);
}


#endif
