#pragma once
#define MAX_VERTEX_NUM 100
#define INFINITY_GRAPH 65537
#define MAX 100
#define m_BTREE 3				//m阶B-树
//#define RADIX_BIT 4				//基数排序位数
#define RADIX  10				//基数排序进制

class Matrix
{
public:
	Matrix();
	Matrix(int l, int r);
	~Matrix();
	void init();
	void init(int I, int J, int N);
	void print();
	Matrix& Add(Matrix &another)throw(int);
	void operator+=(const Matrix &another);
	friend Matrix operator+(Matrix& A, Matrix& B)throw(int);
	Matrix& Subtract(Matrix &another)throw(int);
	Matrix& operator-(Matrix& another)throw(int);
	friend Matrix operator-(Matrix& A, Matrix& B)throw(int);
	Matrix& warshall()throw(int);   //华沙算法 transitive  
	void SymmerticMatric(int* &result)throw(int);  //保存对称矩阵
	int* TriangleMatrix(int up)throw(int);   //保存三角矩阵 up=1:上三角矩阵 up=0:下三角矩阵
	int getNum(int i, int j);
	int getLine();
	int getRow();

private:
	int lines, rows;
	int **num;

};
 
void clean_stdin();		//相当于fflush(stdin)

//栈
typedef struct node
{
	char data;
	int intdata;
	struct node* next;
}Node, *pNode;

typedef struct stack  //链栈
{
	pNode top;
}Stack;

void InitStack(Stack& s);
void push(Stack& s, char c);
void push(Stack& s, int n);
char pop(Stack& s);
int popint(Stack& s);
bool isEmpty(Stack& s);

//队列
typedef struct cqueue  //链队列 
{
	pNode front;
	pNode rear;
}CQueue;

void InitQueue(CQueue& q);
void EnQueue(CQueue& q, char x);
void EnQueue(CQueue& q, int x);
char DeQueue(CQueue& q);
int DeQueueInt(CQueue& q);
bool isEmpty(CQueue& q);

//KMP字符串匹配
void get_next(const char* str, int* next);
void get_nextval(const char* str, int* next);



//三元组存储矩阵
typedef struct tuple
{
	int i, j, v;   //i为行号，j为列号，v为值 
}Tuple;

typedef struct sparmattp
{
	int m, n, t,max;  //m为行数，n为列数，t为非零元素个数,max为最大存储个数
	Tuple data[10];
}Spar;

void InitSpar(Spar& s, int m, int n);  //初始化矩阵 
Spar Multiply(Spar a, Spar b);        //矩阵相乘 
Spar Add(Spar a, Spar b);			 //矩阵相加 
void Print(Spar s);					 //打印矩阵 
Spar Tranpos(Spar a);
int Find(Spar s, int m, int n);		//寻找矩阵第m行第n列的数据 

#if 0
Matrix FromSymmetric(const int* num, int len);		//对称矩阵的恢复
Matrix FromTriangle(const int* num, int len, int up);//三角矩阵的恢复 up=1：还原为上三角 up=0：下三角
#endif

//十字链表存储矩阵
typedef struct OLnode
{
	int i, j;
	int e;
	struct OLnode *right, *down;
}OLNode,*OList;

typedef struct
{
	OList *Rhead, *Chead;
	int mu, nu, tu;
}CrossList;

int CreateCL(CrossList* M);
void DestroyCL(CrossList* M);

typedef struct ChainMatrix
{
	int row, col;
	struct ChainMatrix *right, *down;
	union {
		int val;
		struct ChainMatrix* next;
	};
}CMatrix,*PMatrix;

void InitPM(PMatrix& p);

//二叉树  
typedef struct btnode
{
	int data;
	struct btnode *lc, *rc;
//	struct btnode *parent;
}BiTNode,*BiTree;

BiTree preCreatBiTree();
void preOrderTraverse(BiTree T);
void inOrderTraverse(BiTree T);
void postOrderTraverse(BiTree T);
bool Similar(BiTree p, BiTree q);

//二叉排序树
typedef struct BinarySortTree
{
	int key;
	struct BinarySortTree* lc;
	struct BinarySortTree* rc;
}BST,*BSTree;

BSTree SearchBST(BSTree T, int key);	//查找
BSTree InsertBST(BSTree T, int key);
BSTree CreateBST(BSTree T);				//创建二叉排序树
BSTree DeleteBST(BSTree T, int key);	//删除值为key的结点
void DeleteFullBST(BSTree* T);			//删除整棵树
void printBST(BSTree T);		//中序遍历
void JudgeBST(BiTree T, bool& flag);		//判断二叉树是不是二叉排序树,flag初始值为true

//平衡二叉树AVL
typedef struct avltree
{
	struct avltree* lc;
	struct avltree* rc;
	int key;
	int bf;			//平衡因子
}AVL,*AVLTree;

typedef enum
{
	EH = 0,	LH = 1,	RH = -1
	//EH：两边高度相同 LH:左边高 RH:右边高
}AVL_T;

int insertAVL(AVLTree* T, int key, bool* taller);	
//taller初始值为false
void LeftBalance(AVLTree* T, bool* taller);		//左边高
void RightBalence(AVLTree* T, bool* taller);		//右边高
void R_Rotate(AVLTree* T);		//顺时针转
void L_Rotate(AVLTree* T);		//逆时针转
AVLTree seachAVL(AVLTree T, int key);		//返回所在结点
int deleteAVL_Node(AVLTree* T, int key, bool *shorter);	//删除某个结点,删除前将shorter置为false
void destroyAVL(AVLTree* T);			//free整棵树
void preOrderTraverse(AVLTree T);		//前序遍历
void inOrderTraverse(AVLTree T);		//中序遍历
void printAVL(AVLTree T, int key, int direction);	//输出整棵树
//调用形式为   printAVL(root,root->key,0) 
//direction为0，表示该节点是根节点 为1表示是父结点的右孩子，为-1表示左孩子


//B-树
typedef struct bmtnode
{
	int key;
//	int num;
	struct BMinusTree* child;
}BMTNode[m_BTREE +1];

typedef struct BMinusTree			//B-树
{
	int num;
	struct BMinusTree* parent;
	BMTNode node;
}BMT,*BMTree;

typedef struct BResult		//B树的查找结果类型 
{
	BMTree p;
	int flag;		//1表示查找成功，0表示查找失败
	int i;			//在节点中关键字序号，[1,m]
}BRNode;

void InitBMTree(BMTree &t);				//初始化
void DestroyBMTree(BMTree &t);			//free整棵树
int Search_Insert(BMTree t,int k);		//找到插入的地方
void Insert(BMTree &t, int k, int i, BMTree ap);	//直接插入到i+1的位置
void NewRoot(BMTree &t, int key,BMTree &ap);		//生成新的根结点
void Spilt(BMTree &t, BMTree &ap);		//分割为两个结点，前一半保留，后一半移入新结点ap
void InsertBMTree(BMTree &t, int key, BMTree q,int i);	//在t树的q结点的第i和第i+1个位置插入key


//邻接矩阵
typedef struct    
{
	char vexs[MAX_VERTEX_NUM];
	int arc[MAX_VERTEX_NUM][MAX_VERTEX_NUM];
	int numVertex, numEdge;
}Graph, *pGraph;

int locate(const Graph g, char ch);
void createGraph(Graph &g, int directed=0);	//directed: 1：有向图 0：无向图
void printGraph(const Graph g);
void find_cycle(const Graph g);		//有向图找简单回路
void dfs(const Graph g, int v);		//深度优先遍历找回路
void Print_dfs(const Graph g, int v, int start);
void prim(const Graph g);			//prim算法，求无向图最小生成树
void PrintMST(Graph g);				//打印最小生成树
bool topologicalsort(const Graph g);	//拓扑排序
void Dijkstra(const Graph g,int start,int** result);		//从某个源点到其余各顶点的最短路径
void Floyd(const Graph g, int result[MAX_VERTEX_NUM][MAX_VERTEX_NUM]);			//Floyd算法 求任意两顶点最短路径

typedef struct Edgenode
{
	int adjvex;
	int weight;
	struct 	Edgenode* next;
}EdgeNode;

typedef struct Vertexnode
{
	char data;
	EdgeNode* firstedge;
}VertexNode, AdjList[MAX_VERTEX_NUM];

typedef struct        //邻接表（有向 无向均可）
{
	AdjList adjlist;
	int numVertex, numEdge;
}GraphList;

typedef struct edge
{
	char start, end;
	int weight;
}Edge;

int locate(const GraphList g, char ch);
void createGraph(GraphList &g, int inverse = 0);	//inverse:0邻接表 1:逆邻接表
void printGraph(const GraphList g);
void dfs(const GraphList g, int v);		//深度遍历找连通分量并输出
void CountCP(const GraphList g,int choose);		//求无向图连通分量个数，choose为1调用dfs，为0为调用bfs
void bfs(const GraphList g, int v);		//广度遍历找连通分量并输出
void CountSCP(const GraphList g);		//求有向图强连通分量
void tarjan(const GraphList g, Stack& s, int v);	//tarjan算法
void kruskal(const GraphList g);			//求无向图最小生成树
//获取g中<start,end>的权值，若不连通，返回无穷大
int getWeight(const GraphList g, int start, int end);
void get_edges(const GraphList g, Edge** edges);		//获取所有边
void sorted_edges(Edge** edges, int elen);	//按从小到大排序所有边
int get_end(int vends[],int i);		
bool topologicalsort(const GraphList g);	//传入邻接表：拓扑排序  传入逆邻接表：逆拓扑排序
void Dijkstra(const GraphList g, int start, int** result);		//从某个源点到其余各顶点的最短路径
void CriticalPath(const GraphList g);		//求关键路径


typedef struct cgNode     //弧的十字链表 
{
	int tailvex, headvex, weight;
	struct cgNode *tlink, *hlink;
}*CGNode;

typedef struct//十字链表顶点 
{
	char data;
	struct cgNode *firstin, *firstout;
}CGVertex;

typedef	struct
{
	CGVertex vertex[MAX_VERTEX_NUM];
	int vertexNum, arcNum;
}OrthList;		//图的十字链表 （有向图）

int locate(const OrthList g, char ch);
void createGraph(OrthList &g);
void printGraph(const OrthList g);

typedef struct ebox
{
	int mark;	//标记是否访问  0:unvisited	1:visited
	int weight, ivex, jvex;
	struct ebox *ilink, *jlink;
}EBox;

typedef struct
{
	char data;
	EBox* firstedge;
}VexBox;

typedef struct		//邻接多重图（无向图） 
{
	VexBox	AML[MAX_VERTEX_NUM];
	int vertexNum, edgeNum;
}AMLGraph;

int locate(AMLGraph g, char ch);
void createGraph(AMLGraph &g);
void printGraph(const AMLGraph g);
bool getVex(const AMLGraph g, int v, char& ch);	//v的范围是 0 -- g.vertexNum - 1 
bool insertVex(AMLGraph& g, char ch);		//插入一个新结点，但不插入边
bool delEdge(AMLGraph& g, char v, char w);	//删除边
bool delVex(AMLGraph& g, char v);		//删除顶点和相应的边
void destroyGraph(AMLGraph& g);			//删除图
int FirstAdjVex(const AMLGraph g, char v);	//返回v的第一个邻接顶点序号 
int NextAdjVex(const AMLGraph g, char v, char w);	//返回v相对于w的下一个邻接顶点的序号 
bool insertEdge(AMLGraph& g, char v, char w);	//增加v到w的边 
void MarkUnvisited(AMLGraph& g);		//将所有边的访问标记置为未访问 
void DFS(const AMLGraph g, int v);		//@DFSTraverse 
void DFSTraverse(const AMLGraph g);		//深度优先遍历 
void BFSTraverse(const AMLGraph g);		//广度优先遍历 

void InsertSort(int a[], int n);		//插入排序(非递减)  从1开始 a[0]作为哨兵  n为a的长度
void BInsertSort(int a[], int n);		//折半插入排序(非递减)  a[0]作为哨兵  n为a的长度
void BubbleSort(int a[], int n);		//双向冒泡排序(非递减) 从0开始
void HillSort(int a[], int n);			//希尔排序(非递减) 增量为n/2
void SelectSort(int a[], int n);		//选择排序
void QuickSort(int a[], int low, int high);		//快速排序
int Partition_QSort(int a[], int low, int high);		//快速排序时的划分

typedef struct
{
	int data[100];	//data[0]不存放数据，用来辅助
	int length;		//元素个数
}heap,*Heap;
void HeapSort(Heap h);				//堆排序(非增序)
void HeapSift(Heap h,int k,int m);	//堆排序调整--筛选（小根堆）
void MergeSort(int a[], int n);			//归并排序 a[0]未用 非降序
										//归并操作   r为待排序数组，s为开始元素下标，t为结束元素下标，q为该部分数组最后元素下标+1
void Merge(int r[], int s1,int t1,int s2,int t2,int &q );				
void RadixSort(int *a, int n, int d);		//基数排序 d为最大位数 不受数组a长度限制
int get_value(int n, int d);		//返回n第d位的数字  d: 从1（个位）开始
