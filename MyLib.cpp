#include "pch.h"
#include "MyLib.h"
#include<stdio.h>
#include<iostream>
#include<cmath>  
#include<stdlib.h>
#include<malloc.h>

#define MAX_AB(a,b) ((a)>(b)?(a):(b))


using namespace std;
int NEED_NEW = 1;
//图部分全局变量
int visited[MAX_VERTEX_NUM];	//0：未访问 1：以访问但还有邻接结点没访问 2：全访问完
bool visite[MAX_VERTEX_NUM];	//访问标志数组（全局变量）

int low[MAX_VERTEX_NUM];		//用于tarjan算法 prim算法
int dfn[MAX_VERTEX_NUM];		//用于tarjan算法
int inStack[MAX_VERTEX_NUM];	//tarjan
int tarjan_time = 0;

int parent[MAX_VERTEX_NUM];		//prim算法
int cost[MAX_VERTEX_NUM][MAX_VERTEX_NUM];	//Floyd算法

BiTree pre_judgeBST = NULL;			//用于JudgeBST

Matrix::Matrix()
{

	//	cin>>lines>>rows;
	//    
	//    num = new int* [lines];   
	//    for(int i = 0;i<lines;i++)
	//    {
	//       num[i] = new int[rows];
	//    }
}

Matrix::Matrix(int l, int r)
{
	lines = l;
	rows = r;
	num = new int*[l];
	for (int i = 0; i < l; i++)
	{
		num[i] = new int[r];
	}
	NEED_NEW = 0;
}

void Matrix::init()
{
	int i, j;
	cout << "init" << endl;
	if (NEED_NEW)
	{
		cout << "input lines and rows:";
		cin >> lines >> rows;
		num = new int*[lines];
		for (int i = 0; i < lines; i++)
		{
			num[i] = new int[rows];
		}
	}

	for (i = 0; i < lines; i++)
	{
		cout << "请输入第" << i + 1 << "行的" << rows << "个元素: ";
		for (j = 0; j < rows; j++)
		{
			cin >> (*(num + i))[j];
		}
	}

}
void Matrix::init(int I, int J, int N)
{
	(*(num + I))[J] = N;
}

void Matrix::print()
{
	int i, j;
	cout << "Print" << endl;
	for (i = 0; i < lines; i++)
	{
		for (j = 0; j < rows; j++)
		{
			cout << (*(num + i))[j] << "  ";
		}
		cout << endl;
	}
}

Matrix& Matrix::operator-(Matrix& another) throw(int)
{
	cout << "Subtract" << endl;
	if (!(lines == another.lines && rows == another.rows))  //error 
		throw(-1);
	static Matrix Answer(lines, rows);

	for (int i = 0; i < Answer.lines; i++)
	{
		for (int j = 0; j < Answer.rows; j++)
			Answer.init(i, j, this->num[i][j] - another.num[i][j]);
	}
	return Answer;
}

void Matrix::operator+=(const Matrix &another)throw(int)
{
	cout << "+=" << endl;
	if (!(lines == another.lines && rows == another.rows))  //error 
		throw(-1);

	for (int i = 0; i < lines; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			num[i][j] += another.num[i][j];
		}
	}
}

Matrix& Matrix::Add(Matrix &another)throw(int)
{
	cout << "Add" << endl;
	if (!(lines == another.lines && rows == another.rows))  //error 
		throw(-1);

	static Matrix Answer(lines, rows);
	int i, j;
	for (i = 0; i < another.getLine(); i++)
	{
		for (j = 0; j < another.getRow(); j++)
		{
			Answer.num[i][j] = this->getNum(i, j) + another.getNum(i, j);
		}
	}
	return Answer;
}

Matrix& Matrix::Subtract(Matrix &another)throw(int)
{
	cout << "Subtract" << endl;
	if (!(lines == another.lines && rows == another.rows))  //error 
		throw(-1);
	static Matrix Answer(lines, rows);
	int i, j;
	for (i = 0; i < another.getLine(); i++)
	{
		for (j = 0; j < another.getRow(); j++)
		{
			Answer.num[i][j] = this->getNum(i, j) - another.getNum(i, j);
		}
	}
	return Answer;
}

Matrix& Matrix::warshall()throw(int)
{
	if (lines != rows)
		throw(-1);
	cout << "warshall" << endl;
	static Matrix Answer(lines, lines);
	int i, j, k;
	for (i = 0; i < lines; i++)
		for (j = 0; j < lines; j++)
			Answer.init(i, j, (*(num + i))[j]);
	for (k = 0; k < lines; k++)
		for (i = 0; i < lines; i++)
			for (j = 0; j < lines; j++)
				if (Answer.getNum(i, j) || (Answer.getNum(i, k) && Answer.getNum(k, j)))
					Answer.init(i, j, 1);
	return Answer;
}

void Matrix::SymmerticMatric(int* &result)throw(int)
{
	if (lines != rows)
		throw(-1);
	result = new int[(lines*(lines + 1)) / 2];
	int i, j, n;
	for (i = 0; i < lines; i++)
		for (j = 0; j <= i; j++)
		{
			n = i*(i + 1) / 2 + j;
			result[n] = num[i][j];
		}
}

int* Matrix::TriangleMatrix(int up)throw(int)
{
	if (lines != rows && lines != 1)
		throw(-1);
	int count = lines*(lines + 1) / 2 + 1;
	static int* result = new int[count];
	int i, j, k;
	if (1 == up)
	{
		for (i = 0; i < lines; i++)
			for (j = i; j < lines; j++)
			{
				k = i*(2 * lines - i + 1) / 2 + j - i;
				*(result + k) = num[i][j];
			}
		*(result + count - 1) = num[1][0];
	}
	else
	{
		for (i = 0; i < lines; i++)
			for (j = 0; j <= i; j++)
			{
				k = i*(i + 1) / 2 + j;
				*(result + k) = num[i][j];
			}
		*(result + count - 1) = num[0][1];
	}
	return result;
}

int Matrix::getNum(int i, int j)
{
	return (*(num + i))[j];
}

int Matrix::getLine()
{
	return lines;
}

int Matrix::getRow()
{
	return rows;
}


Matrix::~Matrix()
{
	for (int i = 0; i < lines; i++)
		delete[] num[i];
	delete[] num;
	cout << "deleteMatrix" << endl;
}

Matrix operator+(Matrix& A, Matrix& B)throw(int)
{
	cout << "ADD" << endl;
	if (!(A.lines == B.lines && A.rows == B.rows))  //error 
		throw(-1);

	static Matrix C(A.lines, A.rows);
	for (int i = 0; i < C.lines; i++)
		for (int j = 0; j < C.rows; j++)
			C.init(i, j, A.num[i][j] + B.num[i][j]);
	return C;
}


Matrix operator-(Matrix& A, Matrix& B)throw(int)
{
	cout << "Subtract" << endl;
	if (!(A.lines == B.lines && A.rows == B.rows))  //error 
		throw(-1);

	static Matrix C(A.lines, A.rows);
	for (int i = 0; i < C.lines; i++)
		for (int j = 0; j < C.rows; j++)
			C.init(i, j, A.num[i][j] - B.num[i][j]);
	return C;

}

void clean_stdin()
{
	char c;
	do {
		c = getchar();
	} while (c != '\n'&&c != EOF);
}//clean_stdin

//栈
void InitStack(Stack& s)
{
	s.top = NULL;
}//InitStack

void push(Stack& s, char c)
{
	pNode p;
	p = (pNode)malloc(sizeof(Node));
	if (!p)
	{
		printf("fail to alloc");
		return;
	}
	p->data = c;
	p->next = s.top;
	s.top = p;
}//push 

void push(Stack& s, int n)
{
	pNode p;
	p = (pNode)malloc(sizeof(Node));
	if (!p)
	{
		printf("fail to alloc");
		return;
	}
	p->intdata = n;
	p->next = s.top;
	s.top = p;
}//push

char pop(Stack& s)
{
	if (isEmpty(s))
	{
		printf("isEmpty");
		return NULL;
	}
	char c = s.top->data;
	pNode p = s.top->next;
	free(s.top);
	s.top = p;
	return c;
}//pop 

int popint(Stack& s)
{
	if (isEmpty(s))
	{
		printf("isEmpty");
		return NULL;
	}
	int c = s.top->intdata;
	pNode p = s.top->next;
	free(s.top);
	s.top = p;
	return c;
}//popint 

bool isEmpty(Stack& s)
{
	return s.top == NULL;
}//isEmpty 


//队列
void InitQueue(CQueue& q)
{
	q.front = (pNode)malloc(sizeof(Node));
	if (q.front == NULL) printf("error to malloc");
	q.rear = q.front;
}//InitQueue

void EnQueue(CQueue& q, char x)
{
	pNode p;
	p = (pNode)malloc(sizeof(Node));
	if (!p)
	{
		printf("fail to malloc");
		return;
	}
	p->data = x;            //用新元素p存储数据， 
	p->next = NULL;
	q.rear->next = p;       //队尾的下一个指向p
	q.rear = p;            // 队尾指向p 
}//EnQueue

void EnQueue(CQueue& q, int x)
{
	pNode p;
	p = (pNode)malloc(sizeof(Node));
	if (!p)
	{
		printf("fail to malloc");
		return;
	}
	p->intdata = x;            //用新元素p存储数据， 
	p->next = NULL;
	q.rear->next = p;       //队尾的下一个指向p
	q.rear = p;            // 队尾指向p 
}//EnQueue

char DeQueue(CQueue& q)
{
	if (isEmpty(q))    //先判断不是空队列 
	{
		printf("isEmpty");
		return NULL;
	}
	char data = q.front->next->data;
	pNode p = q.front->next;        //处理掉旧队头，并让q.next指向新队头 
	q.front->next = p->next;
	free(p);
	return data;
}//DeQueue

int DeQueueInt(CQueue& q)
{
	if (isEmpty(q))    //先判断不是空队列 
	{
		printf("isEmpty");
		return NULL;
	}
	int data = q.front->next->intdata;
	pNode p = q.front->next;        //处理掉旧队头，并让q.next指向新队头 
	q.front->next = p->next;
	free(p);
	return data;
}//DeQueueInt

bool isEmpty(CQueue& q)
{
	return q.front->next == NULL;
}//isEmpty

//KMP
void get_next(const char* str, int* next)
{
	int i = 0, j = -1;
	next[0] = -1;
	while (i < strlen(str))
	{
		if (str[i] == str[j] || -1 == j)
			next[++i] = ++j;
		else j = next[j];
	}
}

void get_nextval(const char* str, int* next)
{
	int i = 0, j = -1;
	next[0] = -1;
	while (i < strlen(str))
	{
		if (str[i] == str[j] || -1 == j)
			if (str[++i] == str[++j])
				next[i] = next[j];
			else	next[i] = j;
		else j = next[j];
	}
}

//三元组矩阵
void InitSpar(Spar& s, int m, int n)
{
	int i, j, v;    //i j为计数器 v为读取的矩阵元素值 
	s.m = m;
	s.n = n;
	s.t = 0;
	for (i = 0; i < m; i++)
	{
		fflush(stdin);
		printf("please input %d elements of No.%d row:", n, i + 1);
		for (j = 0; j < n; j++)
		{
			cin >> v;
			if (v)
			{
				if (s.t < s.max)
				{
					s.data[s.t].i = i + 1;
					s.data[s.t].j = j + 1;
					s.data[s.t].v = v;
					s.t++;
				}
				else	printf("The maximum number of stores has been exceeded\n");
			}
		}
	}
}//InitSpar

Spar Multiply(Spar a, Spar b)
{
	Spar c;
	c.m = a.m;
	c.n = b.n;
	int tc = 0;
	if (!(a.t*b.t))
	{
		c.t = 0;
		return c;
	}
	int k, num, A, B;
	for (int i = 1; i <= c.m; i++)
	{
		for (int j = 1; j <= c.n; j++)
		{
			num = 0;
			k = 1;
			while (k <= a.n)
			{
				A = Find(a, i, k);
				B = Find(b, k, j);
				if (A && B)
					num += A*B;
				k++;
			}
			if (num)
			{
				c.data[tc].i = i;
				c.data[tc].j = j;
				c.data[tc++].v = num;
			}
		}
	}
	return c;
}//Multiply

Spar Add(Spar a, Spar b)
{
	Spar c;
	c.m = a.m;
	c.n = a.n;
	c.t = 0;
	int ta = 0, tb = 0;
	for (int i = 1; i <= c.m; i++)
		for (int j = 1; j <= c.n; j++)
		{
			c.data[c.t].v = 0;
			if (a.data[ta].i == i && a.data[ta].j == j)
			{
				c.data[c.t].v += a.data[ta].v;
				ta++;
			}

			if (b.data[tb].i == i && b.data[tb].j == j)
			{
				c.data[c.t].v += b.data[tb].v;
				tb++;
			}
			if (c.data[c.t].v)
			{
				c.data[c.t].i = i;
				c.data[c.t].j = j;
				c.t++;
			}
		}
	return c;
}//Add

void Print(Spar s)
{
	int i, j, k = 0, v;
	for (i = 0; i < s.m; i++)
	{
		for (j = 0; j < s.n; j++)
		{
			if (s.data[k].i == i + 1 && s.data[k].j == j + 1)
			{
				v = s.data[k].v;
				k++;
			}
			else v = 0;
			printf("%d	", v);
		}
		printf("\n");
	}
}//Print

int Find(Spar s, int m, int n)
{
	int i;
	for (i = 0; i < s.t; i++)
	{
		if (s.data[i].i == m && s.data[i].j == n)
			return s.data[i].v;
	}
	return 0;
}

Spar Tranpos(Spar a)
{
	Spar b;
	b.m = a.n;
	b.n = a.m;
	b.t = a.t;
	int i, j, pos;
	int *num = new int[a.m + 1];
	int *cpot = new int[a.m + 1];
	//num存放a中每列非零元素个数，cpot存放a中每列第一个非零元素在b中的位置 
	if (b.t)
	{
		for (i = 1; i <= a.m; i++)
			num[i] = 0;
		for (i = 0; i < a.t; i++)
			num[a.data[i].j]++;
		cpot[1] = 0;
		for (i = 2; i <= a.n; i++)
			cpot[i] = cpot[i - 1] + num[i - 1];
		for (i = 0; i < a.t; i++)
		{
			j = a.data[i].j;
			pos = cpot[j];
			b.data[pos].i = a.data[i].j;
			b.data[pos].j = a.data[i].i;
			b.data[pos].v = a.data[i].v;
			cpot[j]++;
		}
	}
	//	printf("tranpos\n");
	return b;
}

//十字链表存储矩阵
int CreateCL(CrossList* M)
{
	int i, j, m, n, t;
	int k, flag;
	int e;
	OList p, q;
	if (M->Rhead)
		DestroyCL(M);
	do {
		flag = 1;
		printf("input lines and rows and numbers of elements:");
		cin >> m >> n >> t;
		if (m < 0 || n < 0 || t<0 || t>m*n)
			flag = 0;
	} while (!flag);
	M->mu = m;
	M->nu = n;
	M->tu = t;
	M->Rhead = (OList*)malloc(sizeof(OLNode)*(m + 1));
	if (!M->Rhead)	exit(-1);
	M->Chead = (OList*)malloc(sizeof(OLNode)*(n + 1));
	if (!M->Chead)	exit(-1);
	for (k = 1; k <= m; k++)
		M->Rhead[k] = NULL;
	for (k = 1; k <= n; k++)
		M->Chead[k] = NULL;

	for (k = 1; k <= t; k++)
	{
		do {
			flag = 1;
			printf("input No.%d element's line and row and value:", k);
			cin >> i >> j >> e;
			if (i <= 0 || j <= 0)	flag = 0;
		} while (!flag);
		p = (OList)malloc(sizeof(OLNode));
		if (p == NULL) exit(-1);
		p->i = i;
		p->j = j;
		p->e = e;
		//行插入
		if (M->Rhead[i] == NULL || M->Rhead[i]->j > j)
		{
			p->right = M->Rhead[i];
			M->Rhead[i] = p;
		}
		else
		{
			for (q = M->Rhead[i]; q->right&&q->right->j < j; q = q->right);
			p->right = q->right;
			q->right = p;
		}
		//列插入
		if (M->Chead[j] == NULL || M->Chead[j]->i > i)
		{
			p->down = M->Chead[j];
			M->Chead[j] = p;
		}
		else
		{
			for (q = M->Chead[j]; q->down&&q->down->i < i; q = q->down);
			p->down = q->down;
			q->down = p;
		}
	}
	return 1;
}//CreateCL

void DestroyCL(CrossList* M)
{
	int i = 1;
	OList temp, temp2;
	while (M->Rhead[i] != NULL)
	{
		temp = M->Rhead[i];
		while (temp->right)
		{
			temp2 = temp->right;
			free(temp);
			temp = temp2;
		}
		free(temp);
		free(M->Rhead[i]);
		i++;
	}
	i = 1;
	while (M->Chead[i] != NULL)
	{
		free(M->Chead[i]);
		i++;
	}
	free(M->Rhead);
	free(M->Chead);
	free(M);
}//DestroyCL

void InitPM(PMatrix& H)
{
	PMatrix p, q, hd[MAX];
	int i, j, k, m, n, t, v;
	printf("input lines and rows and number of elements:");
	cin >> m >> n >> t;
	H = (PMatrix)malloc(sizeof(CMatrix));
	H->col = m;
	H->row = n;
	hd[0] = H;
	int s = m > n ? m : n;
	for (i = 1; i <= s; i++)
	{
		p = (PMatrix)malloc(sizeof(CMatrix));
		p->row = 0;	p->col = 0;
		p->right = p;	p->down = p;
		hd[i] = p;
		hd[i - 1]->next = p;
	}
	hd[s]->next = H;
	for (k = 1; k <= t; k++)
	{
		printf("input Np.%d elements(row line value):", k);
		cin >> i >> j >> v;
		p = (PMatrix)malloc(sizeof(CMatrix));
		p->row = i; p->col = j; p->val = v;
		q = hd[i];
		while (q->right != hd[i] && q->right->col < j)
			q = q->right;
		p->right = q->right;
		q->right = p;
		q = hd[j];
		while (q->down != hd[j] && q->down->row < i)
			q = q->down;
		p->down = q->down;
		q->down = p;
	}
}//InitPM


BiTree preCreatBiTree()
{
	BiTree T;
	int v;
	char ch;
	cin >> ch;
	if (ch == '#') T = NULL;
	else
	{
		v = ch - '0';
		T = (BiTree)malloc(sizeof(BiTNode));
		T->data = v;
		T->lc = preCreatBiTree();
		T->rc = preCreatBiTree();
	}
	return T;
}// preCreatBiTree

void preOrderTraverse(BiTree T)
{
	if (T)
	{
		cout << T->data << endl;
		if (T->lc != NULL)		preOrderTraverse(T->lc);
		if (T->rc != NULL)		preOrderTraverse(T->rc);
	}
}//preOrderTraverse

void inOrderTraverse(BiTree T)
{
	if (T)
	{
		if (T->lc != NULL)		inOrderTraverse(T->lc);
		cout << T->data << endl;
		if (T->rc != NULL)		inOrderTraverse(T->rc);
	}
}//inOrderTraverse

void postOrderTraverse(BiTree T)
{
	if (T)
	{
		if (T->lc != NULL)		postOrderTraverse(T->lc);
		if (T->rc != NULL)		postOrderTraverse(T->rc);
		cout << T->data << endl;
	}
}//postOrderTraverse

bool Similar(BiTree p, BiTree q)
{
	return ((p == NULL&&q == NULL) || (p != NULL && q != NULL && p->data == q->data
		&& Similar(p->lc, q->lc) && Similar(p->rc, q->rc)));
}//Similar

 //二叉排序树
BSTree SearchBST(BSTree T, int key)
{
	if (T != NULL && T->key == key)
		return T;
	else
	{
		if (T->key > key  && T->lc != NULL)
			return SearchBST(T->lc, key);
		else if (T->key < key && T->rc != NULL)
			return SearchBST(T->rc, key);
	}
	return NULL;
}//SearchBST

BSTree InsertBST(BSTree T, int key)
{
	BSTree S;
	if (T == NULL)
	{
		S = (BSTree)malloc(sizeof(BST));
		S->key = key;
		S->lc = NULL;
		S->rc = NULL;
		return S;
	}
	else
	{
		if (T->key > key)
			T->lc = InsertBST(T->lc, key);
		else if (T->key < key)
			T->rc = InsertBST(T->rc, key);
	}
	return T;
}//InsertBST

BSTree CreateBST(BSTree T)
{
	int e;
	printf("please input the key:(input 0 to exit)");
	cin >> e;
	while (e != 0)
	{
		clean_stdin();
		T = InsertBST(T, e);
		printf("please input the key:(input 0 to exit)");
		cin >> e;
	}
	return T;
}//CreateBST

BSTree DeleteBST(BSTree T, int key)
{
	BSTree f = NULL;	//指向p的父亲结点
	BSTree p = T;
	while (p != NULL)
	{
		if (p->key == key)
			break;
		f = p;
		if (p->key > key)
			p = p->lc;
		else p = p->rc;
	}
	if (p == NULL || p->key != key)		//没找到
	{
		printf("node %d not found\n", key);
		return T;
	}
	BSTree q = p;		//q是key所在结点
	if (p->lc != NULL && p->rc != NULL)		//有左右子树
	{
		BSTree s = p->lc;
		while (s->rc)	//左子树的最右下结点 p的前驱
		{
			q = s;
			s = s->rc;
		}
		p->key = s->key;
		if (p != q)
			q->rc = s->lc;
		if (p == q)
			q->lc = s->lc;
		free(s);
		return T;
	}
	else		//只有一个方向的子树
	{
		if (p->lc == NULL)
			p = p->rc;
		else if (p->rc == NULL)
			p = p->lc;
	}
	if (q == f->lc)
		f->lc = p;
	else f->rc = p;
	free(q);
	return T;
}//DeleteBST

void DeleteFullBST(BSTree* T)
{
	if (NULL != *T)
	{
		DeleteFullBST(&((*T)->lc));
		DeleteFullBST(&((*T)->rc));
		free(*T);
		*T = NULL;
	}
}//DeleteFullBST

void printBST(BSTree T)
{
	if (T == NULL)
		return;
	printBST(T->lc);
	printf("%d\n", T->key);
	printBST(T->rc);
}//printBST

void JudgeBST(BiTree T, bool& flag)		//调用前将pre_judgeBST置为null
{
	if (T != NULL && flag)		//利用中序遍历，为递增序列即为二叉排序树
	{
		JudgeBST(T->lc, flag);
		if (pre_judgeBST == NULL)
			pre_judgeBST = T;
		else if (pre_judgeBST->data < T->data)
			pre_judgeBST = T;
		else flag = false;
		JudgeBST(T->rc, flag);
	}
}//JudgeBST

 //平衡二叉树AVL

int insertAVL(AVLTree* T, int key, bool* taller)
{
	if (NULL == *T)
	{
		(*T) = (AVLTree)malloc(sizeof(AVL));
		(*T)->key = key;
		(*T)->lc = NULL;
		(*T)->rc = NULL;
		(*T)->bf = EH;
		*taller = true;
	}
	else
	{
		if (key == (*T)->key)
		{
			*taller = false;
			return 0;
		}
		if (key < (*T)->key)
		{
			if (!insertAVL(&(*T)->lc, key, taller))
				return 0;
			if (*taller)
			{
				switch ((*T)->bf)
				{
				case EH:
					(*T)->bf = LH;
					break;
				case LH:		
					LeftBalance(T, taller);
					*taller = false;
					break;
				case RH:		
					(*T)->bf = EH;
					*taller = false;
					break;
				default:
					break;
				}
			}
		}
		else
		{
			if (!insertAVL(&(*T)->rc, key, taller))
				return 0;
			if (*taller)
			{
				switch ((*T)->bf)
				{
				case EH:
					(*T)->bf = RH;
					break;
				case LH:
					(*T)->bf = EH;
					*taller = false;
					break;
				case RH:
					RightBalence(T, taller);
					*taller = false;
					break;
				default:
					break;
				}
			}
		}
	}
	return 1;
}//insertAVL

void LeftBalance(AVLTree* T, bool* taller)
{
	AVLTree temp = (*T)->lc;
	AVLTree LR = NULL;
	switch (temp->bf)
	{
	case LH:			//LL型
		(*T)->bf = temp->bf = EH;
		R_Rotate(T);
		break;
	case EH:			//不可能情况
		(*T)->bf = LH;
		*taller = true;
		break;
	case RH:			//LR型
		LR = temp->rc;
		switch (LR->bf)
		{
		case LH:
			(*T)->bf = RH;
			temp->bf = EH;
			break;
		case EH:			//LR是叶子结点
			(*T)->bf = temp->bf = EH; 
			break;
		case RH:
			(*T)->bf = EH;
			temp->bf = LH;
			break;
		default:
			break;
		}
		LR->bf = EH;
		L_Rotate(&temp);
		(*T)->lc = temp;
		R_Rotate(T);
		break;
	default:
		break;
	}
}//LeftBalance

void RightBalence(AVLTree* T, bool* taller)
{
	AVLTree temp = (*T)->rc;
	AVLTree RL = NULL;
	switch (temp->bf)
	{
	case LH:			//RL型
		RL = temp->lc;
		switch (RL->bf)
		{
		case LH:
			(*T)->bf = EH;
			temp->bf = RH;
			break;
		case RH:
			(*T)->bf = LH;
			temp->bf = EH;
			break;
		case EH:
			(*T)->bf = temp->bf = EH;
			break;
		default:
			break;
		}
		RL->bf = EH;
		R_Rotate(&temp);
		(*T)->rc = temp;
		L_Rotate(T);
		break;
	case EH:
		*taller = true;
		(*T)->bf = RH;
		break;
	case RH:				//RR型
		temp->bf = (*T)->bf = EH;
		L_Rotate(T);
		break;
	default:
		break;
	}
}//RightBalence

void R_Rotate(AVLTree* T)
{
	AVLTree temp = (*T)->lc;
	(*T)->lc = temp->rc;
	temp->rc = *T;
	*T = temp;
}//R_Rotate

void L_Rotate(AVLTree* T)
{
	AVLTree temp = (*T)->rc;
	(*T)->rc = temp->lc;
	temp->lc = *T;
	*T = temp;
}//L_Rotate

AVLTree seachAVL(AVLTree T, int key)
{
	if (T != NULL && T->key == key)
		return T;
	else
	{
		if (T->key > key && T->lc != NULL)
			return seachAVL(T->lc, key);
		else if (T->key < key && T->rc != NULL)
			return seachAVL(T->rc, key);
	}
	return NULL;
}//seachAVL

int deleteAVL_Node(AVLTree* T, int key,bool *shorter)
{
	if (NULL == (*T))
	{
		printf("no such node\n");
		*shorter = false;
		return 0;
	}
	if ((*T)->key == key)
	{
		AVLTree temp = *T;
		if ((*T)->lc == NULL)
		{
			*T = (*T)->rc;
			free(temp);
			*shorter = true;
		}
		else if((*T)->rc == NULL)
		{
			*T = (*T)->lc;
			free(temp);
			*shorter = true;
		}
		else		//被删除结点有左右孩子
		{
			if ((*T)->bf >= 0)		//被删结点左子树不比右子树矮
			{
				AVLTree front = temp->lc;
				while(front->rc != NULL)		//找到中序遍历前驱结点
					front = front->rc;
				temp->key = front->key;
				if (!deleteAVL_Node(&(temp->lc), temp->key, shorter))
					return 0;
				if(*shorter)
					switch ((*T)->bf)
					{
					case EH:
						(*T)->bf = RH;
						*shorter = false;
						break;
					case LH:
						(*T)->bf = EH;
						*shorter = true;
						break;
					case RH:
					{
						AVLTree RR = (*T)->rc;
						RR->bf = EH;
						(*T)->bf = EH;
						L_Rotate(T);
						break;
					}
					default:
						break;
					}
				/*switch ((*T)->bf)  //(*T)->bf >= 0
				{
					AVLTree LL = (*T)->lc;
					AVLTree LLL = LL->lc;
				case EH:
				{					
					switch (LL->bf)
					{
					case LH:
						switch (LLL->bf)
						{
						case EH:							
							LLL->bf = RH;
							R_Rotate(&LL);	
							(*T)->lc = LLL;
							break;
						case LH:
							(*T)->bf = RH;
							LLL->bf = EH;
							LL->bf = EH;
							R_Rotate(&LL);
							(*T)->lc = LLL;
							break;
						case RH:
							(*T)->bf = RH;
							LL->bf = EH;
							LLL->bf = EH;
							AVLTree LLLR = LLL->rc;
							L_Rotate(&LLL);
							LL->lc = LLLR;
							R_Rotate(&LL);
							break;
						default:
							break;
						}
						break;	
					case EH:
						LL->bf = LH;
						break;
					case RH:
						LL->bf = EH;
						break;
					default:
						break;
					}
					*shorter = false;
				}
					break;				
				case LH:
					switch (LL->bf)
					{
					case EH:
						LL->bf = LH;
						*shorter = false;
						break;
					case RH:
						LL->bf = EH;
						*shorter = true;
						break;
					case LH:
						LL->bf = EH;
						LLL->bf = EH;
						R_Rotate(&LL);
						(*T)->lc = LLL;
						break;
					default:
						break;
					}
					break;
				default:
					break;
				}*/
			}
			else			//右子树更高
			{
				AVLTree front = temp->rc;
				while (front->lc != NULL)		//找到中序遍历后继结点
					front = front->lc;
				temp->key = front->key;
				if (!deleteAVL_Node(&(temp->rc), temp->key, shorter))
					return 0;		
				if (*shorter)				//(*T)->bf == RH
					switch ((*T)->bf)
					{
					case EH:
						(*T)->bf = LH;
						*shorter = false;
						break;
					case LH:
					{
						AVLTree L = (*T)->lc;
						L->bf = EH;
						(*T)->bf = EH;
						R_Rotate(T);
						break;
					}
					case RH:
						(*T)->bf = EH;
						*shorter = true;
						break;
					default:
						break;
					}
				/*{
					AVLTree RR = (*T)->rc;
					AVLTree RRL = RR->lc;
					switch (RR->bf)
					{
					case EH:
						*shorter = false;
						if (RRL->bf == LH)
						{
							RR->bf = RH;
							RRL->bf = EH;
						}
						break;
					case RH:
						break;
					case LH:
						switch (RRL->bf)
						{
						case EH:
							RRL->bf = RH;
							*shorter = false;
							break;
						case LH:
							RRL->bf = EH;
							*shorter = true;
							break;
						case RH:
						{
							AVLTree RRLR = RRL->rc;					
							switch (RRLR->bf)
							{
							case EH:
								*shorter = false;
								RRLR->bf = LH;
								L_Rotate(&RRL);
								break;
							case LH:
							{
								*shorter = true;
								RRL->bf = EH;
								RRLR->bf = EH;
								AVLTree RRLRL = RRLR->lc;
								R_Rotate(&RRLR);
								RRL->rc = RRLRL;
								L_Rotate(&RRL); 
							}							
								break;
							case RH:
								*shorter = true;
								RRLR->bf = EH;
								RRL->bf = EH;
								L_Rotate(&RRL);
								break;
							default:
								break;
							}
						}
							break;
						default:
							break;
						}
						break;
					default:
						break;
					}
					*shorter = true;
					(*T)->bf = EH;
				}*/
			}
		}
	}
	else if ((*T)->key < key)		//被删结点可能在右子树
	{
		if (!deleteAVL_Node(&(*T)->rc, key, shorter))
			return 0;
		if (*shorter)
		{
			switch ((*T)->bf)
			{
			case EH:
				(*T)->bf = LH;
				*shorter = false;
				break;
			case LH:
			{
				AVLTree L = (*T)->lc;
				L->bf = EH;
				(*T)->bf = EH;
				R_Rotate(T);
				break;
			}				
			case RH:
				(*T)->bf = EH;
				*shorter = true;
				break;
			default:
				break;
			}
		}
	}
	else		//被删结点可能在左子树
	{
		if (!deleteAVL_Node(&(*T)->lc, key, shorter))
			return 0;
		if (*shorter)
		{
			switch ((*T)->bf)
			{
			case EH:
				(*T)->bf = RH;
				*shorter = false;
				break;
			case LH:
				(*T)->bf = EH;
				*shorter = true;
				break;
			case RH:
			{
				AVLTree RR = (*T)->rc;
				RR->bf = EH;
				(*T)->bf = EH;
				L_Rotate(T);
				break;
			}				
			default:
				break;
			}
		}
	}
	return 1;
}//deleteAVL_Node

void destroyAVL(AVLTree* T)
{
	if (NULL != *T)
	{
		destroyAVL(&((*T)->lc));
		destroyAVL(&((*T)->rc));
		free(*T);
		*T = NULL;
	}
}//destroy


void preOrderTraverse(AVLTree T)
{
	if (T)
	{
		cout << T->key << endl;
		if (T->lc != NULL)		preOrderTraverse(T->lc);
		if (T->rc != NULL)		preOrderTraverse(T->rc);
	}
}//preOrderTraverse

void inOrderTraverse(AVLTree T)
{
	if (T)
	{
		if (T->lc != NULL)		inOrderTraverse(T->lc);
		cout << T->key << endl;
		if (T->rc != NULL)		inOrderTraverse(T->rc);
	}
}//inOrderTraverse

void printAVL(AVLTree T, int key, int direction)
{
	if (NULL != T)
	{
		if (!direction)		//根结点
			printf("%d(%d) is root\n", T->key,T->bf);
		else
			if (direction == 1)		//父亲的右孩子
				printf("%d(%d) is %d's right child\n", T->key,T->bf, key);
			else if (direction == -1)	//左孩子
				printf("%d(%d) is %d's left child\n", T->key,T->bf, key);
		printAVL(T->lc, T->key, -1);
		printAVL(T->rc, T->key, 1);
	}
}//printAVL

//B-树
void InitBMTree(BMTree &t)
{
	t = NULL;
}//InitBMTree

void DestroyBMTree(BMTree &t)
{
	int i;
	if (t != NULL)
	{
		for (i = 0; i <= t->num; i++)
			DestroyBMTree(t->node[i].child);
		free(t);
		t = NULL;
	}
}//DestroyBMTree

int Search_Insert(BMTree t, int k)
{
	int i;
	int j = 0;
	for (i = 1; i <= t->num; i++)
		if (t->node[i].key < k)
			j = i;
	return j;
}//Search_Insert

void Insert(BMTree &t, int k, int i, BMTree ap)
{
	int j;
	if (k == t->node[i + 1].key)
		return;
	for (j = t->num; j > i; j--)
		t->node[j + 1] = t->node[j];
	t->node[i + 1].key = k;
	t->node[i + 1].child = ap;
	t->num++;
}//Insert

void NewRoot(BMTree &t, int key, BMTree &ap)
{
	BMTree root;
	root = (BMTree)malloc(sizeof(BMT));
	if (root == NULL)
	{
		printf("fail to malloc   NewRoot\n");
		return;
	}
	root->node[0].child = t;
	t = root;
	if (t->node[0].child != NULL)
		t->node[0].child->parent = root;
	t->parent = NULL;
	t->num = 1;
	t->node[1].key = key;
	t->node[1].child = ap;
	if (t->node[1].child != NULL)
		t->node[1].child->parent = t;
}//NewRoot

void Spilt(BMTree &t, BMTree &ap)
{
	int i, s = (m_BTREE + 1) / 2;
	ap = (BMTree)malloc(sizeof(BMT));
	if (ap == NULL)
	{
		printf("fail to malloc   Spilt\n");
		return;
	}
	ap->node[0].child = t->node[s].child;
	for (i = s + 1; i <= m_BTREE; i++)
	{
		ap->node[i - s] = t->node[i];
		if (ap->node[i - s].child != NULL)
			ap->node[i - s].child->parent = ap;
	}
	ap->num = m_BTREE - s;
	ap->parent = t->parent;
	t->num = s - 1;
}//Spilt

void InsertBMTree(BMTree &t, int key, BMTree q, int i)
{
	BMTree temp = NULL;
	bool finished = false;
	int s;
	int key_temp;
	while (q != NULL && !finished)
	{
		Insert(q, key, i, temp);
		if (q->num < m_BTREE)
			finished = true;
		else
		{
			s = (m_BTREE + 1) / 2;
			key_temp = q->node[s].key;
			Spilt(q, temp);
			q = q->parent;
			if (q != NULL)
				i = Search_Insert(q, key_temp);
		}
	}
	if (!finished)		//t是空树，需要新的根结点
		NewRoot(t, key_temp, temp);
}//InsertBMTree


//邻接矩阵
int locate(const Graph g, char ch)
{
	for (int i = 0; i<g.numVertex; i++)
		if (g.vexs[i] == ch)
			return i;
	return -1;
}//locate

void createGraph(Graph &g, int directed)
{
	int i;
	printf("input number of Vertexs:");
	cin >> g.numVertex;
	clean_stdin();
	printf("input number of Edges:");
	cin >> g.numEdge;
	clean_stdin();
	for (i = 0; i < g.numVertex; i++)
	{
		printf("input the No.%d vex:", i + 1);
		g.vexs[i] = getchar();
		clean_stdin();
		while (g.vexs[i] == '\n')
		{
			printf("input the vex again: ");
			g.vexs[i] = getchar();
		}
	}
	for (i = 0; i < g.numVertex; i++)			//初始化所有边
		for (int j = 0; j < g.numVertex; j++)
			g.arc[i][j] = INFINITY_GRAPH;
	for (i = 0; i < g.numEdge; i++)
	{
		int n, m, w;
		char p, q;
		printf("input i and j of (Vi,Vj) and weight\n");
		printf("now input i:");
		while ((p = getchar()) == '\n')
		{
			clean_stdin();
			printf("input i again\n");
		}
		clean_stdin();
		printf("now input j:");
		while ((q = getchar()) == '\n')
		{
			clean_stdin();
			printf("input j again\n");
		}
		clean_stdin();
		printf("now input weight:");
		cin >> w;
		clean_stdin();
		m = locate(g, p);
		n = locate(g, q);
		if (-1 == m || -1 == n)
		{
			printf("i or j not found\n");
			i--;
		}
		else
		{
			g.arc[m][n] = w;
			if(!directed)  g.arc[n][m] = g.arc[m][n];
		}
	}
}//createGraph

void printGraph(const Graph g)
{
	int i, j;
	for (i = 0; i < g.numVertex; i++)
	{
		for (j = 0; j < g.numVertex; j++)
			printf("%d	", g.arc[i][j]);
		printf("\n");
	}
}//printGraph

void find_cycle(const Graph g)
{
	int i;
	for (i = 0; i < g.numVertex; i++)
		visited[i] = 0;
	for (i = 0; i < g.numVertex; i++)
		if (!visited[i])		//未访问的结点
			dfs(g,i);
}//find_cycle

void dfs(const Graph g, int v)		//深度优先遍历找回路
{
	visited[v] = 1;
	for (int i = 0; i < g.numVertex;i++)
		if(g.arc[v][i])
			if (visited[i] == 1)
				Print_dfs(g, i, i);
			else	if (!visited[i])   dfs(g, i);
	visited[v] = 2;
}//dfs

void Print_dfs(const Graph g, int v, int start)
{
	int i;
	for (i = 0; i < g.numVertex;i++)
		if (g.arc[v][i] && visited[i] == 1)
		{
			printf("%c	", g.vexs[v]);
			if (i == start)
				printf("\n");
			else  Print_dfs(g, i, start);
			break;
		}
}//Print_dfs

void prim(const Graph g)			//prim算法（邻接矩阵），求无向图最小生成树
{
	int i, j, k, min;
	int lowcost[MAX_VERTEX_NUM], closest[MAX_VERTEX_NUM], used[MAX_VERTEX_NUM];
	//lowcost：到某个点的最小权值	closest：到某点权值最小的另一个点	used：该结点是否已找出 0：未找出
	memset(closest, 0, sizeof(closest));
	memset(used, 0, sizeof(used));
	memset(parent, -1, sizeof(parent));
	for (i = 0; i < g.numVertex; i++)
		lowcost[i] = g.arc[0][i];
	used[0] = 1;
	for (i = 0; i < g.numVertex - 1; i++)	// 总共numVertex-1条边
	{
		j = 0;
		min = INFINITY_GRAPH;
		for (k = 1; k < g.numVertex;k++)		//找出权值最小的未used的结点
			if (0 == used[k] && lowcost[k] < min)
			{
				min = lowcost[k];
				j = k;
			}
		parent[j] = closest[j];
		used[j] = 1;
		for (k = 0; k < g.numVertex;k++)				//更新lowcost数组
			if (0 == used[k] && g.arc[k][j] < lowcost[k])
			{
				lowcost[k] = g.arc[k][j];
				closest[k] = j;
			}
	}
}//prim

void PrintMST(Graph g)		
{
	prim(g);
	int weight = 0;
	for (int i = 1; i < g.numVertex; i++)
	{
		printf("%c--%c\n", g.vexs[i], g.vexs[parent[i]]);
		weight += g.arc[i][parent[i]];
	}
	printf("weight in total: %d\n", weight);
}//PrintMST

bool topologicalsort(const Graph g)	//拓扑排序
{
	int i, j, count = 0;
	int indegree[MAX_VERTEX_NUM];		//存储入度
	memset(indegree, 0, sizeof(indegree));
	memset(visite, false, sizeof(visite));
	for (i = 0; i < g.numVertex; i++)
		for (j = 0; j < g.numVertex; j++)
			if (g.arc[i][j] != INFINITY_GRAPH)
				indegree[j]++;
	for (i = 0; i < g.numVertex; i++)
	{
		if (0 == indegree[i] && false == visite[i])
		{
			count++;
			visite[i] = true;
			printf("%d	:%c	", count, g.vexs[i]);
			for (j = 0; j < g.numVertex; j++)
				if (g.arc[i][j] != INFINITY_GRAPH)
					indegree[j]--;
		}
	}
	if (count < g.numVertex)
		return false;
	return true;
}//topologicalsort

void Dijkstra(const Graph g, int start, int** result)		//从某个源点到其余各顶点的最短路径
{
	int i, j, min, pos, temp;
	*result = (int*)malloc(sizeof(int)*g.numVertex);
	if (*result == NULL)
	{
		printf("fail to malloc\n");
		return;
	}
	for (i = 0; i < g.numVertex; i++)
	{
		(*result)[i] = g.arc[start][i];
		if (INFINITY_GRAPH != g.arc[start][i])
			parent[i] = start;
		else  parent[i] = -1;
	}
	memset(visite, false, sizeof(visite));
	visite[start] = true;
	(*result)[start] = 0;
	parent[start] = -1;
	for (i = 1; i < g.numVertex; i++)
	{
		min = INFINITY_GRAPH;
		for (j = 0; j < g.numVertex; j++)
		{
			if (visite[j] == false && (*result)[j] < min)
			{
				min = (*result)[j];
				pos = j;
			}			
		}
		if (visite[pos] == true)	//找不到新的（可能是非连通图）
			break;
		visite[pos] = true;
		for (j = 0; j < g.numVertex; j++)
		{
			temp = (g.arc[pos][j] == INFINITY_GRAPH ? INFINITY_GRAPH : (g.arc[pos][j] + (*result)[pos]));
			//防止溢出
			if (false == visite[j] && (*result)[j] > temp)
			{
				(*result)[j] = temp;
				parent[j] = pos;
			}
		}
	}
	for (i = 0; i < g.numVertex; i++)		//输出结果
	{
		printf("\n%c  :", g.vexs[i]);
		pos = parent[i];
		while (pos != -1)
		{
			printf("%c	", g.vexs[pos]);
			pos = parent[pos];
		}
		printf("\nlength: %d", (*result)[i]);
	}
}//Dijkstra

void Floyd(const Graph g,int result[MAX_VERTEX_NUM][MAX_VERTEX_NUM])			//Floyd算法 求任意两顶点最短路径
{
	int i, j, k, temp;
	for (i = 0; i < g.numVertex; i++)
		for (j = 0; j < g.numVertex; j++)
			if (i == j)
				result[i][j] = 0;
			else	result[i][j] = g.arc[i][j];
			for (k = 0; k < g.numVertex; k++)
				for (i = 0; i < g.numVertex; i++)
					for (j = 0; j < g.numVertex; j++)
					{
						if (result[i][k] == INFINITY_GRAPH || result[k][j] == INFINITY_GRAPH)
							temp = INFINITY_GRAPH;
						else temp = result[i][k] + result[k][j];
						if (temp < result[i][j])
							result[i][j] = temp;
					}
			for (i = 0; i < g.numVertex; i++)
			{
				for (j = 0; j < g.numVertex; j++)
					printf("%d\t", result[i][j]);
				printf("\n");
			}
}//Floyd

//邻接表
int locate(const GraphList g, char ch)
{
	for (int i = 0; i < g.numVertex; i++)
		if (g.adjlist[i].data == ch)
			return i;
	printf("%c not found", ch);
	return -1;
}//locate

void createGraph(GraphList &g, int inverse)
{
	int i, w, temp;
	char p, q;
	EdgeNode *e, *f;
	printf("1: directed graph   0: undirected graph\n");
	printf("your choice:");
	cin >> temp;
	while (temp != 1 && temp != 0)
	{
		printf("input again( 1 or 0):");
		cin >> temp;
	}
	clean_stdin();
	printf("input the number of vertexes:");
	cin >> g.numVertex;
	clean_stdin();
	printf("input the number of edges:");
	cin >> g.numEdge;
	clean_stdin();
	for (i = 0; i < g.numVertex; i++)
	{		
		g.adjlist[i].firstedge = NULL;
		printf("input No.%d vertex:", i + 1);
		while ((g.adjlist[i].data = getchar()) == '\n')
		{
			clean_stdin();
			printf("\ninput again:");
		}
		clean_stdin();
	}
	for (i = 0; i < g.numEdge; i++)
	{
		printf("input i and j of (Vi,Vj) and weight\n");
		printf("now input i:");
		while ((p = getchar()) == '\n')
		{
			printf("\ninput i again:");
			fflush(stdin);
		}
		clean_stdin();
		printf("now input j:");
		while ((q = getchar()) == '\n')
		{
			printf("\ninput j again:");
			fflush(stdin);
		}
		clean_stdin();
		printf("now input weight:");
		cin >> w;
		clean_stdin();
		int m, n;
		m = locate(g, p);
		n = locate(g, q);
		if (-1 == m || -1 == n)
		{
			printf("Vi or Vj not found");
			i--;
		}
		else
		{
			if (inverse == 0)		//邻接表
			{
				e = (EdgeNode*)malloc(sizeof(EdgeNode));
				if (NULL == e)
					printf("fail to malloc edgenode");
				else
				{
					e->adjvex = n;
					e->weight = w;
					e->next = g.adjlist[m].firstedge;
					g.adjlist[m].firstedge = e;
				}
				if (!temp)
				{
					f = (EdgeNode*)malloc(sizeof(EdgeNode));
					if (NULL == f)
						printf("fail to malloc edgenode");
					else
					{
						f->adjvex = m;
						f->weight = w;
						f->next = g.adjlist[n].firstedge;
						g.adjlist[n].firstedge = f;
					}
				}
			}	
			else     //逆邻接表
			{
				e = (EdgeNode*)malloc(sizeof(EdgeNode));
				if (NULL == e)
					printf("fail to malloc edgenode");
				else
				{
					e->adjvex = m;
					e->weight = w;
					e->next = g.adjlist[n].firstedge;
					g.adjlist[n].firstedge = e;
				}
				if (!temp)
				{
					f = (EdgeNode*)malloc(sizeof(EdgeNode));
					if (NULL == f)
						printf("fail to malloc edgenode");
					else
					{
						f->adjvex = n;
						f->weight = w;
						f->next = g.adjlist[m].firstedge;
						g.adjlist[m].firstedge = f;
					}
				}
			}			
		}
	}
}//createGraph

void printGraph(const GraphList g)
{
	int i = 0;
	EdgeNode* e = NULL;
	for (i = 0; i < g.numVertex; i++)
	{
		printf("\nVertex: %c	", g.adjlist[i].data);
		e = g.adjlist[i].firstedge;
		while (e != NULL)
		{
			printf("%c(%d)	", g.adjlist[e->adjvex].data, e->weight);
			e = e->next;
		}
	}
}//printGraph

void CountCP(const GraphList g, int choose)		//求无向图连通分量个数
{
	int k = 0, i;
	for (i = 0; i < g.numVertex; i++)
		visite[i] = false;
	if (choose)
	{
		for (i = 0; i < g.numVertex;i++)
				if (visite[i] == false)
				{
					printf("\nNo.%d connected components:\n", ++k);
					dfs(g, i);
				}
		return;
	}	
	for (i = 0; i < g.numVertex; i++)
		if (visite[i] == false)
		{
			printf("\nNo.%d connected components:\n", ++k);
			bfs(g, i);
		}
}//CountSCP

void dfs(const GraphList g, int v)		//深度遍历找连通分量并输出
{
	visite[v] = true;
	printf("%c	", g.adjlist[v].data);
	EdgeNode* p = g.adjlist[v].firstedge;
	while (p != NULL)
	{
		if (visite[p->adjvex] == false)
			dfs(g, p->adjvex);
		p = p->next;
	}
}//dfs

void bfs(const GraphList g, int v)		//广度遍历找连通分量并输出
{
	visite[v] = true;
	printf("%c	", g.adjlist[v].data);
	EdgeNode* p = NULL;
	CQueue Q;
	InitQueue(Q);
	EnQueue(Q, v);
	while (!isEmpty(Q))
	{
		v = DeQueueInt(Q);
		p = g.adjlist[v].firstedge;
		while (p != NULL)
		{
			if (visite[p->adjvex] == false)
			{
				printf("%c	", g.adjlist[p->adjvex].data);
				visite[p->adjvex] = true;
				EnQueue(Q, p->adjvex);
			}
			p = p->next;
		}
	}
}//bfs

void CountSCP(const GraphList g)		//求有向图强连通分量
{
	int i;
	for (i = 0; i < g.numVertex; i++)
	{
		visited[i] = 0;
		inStack[i] = 0;
	}
	Stack s;
	InitStack(s);
	for (i = 0; i < g.numVertex; i++)
		if (0 == visited[i])
			tarjan(g, s, i);
	tarjan_time = 0;
}//CountSCP

void tarjan(const GraphList g, Stack& s, int v)	//tarjan算法
{
	dfn[v] = low[v] = tarjan_time++;
	push(s, v);
	visited[v] = 1;
	inStack[v] = 1;
	EdgeNode* p = g.adjlist[v].firstedge;
	while (p != NULL)
	{
		int i = p->adjvex;
		if (0 == visited[i])
		{
			tarjan(g, s, i);
			low[v] = low[v] < low[i] ? low[v] : low[i];
		}
		else if(inStack[i])
			low[v] = low[v] < dfn[i] ? low[v] : dfn[i];
		p = p->next;
	}
	if (dfn[v] == low[v])
	{
		int vtx;
		printf("\n strong connected component:	");
		do {
			vtx = popint(s);
			inStack[vtx] = 0;
			printf("%c	", g.adjlist[vtx].data);
		} while (vtx != v);
	}
}//tarjan

void kruskal(const GraphList g)			//求无向图最小生成树
{
	int i, m, n, p1, p2;
	int length, index = 0;
	int vends[MAX_VERTEX_NUM];		//保存已有最小生成树中结点在该树的终点
	Edge rets[MAX_VERTEX_NUM];		//保存最小生成树的边
	Edge* edges;

	memset(vends, 0, sizeof(vends));
	get_edges(g,&edges);
	sorted_edges(&edges, g.numEdge);
	for (i = 0; i < g.numEdge; i++)
	{
		p1 = locate(g, edges[i].start);
		p2 = locate(g, edges[i].end);

		m = get_end(vends, p1);
		n = get_end(vends, p2);
		if (m != n)
		{
			vends[m] = n;
			rets[index++] = edges[i];
		}
		if (index == g.numVertex - 1)
			break;
	}
	free(edges);
	length = 0;
	for (i = 0; i < index; i++)
	{
		printf("%c-->%c\n", rets[i].start, rets[i].end);
		length += rets[i].weight;
	}		
	printf("kruskal : %d\n", length);		
}//kruskal

int getWeight(const GraphList g, int start, int end)
{
	if (start == end)
		return 0;
	EdgeNode* node = g.adjlist[start].firstedge;
	while (node != NULL)
	{
		if (node->adjvex == end)
			return node->weight;
		node = node->next;
	}
	return INFINITY_GRAPH;
}//getWeight

void get_edges(const GraphList g,Edge** edges)
{
	int i;
	int index = 0;
	EdgeNode* node = NULL;
	*edges = (Edge*)malloc(g.numEdge*sizeof(Edge));
	for (i = 0; i < g.numVertex; i++)
	{
		node = g.adjlist[i].firstedge;
		while (node != NULL)
		{
			if (node->adjvex > i)
			{
				(*edges)[index].start = g.adjlist[i].data;
				(*edges)[index].end = g.adjlist[node->adjvex].data;
				(*edges)[index++].weight = node->weight;
			}
			node = node->next;
		}
	}
}//get_edges

void sorted_edges(Edge** edges, int elen)
{
	int i, j, k;
	for (i = 0; i < elen; i++)		//选择排序
	{
		k = i;
		for (j = i + 1; j < elen; j++)
		{
			if ((*edges)[k].weight >(*edges)[j].weight)
				k = j;
		}
		Edge temp = (*edges)[i];
		(*edges)[i] = (*edges)[k];
		(*edges)[k] = temp;
	}
}//sorted_edges

int get_end(int vends[], int i)
{
	while (vends[i] != 0)
		i = vends[i];
	return i;
}//get_end

bool topologicalsort(const GraphList g)	
{
	int i, j, count = 0;
	int degree[MAX_VERTEX_NUM];
	memset(degree, 0, sizeof(degree));
	memset(visite, false, sizeof(visite));
	EdgeNode* p = NULL;
	for (i = 0; i < g.numVertex; i++)
	{
		p = g.adjlist[i].firstedge;
		while (p != NULL)
		{
			j = p->adjvex;
			degree[j]++;
			p = p->next;
		}
	}
	Stack s;
	InitStack(s);
	for (i = 0; i < g.numVertex; i++)		//统计入度
		if (degree[i] == 0)
			push(s, i);
	while (!isEmpty(s))
	{
		j = popint(s);							//将该结点的后续结点的度减一
		printf("%d	%c\n", j, g.adjlist[j].data);
		count++;
		for (p = g.adjlist[j].firstedge; p != NULL; p = p->next)
			if (!(--degree[p->adjvex]))
				push(s, p->adjvex);
	}
	if (count < g.numVertex)
		return false;
	return true;
}//topologicalsort

void Dijkstra(const GraphList g, int start, int** result)		//从某个源点到其余各顶点的最短路径
{
	int i, j, k, min, pos, temp;
	*result = (int*)malloc(sizeof(int)*g.numVertex);
	if (*result == NULL)
	{
		printf("fail to malloc\n");
		return;
	}
	memset(visite, false, sizeof(visite));
	memset(parent, -1, sizeof(parent));
	//memset(*result, INFINITY_GRAPH, sizeof(*result));
	for (i = 0; i < g.numVertex; i++)
		(*result)[i] = INFINITY_GRAPH;
	EdgeNode* p = g.adjlist[start].firstedge;
	while (p != NULL)
	{
		i = p->adjvex;
		(*result)[i] = p->weight;
		parent[i] = start;
		p = p->next;
	}
	visite[start] = true;
	(*result)[start] = 0;
	parent[start] = -1;
	for (i = 0; i < g.numVertex; i++)
	{
		min = INFINITY_GRAPH;
		for (j = 0; j < g.numVertex; j++)
		{
			if (false == visite[j] && (*result)[j] < min)
			{
				min = (*result)[j];
				pos = j;
			}
		}
		if (visite[pos] == true)	//找不到新的（可能是非连通图）
			break;
		visite[pos] = true;
		p = g.adjlist[pos].firstedge;
		while (p != NULL)
		{
			k = p->adjvex;
			temp = (*result)[pos] + p->weight;
			if (false == visite[k] && (*result)[k] > temp)
			{
				(*result)[k] = temp;
				parent[k] = pos;
			}
			p = p->next;
		}
	}
	for (i = 0; i < g.numVertex; i++)		//输出结果
	{
		printf("\n%c  :", g.adjlist[i].data);
		pos = parent[i];
		while (pos != -1)
		{
			printf("%c	", g.adjlist[pos].data);
			pos = parent[pos];
		}
		printf("\nlength: %d", (*result)[i]);
	}
}//Dijkstra

void CriticalPath(const GraphList g)		//求关键路径
{
	int indegree[MAX_VERTEX_NUM];
	memset(indegree, 0, sizeof(indegree));
	int i, j, k;
	EdgeNode* p = NULL;
	for (i = 0; i < g.numVertex; i++)		//统计入度
	{
		p = g.adjlist[i].firstedge;
		while (p != NULL)
		{
			j = p->adjvex;
			indegree[j]++;
			p = p->next;
		}
	}
	int ve[MAX_VERTEX_NUM], vl[MAX_VERTEX_NUM];
	//ve: 事件最早开始时间	vl:事件最迟开始时间
	int e[MAX_VERTEX_NUM], l[MAX_VERTEX_NUM];
	//e:活动最早开始时间	l:活动最迟开始时间	
	memset(ve, 0, sizeof(ve));
	int front = -1, rear = -1;
	int queue[MAX_VERTEX_NUM];
	for (i = 0; i < g.numVertex; i++)
		if (0 == indegree[i])
			queue[++rear] = i;
	int count = 0;	//计数器
	while (front != rear)
	{
		i = queue[++front];
		count++;
		p = g.adjlist[i].firstedge;
		while (p != NULL)
		{
			j = p->adjvex;
			if ((ve[i] + p->weight) > ve[j])
				ve[j] = ve[i] + p->weight;
			if (!(--indegree[j]))
				queue[++rear] = j;
			p = p->next;
		}
	}
	if (count < g.numVertex)
	{
		printf("loop exist\n");
		printf("exit\n");
		return;
	}
	for (i = 0; i < g.numVertex; i++)	//对vl赋初值
		vl[i] = ve[g.numVertex - 1];
	for (i = g.numVertex - 2; i >= 0; i--)
	{
		j = queue[i];
		p = g.adjlist[j].firstedge;
		while (p != NULL)
		{
			k = p->adjvex;
			if (vl[k] - p->weight < vl[j])
				vl[j] = vl[k] - p->weight;
			p = p->next;
		}
	}
	printf("\n开始事件\t结束事件\t最早开始时间\t最迟开始时间\n");
	j = 0;
	for (i = 0; i < g.numVertex; i++)		//输出结果
	{
		p = g.adjlist[i].firstedge;
		while (p != NULL)
		{
			k = p->adjvex;
			e[++j] = ve[i];
			l[j] = vl[k] - p->weight;
			printf("%c\t\t%c\t\t%d\t\t%d", g.adjlist[i].data, g.adjlist[k].data, e[j], l[j]);
			if (e[j] == l[j])
				printf("\t关键路径");
			printf("\n");
			p = p->next;
		}
	}
}//CriticalPath

//图的十字链表 
int locate(const OrthList g, char ch)
{
	int i;
	for (i = 0; i<g.vertexNum; i++)
		if (g.vertex[i].data == ch)
			return i;
	return -1;
}//locate

void createGraph(OrthList &g)
{
	int i, j, k, w;
	char vt, vh;
	printf("input the number of vertexes:");
	cin >> (g.vertexNum);
	fflush(stdin);
	printf("input the number of arcs:");
	cin>> (g.arcNum);
	fflush(stdin);
	for (i = 0; i<g.vertexNum; i++)
	{
		printf("input No.%d vertex:", i + 1);
		cin>> (g.vertex[i].data);
		fflush(stdin);
		g.vertex[i].firstin = NULL;
		g.vertex[i].firstout = NULL;
	}
	for (i = 0; i<g.arcNum; i++)
	{
		printf("input i and j of (Vi,Vj) and weight\n");
		printf("now input i:");
		cin >> vt;
		fflush(stdin);
		printf("now input j:");
		cin >> vh;
		fflush(stdin);
		j = locate(g, vt);
		k = locate(g, vh);
		if (j == -1 || k == -1)
		{
			printf("Vi or Vj not found\n");
			i--;
		}
		else
		{
			printf("now input weight:");
			cin >> w;
			fflush(stdin);
			CGNode p = (CGNode)malloc(sizeof(struct cgNode));
			if (p == NULL)
				printf("fail to malloc CGNode\n");
			else
			{
				p->tailvex = j;
				p->headvex = k;
				p->weight = w;
				p->tlink = g.vertex[j].firstout;
				g.vertex[j].firstout = p;
				p->hlink = g.vertex[k].firstin;
				g.vertex[k].firstin = p;
			}
		}
	}
}//createGraph

void printGraph(const OrthList g)
{
	CGNode p = NULL;
	for (int i = 0; i<g.vertexNum; i++)
	{
		p = g.vertex[i].firstout;
		printf("\n%c:	", g.vertex[i].data);
		while (p)
		{
			printf("(%c,%c)	%d", g.vertex[p->tailvex].data, g.vertex[p->headvex].data, p->weight);
			p = p->tlink;
		}
	}
}//printGraph

//邻接多重表
int locate(AMLGraph g, char ch)
{
	for (int i = 0; i<g.vertexNum; i++)
		if (g.AML[i].data == ch)
			return i;
	return -1;
}//locate

void createGraph(AMLGraph &g)
{
	int i, j, k, w;
	char vt, vh;
	EBox* p = NULL;
	printf("input the number of vertexes:");
	cin >> (g.vertexNum);
	fflush(stdin);
	printf("input the number of edges:");
	cin >> (g.edgeNum);
	fflush(stdin);
	for (i = 0; i<g.vertexNum; i++)
	{
		printf("input No.%d vertex:", i + 1);
		cin >> g.AML[i].data;
		g.AML[i].firstedge = NULL;
		fflush(stdin);
	}
	for (i = 0; i<g.edgeNum; i++)
	{
		printf("input i and j of (Vi,Vj) and weight\n");
		printf("now input i:");
		cin >> vt;
		fflush(stdin);
		printf("now input j:");
		cin >> vh;
		fflush(stdin);
		j = locate(g, vt);
		k = locate(g, vh);
		if (j == -1 || k == -1)
		{
			printf("Vi or Vj not found\n");
			i--;
		}
		else
		{
			printf("now input the weight:");
			cin >> w;
			p = (EBox*)malloc(sizeof(EBox));
			if (p == NULL)
			{
				printf("FAIL TO MALLOC EBOX\n");
			}
			else
			{
				p->mark = 0;
				p->ivex = j;
				p->jvex = k;
				p->weight = w;
				p->ilink = g.AML[j].firstedge;
				g.AML[j].firstedge = p;
				p->jlink = g.AML[k].firstedge;
				g.AML[k].firstedge = p;
			}
		}
	}
}//createGraph

void printGraph(const AMLGraph g)
{
	EBox* p = NULL;
	char data;
	for (int i = 0; i<g.vertexNum; i++)
	{
		printf("%c:\n", g.AML[i].data);
		p = g.AML[i].firstedge;
		while (p)
		{
			if (p->ivex == i)	data = g.AML[p->jvex].data;
			else	data = g.AML[p->ivex].data;
			printf("(%c,%c) %d\n", g.AML[i].data, data, p->weight);
			if (p->ivex == i)	p = p->ilink;
			else 	p = p->jlink;
		}
	}
}//printGraph

bool getVex(const AMLGraph g, int v, char& ch)
{
	if (v >= g.vertexNum || v < 0)
		return false;
	ch = g.AML[v].data;
	return true;
}//getVex

bool insertVex(AMLGraph& g, char ch)
{
	if (g.vertexNum == MAX_VERTEX_NUM)
	{
		printf("this graph is full\n");
		return false;
	}
	if (locate(g, ch) >= 0)
	{
		printf("vertex %c is in the graph already\n", ch);
		return false;
	}
	g.AML[g.vertexNum].data = ch;
	g.AML[g.vertexNum++].firstedge = NULL;
	return true;
}//insertVex

bool delEdge(AMLGraph& g, char v, char w)
{
	int i = locate(g, v);
	int j = locate(g, w);
	if (i == j || i < 0 || j < 0)
	{
//		printf("fail to find  or  v is the same as w\n");
		return false;
	}
	EBox *p, *q=NULL;

	p = g.AML[i].firstedge;
	if (p && p->jvex == j)
		g.AML[i].firstedge = p->ilink;
	else if (p && p->ivex == j)
		g.AML[i].firstedge = p->jlink;
	else
	{
		while (p)
		{
			if (p->ivex == i && p->jvex != j)
			{
				q = p;
				p = p->ilink;
			}
			else if (p->jvex == i && p->ivex != j)
			{
				q = p;
				p = p->jlink;
			}
			else	break;
		}
		if (!p)
		{
	//		printf("fail to find the edge\n");
			return false;
		}
		if (p->ivex == i && p->jvex == j)
			if (q->ivex == i)	q->ilink = p->ilink;
			else	q->jlink = p->ilink;
		else if (p->ivex == j && p->jvex == i)
			if (q->ivex == i)	q->ilink = p->jlink;
			else q->jlink = p->jlink;
	}
	//从另一个顶点删除该边 
	p = g.AML[j].firstedge;
	if (p && p->jvex == i)
		g.AML[j].firstedge = p->ilink;
	else if (p && p->ivex == j)
		g.AML[j].firstedge = p->jlink;
	else
	{
		while (p)
		{
			if (p->ivex == j && p->jvex != i)
			{
				q = p;
				p = p->ilink;
			}
			else if (p->jvex == j && p->ivex != i)
			{
				q = p;
				p = p->jlink;
			}
			else	break;
		}
		if (!p)
		{
	//		printf("fail to find the edge\n");
			return false;
		}
		if (p->ivex == i && p->jvex == j)
			if (q->ivex == j)	q->ilink = p->jlink;
			else	q->jlink = p->jlink;
		else if (p->ivex == j && p->jvex == i)
			if (q->ivex == j)	q->ilink = p->ilink;
			else q->jlink = p->ilink;
	}
	free(p);
	g.edgeNum--;
	return true;
}//delEdge

bool delVex(AMLGraph& g, char v)
{
	int i, j;
	char w;
	EBox* p = NULL;
	i = locate(g, v);
	if (i < 0)
	{
		printf("fail to find vextex %c\n", v);
		return false;
	}
	for (j = 0; j<g.vertexNum; j++)
	{
		if (i == j)	continue;
		getVex(g, j, w);
		delEdge(g, v, w);
	}
	for (j = i + 1; j<g.vertexNum; j++)
		g.AML[j - 1] = g.AML[j];
	g.vertexNum--;
	for (j = i; j<g.vertexNum; j++)
	{
		p = g.AML[j].firstedge;
		if (p)
		{
			if (p->ivex == j + 1)
			{
				p->ivex--;
				p = p->ilink;
			}
			else
			{
				p->jvex--;
				p = p->jlink;
			}
		}
	}
	return true;
}//delVex

void destroyGraph(AMLGraph& g)
{
	for (int i = g.vertexNum - 1; i >= 0; i--)
		delVex(g, g.AML[i].data);
}//destroyGraph

int FirstAdjVex(const AMLGraph g, char v)
{
	int i = locate(g, v);
	if (i<0)	return -1;
	if (g.AML[i].firstedge)
		if (g.AML[i].firstedge->ivex == i)
			return g.AML[i].firstedge->jvex;
		else	return g.AML[i].firstedge->ivex;
		return -1;
}//FirstAdjVex 

int NextAdjVex(const AMLGraph g, char v, char w)
{
	int i, j;
	EBox* p = NULL;
	i = locate(g, v);
	j = locate(g, w);
	if (i<0 || j<0)	return -1;
	p = g.AML[i].firstedge;
	while (p)
		if (p->ivex == i && p->jvex != j)
			p = p->ilink;
		else if (p->jvex == i && p->ivex != j)
			p = p->jlink;
		else break;
		if (p && p->ivex == i && p->jvex == j)
		{
			p = p->ilink;
			if (p && p->ivex == i)	return p->jvex;
			else if (p && p->jvex == i)	return p->ivex;
		}
		if (p && p->ivex == j && p->jvex == i)
		{
			p = p->jlink;
			if (p && p->ivex == i)	return p->jvex;
			else if (p && p->jvex == i)	return p->ivex;
		}
		return -1;
}//NextAdjVex 

bool insertEdge(AMLGraph& g, char v, char w)
{
	int i, j, weight;
	EBox* p;
	i = locate(g, v);
	j = locate(g, w);
	if (i < 0 || j < 0)
	{
		printf("fail to find vertex %c or %c\n", v, w);
		return false;
	}
	p = (EBox*)malloc(sizeof(EBox));
	if (p == NULL)
	{
		printf("fail to malloc edge\n");
		return false;
	}
	printf("now input the weight:");
	cin >> weight;
	p->mark = 0;
	p->ivex = i;
	p->jvex = j;
	p->weight = weight;
	p->ilink = g.AML[i].firstedge;
	g.AML[i].firstedge = p;
	p->jlink = g.AML[j].firstedge;
	g.AML[j].firstedge = p;
	g.edgeNum++;
	return true;
}// insertEdge

void MarkUnvisited(AMLGraph& g)
{
	int i;
	EBox* p;
	for (i = 0; i<g.vertexNum; i++)
	{
		p = g.AML[i].firstedge;
		while (p)
		{
			p->mark = 0;
			if (p->ivex == i)	p = p->ilink;
			else p = p->jlink;
		}
	}
}//MarkUnvisited

 
void DFS(const AMLGraph g, int v)
{
	int i;
	EBox* p;
	printf("%c\n", g.AML[v].data);
	visite[v] = true;
	p = g.AML[v].firstedge;
	while (p)
	{
		i = p->ivex == v ? p->jvex : p->ivex;
		if (!visite[i])
			DFS(g, i);
		p = p->ivex == v ? p->ilink : p->jlink;
	}
}//DFS

void DFSTraverse(const AMLGraph g)
{
	int v;
	for (v = 0; v<g.vertexNum; v++)
		visite[v] = false;
	for (v = 0; v<g.vertexNum; v++)
		if (!visite[v])
			DFS(g, v);
	printf("\n");
}//DFSTraverse

void BFSTraverse(const AMLGraph g)
{
	int v, u, w;
	char w1, u1;
	CQueue q;
	for (v = 0; v<g.vertexNum; v++)
		visite[v] = 0;
	InitQueue(q);
	for (v = 0; v < g.vertexNum;v++)
		if (!visite[v])
		{
			visite[v] = true;
			printf("%c\n", g.AML[v].data);
			EnQueue(q, v);
			while (!isEmpty(q))
			{
				u = DeQueueInt(q);
				u1 = g.AML[u].data;
				for (w = FirstAdjVex(g, u1); w >= 0;)
				{
					if (!visite[w])
					{
						visite[w] = true;
						printf("%c\n", g.AML[w].data);
						EnQueue(q, w);
					}
					w1 = g.AML[w].data;
					w = NextAdjVex(g, u1, w1);
				}
			}
		}
	printf("\n");
}//BFSTraverse

void InsertSort(int a[],int n)
{
	int i, j;
	for (i = 2; i <= n; i++)
	{
		a[0] = a[i];
		j = i - 1;
		while (a[0] < a[j])
		{
			a[j + 1] = a[j];
			j = j - 1;
		}
		a[j + 1] = a[0];
	}
}//InsertSort

void BInsertSort(int a[], int n)
{
	int i,j, low, high, mid;
	for (i = 2; i <= n; i++)
	{
		a[0] = a[i];
		low = 1;
		high = i - 1;
		while (low <= high)
		{
			mid = (low + high) / 2;
			if (a[0] < a[mid])
				high = mid - 1;
			else low = mid + 1;
		}
		for (j = i - 1; j >= high + 1; --j)
			a[j + 1] = a[j];
		a[high + 1] = a[0];
	}
}//BInsertSort

void BubbleSort(int a[], int n)
{
	int change = 1;
	int low = 0;
	int high = n - 1;
	int i;
	while (low <= high && change)
	{
		change = 0;
		for(i = low; i < high;i++)
			if (a[i]>a[i + 1])
			{
				a[i] ^= a[i + 1];
				a[i + 1] ^= a[i];
				a[i] ^= a[i + 1];
			}
		high--;
		for (i = high; i > low;i--)
			if (a[i] < a[i - 1])
			{
				a[i] ^= a[i - 1];
				a[i - 1] ^= a[i];
				a[i] ^= a[i - 1];
				change = 1;
			}
		low--;
	}
}//BubbleSort

void HillSort(int a[], int n)
{
	int k = n / 2;		//增量
	int i, j, t;
	while (k > 0)
	{
		for (i = k; i < n; i++)
		{
			t = a[i];
			j = i - k;
			while (j >= 0 && t < a[j])
			{
				a[j + k] = a[j];
				j = j - k;
			}
			a[j + k] = t;
		}
		k /= 2;
	}
}//HillSort

void SelectSort(int a[], int n)
{
	int i, j, min, pos;
	for (i = 0; i < n; i++)
	{
		min = a[i];
		for (j = i + 1; j < n;j++)
			if (a[j] < min)
			{
				min = a[j];
				pos = j;
			}
		if (min != a[i])
		{
			a[i] ^= a[pos];
			a[pos] ^= a[i];
			a[i] ^= a[pos];
		}		
	}
}//SelectSort

void QuickSort(int a[], int low, int high)
{
	if (low < high)
	{
		int pivotloc = Partition_QSort(a, low, high);
		QuickSort(a, low, pivotloc - 1);
		QuickSort(a, pivotloc + 1, high);
	}
}//QuickSort

int Partition_QSort(int a[], int low, int high)
{
	int key = a[low];
	while (low < high)
	{
		while (low < high && a[high] >= key)
			high--;
		a[low] = a[high];
		while (low < high && a[low] <= key)
			low++;
		a[high] = a[low];
	}
	a[low] = key;
	return low;
}//Partition_QSort

void HeapSort(Heap h)
{
	int i;
	for (i = h->length / 2; i >= 1; i--)
		HeapSift(h, i, h->length);
	for (i = h->length; i >= 2; i--)
	{
		h->data[0] = h->data[i];		//交换i和1
		h->data[i] = h->data[1];
		h->data[1] = h->data[0];
		HeapSift(h, 1, i - 1);			//此时堆顶聚集最小的元素，然后对原来的data[i]排序
	}
}//HeapSort

void HeapSift(Heap h, int k, int m)
{
	int i, j, finished;
	i = k;
	j = 2 * i;
	h->data[0] = h->data[k];
	finished = 0;
	while (j <= m && !finished)
	{
		if (j < m && h->data[j + 1] < h->data[j])		//找两棵子树中 值小的那棵
			j++;
		if (h->data[0] <= h->data[j])				//若要改为大根堆，修改两个if
			finished = 1;
		else
		{
			h->data[i] = h->data[j];
			i = j;
			j = j * 2;
		}
	}
	h->data[i] = h->data[0];
}//HeapSift

void MergeSort(int a[], int n)
{
	int len = 1;
	int s1, s2, t1, t2, q;
	while (len < n)
	{
		s1 = 1;
		q = 1;
		while (s1 + len <= n)
		{
			s2 = s1 + len;
			t1 = s2 - 1;
			t2 = s2 + len - 1;
			if (t2 > n)
				t2 = n;
			Merge(a, s1, t1, s2, t2, q);
			s1 = q;
		}		
		len *=2;
	}
}//MergeSort

void Merge(int r[], int s1, int t1, int s2, int t2, int &q)
{
	int* t = NULL;
	t = (int*)malloc(sizeof(int)*(s2 - s1 + 1));
	if (t == NULL)
	{
		printf("malloc failure\n");
		return;
	}
	int p = 0;
	while (s1 <= t1 && s2 <= t2)
		t[p++] = (r[s1] <= r[s2]) ? r[s1++] : r[s2++];
	while (s1 <= t1)
		t[p++] = r[s1++];
	while (s2 <= t2)
		t[p++] = r[s2++];
	for (int i = 0; i <= p - 1; i++, q++)
		r[q] = t[i];
}//Merge

void RadixSort(int *a, int n, int d)
{
	int i;
	int k[RADIX] = { 0 };
	int *temp = NULL;
	int *b = NULL;
	temp = (int*)malloc(n*sizeof(int));
	b = (int*)malloc(n*sizeof(int));
	if (temp == NULL || b == NULL)
	{
		printf("malloc failure\n");
		return;
	}
	for (int j = 1; j <= d; j++)
	{
		for (i = 0; i < n; i++)
			b[i] = get_value(a[i], j);
		memset(temp, 0, sizeof(temp));
		memset(k, 0, sizeof(k));
		for (i = 0; i < n; i++)
			k[b[i]]++;					//记录与数组下标相等的数值的个数 
		for (i = 1; i < RADIX; i++)		
			k[i] += k[i - 1];			//储存自己数组下标数值在目标数组对应的位置 
		for (i = n - 1; i >= 0; i--)
			temp[--k[b[i]]] = a[i];		//将原数组按大小顺序储存到另一个数组  
		for (i = 0; i < n; i++)
			a[i] = temp[i];
	}
}//RadixSort

int get_value(int n, int d)
{
	int b = n;
	d--;
	while (d-- && b != 0)
		b /= RADIX;
	return	b%RADIX;
}//get_value

#if 0
Matrix FromSymmetric(const int *num, int len)
{
	static Matrix result(len, len);
	int i, j, k;
	for (i = 0; i < len; i++)
	{
		for (j = 0; j < len; j++)
		{
			if (i < j)
				k = j*(j + 1) / 2 + i;
			else
				k = i*(i + 1) / 2 + j;
			result.init(i, j, num[k]);
		}
	}
	return result;
}

Matrix FromTriangle(const int* num, int len, int up)
{
	static Matrix result(len, len);
	int i, j, k;
	k = len*(len + 1) / 2 ;
	if(1 == up)
	{
		for (i = 1; i < len; i++)
			for (j = 0; j < i;j++)
				result.init(i, j, *(num + k));
		for (i = 0; i < len; i++)
			for (j = i; j < len; j++)
			{
				k = i*(2 * len - i + 1) / 2 + j - i;
				result.init(i, j, *(num + k));
			}
	}
	else
	{
		for (i = 0; i < len-1; i++)
			for (j = i+1; j < len; j++)
				result.init(i, j, *(num + k));
		for (i = 0; i < len; i++)
			for (j = 0; j <= i; j++)
			{
				k = i*(i + 1) / 2 + j;
				result.init(i, j, *(num + k));
			}
	}
	return result;
}
#endif