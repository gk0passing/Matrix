#include<iostream>
#include<algorithm>
#include<cmath>
#define EPS 1e-10 
using namespace std;
class Matrix
{
	private:
		int row;
		int col;
		double **p;
	public:
		void init();
		Matrix();//无参构造函数 
		Matrix(int ,int );//两个参数构造函数 
		Matrix(int,int,double);//三个参数构造函数 
		~Matrix();//析构函数 
		Matrix& operator=(const Matrix&);//‘=’运算符重载 ，矩阵赋值 
		Matrix& operator=(double **);// ‘=’运算符重载，指针赋值 
		Matrix& operator+=(const Matrix&);//‘+=’运算符重载 
		Matrix& operator-=(const Matrix&);//‘-=’运算符重载 
		Matrix& operator*=(const Matrix&);//‘*=’运算符重载 
		Matrix operator*(const Matrix &);//‘*’运算符重载 
		double& operator()(int,int);//‘（）’运算符重载 
		static Matrix solve(const Matrix&, Matrix&);//Ax=b求解 
		void show();//输出函数 
		static double det(double**,int);//行列式的值求解 
		int rank1();//矩阵的秩 
		static Matrix eye(int );//生成单位矩阵 
		void swaprows(int a,int b);//行变换函数 
		int rows();//显示行 
		int cols();//显示列 
		static double **ch(Matrix&);//转化成指针类型 
		static Matrix T(const Matrix &);//矩阵的转置 
		friend istream& operator>>(istream&,Matrix &);//输入运算符重载 
		void triu();//转化为上三角矩阵 
		void inv(Matrix&);//矩阵求逆
		static void QR(Matrix*, Matrix*, Matrix*);//QR分解 
		static Matrix evalue(Matrix*, Matrix*, Matrix*);//求矩阵特征值 
		static void Eig(Matrix *A, Matrix *eigenVector, Matrix *eigenValue);//求矩阵特征向量 
		static Matrix q_pow(Matrix& A,int b);
};
int Matrix::rows()//0显示行
{
	return row;
}
int Matrix::cols()//1显示列
{
	return col;
}
void Matrix::swaprows(int a,int b)//2行变换函数
{
	double *t=p[a-1];
	p[a-1]=p[b-1];
	p[b-1]=t;
}
void Matrix::show()//3输出矩阵
{
	if(p){
//	cout<<" row "<<row<<" col "<<col<<endl;
	for(int i=0; i<row; i++)
		{
			for(int j=0; j<col; j++)
				{
					if(fabs(p[i][j])<1e-10)//精度设置 
					p[i][j]=0;
					cout<<p[i][j]<<' ';
				//	if()
				}
			
			cout<<endl;
		}
	cout<<endl;
	cout<<"矩阵输出完成！"<<endl; 
}
}
void Matrix::init()//4初始化函数
{
	p=new double* [row];
	for(int i=0; i<row; i++)
		p[i]=new double[col];
}
Matrix::Matrix()//5无参构造函数 
{
	col=0;
	row=0;
	p=NULL;
}
Matrix::Matrix(int rows, int cols)//6两个参数构造函数
{
	row = rows;
	col = cols;
	init();
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			p[i][j] = 0;
}
Matrix::Matrix(int rows ,int cols ,double value)//7三个参数构造函数
{
	row=rows;
	col=cols;
	init();
	for(int i=0; i<row; i++)
		for(int j=0; j<col; j++)
			p[i][j]=value;
}
Matrix::~Matrix()//8析构函数
{
	for (int i = 0; i < row; ++i)
		delete[] p[i];
	delete[] p;
}
Matrix& Matrix::operator=(const Matrix &M)//9重载运算符=
{
	if(this==&M)
		return *this;
	if(row!=M.row||col!=M.col)
		{
			for(int i=0; i<row; i++)
				delete[] p[i];
			delete[] p;
		}
	row=M.row;
	col=M.col;
	init();
	for(int i=0; i<row; i++)
		for(int j=0; j<col; j++)
			p[i][j]=M.p[i][j];
	return *this;
}
Matrix& Matrix::operator=(double **t)// ‘=’运算符重载，指针赋值 
{
	for(int i=0;i<row;i++)
	for(int j=0;j<col;j++)
	{
		p[i][j]=t[i][j];
	}
	return *this;
}
Matrix& Matrix::operator+=(const Matrix &M)//10重载运算符+=
{
	for(int i=0; i<row; i++)
		for(int j=0; j<col; j++)
			p[i][j]+=M.p[i][j];
	return *this;
}
Matrix& Matrix::operator-=(const Matrix &M)//11重载运算符-=
{
	for(int i=0; i<row; i++)
		for(int j=0; j<col; j++)
			p[i][j]-=M.p[i][j];
	return *this;
}
Matrix& Matrix::operator*=(const Matrix &M)//12重载运算符*=
{
	Matrix t(row,M.col);
	for(int i=0; i<t.row; i++)
		for(int j=0; j<t.col; j++)
			for(int k=0; k<col; k++)
				t.p[i][j]=t.p[i][j]+p[i][k]*M.p[k][j];
	col=M.col;
	init();
	*this=t;
	return *this;
}
Matrix Matrix::operator*(const Matrix &M)//13重载运算符*
{
	Matrix a(row,M.col,0);
	for(int i=0; i<a.row; i++)
		for(int j=0; j<a.col; j++)
			for(int k=0; k<col; k++)
				a.p[i][j]=a.p[i][j]+p[i][k]*M.p[k][j];
	return a;
}
Matrix Matrix::eye(int n)//14产生单位矩阵
{
	Matrix A(n,n,0);
	for(int i=0; i<n; i++)
		A.p[i][i]=1;
	return A;
}
Matrix Matrix::solve(const Matrix &A, Matrix &b)//15解Ax=b
{
	for(int i=0; i<A.row; i++)
		if(A.p[i][i]==0)
			cout<<"error"<<endl;
	double *temp=new double[A.row];
	Matrix x(b.row,1);
	for(int k=0; k<A.row-1; k++)
		{
			for(int i=k+1; i<A.row; i++)
				temp[i]=A.p[i][k]/A.p[k][k];
			for(int i=k+1; i<A.row; i++)
				{
					for(int j=0; j<A.row; j++)
						A.p[i][j]=A.p[i][j]-temp[i]*A.p[k][j];
					b.p[i][0]=b.p[i][0]-temp[i]*b.p[k][0];
				}
		}
	delete temp;
	x.p[x.row-1][0]=b.p[b.row-1][0]/A.p[A.row-1][A.row-1];
	for(int i=A.row-2; i>=0; i--)
		{
			double sum=0;
			for(int j=i+1; j<A.row; j++)
				sum=sum+A.p[i][j]*x.p[j][0];
			x.p[i][0]=(b.p[i][0]-sum)/A.p[i][i];
		}
	return x;
}
Matrix Matrix::T(const Matrix &M)//16矩阵转置
{
	Matrix mt(M.col,M.row);
	for(int i=0; i<mt.row; i++)
		for(int j=0; j<mt.col; j++)
			mt.p[i][j]=M.p[j][i];
	return mt;
}
istream& operator>>(istream &in,Matrix &M)//17矩阵输入
{
	for(int i=0; i<M.row; i++)
		for(int j=0; j<M.col; j++)
			in>>M.p[i][j];
			cout<<"矩阵输入完成！"<<endl;
	return in;
}
double& Matrix::operator()(int r,int c)//18重载运算符（）
{
	return(*(*p+r*col+c));
}
double** Matrix::ch(Matrix &M)//19转化
{
	double **p;
	p=new double* [M.row];
	for(int i=0; i<M.row; i++)
		p[i]=new double[M.col];
	for(int i=0; i<M.row; i++)
		for(int j=0; j<M.col; j++)
			p[i][j]=M.p[i][j];
	return p;
}
double Matrix::det(double **p,int n)//20求行列式的值
{
	if(n==1)
		{
			return p[0][0];
		}
	double ans = 0;
	double **temp;
	temp=new double* [n];
	for(int i=0; i<n; i++)
		temp[i]=new double[n];
	int i,j,k;
	for(i=0; i<n; i++)
		{
			for(j=0; j<n-1; j++)
				{
					for(k=0; k<n-1; k++)
						{
							temp[j][k] = p[j+1][(k>=i)?k+1:k];

						}
				}
			double t = det(temp,n-1);
			if(i%2==0)
				{
					ans += p[0][i]*t;
				}
			else
				{
					ans -=  p[0][i]*t;
				}
		}
	return ans;
}
void Matrix::triu()//21化为上三角矩阵 
{
	int times = row < col ? row : col;
	for (int i = 0; i < times-1; i++)
	{
		
		int k;
		for (k = i; k < row; k++)
		{
			if (fabs(p[k][i]) > 1e-10) 
				break;
		}
		if (k <row)
		{
			for (int j = i; j < col; j++)
			{
				double temp;
				temp=p[i][j];
				p[i][j]=p[k][j];
				p[k][j]=temp;
			}
			double a;
			for (int j = i+1; j <row; j++)
			{
				a = -p[j][i] / p[i][i];
				for (k = i; k<col; k++)
				{
					p[j][k] += a * p[i][k];
				}
			}
		}
	}
}
int Matrix::rank1()//22矩阵求秩 
{
	Matrix a(row, col);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			a.p[i][j] = p[i][j];
		}
	}
	a.triu(); 
	int amount = 0;  
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j <col; j++)
		{
			if (fabs(a.p[i][j]) > 1e-10)
			{
				amount++; break;
			}
		}
	}
	return amount;
}
void Matrix::inv(Matrix &b)//23矩阵求逆 
{
	
	int i, j, k;
	float max, t;
	Matrix temp(row,col);          //临时矩阵
	//将A矩阵存放在临时矩阵t[n][n]中
	temp=p;
	if(temp.rank1()!=row)
	{
		cout << "There is no inverse matrix!";
		return;
	}
	
	//初始化B矩阵为单位阵
	Matrix B=eye(col);
	for (i = 0; i < col; i++)
		{
			//寻找主元
			max = temp.p[i][i];
			k = i;
			for (j = i+1; j < row; j++)
				{
					if (fabs(temp.p[j][i]) > fabs(max))
						{
							max = temp.p[j][i];
							k = j;
						}
				}
			//如果主元所在行不是第i行，进行行交换
			if (k != i)
				{
					for (j = 0; j < row; j++)
						{
							t = temp.p[i][j];
							temp.p[i][j] = temp.p[k][j];
							temp.p[k][j] = t;
							//B伴随交换
							t = B.p[i][j];
							B.p[i][j] = B.p[k][j];
							B.p[k][j] = t;
						}
				}
			
			//判断主元是否为0, 若是, 则矩阵A不是满秩矩阵,不存在逆矩阵
//			if (temp.p[i][i] == 0)
//				{
//					cout << "There is no inverse matrix!";
//					return;
//				}
			//消去A的第i列除去i行以外的各行元素
			t = temp.p[i][i];
			for (j = 0; j <row; j++)
				{
					temp.p[i][j] = temp.p[i][j] / t;        //主对角线上的元素变为1
					B.p[i][j] = B.p[i][j] / t;        //伴随计算
				}
			for (j = 0; j <row; j++)        //第0行.第n行
				{
					if (j != i)                //不是第i行
						{
							t = temp.p[j][i];
							for (k = 0; k <row; k++)        //第j行元素 - i行元素*j列i行元素
								{
									temp.p[j][k] = temp.p[j][k] - temp.p[i][k]*t;
									B.p[j][k] = B.p[j][k] - B.p[i][k]*t;
								}
						}
				}
				//	cout<<1<<endl;
		}
	//	B.show();
	b=B;
}
Matrix Matrix::evalue(Matrix *A,Matrix *Q,Matrix *R)//24矩阵求特征值 
{
	Matrix temp(A->row,1,0);
	for (int k=0; k<5;k++)
		{
			QR(A,Q,R);
			*A=*R*(*Q);
		}
	for(int k=0;k<A->col;++k)
	{
		temp.p[k][0]=A->p[k][k];
	}
	return temp;
}
void Matrix::QR(Matrix *A, Matrix *Q, Matrix *R)//25矩阵QR分解 
{
	Matrix col_A=Matrix(A->row,1,0);//用来存A的每一列
	Matrix col_Q=Matrix(A->row,1,0);//用来存Q的每一列
	*Q=Matrix(A->row,A->col,0);
	*R=Matrix(A->row,A->col,0);
	//施密特正交化
	for (int j = 0; j < A->col; j++)
		{
			for (int i = 0; i < A->col; i++)//把A的第j列存入col_A中
				{
					col_A.p[i][0] = A->p[i][j];
					col_Q.p[i][0] = A->p[i][j];
				}
			for (int k = 0; k < j; k++)//计算第j列以前
				{
					R->p[k][j] = 0;
					for (int i1 = 0; i1 < col_A.row; i1++)
						{
							//R=Q'A(Q'即Q的转置) 即Q的第k列和A的第j列做内积
							R->p[k][j] += col_A.p[i1][0] * Q->p[i1][k];//Q的第k列
						}
					for (int i2 = 0; i2 < A->col; i2++)
						{
							col_Q.p[i2][0] -= R->p[k][j] * Q->p[i2][k];
						}
				}
			double ans=0.0;
			for(int i=0; i<col_Q.row*col_Q.col; i++)
				{
					ans+=col_Q.p[i][0]*col_Q.p[i][0];
				}
			ans=(double)sqrt(ans);
			R->p[j][j] = ans;
			for (int i3 = 0; i3 < Q->col; i3++)
				{
					//单位化Q
					Q->p[i3][j] = col_Q.p[i3][0] / ans;
				}
		}
	
} 
void Matrix::Eig(Matrix *A, Matrix *eigenVector, Matrix *eigenValue)//26矩阵求特征向量 
{
	int num = A->col;
	double eValue;
	Matrix temp(A->row,A->col);
//	InitMatrix(&temp, A->row, A->col);
	//CopyMatrix(A, &temp);
	for (int count = 0; count < num; count++)
	{
		eValue = eigenValue->p[count][0];//当前的特征值
		temp=*A;
	//	CopyMatrix(A, &temp);//这个每次都要重新复制，因为后面会破坏原矩阵(刚开始没注意到这个找bug找了好久。。)
		for (int i = 0; i < temp.row; i++)
		{
			temp.p[i][i] -= eValue;
		}
		//将temp化为阶梯型矩阵(归一性)对角线值为一
		for (int i = 0; i < temp.row - 1; i++)
		{
			double coe = temp.p[i][i];
			for (int j = i; j < temp.col; j++)
			{
				temp.p[i][j] /= coe;//让对角线值为一
			}
			for (int i1 = i + 1; i1 < temp.row; i1++)
			{
				coe = temp.p[i1][i];
				for (int j1 = i; j1 < temp.col; j1++)
				{
					temp.p[i1][j1] -= coe * temp.p[i][j1];
				}
			}
		}
		//让最后一行为1
		double sum1 = eigenVector->p[eigenVector->row - 1][count] = 1;
		for (int i2 = temp.row - 2; i2 >= 0; i2--)
		{
			double sum2 = 0;
			for (int j2 = i2 + 1; j2 < temp.col; j2++)
			{
				sum2 += temp.p[i2][j2] * eigenVector->p[j2][count];
			}
			sum2 = -sum2 / temp.p[i2][i2];
			sum1 += sum2 * sum2;
			eigenVector->p[i2][count] = sum2;
		}
		sum1 = sqrt(sum1);//当前列的模
		for (int i = 0; i < eigenVector->row; i++)
		{
			//单位化
			eigenVector->p[i][count] /= sum1;
		}
	}
}
Matrix Matrix::q_pow(Matrix& A,int b)//矩阵快速幂 
{  
    Matrix C=eye(A.row);
    while(b){  
        if(b&1) C = C*A;  
        A = A*A;  
        b >>= 1;  
    }  
    return C;  
} 
int main() //主函数 
{
	//Ax=b
//	 Matrix a(3,3);
//	 cin>>a;
//	 Matrix b(3,1);
//	 cin>>b;
//	 Matrix c=Matrix::solve(a,b);
//	 c.show();

	 //矩阵快速幂 
	 //1 2 3
	 //4 5 6
	 //7 8 9 
//	 Matrix a(3,3);
//	 cin>>a;
//	 Matrix C=Matrix::q_pow(a,10);
//	 C.show();

	 //矩阵求逆 
//		3 -1 0
//		-2 1 1
//		2 -1 4
//		不可行 
//	    1 7 6
//      2 3 7
//     2 14 12
//	 Matrix a(3,3);
//	 cin>>a;
//	 Matrix b;
//	 a.inv(b);
//	 b.show();

	//矩阵特征值和特征向量
//	Matrix a(3,3);
//	cin>>a;
//	Matrix A=a;
//	Matrix Q;
//	Matrix R;
//	Matrix t=Matrix::evalue(&A,&Q,&R);
//	t.show();
//	Matrix::Eig(&a,&Q,&t);
//	Q.show();
}







