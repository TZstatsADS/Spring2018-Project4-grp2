#include<iostream>
#include<cstring>
#include<cstdlib>
#include<algorithm>
#include<ctime>
#include<cmath> 
using namespace std;


//parameter you can change

const int C = 9;
const double ErrorFix = 1.0;
const int MaxIter =20;

//parameter you can't change

const int MaxTestNum = 199876;
const int MaxUser = 74418;
const int MaxMovie = 1648;
const int MaxN = 5055;
const int MaxM = 1648;
const int MaxRate = 6;
const int MaxEM = 50;
int UserToN[MaxUser + 1];
int MovieToM[MaxMovie + 1];
int RateTrain[MaxN + 1][MaxM + 1];
double miu[MaxEM][C + 1];
double gamma[MaxEM][C + 1][MaxM + 1][MaxRate + 1];
double pi[MaxEM][MaxN + 1][C + 1];
int clock_start;
int Score[MaxTestNum + 1];
double Predict[MaxTestNum + 1];
double PredictMAE[MaxIter + 1];

void init()
{
	srand(time(0));
	int x,y,z,w;
	
	freopen("UserList.txt", "r", stdin);
	for (int i = 0; i < MaxUser; i ++)
	{
		scanf("%d%d%d", &x, &y, &z);
		UserToN[y] = z;
	}
	
	memset(RateTrain, -1, sizeof(RateTrain));
	freopen("data_train.txt", "r", stdin);
	int Column = 809526; 
	for (int i = 0; i < Column; i ++)
	{
		scanf("%d%d%d%d", &x, &y, &z, &w);
		RateTrain[UserToN[z]][y] = w;
	}
	
	//printf("initial done. ");
	//int clock_now = clock();
	//printf("%d ms\n", clock_now - clock_start);
	
}

void printMiu(int cnt)
{
	printf("Miu: ");
	for (int i = 1; i <= C; i ++)
		printf("%.5f ", miu[cnt][i]);
	printf("\n");
}

void EM1()
{
	//step 1: initial
	double sum = 0;
	
	memset(Score, 0, sizeof(Score));
	memset(Predict, 0, sizeof(Predict));
	for (int i = 1; i <= C; i ++)
	{
		miu[0][i] = double(rand() + ErrorFix);
		sum +=miu[0][i];
	}
	for (int i = 1; i <= C; i ++)
		miu[0][i] /= sum;
	for (int i = 1; i <= C; i ++)
		for (int j = 1; j <= MaxM; j ++)
		{
			sum = 0;
			for (int k = 1; k <= MaxRate; k ++)
			{
				gamma[0][i][j][k] = double(rand() + ErrorFix);
				sum += gamma[0][i][j][k];
			}
			for (int k = 1; k <= MaxRate; k ++)
				gamma[0][i][j][k] /= sum;
		}
	//printf("EM Step1 done! ");
	//int clock_now = clock();
	//printf("%d ms ", clock_now - clock_start);
	//printMiu(0);
}
double GetLogPhi(int Class, int UserID, double coef, int cnt)// -> log(coef) + sigma(log(Di,c))
{
	double tmp = log(coef);
	for (int j = 1; j <= MaxM; j ++)
		if (RateTrain[UserID][j] > -1)
			tmp += log(gamma[cnt][Class][j][RateTrain[UserID][j]]);
	//if (rand()%10 == 0)
	//	cout<<tmp<<" ";
	return tmp;
}

void EM2(int cnt)
{
	//step 2: E
	for (int i = 1; i <= MaxN; i ++)
	{
		double logPhi[C + 1];
		double logPhiMin = 0;
		double logPhiSum = 0;
		for (int j = 1; j <= C; j ++)
		{
			logPhi[j] = GetLogPhi(j,i, miu[cnt][j], cnt);
			if (logPhiMin > logPhi[j])
				logPhiMin = logPhi[j];
		}
		for (int j = 1; j <= C; j ++)
		{
			logPhi[j] -= logPhiMin;
			if (logPhi[j] > 500.0) // avoid INF
				logPhi[j] = 500.0;
			logPhi[j] = exp(logPhi[j]);
			logPhiSum += logPhi[j];			
		}
		for (int j = 1; j <= C; j ++)
			pi[cnt][i][j] = logPhi[j] / logPhiSum;		
	}
	
	//printf("EM Step2 done!! ");
	//int clock_now = clock();
	//printf("%d ms\n", clock_now - clock_start);
	
}

void EM3(int cnt)
{
	for (int j = 1; j <= C; j ++)
		for (int i = 1; i <= MaxN; i ++)
			miu[cnt + 1][j] += pi[cnt][i][j] / double(MaxN);
	
	for (int j = 1; j <= MaxM; j ++)
		for (int k = 1; k <= MaxRate; k ++)
			for (int c = 1; c <= C; c ++)
			{
				double Num = 0;
				double Den = 0;
				for (int i = 1; i <= MaxN; i ++)
				{
					if (RateTrain[i][j] == k)
						Num += pi[cnt][i][c];
					if (RateTrain[i][j] > -1)
						Den += pi[cnt][i][c];
				}
				gamma[cnt + 1][c][j][k] = Num / Den;
			}
			
	//printf("EM Step3 done!!! ");
	//int clock_now = clock();
	//printf("%d ms ", clock_now - clock_start);
	//printMiu(cnt + 1);
}

bool EM4(int cnt)
{
	double Delta[C + 1];
	double sum = 0;
	for (int i = 1; i <= C; i++)
		Delta[i] = miu[cnt + 1][i] - miu[cnt][i];
	for (int i = 1; i <= C; i++)
		sum += Delta[i] * Delta[i];
		
	//printf("EM Step4 done!!!! \n");
	if (sum < 1e-4)
	//if (cnt >= 3)
		return true;
	return false;	
}

int EM()
{
	EM1();
	int EMcnt = 0;
	while (1)
	{
		EM2(EMcnt);
		EM3(EMcnt);
		if (EM4(EMcnt))
			return EMcnt;
		EMcnt ++;
		printf("No break. EM's Cnt: %d\n", EMcnt);
	}
}

double PredictScore(int cnt)
{
	int x,y,z,w;
	double MAE = 0;
	freopen("data_test.txt", "r", stdin);
	for (int now = 0; now < MaxTestNum; now ++)
	{
		scanf("%d%d%d%d", &x, &y, &z, &w);
		int i = UserToN[z];
		int b = y;
		Score[now] = w;
		
		double num[C + 1];
		double den[C + 1];
		double denMin = 0;
		double Sumnum = 0;
		double Sumden = 0;
		for (int c = 1; c <= C; c++)
		{
			den[c] = log(miu[cnt][c]);
			for (int j = 1; j <= MaxM; j ++)
				if (RateTrain[i][j] > -1)
					den[c] += log(gamma[cnt][c][j][RateTrain[i][j]]);
		}
		for (int c = 1; c <= C; c++)
			denMin = min(denMin, den[c]);
		for (int c = 1; c <= C; c++)
		{
			den[c] -= denMin;
			if (den[c] >500)
				den[c] = 500;
		}
			
		for (int c = 1; c <= C; c++)
			Sumden += exp(den[c]);
		
		for (int k = 1; k <= MaxRate; k ++)
		{
			Sumnum = 0;
			for (int c = 1; c <= C; c++)
				num[c] = den[c] + log(gamma[cnt][c][b][k]);
			for (int c = 1; c <= C; c++)
				Sumnum += exp(num[c]);
			Predict[now] += double(k) *Sumnum / Sumden;
		}
		MAE += abs(double(Score[now])- Predict[now]);
		if (now%10000 == 0)
			printf("Predict-> Score: %d, Pred: %.3f, MAE=%.4f\n", Score[now], Predict[now], MAE / (now+1));
		
	}
	return MAE / MaxTestNum;

}

int main()
{
	init();
	
	for (int iter = 0; iter <= MaxIter; iter ++)
	{
		clock_start = clock();
		//cout<<(2e-100 / 5e133)<<endl;
		PredictMAE[iter] = PredictScore(EM() + 1);
		int clock_now = clock();
		printf("%d ms, MAE= ", clock_now - clock_start);
		printf("%.4f\n", PredictMAE[iter]);
	} 
	freopen("Output.txt","w",stdout);
	for (int iter = 0; iter < MaxIter; iter ++)
		printf("%.4f\n", PredictMAE[iter]);
	return 0;
}
