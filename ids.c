#if 1
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#define N 10
#define P 1
#define MAX_LINE_LENGTH 5000
#define MAX_INVAILD_KICK 5

void load_data(int run_num,int itr,int dr,int dr_type,int oprtct,char *fn,int mui);
void init_population_solution();
void init_solution(int x[][N],int i);
void sms_iterated_local_search(int itr,int dr,int dr_type,int oprtct,int mui);
void two_swap();
void one_move();//foreInsertion
void operation_tactics(int oprtct);
int dynamic_programming_ts();
int dynamic_programming_om();
int dynamic_programming_om_ts();
int dynamic_programming(int oprtct);
void backtracking();//find the corresponding sequence
int perturbation(int type,int f);

int s[N][3];//问题描述0;processing time,1:weights;2:due dates
int x[P][N];//种群,p为种群个数
int sm[N],cm[N],tm[N],sb[N];//初始解序列
int csm,csm_pre,csb;//cs0:初始解的代价;csm:邻域解中最优解的代价
int ts[N][N];//two swap产生的惩罚差值
int om[N][N];//one move (foreinsertion)产生的惩罚差值
int f[N+1];//当前元素所能得到的最大惩罚减少量，动态规划中用到的迭代值
int list_move[3*N];//i,j,type记录move动作
int mov_cnt;
int min_index,max_index;

FILE *fpr,*fpw;

int main()
{
	srand((unsigned int)time(0));

	load_data(20,1000,6,1,1,"wtt10",1000);

	system("pause");
	return EXIT_SUCCESS;
}
//itr迭代次数，dr扰动个数,dr_type=0插入，1交换;oprtct=0插入，1交换,2混合;
void load_data(int run_num,int itr,int dr,int dr_type,int oprtct,char *fn,int mui)
{
	char buffer[MAX_LINE_LENGTH];
	int i,j,k,cnt;
	char *token;
	static char whitespace[]=" \t\f\r\v\n";
	char fnr[100]="WTT_TestData\\";
	char fnw[100]="test1\\";
	char fnw_temp[20];
	strcat(fnr,fn);
	strcat(fnr,".txt");
	fpr=fopen(fnr,"r");

	strcat(fnw,fn);
	strcat(fnw,"_");
	itoa(run_num,fnw_temp,10);
	strcat(fnw,fnw_temp);
	strcat(fnw,"_");
	itoa(itr,fnw_temp,10);
	strcat(fnw,fnw_temp);
	strcat(fnw,"_");
	itoa(dr,fnw_temp,10);
	strcat(fnw,fnw_temp);
	strcat(fnw,"_");
	itoa(dr_type,fnw_temp,10);
	strcat(fnw,fnw_temp);
	strcat(fnw,"_");
	itoa(oprtct,fnw_temp,10);
	strcat(fnw,fnw_temp);
	strcat(fnw,"_cm");
	itoa(mui,fnw_temp,10);
	strcat(fnw,fnw_temp);
	strcat(fnw,".txt");
	fpw=fopen(fnw,"a");

	cnt=0;
	if(fpr!=NULL)
	{
		i=j=0;
		while(fgets(buffer,MAX_LINE_LENGTH,fpr)!=NULL)
		{
			//		fputs(buffer,stdout);
			for(token=strtok(buffer,whitespace);token!=NULL;token=strtok(NULL,whitespace))
			{
				s[j][i]=atoi(token);
				j+=1;
				if(j==N)
				{
					i+=1;
					j=0;
				}
				if(i==3)
				{
					cnt+=1;
					printf("this is NO.%d:\n",cnt);
					fprintf(fpw,"this is NO.%d\n",cnt);
					//for(i=0;i<3;i++)
					//	for(j=0;j<N;j++)
					//		printf("%d  ",s[j][i]);
					if(cnt>=121&&cnt<=125)
					{
						for(k=0;k<run_num;k++)
						{
							init_population_solution();
							sms_iterated_local_search(itr,dr,dr_type,oprtct,mui);
							fflush(fpw);
						}
					}
					i=j=0;
				}
			}
		}
	}else
	{
		perror("can't open file");
		exit(EXIT_FAILURE);
	}
	if(fclose(fpr)!=0)
	{
		perror("fclose,fpr");
		exit(EXIT_FAILURE);
	}	
	if(fclose(fpw)!=0)
	{
		perror("fclose,fpw");
		exit(EXIT_FAILURE);
	}
}

void init_population_solution()
{
	int i,j,r,temp;
	for(i=0;i<P;i++)
	{
		for(j=0;j<N;j++)
			x[i][j]=j;
		for(j=0;j<N;j++)
		{
			r=rand()%(N-j)+j;
			if(r!=j)
			{
				temp=x[i][j];
				x[i][j]=x[i][r];
				x[i][r]=temp;
			}
		}
	}
}

void init_solution(int x[][N],int p)
{
	int i;
	for(i=0;i<N;i++)
		sm[i]=x[p][i];
	csm=0;
	cm[0] = s[sm[0]][0];
	tm[0] = cm[0] - s[sm[0]][2];
	csm += tm[0] > 0 ? tm[0]*tm[0] * s[sm[0]][1]: 0;
	for (i = 1; i < N; i++)
	{
		cm[i] = cm[i - 1] + s[sm[i]][0];
		tm[i] = cm[i] - s[sm[i]][2];
		csm += tm[i] > 0 ? tm[i] * tm[i] *s[sm[i]][1] : 0;
	}
	//csm_pre=csm;
	/**	 printf("\nhere is the initial solution:\ns0:");
	for(i=0;i<N;i++)
	printf("%d ",s0[i]);
	printf("\nc0:");
	for(i=0;i<N;i++)
	printf("%d ",c0[i]);
	printf("\nt0:");
	for(i=0;i<N;i++)
	printf("%d ",t0[i]);
	printf("\t\tcs0:%d",cs0);
	**/
	//printf("csm0:%d\n",csm);
}

//求解交换的ts矩阵
void two_swap()
{
	int i,j,k,temp,dti,dtj,dtm,dt,tempij;//dt i,j,middle的变化量
	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			if(j>=min_index&&i<=max_index)
			{
				ts[i][j]=0;
				if((cm[j]>=s[sm[j]][2])&&(s[sm[i]][2]>s[sm[j]][2]||
					s[sm[i]][2]>cm[j]||
					s[sm[i]][1]<s[sm[j]][1]||
					s[sm[i]][0]>s[sm[j]][0]
				))
				{
					dti = cm[j] - cm[i];//正为推迟，向右，
					dtj = cm[i] + s[sm[j]][0] - cm[j] - s[sm[i]][0];
					dtm = s[sm[j]][0] - s[sm[i]][0];//正为推迟，负为提前

					// for job i
					if (tm[i] <= 0)//non-late job
					{
						temp = dti + tm[i];
						if(temp>0)
							ts[i][j] += temp *temp* s[sm[i]][1];
					}
					else//late job
					{
						ts[i][j] += dti *(dti+2*tm[i]) * s[sm[i]][1];
					}
					// for job j, dtj<0
					if (tm[j] <= 0)//non-late job
					{
						temp = dtj + tm[j];
						if(temp>0)
							ts[i][j] += temp *temp * s[sm[j]][1];
					}
					else//late job
					{
						if (dtj + tm[j] <= 0)
							ts[i][j] -= tm[j] *tm[j] * s[sm[j]][1];
						else
							ts[i][j] += dtj *(dtj+2*tm[j]) * s[sm[j]][1];//change for wtt
					}

					if(dtm<0)//dtm<0
					{
						for (k = i+1; k < j; k++)
						{
							if (tm[k] > 0)//late job
							{
								if(dtm + tm[k] <= 0)
									ts[i][j] -= tm[k] *tm[k] * s[sm[k]][1];
								else
									ts[i][j] += dtm *(dtm+2*tm[k]) * s[sm[k]][1];
							}
						}
					}else if(dtm>0)//此种情况，ts只可能增加；若相等，则只与i，j有关
					{	 
						if(ts[i][j]<0)
						{
							for (k = i+1; k < j; k++)
							{
								if (tm[k] <= 0)//non-late job
								{
									temp = dtm + tm[k];
									if(temp>0)
									{
										ts[i][j] += temp *temp * s[sm[k]][1];
										if(ts[i][j]>=0)
											break;
									}
								}
								else//late job
								{
									ts[i][j] += dtm *(dtm+2*tm[k]) * s[sm[k]][1];
									if(ts[i][j]>=0)
										break;
								}
							}
						}
					}

				}
			}
		}
	}/*
	 for(i=0;i<N;i++)
	 {
	 printf("\n");
	 for(j=0;j<N;j++)
	 {
	 printf("%d  ",ts[i][j]);
	 }
	 }*/
}
//求解向前插入的om矩阵
void one_move()
{   
	int i,j,k,temp,dti,dtj,dtm,dt;//dt i,j,middle的变化量
	for(j=N-1;j>0;j--)
	{
		dtj=0;
		for(i=j-1;i>=0;i--)//j is set before i
		{
			dtj+=s[sm[i]][0];//提前
			/*if(j>=min_index&&i<=max_index)
			{*/
			if(i==j-1)
			{
				om[i][j]=0;
				temp=tm[j]-dtj;//对于job j,提前;
				if(tm[j]>0)//延迟
				{
					if(temp>0)
						om[i][j]-=dtj*s[sm[j]][1];
					else
						om[i][j]-=tm[j]*s[sm[j]][1];
				}
				temp=tm[i]+s[sm[j]][0];//对于job i,推迟
				if(tm[i]>0)//延迟
				{
					om[i][j]+=s[sm[j]][0]*s[sm[i]][1];
				}else
				{
					if(temp>0)
						om[i][j]+=temp*s[sm[i]][1];
				}
			}else
			{
				om[i][j]=om[i+1][j];
				//对于job j,提前;
				temp=tm[j]-dtj;
				if(temp>0)//当前一步延迟,则前一步更延迟
				{
					om[i][j]-=s[sm[i]][0]*s[sm[j]][1];
				}else//当前一步不延迟
				{
					if(temp+s[sm[i]][0]>0)//前一步延迟
						om[i][j]-=(temp+s[sm[i]][0])*s[sm[j]][1];
				}
				//对于job i,推后，直接增加
				temp=tm[i]+s[sm[j]][0];
				if(temp>0)//延迟
				{
					om[i][j]+=temp*s[sm[i]][1];
				}
				//}
			}
		}
	}
	/**
	for (j = 0; j < N; j++)
	{
	dti=0;
	for (i = j+1; i < N; i++)
	{
	dti+=s[sm[i-1]][0];
	if(i>=min_index&&j<=max_index)
	{
	om[j][i]=0;
	if ((cm[i]>s[sm[i]][2])&&(s[sm[j]][2]>s[sm[i]][2]||
	s[sm[j]][2]>cm[i]||
	s[sm[j]][1]<s[sm[i]][1]||
	s[sm[j]][0]>s[sm[i]][0]
	))//i is set before j
	{
	if(i==j+1)//首次情况
	{
	om[j][i]=0;
	temp=tm[i]-dti;//对于job i,提前;
	if(tm[i]>0)//延迟
	{
	if(temp>0)
	om[j][i]-=dti*s[sm[i]][1];
	else
	om[j][i]-=tm[i]*s[sm[i]][1];
	}
	temp=tm[j]+s[sm[i]][0];//对于job j,推迟
	if(tm[j]>0)//延迟
	{
	om[j][i]+=s[sm[i]][0]*s[sm[j]][1];
	}else
	{
	if(temp>0)
	om[j][i]+=temp*s[sm[j]][1];
	}

	}else//非首次情况
	{

	}

	}
	}
	}
	}
	}
	**/
}
//动态规划求解最大惩罚减小量的移动集合，返回最大惩罚减小量
int dynamic_programming_ts()
{
	int i,j;
	mov_cnt=0;
	//printf("\n***\n");
	for(i=0;i<N+1;i++)
		f[i]=0;
	f[0]=f[1]=0;
	if(ts[0][1]<0)
	{
		f[2]=ts[0][1];
		//记录变动
		list_move[mov_cnt]=0;
		list_move[mov_cnt+1]=1;
		list_move[mov_cnt+2]=1;//1 表示two_swap
		mov_cnt+=3;
	}
	for(j=3;j<=N;j++)
	{
		int i_index,temp;
		f[j]=f[j-1];//i=j-1,j is simply appended to j-1][no swap
		for(i=1;i<=j-1;i++)
		{
			temp=f[i-1]+ts[i-1][j-1];
			if(f[j]>temp)
			{
				f[j]=temp;
				i_index=i-1;
			}
		}
		if(f[j]<f[j-1])
		{
			list_move[mov_cnt]=i_index;
			list_move[mov_cnt+1]=j-1;
			list_move[mov_cnt+2]=1;//1 表示two_swap
			mov_cnt+=3;
			//printf("%d\t",i_index);
		}
	}
	//printf("\nf[n]:%d\n",f[N]);
	//if(f[N]==0)
	//	system("pause");
	//for(i=0;i<mov_cnt;i+=3)
	//{
	//	printf("(%d,%d,%d) ",list_move[i],list_move[i+1],list_move[i+2]);
	//}
	return f[N];
}
int dynamic_programming_om()
{
	int i,j;
	mov_cnt=0;
	for(i=0;i<N+1;i++)
		f[i]=0;
	f[0]=f[1]=0;
	if(om[0][1]<0)
	{
		f[2]=om[0][1];
		//记录变动
		list_move[mov_cnt]=0;
		list_move[mov_cnt+1]=1;
		list_move[mov_cnt+2]=0;//1 表示two_swap
		mov_cnt+=3;
	}
	for(j=3;j<=N;j++)
	{
		int i_index,temp;
		f[j]=f[j-1];//i=j-1,j is simply appended to j-1][no swap
		for(i=1;i<=j-1;i++)
		{
			temp=f[i-1]+om[i-1][j-1];
			if(f[j]>temp)
			{
				f[j]=temp;
				i_index=i-1;
			}
		}
		if(f[j]<f[j-1])
		{
			list_move[mov_cnt]=i_index;
			list_move[mov_cnt+1]=j-1;
			list_move[mov_cnt+2]=0;//1 表示two_swap
			mov_cnt+=3;
		}
	}
	//printf("\nf[n]:%d\n",f[N]);
	//if(f[N]==0)
	//	system("pause");
	//for(i=0;i<mov_cnt;i+=3)
	//{
	//	printf("(%d,%d,%d) ",list_move[i],list_move[i+1],list_move[i+2]);
	//}
	return f[N];
}

int dynamic_programming_om_ts()
{
	int i,j,min,min_type;
	mov_cnt=0;
	for(i=0;i<N+1;i++)
		f[i]=0;
	f[0]=f[1]=0;
	if(om[0][1]<ts[0][1])
	{
		min=om[0][1];
		min_type=0;
	}else
	{
		min=ts[0][1];
		min_type=1;
	}
	if(min<0)
	{
		f[2]=ts[0][1];
		//记录变动
		list_move[mov_cnt]=0;
		list_move[mov_cnt+1]=1;
		list_move[mov_cnt+2]=min_type;//1 表示two_swap
		mov_cnt+=3;
	}
	for(j=3;j<=N;j++)
	{
		int i_index,temp;
		f[j]=f[j-1];//i=j-1,j is simply appended to j-1,no swap
		for(i=1;i<=j-1;i++)
		{
			temp=f[i-1]+ts[i-1][j-1];
			if(f[j]>temp)
			{
				f[j]=temp;
				i_index=i-1;
				min_type=1;
			}
		}
		for(i=1;i<=j-1;i++)
		{
			temp=f[i-1]+om[i-1][j-1];
			if(f[j]>temp)
			{
				f[j]=temp;
				i_index=i-1;
				min_type=0;
			}
		}
		if(f[j]<f[j-1])
		{
			list_move[mov_cnt]=i_index;
			list_move[mov_cnt+1]=j-1;
			list_move[mov_cnt+2]=min_type;//1 表示two_swap
			mov_cnt+=3;
		}
	}
	//printf("\nf[n]:%d\n",f[N]);
	//if(f[N]==0)
	//	system("pause");
	//for(i=0;i<mov_cnt;i+=3)
	//{
	//	printf("(%d,%d,%d) ",list_move[i],list_move[i+1],list_move[i+2]);
	//}
	return f[N];
}
//回溯，找最优序列，并重置sm,cm,tm,csm
void backtracking()
{
	int i,temp,j,is_first,max_delta;
	int back_j=N;
	is_first=1;
	//max_delta=0;
	//printf("\nthe backtracking order:\n");
	////printf("\n");
	//for(i=0;i<N;i++)
	//	printf("%d ",sm[i]);
	////	printf("\n");
	for(i=mov_cnt-1;i>=0;i-=3)
	{
		if(list_move[i-1]<back_j)
		{
			//printf("^(%d,%d,%d) ",list_move[i-2],list_move[i-1],list_move[i]);
			/*if(list_move[i-1]-list_move[i-2]>max_delta)
			max_delta=list_move[i-1]-list_move[i-2];*/
			back_j=list_move[i-2];
			if(is_first)
			{
				min_index=list_move[i-2];
				max_index=list_move[i-1];
				is_first=0;
			}else
			{
				if(min_index>list_move[i-2])
					min_index=list_move[i-2];
			}

			//make move
			if(list_move[i]==0)
			{
				//om[j][i]=0; i is set before j,j<i,i=list_move[i-2]
				temp=sm[list_move[i-1]];
				for(j=list_move[i-1];j>list_move[i-2];j--)
					sm[j]=sm[j-1];
				sm[j]=temp;
			}else if(list_move[i]==1)
			{
				temp=sm[list_move[i-2]];
				sm[list_move[i-2]]=sm[list_move[i-1]];
				sm[list_move[i-1]]=temp;
			}
		}
	}
	//printf("max_delta: %d\t",max_delta);
	//printf("\n");
	//for(i=0;i<N;i++)
	//	printf("%d ",sm[i]);
	//printf("\n");
	//update cm,tm,csm
	//csm=0;// reset to 0
	if(mov_cnt>0)
	{	
		cm[0] = s[sm[0]][0];
		tm[0] = cm[0] - s[sm[0]][2];
		//csm += tm[0] > 0 ? tm[0] *tm[0] * s[sm[0]][1] : 0;
		for (i = 1; i < N; i++)
		{
			cm[i] = cm[i - 1] + s[sm[i]][0];
			tm[i] = cm[i] - s[sm[i]][2];
			//csm += tm[i] > 0 ? tm[i] * tm[i] * s[sm[i]][1] : 0;
		}
	}//正确后考虑删去csm and csm_pre
	/*printf("\nsm:");
	for(i=0;i<N;i++)
	{
	printf("%d ",sm[i]);
	}
	printf("\ncm:");
	for(i=0;i<N;i++)
	printf("%d ",cm[i]);
	printf("\ntm:");
	for(i=0;i<N;i++)
	printf("%d ",tm[i]);*/
	//for(i=0;i<N;i++)
	//{
	//	for(j=0;j<N;j++)
	//		printf("%d ",om[i][j]);
	//	printf("\n");
	//}
	//printf("\ncsm:%d,f[N]:%d\n",csm,f[N]);
	//if(csm-csm_pre!=f[N])
	//{
	//	printf("csm:%d,csm_pre:%d,f[N]:%d\n",csm,csm_pre,f[N]);
	//	printf("\nthere appears a error!!\n");	
	//	system("pause");
	//}
	//csm_pre=csm;
	csm+=f[N];
}
//扰动,type为扰动类型，f为扰动数量
int perturbation(int type,int f)
{
	int i,j,k,r1,r2,temp,is_first;
	is_first=1;
	//f-=rand()%2;
	for(i=0;i<f;i++)
	{
		r1=rand()%N;
		r2=rand()%N;
		while(r1==r2)
			r2=rand()%N;
		if(type==1)
		{
			temp=sm[r1];
			sm[r1]=sm[r2];
			sm[r2]=temp;
			if(is_first)
			{
				if(r1<r2)
				{
					min_index=r1;
					max_index=r2;
				}else
				{
					min_index=r2;
					max_index=r1;
				}
				is_first=0;
			}else
			{
				if(r1<r2)
				{
					if(min_index>r1)
						min_index=r1;
					if(max_index<r2)
						max_index=r2;
				}else
				{
					if(min_index>r2)
						min_index=r2;
					if(max_index<r1)
						max_index=r1;
				}
			}
		}else if(type==0){
			if(r1<r2)//r2插入到r1前面
			{
				temp=sm[r2];
				for(j=r2;j>r1;j--)
					sm[j]=sm[j-1];
				sm[r1]=temp;
			}else//r2插入r1后面
			{
				temp=sm[r2];
				for(j=r2;j<r1;j++)
					sm[j]=sm[j+1];
				sm[r1]=temp;
			}
			if(is_first)
			{
				if(r1<r2)
				{
					min_index=r1;
					max_index=r2;
				}else
				{
					min_index=r2;
					max_index=r1;
				}
				is_first=0;
			}else
			{
				if(r1<r2)
				{
					if(min_index>r1)
						min_index=r1;
					if(max_index<r2)
						max_index=r2;
				}else
				{
					if(min_index>r2)
						min_index=r2;
					if(max_index<r1)
						max_index=r1;
				}
			}
		}
	}

	csm=0;// reset to 0
	cm[0] = s[sm[0]][0];
	tm[0] = cm[0] - s[sm[0]][2];
	csm += tm[0] > 0 ? tm[0] *tm[0] * s[sm[0]][1] : 0;
	for (i = 1; i < N; i++)
	{
		cm[i] = cm[i - 1] + s[sm[i]][0];
		tm[i] = cm[i] - s[sm[i]][2];
		csm += tm[i] > 0 ? tm[i] *tm[i] * s[sm[i]][1] : 0;
	}
	//csm_pre=csm;
	//printf("\tcsm:%d\n",csm);
	return csm;
}
//操作策略,0表示one_move,2表示two_swap
void operation_tactics(int oprtct)
{
	if(oprtct==0)
		one_move();
	else if(oprtct==1)
		two_swap();
	else if(oprtct==2)
	{
		one_move();
		two_swap();
	}
}
int dynamic_programming(int oprtct)
{
	if(oprtct==0)
		return dynamic_programming_om();
	else if(oprtct==1)
		return dynamic_programming_ts();
	else if(oprtct==2)
	{
		return dynamic_programming_om_ts();
	}
}
void sms_iterated_local_search(int itr,int dr,int dr_type,int oprtct,int mui)
{
	clock_t start,finish,mid_point;
	int duration,itr_cnt,itr_cnt_pre,i,j,res_pert,res_div,unimp_cnt,duration1,duration2;
	char buffer[50];
	itr_cnt=itr_cnt_pre=unimp_cnt=0;
	min_index=0;
	max_index=N-1;
	init_solution(x,0);
	start=clock();
	operation_tactics(oprtct);
	while(dynamic_programming(oprtct)<0)
	{
		itr_cnt+=1;
		//printf("\niteration %d:",itr_cnt);
		backtracking();
		operation_tactics(oprtct);
	}
	itr_cnt_pre=itr_cnt;
	//	printf("\ncsm:%d\titer:%d",csm,itr_cnt);
	csb=csm;
	for(j=0;j<N;j++)
		sb[j]=sm[j];
	finish=clock();
	duration=finish-start;
	//perturb the local optimal solution
	for(i=0;i<itr;i++)
	{
		//printf("\nptb:%d",i+1);
		perturbation(dr_type,dr);
		operation_tactics(oprtct);
		while(dynamic_programming(oprtct)<0)
		{
			itr_cnt+=1;
			//printf("\niteration %d:",itr_cnt);
			backtracking();
			operation_tactics(oprtct);
		}
		//		printf("\ncsm:%d\titer:%d",csm,itr_cnt-itr_cnt_pre);
		itr_cnt_pre=itr_cnt;
		if(csm<csb)
		{
			csb=csm;
			for(j=0;j<N;j++)
				sb[j]=sm[j];
			unimp_cnt=0;
		finish=clock();
		duration=finish-start;
		}
		else
			unimp_cnt+=1;
		if(unimp_cnt<mui)
		{
			for(j=0;j<N;j++)//每次都从历史最优解扰动
				sm[j]=sb[j];
		}
	}
	//printf("%d\t%d\t%d\n",csb,duration,itr_cnt);
	printf("%d\t%d\n",csb,duration);
	fprintf(fpw,"%d\t%d\n",csb,duration);
	////add the best solution found so far
	//for(j=0;j<N;j++)
	//	fprintf(fpw,"%d\t",sb[j]);
	//fprintf(fpw,"\n\n");
}
#endif
