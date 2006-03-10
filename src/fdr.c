#include <math.h>
#include <stdlib.h>

void compute_t_statistic(double *m_array, int *columns_a,int *columns_b,int *n_a,int *n_b,int *n_total,int *genes_num,double *t)
{
	int i,j=0;
	double sum_a,sum_b,mean_a,mean_b,var_a,var_b,s_pooled;
	
	for (j=0;j<*genes_num;j++)
	{
		sum_a=0;
		sum_b=0;
		var_a=0;
		var_b=0;
		
		for (i=0;i<*n_a;i++)
			sum_a+=m_array[j+columns_a[i]*(*genes_num)];	
		
		for (i=0;i<*n_b;i++)
			sum_b+=m_array[j+columns_b[i]*(*genes_num)];	
		
		mean_a=sum_a/(*n_a);
		mean_b=sum_b/(*n_b);
		
		for (i=0;i<*n_a;i++)
			var_a+= (m_array[j+columns_a[i]*(*genes_num)]-mean_a)*(m_array[j+columns_a[i]*(*genes_num)]-mean_a)/(*n_a-1);
		
		for (i=0;i<*n_b;i++)
			var_b+= (m_array[j+columns_b[i]*(*genes_num)]-mean_b)*(m_array[j+columns_b[i]*(*genes_num)]-mean_b)/(*n_b-1);
		s_pooled=((var_a*(*n_a-1))+(var_b*(*n_b-1)))/(*n_a+*n_b-2);
		t[j]=(mean_a-mean_b)/sqrt(s_pooled*((1/(*n_a))+(1/(*n_b))));
	}

}

void permute(int *a,int n)
{
	int i,temp,r;
	for (i=0;i<n;i++)
	{
		r=i+(int)(rand()%(n-i));
		temp=a[i];
		a[i]=a[r];
		a[r]=temp;

	}
}




void compute_resampling_t_stat(double *m_array, int *n_a, int *n_b, int *genes_num, int *perms_num,double *statistic_array)
{
	int i,per;
	int n_total;
	n_total=(*n_a+*n_b);
	int a[n_total];
	
	for (i=0;i<n_total;i++)
	{
		a[i]=i;
	}

	for(per=0;per<*perms_num;per++)
	{
		permute(a,n_total);
		//compute_t_statistic(m_array,a,(a+*n_a),n_a,n_b,n_total,genes_num,statistic_array+(*genes_num)*(per));
		
		compute_t_statistic(m_array,a,(a+*n_a),n_a,n_b,&n_total,genes_num,(statistic_array+(*genes_num)*(per)));
	}
}

///////////////////////////////////////////////////////
//			adaptive		     //
///////////////////////////////////////////////////////

double q(int i,int j,int *m,double pvaluei)
{
	double a;
	a=(pvaluei/i)*((*m)-(j));
	if (a<1) a=(a*(1/(1-pvaluei*((*m)-j)/i)));
	else a=1;	
	if (a<1) return(a);
	else return 1;
}

		

void adaptive(double *pvalues,int *m,double *adapk)
{
	int j,k;
	double qij,mini=0;
	double a,b;
	double first_stage[*m];
	for (j=*m;j>=1;j--)
	{
		first_stage[j-1]=(pvalues[j-1]/j)*(*m);
		if (j<*m)
			if (first_stage[j-1]>first_stage[j])
				first_stage[j-1]=first_stage[j];
	}
		
	for (k=*m;k>0;k--)
	{
		for (j=1;j<=k;j++)
		{		
			a=q(k,j,m,pvalues[k-1]);
			b=first_stage[(j-1)];
		
			if (a>b) qij=a;
			else qij=b;
			if ((j==1)||(qij<mini)) mini=qij;
		}
		if (k==(*m))
			adapk[k-1]=mini;
		else if (k<(*m))
		{
			if (adapk[k]<mini) adapk[k-1]=adapk[k];
			else	adapk[k-1]=mini;
	
		}
		if (adapk[k-1]>first_stage[k-1])
			adapk[k-1]=first_stage[k-1];
	}

}

void check(int *a,int *n)
{
	int i,temp,r;
	for (i=0;i<*n;i++)
	{
		r=i+(int)(rand()%(*n-i));
		temp=a[i];
		a[i]=a[r];
		a[r]=temp;
		
	}
	
}


