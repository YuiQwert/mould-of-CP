# 补充 plus

## $dfs$​​ 判负环

例题：最小圈

```cpp
inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=3030, M=20020;
const double eps=1e-10;

int n, m;
struct Edge{
	int to, next;
	double w;
}e[M];

int h[N], tot;

void add(int u, int v, double w){
	e[tot].to=v, e[tot].w=w, e[tot].next=h[u], h[u]=tot++;
}

struct Buf{
	int u, v;
	double w;
}buf[M];

double d[N];
bool vis[N], st[N];

bool spfa(int u, double x){
	vis[u]=true;
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(d[go]>d[u]+e[i].w-x){
			d[go]=d[u]+e[i].w-x;
			if(vis[go] || spfa(go, x)) return true;
		}
	}
	vis[u]=false;
	return false;
}

bool ok(double x){
	rep(i,1,n) vis[i]=0, d[i]=0;
	rep(i,1,n) if(spfa(i, x)) return true;
	return false;
}

int main(){
	memset(h, -1, sizeof h);
	read(n), read(m);
	
	rep(i,1,m){
		int u, v; read(u), read(v);
		double w; scanf("%lf", &w);
		add(u, v, w);
	}
			
	double l=-1e7, r=1e7;
	while(l+eps<r){
		double mid=(l+r)/2;
		if(ok(mid)) r=mid;
		else l=mid;
	}

	printf("%.8lf\n", l);

	return 0;
}
```



## cdq 分治

给定 $n$ 个元素（编号 $1∼n$），其中第 $i$ 个元素具有 $ai,bi,ci$ 三种属性。

设 $f(i)$ 表示满足以下 $4$ 个条件：

1. $aj≤ai$
2. $bj≤bi$
3. $cj≤ci$
4. $j≠i$

的 $j$ 的数量。

对于 $d∈[0,n)$ 求满足 $f(i)=d$ 的 $i$ 的数量。

```cpp
#pragma GCC optimize("O3")
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;

const int N=1e5+5, M=2e5+5;
int n,m;
struct data{
    int a,b,c,cnt,res;

    bool operator<(const data &o)const{
        if(a!=o.a) return a<o.a;
        if(b!=o.b) return b<o.b;
        return c<o.c;
    }

    bool operator==(const data &o)const{
        return a==o.a && b==o.b && c==o.c;
    }
}e[N], tmp[N];
int tot;
int tr[M];

int lowbit(int x){return x&-x;}

void add(int p,int k){
    for(;p<M;p+=lowbit(p)) tr[p]+=k;
}

int query(int p){
    int res=0;
    for(;p;p-=lowbit(p)) res+=tr[p];
    return res;
}

void cdq(int l,int r){
    if(l>=r) return;
    int mid=l+r>>1;
    cdq(l,mid), cdq(mid+1,r);
    // 这里定义左区间对应的指针为 j, 右区间对应的指针为 i.
    for(int j=l, i=mid+1, k=l;k<=r;k++)
        if(i>r || j<=mid && e[j].b<=e[i].b) add(e[j].c,e[j].cnt), tmp[k]=e[j++]; // 如果说右区间的指针已经走到边界，或者左区间的b值比较小。
        else e[i].res+=query(e[i].c), tmp[k]=e[i++];
    for(int j=l;j<=mid;j++) add(e[j].c,-e[j].cnt); // 恢复桶
    for(int k=l;k<=r;k++) e[k]=tmp[k];  // 完成排序，复制
}

int ans[M];

int main(){
    cin>>n>>m;
    for(int i=1;i<=n;i++){
        int a,b,c; cin>>a>>b>>c;
        e[i]={a,b,c,1};
    }

    sort(e+1,e+1+n);

    for(int i=1;i<=n;i++)
        if(e[i]==e[tot]) e[tot].cnt++;
        else e[++tot]=e[i];

    cdq(1,tot);

    for(int i=1;i<=tot;i++) ans[e[i].res+e[i].cnt-1]+=e[i].cnt;
    for(int i=0;i<n;i++) cout<<ans[i]<<'\n';

    return 0;
}
```



## DLX

```cpp
// Problem: 精确覆盖问题
// Contest: AcWing
// URL: https://www.acwing.com/problem/content/1069/
// Memory Limit: 64 MB
// Time Limit: 1000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

#include<bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << ": " << (x) << endl
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)

using pii = pair<int, int>;
using ll = long long;

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int M=505;
const int N=M*(M+1); // 1 的个数

int n, m;
int w[M][M];

struct DLX{
	int l[N], r[N], u[N], d[N]; // 四个指针
	int s[N]; // s[j] 表示第 j 列 1 的数量
	int row[N], col[N]; // 相应标号所在行、列
	int idx; // 点的标号
	
	int ans[N], top;
	
	void init(){ // 一开始有 m+1 个点
		rep(i,0,m) l[i]=i-1, r[i]=i+1, u[i]=d[i]=i;
		l[0]=m, r[m]=0;
		idx=m+1;
	}
	
	void add(int &hh, int &tt, int x, int y){
		row[idx]=x, col[idx]=y, s[y]++;
		u[idx]=y, d[idx]=d[y], u[d[y]]=idx, d[y]=idx;
		r[hh]=l[tt]=idx, r[idx]=tt, l[idx]=hh;
		tt=idx++;
	}
	
	void build(){
		rep(i,1,n){
			int hh=idx, tt=idx;
			rep(j,1,m){
				int x=w[i][j];
				if(x) add(hh, tt, i, j);
			}
		}
	}
	
	void remove(int p){
		r[l[p]]=r[p], l[r[p]]=l[p];
		for(int i=d[p]; i!=p; i=d[i]) for(int j=r[i]; j!=i; j=r[j])
			s[col[j]]--, u[d[j]]=u[j], d[u[j]]=d[j];
	}
	
	void resume(int p){
		for(int i=u[p]; i!=p; i=u[i]) for(int j=l[i]; j!=i; j=l[j])
			u[d[j]]=j, d[u[j]]=j, s[col[j]]++;
		r[l[p]]=p, l[r[p]]=p;
	}
	
	bool dfs(){
		if(!r[0]) return true;
		int p=r[0];
		for(int i=r[0]; i; i=r[i]) if(s[i]<s[p]) p=i; // 找到点最少的列
		remove(p);
		for(int i=d[p]; i!=p; i=d[i]){
			ans[++top]=row[i];
			for(int j=r[i]; j!=i; j=r[j]) remove(col[j]); 
			if(dfs()) return true;
			for(int j=l[i]; j!=i; j=l[j]) resume(col[j]);
			top--;
		}
		resume(p);
		return false;
	}
	
	void solve(){
		if(dfs()){
			rep(i,1,top) cout<<ans[i]<<' ';
			puts("");
		}
		else puts("No Solution!");
	}
}dlx;

int main(){
	cin>>n>>m;
	rep(i,1,n) rep(j,1,m) read(w[i][j]);
	
	dlx.init();
	dlx.build();
	dlx.solve();
		
	return 0;
}
```



## NTT（三模数）

```cpp
// Problem: P4245 【模板】任意模数多项式乘法
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P4245
// Memory Limit: 500 MB
// Time Limit: 2000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

#include<bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << ": " << (x) << endl
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)

using pii = pair<int, int>;
using ll = long long;

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=3e5+5;
const ll m1=998244353, m2=1004535809, m3=469762049, M=m1*m2, rt=3;

int n, m, P;
ll a[3][N], b[3][N], ans[N];

int rev[N], tot=1, bit;

ll fpow(ll x, int p, ll mod){
	int res=1;
	for(; p; p>>=1, x=x*x%mod) if(p&1) res=res*x%mod;
	return res;
}

ll inv(ll x, ll mod){
	return fpow(x, mod-2, mod);
}

ll mul(ll x, int p, ll mod){
	ll res=0;
	for(; p; p>>=1, x=(x+x)%mod) if(p&1) res=(res+x)%mod;
	return res;
}

void NTT(ll *a, int type, int mod){
	for(int i=0; i<tot; i++){
		a[i]%=mod;
		if(i<rev[i]) swap(a[i], a[rev[i]]);
	}
	
	for(int mid=1; mid<tot; mid<<=1){
		ll w1=fpow(rt, (type==1? (mod-1)/(mid<<1): mod-1-(mod-1)/(mid<<1)), mod);
		for(int i=0; i<tot; i+=mid*2){
			ll wk=1;
			for(int j=0; j<mid; j++, wk=wk*w1%mod){
				auto x=a[i+j], y=wk*a[i+j+mid]%mod;
				a[i+j]=(x+y)%mod, a[i+j+mid]=(x-y+mod)%mod;
			}
		}
	}
	
	if(type==-1){
		for(int i=0; i<tot; i++) a[i]=a[i]*inv(tot, mod)%mod;
	}
}

void CRT(){
	for(int i=0; i<tot; i++){
		ll res=0;
		(res+=mul(a[0][i]*m2%M, inv(m2, m1), M))%=M;
		(res+=mul(a[1][i]*m1%M, inv(m1, m2), M))%=M;
		a[1][i]=res;
	}
	for(int i=0; i<tot; i++){
		ll res=(a[2][i]-a[1][i]%m3+m3)%m3*inv(M%m3, m3)%m3;
		ans[i]=(M%P*res%P+a[1][i]%P)%P;
	}
}

void solve(int k, int mod){
	NTT(a[k], 1, mod), NTT(b[k], 1, mod);
	for(int i=0; i<tot; i++) a[k][i]=a[k][i]*b[k][i]%mod;
	NTT(a[k], -1, mod);
}

int main(){
	cin>>n>>m>>P;
	rep(i,0,n){
		int t; read(t);
		rep(j,0,2) a[j][i]=t%P;
	}
	rep(i,0,m){
		int t; read(t);
		rep(j,0,2) b[j][i]=t%P;
	}
	
	while(tot<=n+m) bit++, tot<<=1;
	for(int i=0; i<tot; i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
	
	solve(0, m1), solve(1, m2), solve(2, m3);
	CRT();
	
	rep(i,0,n+m) cout<<ans[i]<<' ';
	cout<<endl;
	
	return 0;
}

```





## 分治 FFT

$f$ 满足递推式 $f(n) = \sum_{i=1}^n f(n-i)g(i)$，现在给你 $n$ 还有 $g(1),g(2)\dots g(n-1)$，求出 $f(0),f(1)\dots f(n-1)$，其中 $f(0) = 1$（首项）。

```cpp
// Problem: P4721 【模板】分治 FFT
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P4721
// Memory Limit: 125 MB
// Time Limit: 5000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

#include<bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << ": " << (x) << endl
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)

using pii = pair<int, int>;
using ll = long long;

#define int long long

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

/////////////////////////////////////////////////////////////////////

const int N=3e5+5, mod=998244353, rt=3;

ll fpow(ll x, int p, ll mod){
	int res=1;
	for(; p; p>>=1, x=x*x%mod) if(p&1) res=res*x%mod;
	return res;
}

ll inv(ll x, ll mod){
	return fpow(x, mod-2, mod);
}

int rev[N], tot=1, bit;

void NTT(ll *a, int type, int mod){
	for(int i=0; i<tot; i++){
		a[i]%=mod;
		if(i<rev[i]) swap(a[i], a[rev[i]]);
	}
	
	for(int mid=1; mid<tot; mid<<=1){
		ll w1=fpow(rt, (type==1? (mod-1)/(mid<<1): mod-1-(mod-1)/(mid<<1)), mod);
		for(int i=0; i<tot; i+=mid*2){
			ll wk=1;
			for(int j=0; j<mid; j++, wk=wk*w1%mod){
				auto x=a[i+j], y=wk*a[i+j+mid]%mod;
				a[i+j]=(x+y)%mod, a[i+j+mid]=(x-y+mod)%mod;
			}
		}
	}
	
	if(type==-1){
		for(int i=0; i<tot; i++) a[i]=a[i]*inv(tot, mod)%mod;
	}
}

int n;
int f[N], g[N];

int A[N], B[N];

void conv(int *A, int *B){
	NTT(A, 1, mod), NTT(B, 1, mod);
	rep(i,0,tot-1) (A[i]*=B[i])%=mod;
	NTT(A, -1, mod);
}

void divi(int l, int r){
	if(l>=r) return;
	int mid=l+r>>1;
	divi(l, mid);
	
	// init
	tot=1, bit=0;
	int len=r-l+1;
	while(tot<=len-2) bit++, tot<<=1;
	for(int i=0; i<tot; i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
	rep(i,0,tot-1) A[i]=B[i]=0;
	rep(i,l,mid) A[i-l]=f[i];
	rep(i,1,r-l) B[i-1]=g[i];
	
	conv(A, B);
	
	rep(i,mid+1,r) (f[i]+=A[i-l-1])%=mod;
	divi(mid+1, r);
}

signed main(){
	cin>>n; n--;
	f[0]=1;
	rep(i,1,n) read(g[i]);
	
	divi(0, n);
	
	rep(i,0,n) cout<<f[i]<<' ';
	cout<<endl; 
	return 0;
}
```



## 多项式求逆

```cpp
// Problem: P4238 【模板】多项式乘法逆
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P4238
// Memory Limit: 125 MB
// Time Limit: 1000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

#include<bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << ": " << (x) << endl
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)

using pii = pair<int, int>;
using ll = long long;

#define int long long

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=3e5+5, rt=3, mod=998244353;

int rev[N], tot=1, bit;

ll fpow(ll x, int p, ll mod){
	int res=1;
	for(; p; p>>=1, x=x*x%mod) if(p&1) res=res*x%mod;
	return res;
}

ll inv(ll x, ll mod){
	return fpow(x, mod-2, mod);
}

ll mul(ll x, int p, ll mod){
	ll res=0;
	for(; p; p>>=1, x=(x+x)%mod) if(p&1) res=(res+x)%mod;
	return res;
}

void NTT(ll *a, int type, int mod){
	for(int i=0; i<tot; i++){
		a[i]%=mod;
		if(i<rev[i]) swap(a[i], a[rev[i]]);
	}
	
	for(int mid=1; mid<tot; mid<<=1){
		ll w1=fpow(rt, (type==1? (mod-1)/(mid<<1): mod-1-(mod-1)/(mid<<1)), mod);
		for(int i=0; i<tot; i+=mid*2){
			ll wk=1;
			for(int j=0; j<mid; j++, wk=wk*w1%mod){
				auto x=a[i+j], y=wk*a[i+j+mid]%mod;
				a[i+j]=(x+y)%mod, a[i+j+mid]=(x-y+mod)%mod;
			}
		}
	}
	
	if(type==-1){
		for(int i=0; i<tot; i++) a[i]=a[i]*inv(tot, mod)%mod;
	}
}

int n;
int A[N], B[N], C[N];

void poly_inv(int sz, int *a, int *b){
	if(sz==1) return b[0]=inv(a[0], mod), void();
	poly_inv(sz+1>>1, a, b);
	
	// init
	bit=0, tot=1;
	while(tot<=(sz-1<<1)) tot<<=1, bit++;
	for(int i=0; i<tot; i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
	
	rep(i,0,sz-1) C[i]=a[i];
	rep(i,sz,tot-1) C[i]=0;
	
	NTT(C, 1, mod), NTT(b, 1, mod);
	rep(i,0,tot-1) b[i]=(2-C[i]*b[i]%mod+mod)%mod*b[i]%mod;
	NTT(b, -1, mod);
	
	rep(i,sz,tot-1) b[i]=0;
}

signed main(){
	cin>>n;
	rep(i,0,n-1) read(A[i]);
	poly_inv(n, A, B);
	rep(i,0,n-1) cout<<B[i]<<' ';
	cout<<endl;
	
	return 0;
}
```



## 矩阵快速幂

```cpp
#include<bits/stdc++.h>
using namespace std;

#define int long long
int n, m;
const int N=3;

int v[N]={1, 1, 1};
int A[N][N]={
    {0, 1, 0},
    {1, 1, 1},
    {0, 0, 1}
};

void mul(int a[], int b[][N]){
    int t[N]={0};
    for(int c=0; c<N; c++) for(int r=0; r<N; r++) t[c]=(t[c]+a[r]*b[r][c])%m;
    memcpy(a, t, sizeof t);
}

void mul(int a[][N], int b[][N]){
    int t[N][N]={0};
    for(int c=0; c<N; c++) for(int r=0; r<N; r++) for(int k=0; k<N; k++)
        t[r][c]=(t[r][c]+a[r][k]*b[k][c])%m;
    memcpy(a, t, sizeof t);
}

signed main(){
    cin>>n>>m;

    n--;
    while(n){
        if(n&1) mul(v, A);
        mul(A, A);
        n>>=1;
    }   
    cout<<v[2]<<endl;

    return 0;
} 
```



## 拉格朗日插值

```cpp
// Problem: P4781 【模板】拉格朗日插值
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P4781
// Memory Limit: 125 MB
// Time Limit: 1000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

#include<bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << ": " << (x) << endl
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)

using pii = pair<int, int>;
using ll = long long;

#define int long long

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=2020, mod=998244353;

int n, x[N], y[N], k;

int fpow(int x, int p){
	int res=1;
	for(; p; p>>=1, x=x*x%mod) if(p&1) res=res*x%mod;
	return res;
}

int inv(int x){
	return fpow(x, mod-2);
}

int lagr(int n, int *x, int *y, int k){
	int res=0;
	rep(i,1,n){
		int nu=y[i], de=1;
		rep(j,1,n) if(j!=i) (nu*=(k-x[j])%mod+mod)%=mod, (de*=(x[i]-x[j])%mod+mod)%=mod;
		(res+=nu*inv(de)%mod)%=mod;
	}
	return res;
}

signed main(){
	cin>>n>>k;
	rep(i,1,n) read(x[i]), read(y[i]);
	
	cout<<lagr(n, x, y, k)<<endl;
	
	return 0;
}
```

