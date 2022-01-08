# dp

## 背包 dp

### 01 背包

```cpp
for(int i=1;i<=n;i++){
    int v,w;
    cin>>v>>w;
    for(int j=t;j>=v;j--){
        f[j]=max(f[j],f[j-v]+w);
    }
}
```

### 完全背包

```cpp
for(int i=1;i<=n;i++)
        for(int j=c[i];j<=t;j++)
            f[j]=max(f[j],f[j-c[i]]+w[i]);
```

### 二维费用 01 背包

```cpp
for(int i=1;i<=n;i++){
    int v1,v2,w;
    cin>>v1>>v2>>w;
    for(int j=t1;j>=v1;j--)
        for(int k=t2;k>=v2;k--)
            f[j][k]=max(f[j][k],f[j-v1][k-v2]+w);
}
```

### 二进制分组优化多重背包

```cpp
for(int i=1;i<=n;i++){
    for(int k=1;k<=s[i];k<<=1){
        s[i]-=k;
        for(int j=t;j>=k*c[i];j--){
            f[j]=max(f[j],f[j-k*c[i]]+w[i]*k);
        }
    }
    if(s[i]){
        for(int j=t;j>=s[i]*c[i];j--)
            f[j]=max(f[j],f[j-s[i]*c[i]]+w[i]*s[i]);
    }
}
```

### 有依赖的背包

$N$ 件物品，背包体积为 $V$，选取一个物品的前提是其父节点对应的物品被选取，最大化价值。

状态表示 $f(u, vol)$ 代表子树 $u$ 在不超过 $vol$ 的背包中所能得到的最大贡献。

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=105;
int f[N][N];
int v[N],w[N];

int head[N],tot;
struct node{int to,next;}e[N];
void add(int u,int v){e[tot].to=v;e[tot].next=head[u];head[u]=tot++;}
int n,t;

void dfs(int cur){
    for(int i=v[cur];i<=t;i++) f[cur][i]=w[cur];
    for(int i=head[cur];~i;i=e[i].next){
        int son=e[i].to;
        dfs(son);
        for(int j=t;j>=v[cur];j--)
            for(int k=0;k<=j-v[cur];k++)
                f[cur][j]=max(f[cur][j],f[cur][j-k]+f[son][k]);
    }
}
int main(){
    memset(head,-1,sizeof head);
    cin>>n>>t;
    
    int root;
    for(int i=1;i<=n;i++){
        cin>>v[i]>>w[i];
        int fa; cin>>fa;
        if(~fa) add(fa,i);
        else root=i;
    }
    
    dfs(root);
    cout<<f[root][t];
    return 0;
}
```



## 状压 dp

```cpp
在 n×n 的棋盘上放 k 个国王，国王可攻击相邻的 8 个格子，求使它们无法互相攻击的方案总数。
```

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N = 12, K = 110, M = 1<<10;

int n,k;
long long f[N][K][M];
int cnt[M];
vector<int> state;
vector<int> head[M];

int lowbit(int x){return x&-x;}
int cal(int x){
    int res=0;
    while(x){x-=lowbit(x);res++;}
    return res;
}

bool check(int state){
    if(!(state&(state<<1)) && !(state&(state>>1))) return true;
    return false;
}
int main(){
    cin>>n>>k;
    for(int i=0;i<1<<n;i++){
        if(check(i)){
            state.push_back(i);
            cnt[i]=cal(i);
        }
    }
    
    for(int i=0;i<state.size();i++)
        for(int j=0;j<state.size();j++)
            {
                int a=state[i]; int b=state[j];
                if(check(a|b) && !(a&b)) head[i].push_back(j);
            }
            
    f[0][0][0]=1;
    for(int i=1;i<=n+1;i++)
        for(int j=0;j<=k;j++)
            for(int a=0;a<state.size();a++){
                int c=cnt[state[a]];
                for(int b:head[a])
                {
                    if(j>=c)
                        f[i][j][a]+=f[i-1][j-c][b];
                }
            }
    cout<<f[n+1][k][0];
    return 0;
}
```



## 区间 dp

> 环形石子合并为例。

```cpp
#include<bits/stdc++.h>
using namespace std;

#define INF 0x3f3f3f3f

const int N = 410;

int f[N][N],g[N][N];
int s[N],w[N];
int n;

int main(){
    cin>>n;
    for(int i=1;i<=n;i++){
        cin>>w[i];
        w[i+n]=w[i];
    }

    memset(f,0x3f,sizeof f);
    memset(g,0xcf,sizeof g);

    for(int i=1;i<=2*n;i++) s[i]=s[i-1]+w[i];

    for(int len=1;len<=n;len++){
        for(int l=1;l+len-1<=n*2;l++){
            int r=l+len-1;

            if(len==1) f[l][r]=g[l][r]=0;
            else{
                for(int k=l;k<r;k++){
                    f[l][r]=min(f[l][r],f[l][k]+f[k+1][r]+s[r]-s[l-1]);
                    g[l][r]=max(g[l][r],g[l][k]+g[k+1][r]+s[r]-s[l-1]);
                }
            }

        }
    }

    int maxv=-INF,minv=INF;
    for(int i=1;i<=n;i++){
        maxv=max(maxv,g[i][i+n-1]);
        minv=min(minv,f[i][i+n-1]);
    }

    cout<<minv<<endl<<maxv<<endl;

    return 0;
}
```



## 树形 dp

> 更具体地说，例题是换根 dp。

给定一个 $n$ 个点的树，请求出一个结点，使得以这个结点为根时，所有结点的深度之和最大。

一个结点的深度之定义为该节点到根的简单路径上边的数量。

```cpp
#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
const int N=1e6+5;

struct node{
	int to, next;
}e[N<<1];
int h[N], tot;
void add(int u, int v){
    e[tot].to=v, e[tot].next=h[u], h[u]=tot++;
}

int n;
ll sz[N], f[N];

int dfs1(int u, int fa, int depth){
	int res=depth;
	sz[u]=1;
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa) continue;
		res+=dfs1(go, u, depth+1);
		sz[u]+=sz[go];
	}
	return res;
}

void dfs2(int u, int fa){
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa) continue;
		f[go]=f[u]+sz[1]-2*sz[go];
		dfs2(go, u);
	}
}

int main(){
	memset(h, -1, sizeof h);
	cin>>n;
	for(int i=1; i<n; i++){
		int u, v; cin>>u>>v;
		add(u, v), add(v, u);
	}
	
	f[1]=dfs1(1, -1, 0);
	dfs2(1, -1);
	
	ll res=0, id;
	
	for(int i=1; i<=n; i++)
		if(f[i]>res) res=f[i], id=i;
	cout<<id<<endl;
	
	return 0;
}
```



## 数位 dp

Windy 定义了一种 Windy 数：不含前导零且相邻两个数字之差至少为 $2$ 的正整数被称为 Windy 数。

Windy 想知道，在 A 和 B 之间，包括 A 和 B，总共有多少个 Windy 数？

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=15;

int f[N][N];

void init(){
    for(int i=0;i<=9;i++) f[1][i]=1;

    for(int i=2;i<=N;i++)
        for(int j=0;j<=9;j++)
            for(int k=0;k<=9;k++)
                if(abs(j-k)>=2) f[i][j]+=f[i-1][k];


}

int dp(int n){
    if(!n) return 0;

    vector<int> v;
    while(n) v.push_back(n%10),n/=10;

    int res=0,last=-2;
    for(int i=v.size()-1;~i;i--){
        int x=v[i];

        for(int j= i==v.size()-1 ;j<x;j++)
            if(abs(j-last)>=2) res+=f[i+1][j];

        if(abs(x-last)>=2) 
            last=x;
        else break;

        if(!i) res++;
    }

    for(int i=1;i<v.size();i++)
        for(int j=1;j<=9;j++)
            res+=f[i][j];

    return res;
}

int main(){
    init();

    int l,r;
    cin>>l>>r;
    cout<<dp(r)-dp(l-1)<<endl;

    return 0;
}
```

