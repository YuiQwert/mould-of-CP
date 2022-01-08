数据结构 pp

## dsu on tree

例题：节点有颜色的树，需要统计每个棵子树内出现次数最多的颜色之和。

```cpp
// Problem: E. Lomsat gelral
// Contest: Codeforces - Educational Codeforces Round 2
// URL: https://codeforces.com/problemset/problem/600/E
// Memory Limit: 256 MB
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

#define int ll

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e5+5, M=N<<1;

int n, w[N];

struct Edge{
	int to, next;
}e[M];

int h[N], tot;

void add(int u, int v){
	e[tot].to=v, e[tot].next=h[u], h[u]=tot++;
}

int sz[N], son[N];

void dfs(int u, int fa){
	sz[u]=1;
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa) continue;
		dfs(go, u);
		sz[u]+=sz[go];
		if(sz[go]>sz[son[u]]) son[u]=go;
	}
}

int res[N];
int cnt[N], ans, mx;

void upd(int u, int fa, int d){
	cnt[w[u]]+=d;
	if(d==1){
		if(cnt[w[u]]>mx) mx=cnt[w[u]], ans=w[u];
		else if(cnt[w[u]]==mx) ans+=w[u];
	}
	
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa) continue;
		upd(go, u, d);
	}
}

void dfs2(int u, int fa, bool del){
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa || go==son[u]) continue;
		dfs2(go, u, 1);
	}
	if(son[u]) dfs2(son[u], u, 0);
	
	cnt[w[u]]++;
	if(cnt[w[u]]>mx) mx=cnt[w[u]], ans=w[u];
	else if(cnt[w[u]]==mx) ans+=w[u];
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa || go==son[u]) continue;
		upd(go, u, 1);
	}
	res[u]=ans;
	
	if(del) upd(u, fa, -1), mx=ans=0;
}

signed main(){
	memset(h, -1, sizeof h);
	cin>>n;
	rep(i,1,n) read(w[i]);
	
	rep(i,1,n-1){
		int u, v; read(u), read(v);
		add(u, v), add(v, u);
	}
	
	dfs(1, -1);
	dfs2(1, -1, 0);

	rep(i,1,n) cout<<res[i]<<' ';
	cout<<endl;
	
	return 0;
}
```



## 线段树优化建图



```cpp
// Problem: B. Legacy
// Contest: Codeforces - Codeforces Round #406 (Div. 1)
// URL: https://codeforces.com/problemset/problem/786/B
// Memory Limit: 256 MB
// Time Limit: 2000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

#include<bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << ": " << (x) << endl
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)

using pli = pair<long long, int>;
using ll = long long;

const int N=2e6+5, M=N+1e5*30;
const ll INF=1e18;

int n, q, S;

struct Edge{
	int to, next;
	ll w;
}e[M];

int h[N], tot;

void add(int u, int v, ll w){
	e[tot].to=v, e[tot].w=w, e[tot].next=h[u], h[u]=tot++;
}

#define ls u<<1
#define rs u<<1|1

int id[N];
int idx, tmp;
int rtL, rtR;
void build(int u, int l, int r, int op){
	if(op==-1){
		id[u]=++idx; 
		int mid=l+r>>1;
		build(ls, l, mid, 0), build(rs, mid+1, r, 1);
		rtL=ls, rtR=rs;
		return;
	}
	if(l==r){
		id[u]=++tmp;
		return;
	}

	int mid=l+r>>1;
	id[u]=++idx; 
	build(ls, l, mid, op), build(rs, mid+1, r, op);
	if(!op) add(id[u], id[ls], 0), add(id[u], id[rs], 0);
	else add(id[ls], id[u], 0), add(id[rs], id[u], 0);
}

void link(int u, int l, int r, int cur, int nl, int nr, ll w, int op){
	if(nl<=l && r<=nr){
		if(!op) add(cur, id[u], w);
		else add(id[u], cur, w);
		return;
	}
	int mid=l+r>>1;
	if(nl<=mid) link(ls, l, mid, cur, nl, nr, w, op);
	if(mid<nr) link(rs, mid+1, r, cur, nl, nr, w, op);
}

bool vis[N];
ll d[N];
priority_queue<pli, vector<pli>, greater<pli>> pq;

#define x first
#define y second

void dijk(){
	rep(i,1,idx) d[i]=1e18;
	d[S]=0;
	pq.push({d[S], S});
	
	while(pq.size()){
		auto t=pq.top(); pq.pop();
		int u=t.y; 
		if(vis[u]) continue;
		vis[u]=true;
		
		for(int i=h[u]; ~i; i=e[i].next){
			int go=e[i].to;
			if(d[go]>e[i].w+d[u]){
				d[go]=e[i].w+d[u];
				pq.push({d[go], go});
			}
		}
	}
	
	rep(i,n+1,n<<1) cout<<(d[i]==INF? -1: d[i])<<' ';
	cout<<endl;
}

int main(){
	memset(h, -1, sizeof h);
	cin>>n>>q>>S;
	idx=n<<1;
	build(1, 1, n<<1, -1);
	rep(i,1,n) add(i, i+n, 0);
	// rep(i,1,idx) debug(id[i]);
	
	while(q--){
		int op; scanf("%d", &op);
		int u, v;
		ll w;
		if(op==1){
			scanf("%d%d%lld", &u, &v, &w);
			add(u+n, v, w);
		}
		else if(op==2){
			int l, r; scanf("%d%d%d%lld", &u, &l, &r, &w);
			link(rtL, 1, n, u, l, r, w, 0);
		}
		else if(op==3){
			int l, r; scanf("%d%d%d%lld", &u, &l, &r, &w);
			link(rtR, n+1, n<<1, u, l+n, r+n, w, 1);
		}
	}
	
	
	dijk();

	return 0;
}
```



## 线段树毒瘤题（序列操作）

lxhgww 最近收到了一个 01 序列，序列里面包含了 n 个数，下标从 0 开始。这些数要么是 0，要么是 1，现在对于这个序列有五种变换操作和询问操作：

- `0 l r` 把 [l, r][*l*,*r*] 区间内的所有数全变成 0
- `1 l r` 把 [l, r][*l*,*r*] 区间内的所有数全变成 1
- `2 l r` 把 [l,r][*l*,*r*] 区间内的所有数全部取反，也就是说把所有的 0 变成 1，把所有的 1 变成 0
- `3 l r` 询问 [l, r][*l*,*r*] 区间内总共有多少个 1
- `4 l r` 询问 [l, r][*l*,*r*] 区间内最多有多少个连续的 1

```cpp
// #pragma GCC optimize("O3")
#include<bits/stdc++.h>
using namespace std;

#define endl '\n'
#define debug(x) cerr << #x << ": " << x << endl
#define pb push_back
#define eb emplace_back
#define set0(a) memset(a,0,sizeof(a))
#define rep(i,a,b) for(int i=(a);i<=(b);i++)

inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e5+5;

int n, m;
int w[N];
int tag1[N<<2], tag2[N<<2];

struct Tree{
	int l, r, len;
	int sum, con[2], pre[2], suf[2]; // longest connected, prefix, suffix
	
	#define ls(u) u<<1
	#define rs(u) u<<1|1
}tr[N<<2];

Tree merge(Tree a, Tree b){
	Tree res={
		a.l, b.r, a.len+b.len,
		a.sum+b.sum
	};
	if(!a.sum) res.pre[0]=a.len+b.pre[0];
	else res.pre[0]=a.pre[0];
	if(a.sum==a.len) res.pre[1]=a.len+b.pre[1];
	else res.pre[1]=a.pre[1];
	
	if(!b.sum) res.suf[0]=b.len+a.suf[0];
	else res.suf[0]=b.suf[0];
	if(b.sum==b.len) res.suf[1]=b.len+a.suf[1];
	else res.suf[1]=b.suf[1];
	
	rep(i,0,1) res.con[i]=max(max(a.con[i], b.con[i]), a.suf[i]+b.pre[i]);
	
	return res;
}

void F(int u, int op){
	if(op!=2){
		int t=op;
		tr[u].sum=t*tr[u].len;
		tr[u].con[t]=tr[u].pre[t]=tr[u].suf[t]=tr[u].len;
		tr[u].con[t^1]=tr[u].pre[t^1]=tr[u].suf[t^1]=0;
		tag1[u]=op, tag2[u]=0;
	}
	else{
		tag2[u]^=1;
		tr[u].sum=tr[u].len-tr[u].sum;
		swap(tr[u].con[0], tr[u].con[1]);			 
		swap(tr[u].pre[0], tr[u].pre[1]);
		swap(tr[u].suf[0], tr[u].suf[1]);
	}
}

void pushdown(int u){
	if(~tag1[u]){
		F(ls(u), tag1[u]), F(rs(u), tag1[u]);
		tag1[u]=-1;
	}
	if(tag2[u]){
		F(ls(u), 2), F(rs(u), 2);
		tag2[u]=0;
	}
}

void build(int u, int l, int r){
	tag1[u]=-1;
	if(l==r){
		int v=w[l];
		tr[u]={l, r, 1, v};
		tr[u].con[0]=tr[u].pre[0]=tr[u].suf[0]=v^1;
		tr[u].con[1]=tr[u].pre[1]=tr[u].suf[1]=v;
		return;
	}
	int mid=l+r>>1;
	build(ls(u), l, mid), build(rs(u), mid+1, r);
	tr[u]=merge(tr[ls(u)], tr[rs(u)]);
}

void modify(int u, int l, int r, int op){
	if(tr[u].l>=l && tr[u].r<=r){
		F(u, op);
		return;
	}
	pushdown(u);
	int mid=tr[u].l+tr[u].r>>1;
	if(l<=mid) modify(ls(u), l, r, op); 
	if(r>mid) modify(rs(u), l, r, op);
	tr[u]=merge(tr[ls(u)], tr[rs(u)]);
}

Tree query(int u, int l, int r){
	if(tr[u].l>=l && tr[u].r<=r) return tr[u];
	pushdown(u);
	int mid=tr[u].l+tr[u].r>>1;
	Tree L={-1}, R={-1}; //
	if(l<=mid) L=query(ls(u), l, r);
	if(r>mid) R=query(rs(u), l, r);
	if(L.l==-1) return R; if(R.l==-1) return L;
	return merge(L, R);
}

int main(){
	read(n), read(m);
	rep(i,1,n) read(w[i]);
	
	build(1, 1, n);
	
	rep(i,1,m){
		int op, l, r; read(op), read(l), read(r); l++, r++;
		if(op<=2) modify(1, l, r, op);
		else{
			Tree t=query(1, l, r);
			if(op==3) cout<<t.sum<<endl;
			else if(op==4) cout<<t.con[1]<<endl;
		}
		
		if(i==5){
			Tree t=query(1, l, r);
			debug(t.con[1]);
		}
	}
    return 0;
}
```



## 势能线段树

例题：

第一行一个整数 n，代表数列中数的个数。

第二行 n 个正整数，表示初始状态下数列中的数。

第三行一个整数 m，表示有 m 次操作。

接下来 m 行每行三个整数 `k l r`。

- k=0 表示给 [l,r][*l*,*r*] 中的每个数开平方（下取整）。
- k=1 表示询问 [l,r][*l*,*r*] 中各个数的和。



```cpp
// #pragma GCC optimize("O3")
#include<bits/stdc++.h>
using namespace std;

#define endl '\n'
#define debug(x) cerr << #x << ": " << x << endl
#define pb push_back
#define eb emplace_back
#define set0(a) memset(a,0,sizeof(a))
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)
#define ceil(a,b) (a+(b-1))/b

#define all(x) (x).begin(), (x).end()
#define SUM(a) accumulate(all(a), 0LL)
#define MIN(a) (*min_element(all(a)))
#define MAX(a) (*max_element(all(a)))
#define lb(a, x) distance(begin(a), lower_bound(all(a), (x)))
#define ub(a, x) distance(begin(a), upper_bound(all(a), (x)))

#define INF 0x3f3f3f3f
#define ll_INF 0x7f7f7f7f7f7f7f7f

using pii = pair<int, int>;
using pdd = pair<double, double>;
using vi = vector<int>;
using vvi = vector<vi>;
using vb = vector<bool>;
using vpii = vector<pii>;
using ll = long long;
using ull = unsigned long long;

#define int ll

inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e5+5;

int n, q, w[N];

struct Node{
	int l, r;
	int mx, sum;
	
	#define ls(u) u<<1
	#define rs(u) u<<1|1
}tr[N<<2];

void pushup(int u){
	tr[u].mx=max(tr[ls(u)].mx, tr[rs(u)].mx);
	tr[u].sum=tr[ls(u)].sum+tr[rs(u)].sum;
}

void build(int u, int l, int r){
	tr[u]={l, r, w[l], w[r]};
	if(l==r) return;
	int mid=l+r>>1;
	build(ls(u), l, mid), build(rs(u), mid+1, r);
	pushup(u);
}

void upd(int u, int l, int r){
	if(tr[u].l==tr[u].r) return tr[u].mx=tr[u].sum=sqrt(tr[u].sum), void();
	int mid=tr[u].l+tr[u].r>>1;
	if(l<=mid && tr[ls(u)].mx>1) upd(ls(u), l, r);
	if(mid<r && tr[rs(u)].mx>1) upd(rs(u), l, r);
	pushup(u); 
}

int query(int u, int l, int r){
	if(l<=tr[u].l && tr[u].r<=r) return tr[u].sum;
	int mid=tr[u].l+tr[u].r>>1, res=0;
	if(l<=mid) res+=query(ls(u), l, r);
	if(mid<r) res+=query(rs(u), l, r);
	return res;
}

signed main(){
	read(n);
	rep(i,1,n) read(w[i]);
	read(q);
	
	build(1, 1, n);
	
	while(q--){
		int op, l, r; read(op), read(l), read(r);
		if(l>r) swap(l, r);
		if(!op) upd(1, l, r);
		else cout<<query(1, l, r)<<endl;
	}
    return 0;
}
```



## 分块

```cpp
// Problem: P2801 教主的魔法
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P2801
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

inline void read(int &x){
    int s=0; x=1;
    char ch=getchar();
    while(ch<'0' || ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0' && ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e6+6, M=1010;

int n, q;
int len;

struct Node{
	int x, w;
	bool operator < (const Node &o)const{
		return w<o.w;
	}
}e[N];

int get(int x){
	return (x+len-1)/len;
}

int getL(int id){
	return (id-1)*len+1;
}

int getR(int id){
	return min(id*len, n);
}

int add[M];

int query(int l, int r, int k){
	int res=0;
	if(get(l)==get(r)){
		int id=get(l);
		if(add[id]){
			rep(i,getL(id),getR(id)) e[i].w+=k;
			add[id]=0;
		} 
		rep(i,getL(id),getR(id)) if(e[i].x>=l && e[i].x<=r && e[i].w>=k) res++;
		return res; 
	}
	int id=get(l);
	if(add[id]){
		rep(i,getL(id),getR(id)) e[i].w+=k;
		add[id]=0;
	} 
	rep(i,getL(id),getR(id)) if(e[i].x>=l && e[i].x<=getR(id) && e[i].w>=k) res++;
	// debug(res);
	id=get(r);
	if(add[id]){
		rep(i,getL(id),getR(id)) e[i].w+=k;
		add[id]=0;
	} 
	rep(i,getL(id),getR(id)) if(e[i].x>=getL(id) && e[i].x<=r && e[i].w>=k) res++;
	// debug(res);
	int L=get(l)+1, R=get(r)-1;
	rep(id,L,R){
		int l=getL(id), r=getR(id);
		while(l<r){
			int mid=l+r>>1;
			if(e[mid].w>=k-add[id]) r=mid;
			else l=mid+1;
		}
		res+=(e[l].w>=k-add[id]? getR(id)-l+1: 0); 
	}
	
	return res;
}

void upd(int l, int r, int k){
	if(get(l)==get(r)){
		int id=get(l);
		rep(i,getL(id),getR(id)) if(e[i].x>=l && e[i].x<=r) e[i].w+=k;
		sort(e+getL(id), e+getR(id)+1); 
		return;
	}
	int id=get(l);
	rep(i,getL(id),getR(id)) if(e[i].x>=l && e[i].x<=getR(id)) e[i].w+=k;
	sort(e+getL(id), e+getR(id)+1); 
	
	id=get(r);
	rep(i,getL(id),getR(id)) if(e[i].x>=getL(id) && e[i].x<=r) e[i].w+=k;
	sort(e+getL(id), e+getR(id)+1); 
	
	int L=get(l)+1, R=get(r)-1; 
	rep(id,L,R) add[id]+=k;
}

void out(){
	rep(i,1,n) cerr<<e[i].w<<' ';
	cerr<<endl; 
}

int main(){
	cin>>n>>q;
	len=sqrt(n+1);
	rep(i,1,n) read(e[i].w), e[i].x=i;

	rep(id,1,get(n)) sort(e+getL(id), e+getR(id)+1);

	while(q--){
		char ch; cin>>ch;
		int l, r, v; read(l), read(r), read(v);
		if(ch=='A') cout<<query(l, r, v)<<endl;
		else upd(l, r, v);
	}
	return 0;
}
```



##  线段树合并

例题

给出一颗树，每个点都有一个权值，最后对于每个点，输出在它的子树中，有多少个点的权值比它大。

```cpp
// #pragma GCC optimize("O3")
#include<iostream>
#include<algorithm>
#include<cstring>
#include<vector>
using namespace std;

#define endl '\n'
#define debug(x) cerr << #x << ": " << x << endl
#define pb push_back
#define eb emplace_back
#define set0(a) memset(a,0,sizeof(a))
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)
#define INF 0x3f3f3f3f
#define ll_INF 0x7f7f7f7f7f7f7f7f

#define lb(a, x) distance(begin(a), lower_bound(all(a), (x)))
#define all(x) (x).begin(), (x).end()

using pii = pair<int, int>;
using pdd = pair<double, double>;
using vi = vector<int>;
using vvi = vector<vi>;
using vb = vector<bool>;
using vpii = vector<pii>;
using ll = long long;
using ull = unsigned long long;

#define x first
#define y second

inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e5+5, M=N<<1;

vi nums;

int find(int x){
	return lb(nums, x)+1;
}

int w[N], n;
int res[N];

struct Edge{
	int to, next;
}e[M];

int h[N], tot;

void add(int u, int v){
	e[tot].to=v, e[tot].next=h[u], h[u]=tot++;
}

struct Node{
	int l, r;
	int cnt;
	
	#define ls tr[u].l
	#define rs tr[u].r
}tr[N*22];

int root[N], idx;

void pushup(int u){
	tr[u].cnt=tr[ls].cnt+tr[rs].cnt;
}

void insert(int &u, int l, int r, int x){
	if(!u) u=++idx;
	if(l==r) return tr[u].cnt++, void();
	int mid=l+r>>1;
	if(x<=mid) insert(ls, l, mid, x);
	else insert(rs, mid+1, r, x);
	pushup(u);
}

int merge(int u, int v){
	if(!u) return v;
	if(!v) return u;
	ls=merge(ls, tr[v].l);
	rs=merge(rs, tr[v].r);
	pushup(u);
	return u;
}

int query(int u, int l, int r, int nl, int nr){
	if(!u) return 0;
	if(nl<=l && r<=nr) return tr[u].cnt;
	int mid=l+r>>1, res=0;
	if(nl<=mid) res+=query(ls, l, mid, nl, nr);
	if(mid<nr) res+=query(rs, mid+1, r, nl, nr);
	return res;
}

void dfs(int u, int fa){
	insert(root[u], 1, n, find(w[u]));
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(go==fa) continue;
		dfs(go, u);
		root[u]=merge(root[u], root[go]);
	}
	if(find(w[u])+1<=n) res[u]=query(root[u], 1, n, find(w[u])+1, n);
}

int main(){
	read(n);
	rep(i,1,n) read(w[i]), nums.pb(w[i]);
	
	sort(all(nums));
	nums.erase(unique(all(nums)), nums.end());
	
	memset(h, -1, sizeof h);
	rep(i,2,n){
		int fa; read(fa);
		add(fa, i), add(i, fa);
	}
	
	dfs(1, -1);
	
	rep(i,1,n) cout<<res[i]<<endl;
	return 0;
}
```

