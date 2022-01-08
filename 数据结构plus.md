ac 自动机

```cpp
#include<bits/stdc++.h>
using namespace std;

#define set0(a) memset(a,0,sizeof(a))
#define rep(i,a,b) for(int i=(a);i<=(b);i++)

const int N=2e5+5;

struct AcAutomaton{
	int tr[N][26], idx;
	int cnt[N], fail[N], id[N];
	int q[N], tt, hh;
	
	void clear(int u){
	    memset(tr[u], 0, sizeof tr[u]), id[u]=cnt[u]=fail[u]=0;
	}
	
	void init(){
		clear(0);
		idx=0;
		tt=-1, hh=0;
	}
	
	void insert(char *s, int index){
		int u=0;
		for(int i=0; s[i]; i++){
			int v=s[i]-'a';
			if(!tr[u][v]){
			     tr[u][v]=++idx;
			     clear(idx);
			}
			u=tr[u][v];
		}
		id[u]=index;
	}
	
	void build(){
		for(int i=0; i<26; i++) if(tr[0][i]){
			fail[tr[0][i]]=0;
			q[++tt]=tr[0][i];
		}
		
		while(tt>=hh){
			int u=q[hh++];
			for(int i=0; i<26; i++){
				int &p=tr[u][i];
				if(p) fail[p]=tr[fail[u]][i], q[++tt]=p;
				else p=tr[fail[u]][i];
			}
		}
	}
	
	void query(char *s){
		int u=0;
		for(int i=0; s[i]; i++){
			u=tr[u][s[i]-'a'];
			for(int p=u; p; p=fail[p]) cnt[id[p]]++;
		}
	}
}ac;

const int M=1e6+5;
char str[155][75], tmp[M];

int main(){
	int n; 
	while(cin>>n, n){
		ac.init();
		
		for(int i=1; i<=n; i++){
			scanf("%s", str[i]);
			ac.insert(str[i], i);
		}
		
		ac.build();
		scanf("%s", tmp);
		
		ac.query(tmp);
		
		int mx=0;
		rep(i,1,n) mx=max(mx, ac.cnt[i]);

		cout<<mx<<endl;
		rep(i,1,n) if(ac.cnt[i]==mx) cout<<str[i]<<endl;
	}
	
	return 0;
}
```

后缀自动机

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=2e6+5, M=N<<1;

#define int long long

struct Edge{
	int to, next;
}e[M];

int h[N], idx;

void add(int u, int v){
	e[idx].to=v, e[idx].next=h[u], h[u]=idx++;
}

char str[N];

struct Node{
	int len, fa; // longest length of state, link
	int ch[26];
}node[N];
int tot=1, last=1;
int f[N];

void ins(int c){
	int p=last, np=last=++tot;
	f[tot]=1;
	node[np].len=node[p].len+1;
	for(; p && !node[p].ch[c]; p=node[p].fa) node[p].ch[c]=np;
	if(!p) node[np].fa=1;
	else{
		int q=node[p].ch[c];
		if(node[q].len==node[p].len+1) node[np].fa=q;
		else{
			int nq=++tot;
			node[nq]=node[q], node[nq].len=node[p].len+1;
			node[q].fa=node[np].fa=nq;
			for(; p && node[p].ch[c]==q; p=node[p].fa) node[p].ch[c]=nq;
		}
	}
}

int res=0;

void dfs(int u){
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		dfs(go);
		f[u]+=f[go];
	}
	if(f[u]>1) res=max(res, f[u]*node[u].len);
}

signed main(){
	scanf("%s", str);
	for(int i=0; str[i]; i++) ins(str[i]-'a');
	memset(h, -1, sizeof h);
	for(int i=2; i<=tot; i++) add(node[i].fa, i);
	dfs(1);
	
	cout<<res<<endl;
		
	return 0;
}
```

ext/rope 块状链表：

```cpp
#include<bits/stdc++.h>
#include<ext/rope>
using namespace std;
using namespace __gnu_cxx;

const int N=25e5;
char s[N];

inline void reads(char *s, int len) {
    s[len]='\0'; len--;
    for(int i=0;i<=len;i++) {
        s[i]='\0';
        while(s[i]<32 || 126<s[i]) s[i]=getchar();
    }
}

rope<char> res;

int main(){
    int q; cin>>q;
    int cur=0;

    while(q--){
        int t;
        string op; cin>>op;
        if(op=="Move"){
            cin>>t; cur=t;
        }
        else if(op=="Insert"){
            cin>>t; reads(s, t);
            res.insert(cur, s);
        }
        else if(op=="Delete"){
            cin>>t;
            res.erase(cur, t);
        }
        else if(op=="Get"){
            cin>>t;
            cout<<res.substr(cur, t)<<endl;
        }
        else if(op=="Prev") cur--;
        else if(op=="Next") cur++;
    }

    return 0;
}
```

莫队

```cpp
#pragma GCC optimize("O3")
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cmath>
using namespace std;

#define endl '\n'

inline int read(){
   int s=0,w=1;
   char ch=getchar();
   while(ch<'0'||ch>'9'){if(ch=='-')w=-1;ch=getchar();}
   while(ch>='0'&&ch<='9') s=s*10+ch-'0',ch=getchar();
   return s*w;
}

const int N=50005, M=2e5+5, Range=1e6+5;
struct query{
    int id, l, r;
}q[M];

int w[N];
int n;
int len;
int cnt[Range], ans[M];

int get(int x){return x/len;}

bool cmp(const query& x, const query& y){
    if(get(x.l)!=get(y.l)) return get(x.l)<get(y.l);
    return x.r<y.r;
}

void add(int v, int& res){
    cnt[v]++;
    if(cnt[v]==1) res++;
}

void del(int v, int& res){
    cnt[v]--;
    if(!cnt[v]) res--;
}

int main(){
    n=read();
    for(int i=1; i<=n; i++) w[i]=read();

    int m; m=read();
    len=sqrt(n);

    for(int i=1; i<=m; i++){
        int l, r; l=read(), r=read();
        q[i]={i,l,r};
    }

    sort(q+1, q+1+m, cmp);

    for(int i=0, j=1, k=1, res=0; k<=m; k++){
        int id=q[k].id, l=q[k].l, r=q[k].r;
        while(i<r) add(w[++i], res);
        while(i>r) del(w[i--], res);
        while(j<l) del(w[j++], res);
        while(j>l) add(w[--j], res);
        ans[id]=res;
    }

    for(int i=1; i<=m; i++) cout<<ans[i]<<endl;

    return 0;
}
```

带修莫队

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=1e4+5, S=1e6+5;

struct Query{
    int id, l, r, t;
}q[N];

struct Modify{
    int p, c;
}c[N];

int cnt[S];
int w[N], ans[N];
int n, m;
int mq, mc;
int len;

int get(int x){
    return x/len;
}

bool cmp(Query &a, Query &b){
    int al=get(a.l), ar=get(a.r), bl=get(b.l), br=get(b.r);
    if(al!=bl) return al<bl;
    if(ar!=br) return ar<br;
    return a.t<b.t; 
}

void add(int v, int &res){
    if(!cnt[v]) res++;
    cnt[v]++;
}

void del(int v, int &res){
    if(cnt[v]==1) res--;
    cnt[v]--;
}

int main(){
    cin>>n>>m;
    for(int i=1; i<=n; i++) cin>>w[i];

    for(int i=0; i<m; i++){
        char op; int l, r; cin>>op>>l>>r;
        if(op=='Q') mq++, q[mq]={mq, l, r, mc};
        else mc++, c[mc]={l, r};
    } 

    len=cbrt((double)n*mc)+1;

    sort(q+1, q+1+mq, cmp);

    for(int i=0, j=1, res=0, t=0, k=1; k<=mq; k++){
        int id=q[k].id, l=q[k].l, r=q[k].r, tm=q[k].t;
        while(i<r) add(w[++i], res);
        while(i>r) del(w[i--], res);
        while(j<l) del(w[j++], res);
        while(j>l) add(w[--j], res);

        while(t<tm){
            t++;
            if(c[t].p>=j && c[t].p<=i){
                del(w[c[t].p], res);
                add(c[t].c, res);
            }
            swap(w[c[t].p], c[t].c);
        }

        while(t>tm){
            if(c[t].p>=j && c[t].p<=i){
                del(w[c[t].p], res);
                add(c[t].c, res);
            }
            swap(w[c[t].p], c[t].c);
            t--;
        }

        ans[id]=res;
    }

    for(int i=1; i<=mq; i++) cout<<ans[i]<<endl;

    return 0;
}
```

点分治

```cpp
#pragma GCC optimize("O3")
#include<bits/stdc++.h>
using namespace std;

#define endl '\n'
#define pb push_back
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define dwn(i,a,b) for(int i=(a);i>=(b);i--)
#define ceil(a,b) (a+(b-1))/b
#define all(x) (x).begin(), (x).end()
using vi = vector<int>;

inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e4+5, M=N<<1;

int n, m;

struct node{
    int to, next, w;
}e[M];

int h[N], tot;

void add(int u, int v, int w){
    e[tot].to=v, e[tot].w=w, e[tot].next=h[u], h[u]=tot++;
}

bool vis[N];

int get_sz(int u, int fa){
    int res=1;
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(go==fa || vis[go]) continue;

        res+=get_sz(go, u);
    }

    return res;
}

int get_wc(int u, int fa, int tot, int &wc){
    if(vis[u]) return 0;
    int sum=1, ms=0;
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(go==fa) continue;

        int t=get_wc(go, u, tot, wc);
        ms=max(ms, t);
        sum+=t;
    }
    ms=max(ms, tot-sum);
    if(ms<=tot/2) wc=u;
    return sum;
}

void get_dis(int u, int fa, int dis, vi &v){
    if(vis[u]) return;
    v.pb(dis);
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(go==fa) continue;

        get_dis(go, u, dis+e[i].w, v);
    }
}

int get(vi v){
    if(v.size()<2) return 0;
    sort(all(v));
    int res=0;
    rep(i,0,v.size()-1) res+=upper_bound(v.begin()+i+1, v.end(), m-v[i])-v.begin()-i-1;
    return res;
}

int calc(int u){
    if(vis[u]) return 0;
    get_wc(u, -1, get_sz(u, -1), u);
    vis[u]=true;

    int res=0;
    vi buf; // 存储重心到其它点的距离。
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        vi tmp;

        get_dis(go, u, e[i].w, tmp);
        res-=get(tmp);

        buf.insert(buf.end(), all(tmp));
    }

    for(int i: buf) if(i<=m) res++;
    res+=get(buf);

    for(int i=h[u]; ~i; i=e[i].next) res+=calc(e[i].to);
    return res;
}

int main(){
    while(cin>>n>>m, n || m){
        rep(i,1,n) h[i]=-1, vis[i]=false;
        tot=0;

        rep(i,1,n-1){
            int u, v, w; read(u), read(v), read(w); u++, v++;
            add(u, v, w), add(v, u, w);
        }   

        cout<<calc(1)<<endl;
    }

    return 0;
}
```

点分树

```cpp
#pragma GCC optimize("O3")
#include<bits/stdc++.h>
using namespace std;

#define int long long
#define endl '\n'
#define debug(x) cerr << #x << ": " << x << endl
#define pb push_back
#define rep(i,a,b) for(int i=(a);i<=(b);i++)

#define all(x) (x).begin(), (x).end()
#define lb(a, x) distance(begin(a), lower_bound(all(a), (x)))
#define ub(a, x) distance(begin(a), upper_bound(all(a), (x)))

inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=15e4+5, M=N<<1;

int n, m, A, age[N];

struct node{
    int to, w, next;
}e[M];

int h[N], tot;

void add(int u, int v, int w){
    e[tot].to=v, e[tot].w=w, e[tot].next=h[u], h[u]=tot++;
}

bool vis[N];

struct Father{
    int u, num, dis;
};

struct Son{
    int age, dis;
    bool operator < (const Son &o)const{
        return age<o.age;
    }
};

vector<Father> f[N];
vector<Son> s[N][3];

int get_sz(int u, int fa){
    if(vis[u]) return 0;
    int res=1;
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(go==fa) continue;
        res+=get_sz(go, u);
    }
    return res;
}

int get_wc(int u, int fa, int tot, int &wc){
    if(vis[u]) return 0;
    int sum=1, ms=0;
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(vis[go] || go==fa) continue;

        int t=get_wc(go, u, tot, wc);
        ms=max(ms, t);
        sum+=t;
    }
    ms=max(ms, tot-sum);
    if(ms<=tot/2) wc=u;
    return sum;
}

void get_dis(int u, int fa, int dis, int wc, int k, vector<Son> &p){
    if(vis[u]) return;
    f[u].pb({wc, k, dis});
    p.pb({age[u], dis});
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(go==fa || vis[go]) continue;
        get_dis(go, u, dis+e[i].w, wc, k, p);
    }
}

void calc(int u){
    if(vis[u]) return;
    get_wc(u, -1, get_sz(u, -1), u);
    vis[u]=true;

    for(int i=h[u], k=0; ~i; i=e[i].next, k++){ // 
        int go=e[i].to;
        if(vis[go]) continue;

        vector<Son> &p=s[u][k];
        p.pb({-1, 0}), p.pb({A+1, 0}); // 哨兵
        get_dis(go, u, e[i].w, u, k, p);

        sort(all(p));

        rep(i,1,p.size()-1) p[i].dis+=p[i-1].dis;
    }
    for(int i=h[u]; ~i; i=e[i].next) calc(e[i].to);
}

int query(int u, int l, int r){
    int res=0;
    for(auto fa: f[u]){ // 遍历与上层的重心的相关信息
        int g=age[fa.u]; // ① 展开和重心节点的讨论
        if(l<=g && g<=r) res+=fa.dis;
        rep(i,0,2){ // 展开和兄弟子树的节点的讨论
            if(i==fa.num) continue;
            vector<Son> &p=s[fa.u][i];
            if(!p.size()) continue; // 空树跳过
            int L=lb(p, Son({l, -1}));
            int R=lb(p, Son({r+1, -1}));
            res+=(R-L)*fa.dis+p[R-1].dis-p[L-1].dis;
        }
    }

    rep(i,0,2){ // 展开和子节点的讨论
        vector<Son> &p=s[u][i];
        if(!p.size()) continue; // 空树跳过
        int L=lb(p, Son({l, -1}));
        int R=lb(p, Son({r+1, -1}));
        res+=p[R-1].dis-p[L-1].dis;
    }

    return res;
}

signed main(){
    memset(h, -1, sizeof h);
    cin>>n>>m>>A;
    rep(i,1,n) read(age[i]);

    rep(i,1,n-1){
        int u, v, w; read(u), read(v), read(w);
        add(u, v, w), add(v, u, w);
    }

    calc(1); 

    int res=0;
    rep(i,1,m){
        int u, l, r; read(u), read(l), read(r);
        l=(l+res)%A, r=(r+res)%A;
        if(l>r) swap(l, r);
        res=query(u, l, r);
        cout<<res<<endl;
    }

    return 0;
}
```

LCT

```cpp
#include<bits/stdc++.h>
using namespace std;

inline void read(int &x) {
    int s=0;x=1;
    char ch=getchar();
    while(ch<'0'||ch>'9') {if(ch=='-')x=-1;ch=getchar();}
    while(ch>='0'&&ch<='9') s=(s<<3)+(s<<1)+ch-'0',ch=getchar();
    x*=s;
}

const int N=1e5+5;

int n, m;

struct Node{
    int s[2], p, v;
    int sum, rev;
}tr[N];
int stk[N];

void pushrev(int x){
    swap(tr[x].s[0], tr[x].s[1]);
    tr[x].rev^=1;
}

void pushup(int x){
    tr[x].sum=tr[tr[x].s[0]].sum^tr[x].v^tr[tr[x].s[1]].sum;
}

void pushdown(int x){
    if(tr[x].rev){
        pushrev(tr[x].s[0]), pushrev(tr[x].s[1]);
        tr[x].rev=0;
    }
}

bool isroot(int x){
    return tr[tr[x].p].s[0]!=x && tr[tr[x].p].s[1]!=x;
}

void rotate(int x){
    int y=tr[x].p, z=tr[y].p;
    int k=tr[y].s[1]==x;
    if(!isroot(y)) tr[z].s[tr[z].s[1]==y]=x;
    tr[x].p=z;
    tr[tr[x].s[k^1]].p=y, tr[y].s[k]=tr[x].s[k^1];
    tr[y].p=x, tr[x].s[k^1]=y;
    pushup(y), pushup(x);
}

void splay(int x){
    int top=0, r=x;
    stk[++top]=r;
    while(!isroot(r)) stk[++top]=r=tr[r].p;
    while(top) pushdown(stk[top--]);

    while(!isroot(x)){
        int y=tr[x].p, z=tr[y].p;
        if(!isroot(y))
            if((tr[y].s[1]==x) ^ (tr[z].s[1]==y)) rotate(x);
            else rotate(y);
        rotate(x);
    }
}

void access(int x){ // 建立从原树的根到 x 的路径，同时 x 成为（splay）根节点。
    int z=x;
    for(int y=0; x; y=x, x=tr[x].p){
        splay(x);
        tr[x].s[1]=y, pushup(x);
    }
    splay(z);
}

// !
void makeroot(int x){ // x 成为原树的根节点
    access(x);
    pushrev(x);
}

int findroot(int x){ // 找到 x 所在原树的根节点，再将原树根节点转到 splay 的根节点
    access(x);
    while(tr[x].s[0]) pushdown(x), x=tr[x].s[0];
    splay(x);
    return x;
}

void split(int x, int y){ // x, y 路径用 splay 维护起来，splay 根为 y
    makeroot(x);
    access(y);
}

void link(int x, int y){ // 若 x, y 不连通，连边
    makeroot(x);
    if(findroot(y)!=x) tr[x].p=y;
}

void cut(int x, int y){ // x, y 间若直接连边，删除之
    makeroot(x); // 让 y 一定是 x 的后继
    if(findroot(y)==x && tr[y].p==x && !tr[y].s[0]){
        tr[x].s[1]=tr[y].p=0;
        pushup(x);
    }
}

int main(){
    cin>>n>>m;
    for(int i=1; i<=n; i++) read(tr[i].v);

    while(m--){
        int op, x, y; read(op), read(x), read(y);
        if(op==0){
            split(x, y);
            cout<<tr[y].sum<<'\n';
        }
        else if(op==1) link(x, y);
        else if(op==2) cut(x, y);
        else if(op==3){
            splay(x);
            tr[x].v=y;
            pushup(x);
        }
    }
    return 0;
}
```

