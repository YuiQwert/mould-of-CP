spfa 判负环

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=2005, M=N<<2;

struct node{
	int to,next,w;
}e[M];
int h[N],tot;
void add(int u,int v,int w){e[tot].to=v, e[tot].w=w, e[tot].next=h[u], h[u]=tot++;}

int n,m;
int d[N];
bool vis[N];
int cnt[N];

bool spfa(){
	memset(vis,false,sizeof vis);
	memset(cnt,0,sizeof cnt);
	memset(d,0x3f,sizeof d);
	queue<int> q;
	q.push(1);
	vis[1]=true;
	d[1]=0;
	
	while(q.size()){
		int hd=q.front();
		q.pop();
		vis[hd]=false;
		for(int i=h[hd];~i;i=e[i].next){
			int go=e[i].to;
			if(d[go]>d[hd]+e[i].w){
				d[go]=d[hd]+e[i].w;
				cnt[go]=cnt[hd]+1;
				if(cnt[go]>=n) return true;
				
				if(!vis[go]){
					vis[go]=true;
					q.push(go);
				}
			}
		}
	}
	
	return false;
}

int main(){
	int T; cin>>T;
	while(T--){
		memset(h,-1,sizeof h);
		tot=0;
		cin>>n>>m;
		while(m--){
			int u,v,w; cin>>u>>v>>w;
			if(w>=0) add(u,v,w), add(v,u,w);
			else add(u,v,w);
		}
		
		puts(spfa()?"YES":"NO");
	}
	return 0;
}
```

欧拉回路

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=1e5+5, M=4e5+5;

struct node{
    int to, next;
}e[M];

int h[N], tot;
int din[N], dout[N];

void add(int u, int v){
    e[tot].to=v, e[tot].next=h[u], h[u]=tot++;
}

int op, n, m;

int used[M];
int path[M], cnt;

void dfs(int u){
    for(int &i=h[u]; ~i; ){
        if(used[i]){
            i=e[i].next;
            continue;
        }

        used[i]=true;
        if(op==1) used[i^1]=true;

        int id;
        if(op==1){
            id=i/2+1;
            if(i&1) id=-id;
        }else id=i+1;

        int go=e[i].to;
        i=e[i].next;
        dfs(go);

        path[++cnt]=id;
    }
}

int main(){
    memset(h, -1, sizeof h);
    cin>>op>>n>>m;

    for(int i=0; i<m; i++){
        int u, v; cin>>u>>v;
        add(u, v);
        if(op==1) add(v, u);
        din[v]++, dout[u]++;
    }   

    // 判断是否每个点出度=入度
    if(op==2){
        for(int i=1; i<=n; i++) if(din[i]!=dout[i]){
            puts("NO");
            return 0;
        } 
    }
    else{
        for(int i=1; i<=n; i++) if(din[i]+dout[i]&1){
            puts("NO");
            return 0;
        } 
    }

    for(int u=1; u<=n; u++) if(~h[u]){
        dfs(u);
        break;
    }

    if(cnt<m){ // 未能画完所有边
        puts("NO");
        return 0;
    }

    puts("YES");
    for(int i=cnt; i; i--) cout<<path[i]<<' ';
    cout<<endl;

    return 0;
}
```

边双

```cpp
int n, m;

struct node{
	int to, next;
}e[M];

int h[N], tot;

void add(int u, int v){
	e[tot].to=v, e[tot].next=h[u], h[u]=tot++;
}

int dfn[N], low[N], ts;
int stk[N], top;
int id[N], bcc_cnt;
bool is_bridge[M];

void tarjan(int u, int from){  // 起始点和从前而来的边
	dfn[u]=low[u]=++ts;
	stk[++top]=u;

	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(!dfn[go]){
			tarjan(go, i);
			low[u]=min(low[u], low[go]);
			if(dfn[u]<low[go]) is_bridge[i]=is_bridge[i^1]=true;
		}
		else if(i!=(from^1)) // 非反向边
			low[u]=min(low[u], dfn[go]);
	}
	
	if(dfn[u]==low[u]){
		++bcc_cnt;
		int y;
		do{
			y=stk[top--];
			id[y]=bcc_cnt;
		}while(y!=u);
	}
}
```

点双

```cpp
int dfn[N], low[N], ts;
int stk[N], top;
int bcc_cnt; // v-bcc cnt
vector<int> bcc[N];
bool cut[N];
int root;

void tarjan(int u){
	dfn[u]=low[u]=++ts;
	stk[++top]=u;
	
	if(u==root && h[u]==-1){
		bcc[++bcc_cnt].push_back(u);
		return;
	}
	
	int cnt=0;
	for(int i=h[u]; ~i; i=e[i].next){
		int go=e[i].to;
		if(!dfn[go]){
			tarjan(go);
			low[u]=min(low[u], low[go]);
			
			if(dfn[u]<=low[go]){
				cnt++;
				if(u!=root || cnt>1) cut[u]=true;
				
				++bcc_cnt;
				int y;
				do{
					y=stk[top--];
					bcc[bcc_cnt].push_back(y);
				}while(y!=go);
				bcc[bcc_cnt].push_back(u);
			}
		}
		else low[u]=min(low[u], dfn[go]);
	}
}

```

2-SAT

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

const int N=2e6+5, M=2e6+5;
int n, m;

struct node{
    int to, next;   
}e[M];

int h[N], tot;

void add(int u, int v){
    e[tot].to=v, e[tot].next=h[u], h[u]=tot++;
}

int ts, dfn[N], low[N];
int stk[N], top;
int id[N], cnt;
bool ins[N];

void tarjan(int u){
    dfn[u]=low[u]=++ts;
    stk[++top]=u, ins[u]=true;
    for(int i=h[u]; ~i; i=e[i].next){
        int go=e[i].to;
        if(!dfn[go]){
            tarjan(go);
            low[u]=min(low[u], low[go]);
        }else if(ins[go]) low[u]=min(low[u], dfn[go]);
    }

    if(dfn[u]==low[u]){
        int y;
        cnt++;
        do{
            y=stk[top--], ins[y]=false, id[y]=cnt; 
        }while(y!=u);
    }
}

int res[N];

int main(){
    memset(h, -1, sizeof h);
    read(n), read(m);

    while(m--){
        int i, a, j, b; read(i), read(a), read(j), read(b);
        i--, j--;
        add(2*i+!a, 2*j+b), add(2*j+!b, 2*i+a);
    }

    for(int i=0; i<2*n; i++) if(!dfn[i]) tarjan(i);

    for(int i=0; i<n; i++) if(id[2*i]==id[2*i+1]){
        puts("IMPOSSIBLE");
        return 0;
    }else res[i]= id[2*i]>id[2*i+1];

    puts("POSSIBLE");
    for(int i=0; i<n; i++) printf("%d ", res[i]);
    cout<<endl;

    return 0;
}
```



prufer 序列

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=1e5+5;

int n, m;
int f[N], d[N], p[N]; // father degree(out) prufer

void tree2prufer(){
    for(int i=1; i<=n-1; i++) cin>>f[i], d[f[i]]++;

    for(int i=0, j=1; i<n-2; j++){
        while(d[j]) j++;
        p[i++]=f[j];
        while(i<n-2 && --d[p[i-1]]==0 && p[i-1]<j) p[i++]=f[p[i-1]];
    }   
    for(int i=0; i<n-2; i++) cout<<p[i]<<' ';
    cout<<endl;
}

void prufer2tree(){
    for(int i=1; i<=n-2; i++) cin>>p[i], d[p[i]]++;
    p[n-1]=n;

    for(int i=1, j=1; i<n; i++, j++){
        while(d[j]) j++;
        f[j]=p[i];
        while(i<n-1 && --d[p[i]]==0 && p[i]<j) f[p[i]]=p[i+1], i++;
    }

    for(int i=1; i<=n-1; i++) cout<<f[i]<<' ';
}

int main(){
    cin>>n>>m;
    if(m&1) tree2prufer();
    else prufer2tree();

    return 0;
}
```

