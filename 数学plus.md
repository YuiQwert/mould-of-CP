线性基

```cpp
#include<bits/stdc++.h>
using namespace std;

using ll=long long;

const int N=3e5+5;

ll a[N], n;

int main(){
    cin>>n;
    for(int i=1; i<=n; i++) cin>>a[i];

    int k=1;
    for(int i=62; ~i; i--){
        for(int j=k; j<=n; j++){
            if(a[j]>>i&1){
                swap(a[j], a[k]);
                break;
            }
        }
        if(!(a[k]>>i&1)) continue;
        for(int j=1; j<=n; j++) if(j!=k && (a[j]>>i&1)) a[j]^=a[k];
        k++;
        if(k==n+1) break;
    }

    ll res=0;
    for(int i=1; i<k; i++) res^=a[i];
    cout<<res;
    return 0;
}
```

FFT

```cpp
#include<bits/stdc++.h>
using namespace std;

const int N=3e5+5;
const double pi=acos(-1);

int n, m;

struct Complex{
    double x, y;
    Complex operator + (const Complex &o)const { return {x+o.x, y+o.y}; }
    Complex operator - (const Complex &o)const { return {x-o.x, y-o.y}; }
    Complex operator * (const Complex &o)const { return {x*o.x-y*o.y, x*o.y+y*o.x}; }
};

Complex a[N], b[N];
int res[N];

int rev[N], bit, tot;

void fft(Complex a[], int inv){
    for(int i=0; i<tot; i++) if(i<rev[i]) swap(a[i], a[rev[i]]);

    for(int mid=1; mid<tot; mid<<=1){
        auto w1=Complex({cos(pi/mid), inv*sin(pi/mid)});
        for(int i=0; i<tot; i+=mid*2){
            auto wk=Complex({1, 0});
            for(int j=0; j<mid; j++, wk=wk*w1){
                auto x=a[i+j], y=wk*a[i+j+mid];
                a[i+j]=x+y, a[i+j+mid]=x-y;
            }
        }
    }
}

int main(){
    cin>>n>>m;
    for(int i=0; i<=n; i++) cin>>a[i].x;
    for(int i=0; i<=m; i++) cin>>b[i].x;

    while((1<<bit)<n+m+1) bit++;

    tot=1<<bit;

    for(int i=0; i<tot; i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));

    fft(a, 1), fft(b, 1);

    for(int i=0; i<tot; i++) a[i]=a[i]*b[i];

    fft(a, -1);

    for(int i=0; i<=n+m; i++) res[i]=(int)(a[i].x/tot+0.5), printf("%d ", res[i]);

    return 0;
}
```

BSGS

```cpp
#include<bits/stdc++.h>
using namespace std;

#define int long long

int bsgs(int a, int p, int b){
    if(1%p==b%p) return 0;
    int k=sqrt(p)+1;
    unordered_map<int, int> hash;

    for(int i=0, j=b%p; i<k; i++){
        hash[j]=i;
        j=j*a%p;
    }   
    int ak=1;
    for(int i=0; i<k; i++) ak=ak*a%p; // get a^k 
    for(int i=1, j=ak; i<=k; i++){
        if(hash.count(j)) return i*k-hash[j];
        j=j*ak%p;
    }
    return -1;
}

signed main(){
    int a, p, b;
    while(cin>>a>>p>>b, a || p || b){
        int res=bsgs(a, p, b);
        if(res==-1) puts("No Solution");
        else cout<<res<<endl;
    }
    return 0;
}
```

