#include <bits/stdc++.h>
#define pb push_back
#define Int int64_t
using namespace std;

using D=long double;
Int N,M,K;
Int X,Y,Z;
Int h,w;

const int MAXN=3e5+5;
const int NAX=105;

vector<int> adj[MAXN];
//vector<pair<int,int>> adj[MAXN];

Int ceilS(Int x,Int y){
    return (x+y-1)/y;
}
const int MOD = 998244353;
//const int MOD=1e9+7;

Int perform(Int x){
    return(x*(x+1))/2LL;
}

long long ModExpo(long long x,unsigned long long y,long long M){
    Int ans=1;
    ans=(long long)ans;
    
    while(y>0){
        if(y&1) ans=((ans%M)*(x%M))%M;
        y>>=1LL;
        x=((x%M)*(x%M))%M;
    }
    return ans;
}

vector<pair<int64_t,int64_t>> factorize(int64_t x){

    vector<pair<int64_t,int64_t>> factors;

    for(int64_t i=2;i*i<=x;++i){
        if(x % i) continue;
        int cnt=0;

        while(x % i ==0){
            x/=i;
            ++cnt;
        }
        factors.pb({i,cnt});
    }
    if(x > 1){
        factors.pb({x,1});
    }
    return factors;
}


long long ModInv(long long x){
    return ModExpo(x,MOD-2,MOD);
}
vector<long long> fac(MAXN);

vector<Int> A;
const int INF=1e9+5;
Int madd(Int x,Int y){
    return (x+y)%MOD;
}

Int mmul(Int x,Int y){
    return ((x%MOD)*(y%MOD))%MOD;
}

Int performMod(Int x){
    return mmul(mmul(x,x+1),ModInv(2))%MOD;
}


void fill(){
    fac[0]=1;
    for(int i=1;i<MAXN;++i){
        fac[i]=mmul(fac[i-1],i)%MOD;
    }
}

long long choose(long long n,long long r){
    if(r==0||r==n) return 1;
    if(r>n) return 0;

    return mmul(mmul(fac[n],ModInv(fac[n-r])),ModInv(fac[r]))%MOD;
}


const int nax=4e6+5;
Int choose3(Int x,Int y){
    return mmul(mmul(mmul(x,x-1)%MOD,x-2)%MOD,ModInv(6)%MOD)%MOD;
}

Int sqSum(Int x){
    return mmul(mmul(x,mmul((2*x+1)%MOD,x+1)),ModInv(6));
}

#define pi 3.1415926535897932385
struct FastFourierTransform{

using cd = complex<double>;
const double PI = acos(-1);

    int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

        void fft(vector<cd> & a, bool invert) {
            int n = a.size();
            int lg_n = 0;
            while ((1 << lg_n) < n)
                lg_n++;

            for (int i = 0; i < n; i++) {
                if (i < reverse(i, lg_n))
                    swap(a[i], a[reverse(i, lg_n)]);
            }

            for (int len = 2; len <= n; len <<= 1) {
                double ang = 2 * PI / len * (invert ? -1 : 1);
                cd wlen(cos(ang), sin(ang));
                for (int i = 0; i < n; i += len) {
                    cd w(1);
                    for (int j = 0; j < len / 2; j++) {
                        cd u = a[i+j], v = a[i+j+len/2] * w;
                        a[i+j] = u + v;
                        a[i+j+len/2] = u - v;
                        w *= wlen;
                    }
                }
            }

            if (invert) {
                for (cd & x : a)
                    x /= n;
            }
        }

        vector<int64_t> multiply(vector<int64_t> const& a, vector<int64_t> const& b) {
            vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
            int n = 1;
            while (n < a.size() + b.size()) 
                n <<= 1;
            fa.resize(n);
            fb.resize(n);

            fft(fa, false);
            fft(fb, false);
            for (int i = 0; i < n; i++)
                fa[i] *=  fb[i];

            fft(fa, true);

            vector<int64_t> result(n);

            for (int i = 0; i < n; i++)
                result[i] = (int64_t)llround(fa[i].real());

            return result;
    };
};

struct MinMax_SegmentTree{
    vector<int64_t> ST;
    int n;
    vector<int64_t> arr;
        MinMax_SegmentTree(int n,vector<Int> arr){
            ST.resize(4*n);
            this->arr=arr;
            this->n=n;
        }

        void build(int l,int r,int pos){
            if(l==r){
                ST[pos]=arr[l];
                return;
            }

            int mid=(l+r)/2;
            build(l,mid,2*pos+1);
            build(mid+1,r,2*pos+2);
            ST[pos]=max(ST[2*pos+1],ST[2*pos+2]);

        }

        int RMQ(int l, int r,int x,int y,int pos){
            if(x>=l && y<=r) return ST[pos];
            if(x>r || y<l) return -1e9;

            int mid=(x+y)/2;
            return max(RMQ(l,r,x,mid,2*pos+1),RMQ(l,r,mid+1,y,2*pos+2));
        }

        void update(int l,int r,int pos,int value,int idx){
            if(idx < l || idx>r) return;
            if(l==r){
                arr[idx]=value;
                ST[pos]=value;
            }
            else{
                int mid=(l+r)/2;
                if(idx >=l && idx<=mid){
                    update(l,mid,2*pos+1,value,idx);
                }
                else update(mid+1,r,2*pos+2,value,idx);
                ST[pos]=max(ST[2*pos+1],ST[2*pos+2]);
            }
            return;
        }


};

struct Fenwick{
    vector<long long> BIT;
    int n;

        Fenwick(int n){
            BIT.resize(n);
            BIT.assign(n,0);
            this->n=n;
        }
        void add(int x,int value){
            while(x<n){
                BIT[x]+=value;
                x+=(x&(-x));
            }
        }
        long long sum(int x){
            long long answer=0;
            while(x){
                answer+=BIT[x];
                x-=(x&(-x));
            }
            return answer;
        }
        
        long long rangeSum(int L,int R){
            return sum(R)-sum(L-1);
        }
};

Int countInversions(vector<Int> arr){
    Int inv=0;
    int n=arr.size();
    Fenwick BIT(MAXN);
    for(int i=(int)arr.size()-1;i>=0;--i){
        inv+=BIT.sum(arr[i]-1);
        BIT.add(arr[i],1);
    }
    
    return inv;
}

struct Trie
{
    int value;
    Trie* children[2];

    Trie(){
        value=-1;
        for(int i=0;i<2;++i){
            children[i]=nullptr;
        }
    }

};

Trie* root;
const int NN=2e5+5;
struct DSU{
    vector<int> SZ,ID;
    int n;
    int cc;

    DSU(int n){
        SZ.resize(n);
        ID.resize(n);
        cc=n;
        this->n=n;
        for(int i=0;i<n;++i){
            SZ[i]=1;
            ID[i]=i;
        }
    }

    int find(int x){
        return (ID[x]==x ? x:ID[x]=find(ID[x]));
    }

    bool same(int x,int y){
        return find(x)==find(y);
    }

    void merge(int x,int y){
        x=find(x);y=find(y);
        if(x==y) return;
        if(SZ[x]<SZ[y]) {
            swap(x,y);
        }

        SZ[x]+=SZ[y];
        ID[y]=x;
        --cc;
    }

};

struct FlowEdge {
    int v, u;
    long long cap, flow = 0;
    FlowEdge(int v, int u, long long cap) : v(v), u(u), cap(cap) {}
};

struct Dinic {
    const long long flow_inf = 1e18;
    vector<FlowEdge> edges;
    vector<vector<int>> adj;
    int n, m = 0;
    int s, t;
    vector<int> level, ptr;
    queue<int> q;

    Dinic(int n, int s, int t) : n(n), s(s), t(t) {
        adj.resize(n);
        level.resize(n);
        ptr.resize(n);
    }

    void add_edge(int v, int u, long long cap) {
        edges.emplace_back(v, u, cap);
        edges.emplace_back(u, v, 0);
        adj[v].push_back(m);
        adj[u].push_back(m + 1);
        m += 2;
    }

    bool bfs() {
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int id : adj[v]) {
                if (edges[id].cap - edges[id].flow < 1)
                    continue;
                if (level[edges[id].u] != -1)
                    continue;
                level[edges[id].u] = level[v] + 1;
                q.push(edges[id].u);
            }
        }
        return level[t] != -1;
    }

    long long dfs(int v, long long pushed) {
        if (pushed == 0)
            return 0;
        if (v == t)
            return pushed;
        for (int& cid = ptr[v]; cid < (int)adj[v].size(); cid++) {
            int id = adj[v][cid];
            int u = edges[id].u;
            if (level[v] + 1 != level[u] || edges[id].cap - edges[id].flow < 1)
                continue;
            long long tr = dfs(u, min(pushed, edges[id].cap - edges[id].flow));
            if (tr == 0)
                continue;
            edges[id].flow += tr;
            edges[id ^ 1].flow -= tr;
            return tr;
        }
        return 0;
    }

    long long flow() {
        long long f = 0;
        while (true) {
            fill(level.begin(), level.end(), -1);
            level[s] = 0;
            q.push(s);
            if (!bfs())
                break;
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, flow_inf)) {
                f += pushed;
            }
        }
        return f;
    }
};

template <typename T, class F = function<T(const T&, const T&)>>
class SparseTable {
 public:
  int n;
  vector<vector<T>> mat;
  F func;

  SparseTable(const vector<T>& a, const F& f) : func(f) {
    n = static_cast<int>(a.size());
    int max_log = 32 - __builtin_clz(n);
    mat.resize(max_log);
    mat[0] = a;
    for (int j = 1; j < max_log; j++) {
      mat[j].resize(n - (1 << j) + 1);
      for (int i = 0; i <= n - (1 << j); i++) {
        mat[j][i] = func(mat[j - 1][i], mat[j - 1][i + (1 << (j - 1))]);
      }
    }
  }

  T get(int from, int to) const {
    assert(0 <= from && from <= to && to <= n - 1);
    int lg = 32 - __builtin_clz(to - from + 1) - 1;
    return func(mat[lg][from], mat[lg][to - (1 << lg) + 1]);
  }
};

Int binomialCoeff(int n, int k) 
{ 
    if(k>n) return 0;
    Int res = 1; 
    if (k > n - k) 
        k = n - k; 
  
    for (int i = 0; i < k; ++i) { 
        res =mmul(res,(n - i)); 
        res =mmul(res,ModInv(i + 1)); 
    } 
  
    return res%MOD; 
} 


Int dfs(vector<Int> &arr,int idx,vector<Int> &s){
    if(idx==arr.size()) {
        return accumulate(s.begin(),s.end(),0LL);
    }
    Int ans=1e9;
        for(int i=0;i<s.size();++i){
            s[i]*=arr[idx];
            ans=min(ans,dfs(arr,idx+1,s));
            s[i]/=arr[idx];
        }

        return ans;
}


Int dp[NAX][NAX][NAX];

int main()
{
    
    #ifndef ONLINE_JUDGE
    freopen("inputf.txt","r",stdin);
    freopen("outputf.txt","w",stdout);
    #endif

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);cout.tie(0);

    cout <<"This is just a default template" <<'\n';

    
        

}