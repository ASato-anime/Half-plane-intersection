#include <bits/stdc++.h>
#define all(x) x.begin(), x.end()
#define sz(x) int(x.size())
#define pb push_back
#define rep(i, a, b) for(int i = a; i < b; i++)
using namespace std;
using ld = long double;//change to double if Time Limit

template<class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T> struct Point{
	T x, y;
	Point(T x = 0, T y = 0): x(x), y(y) {}
	using P = Point;
	P operator + (P p) const { return P(x + p.x, y + p.y); }
	P operator - (P p) const { return P(x - p.x, y - p.y); }
	P operator * (T d) const { return P(x * d, y * d); }
	P operator / (T d) const { return P(x / d, y / d); }
	bool operator == (P p) const { return tie(x, y) == tie(p.x, p.y); }
	bool operator < (P p) const { return tie(x, y) < tie(p.x, p.y); }
	T cross(P p) const { return x * p.y - y * p.x; }
	T cross(P a, P b) const { return (a - *this).cross(b - *this); }
	T dot(P p) const { return x * p.x + y * p.y; }
	T dot(P a, P b) const { return (a - *this).cross(b - *this); }
	T dist2() const { return dot(*this); }
	ld dist() const { return sqrtl(dist2()); }
	ld angle() const { return atan2l(y, x); }
	P perp() const { return P(-y, x); }
	friend ostream& operator << (ostream& os, P p) {
		return os << '(' << p.x << ", " << p.y << ')';
	}
};
template<class P> ld segDist(P& s, P& e, P& p) {
	if (s==e) return (p-s).dist();
	auto d = (e-s).dist2(), t = min(d,max(ld(0),(p-s).dot(e-s)));
	return ((p-s)*d-(e-s)*t).dist()/d;
}
const ld eps = 1e-8;
template<class P> bool onSegment(P s, P e, P p) { return segDist(s,e,p) <= eps; }
template<class P> bool inPolygon(vector<P> &p, P a, bool strict = true) {
	int cnt = 0, n = sz(p);
	rep(i,0,n) {
		P q = p[(i + 1) % n];
		if (onSegment(p[i], q, a)) return !strict;
		cnt ^= ((a.y<p[i].y) - (a.y<q.y)) * a.cross(p[i], q) > 0;
	}
	return cnt;
}
template<class P> vector<P> segInter(P a, P b, P c, P d) {
	auto oa = c.cross(d, a), ob = c.cross(d, b),
	oc = a.cross(b, c), od = a.cross(b, d);
	if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
		return {(a * ob - b * oa) / (ob - oa)};
	set<P> s;
	if (onSegment(c, d, a)) s.insert(a);
	if (onSegment(c, d, b)) s.insert(b);
	if (onSegment(a, b, c)) s.insert(c);
	if (onSegment(a, b, d)) s.insert(d);
	return {all(s)};
}
using P = Point<ld>;

struct HP {
	P p, pq;
	ld angle;
	HP() {}
	HP(P a, P b) : p(a), pq(b - a) { angle = pq.angle(); }
	bool out(P r) { return pq.cross(r - p) < -eps; }
	bool operator<(HP e)const{return angle<e.angle;}
	P operator ^ (HP t){
		ld alpha = (t.p - p).cross(t.pq) / pq.cross(t.pq);
		return p + (pq * alpha);
	}
};
vector<P> intersect(vector<HP>& H){//with box
	sort(all(H));
	deque<HP> dq;
	rep(i, 0, sz(H)) {
		while(sz(dq) > 1 && H[i].out(dq.back() ^ dq[sz(dq)-2]))
			dq.pop_back();
		while(sz(dq) > 1 && H[i].out(dq[0] ^ dq[1]))
			dq.pop_front();
		if(sz(dq) && fabsl(H[i].pq.cross(dq.back().pq)) < eps){
			if(H[i].pq.dot(dq.back().pq) < 0)return{};
			if(H[i].out(dq.back().p)) dq.pop_back();
			else continue;
		}
		dq.pb(H[i]);
	}
	while(sz(dq) > 2 && dq[0].out(dq.back() ^ dq[sz(dq) - 2]))
		dq.pop_back();
	while(sz(dq) > 2 && dq.back().out(dq[0] ^ dq[1]))
		dq.pop_front();
	if (sz(dq) < 3) return {};
	vector<P> cht(sz(dq));
	rep(i, 0, sz(dq)) cht[i] = dq[i] ^ dq[(i + 1) % sz(dq)];
	return cht;
}

main(){
	ios_base :: sync_with_stdio(0);
	cin.tie(0); cout.tie(0);
	int n;
	cin >> n;
	vector<P> a(n);
	for(int i = 0; i < n; i++) cin >> a[i].x >> a[i].y;
	ld ans = 0;
	for(int i = 0; i < n; i++){
		P box[] = {
			P(-1e5, -1e5),
			P(1e5, -1e5),
			P(1e5, 1e5),
			P(-1e5, 1e5)
		};
		vector<HP> vec;
		for(int j = 0; j < 4; j++) vec.pb(HP(box[j], box[(j + 1) % 4]));
		for(int j = 0; j < n; j++){
			if(i == j) continue;
			P p = (a[i] + a[j]) / 2;
			vec.pb(HP(p, p + (a[j] - a[i]).perp()));
		}
		vector<P> ch = intersect(vec);
		for(int j = 0; j < sz(ch); j++){
			if(inPolygon(a, ch[j])) ans = max(ans, (ch[j] - a[i]).dist());
			for(int k = 0; k < n; k++){
				vector<P> inter = segInter(a[k], a[(k + 1) % n], ch[j], ch[(j + 1) % sz(ch)]);
				if(sz(inter) == 1) ans = max(ans, (a[i] - inter[0]).dist());
			}
		}
	}
	cout << fixed << setprecision(9) << ans;
}
