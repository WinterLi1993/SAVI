#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include <signal.h>
#include "fasd/dcdflib.h"

using namespace std;

unsigned line;

double sum_log(double a, double b) { // use this for logs of probs, i.e. a, b<=0
	if (a == -numeric_limits<double>::infinity()) return b;
	else if (a<b) return b + log((double)1+exp(a-b));
	else return a + log((double)1+exp(b-a));
}

double diff_log(double a, double b) { //
	if (b == -numeric_limits<double>::infinity()) return a;
	else return a + log((double)1-exp(b-a));
}

double and_log(double a, double b) {
	return a + b;
}

double or_log(double a, double b) {
	double s;
	
	s = diff_log(sum_log(a, b), a + b);	
	if (s<=0) return s;
	else return 0;
}

double not_log(double a) {
	return diff_log(0, a);
}

double lognchoosek(double n, double k, double l) {
	double l1;
	
	if (l == 0)
		return 0;
	else if (l == 1)
		return log(n) - log(k);
	else {
		l1 = floor(l/2);
		return lognchoosek(n, k, l1) + lognchoosek(n-l1, k-l1, l-l1);
	}
}

double lognchoosek(vector<double> &f, double n, double k) {
	return f[n] - f[k] - f[n-k];
}

double logbinocdf(vector<double> &f, double n, double k1, double k2, double x) {
	double l1, r1, r2;

	if (x == 0)
		if (k1 == 0) return 0;
		else return -numeric_limits<double>::infinity();
	else if (x == 1) 
		if (k2 == n) return 0;
		else return -numeric_limits<double>::infinity();
		
	if (k1 < 0 || k1 > k2 || k2 > n) return -numeric_limits<double>::infinity();

	if (k1 == k2) {
		return lognchoosek(f, n, k1)  + k1*log(x) + (n-k1)*log(1-x);
	} else {
		l1 = floor((k2-k1+1)/2);
		r1 = logbinocdf(f, n, k1, k1+l1-1, x);
		r2 = logbinocdf(f, n, k1+l1, k2, x);
		return sum_log(r1, r2);
	}
}

double logbinocdf2t(vector<double> &f, double n, double k1, double x) {
	double k2, r1, r2, r3;

	r3 = logbinocdf(f, n, k1, k1, x);
	r1 =  -numeric_limits<double>::infinity();
	for(k2 = 0; k2 <= n; k2++) {
		r2 = logbinocdf(f, n, k2, k2, x);
		if (r2 <= r3) r1 = sum_log(r1, r2);
	}
	return r1;
}


double loghygepdf(vector<double> &f, double a, double n, double b, double t) {
	return lognchoosek(f, a, n) + lognchoosek(f, b, t) - lognchoosek(f, a+b, n+t);
}

double loghygecdf2d(vector<double> &f, double x, double a, double b) {
	double r1, r2;
	unsigned n, t;
	
	r1 = -numeric_limits<double>::infinity();
	for(n = 0; n <= a; n++)
		for(t = 0; t <= b; t++)
			if (abs(b*n - a*t) >= x) {
				r2 = loghygepdf(f, a, n, b, t);
				if (r1 == -numeric_limits<double>::infinity()) r1 = r2;
				else r1 = sum_log(r1, r2);
			}
	return r1 - log(a + b + 1);
}

double loghygecdf1d2t(vector<double> &f, double x, double a, double b, double t) {
	double r1, r2, r3, k;

	r3 = loghygepdf(f, a, x, b, t);
	r1 = -numeric_limits<double>::infinity();	
	for(k = 0; k <= a; k++) {
		r2 = loghygepdf(f, a, k, b, t);
		if (r2 <= r3) r1 = sum_log(r1, r2);
	}
	return r1 + log(b + 1) - log(a + b + 1) ;
	
}

double betafunc(double xx, double xm, double xn, double &w, double &w1) {
	double a, b, x, y;
	int err;
	
	a = xn + 1;
	b = xm - xn + 1;
	x = xx;
	y = 1 - xx; 
	beta_inc(&a, &b, &x, &y, &w, &w1, &err);
}


double intbetafunc1(double xx, double zz, double xm, double xn, double &w1, double &w2) {
	double wx1, wz1, wx2, wz2, w;
	double a, b, q;
	double x, y;
	int err;
	
	a = xn + 1;
	b = xm - xn + 1;

	x = xx;
	y = 1 - xx;	
	beta_inc(&a, &b, &x, &y, &wx1, &w, &err);
	beta_inc(&a, &b, &y, &x, &w, &wx2, &err);
	
	x = zz;
	y = 1 - zz;	
	beta_inc(&a, &b, &x, &y, &wz1, &w, &err);
	beta_inc(&a, &b, &y, &x, &w, &wz2, &err);
	
	w1 = wz1 - wx1;
	w2 = wz2 - wx2;
}

double intbetafunc2(double xx, double zz, double a, double b, double &w1) {
	double wx1, wz1, w;
	double q;
	double x, y;
	int err;
	
	x = xx;
	y = 1 - xx;	
	beta_inc(&a, &b, &x, &y, &wx1, &w, &err);
	
	x = zz;
	y = 1 - zz;	
	beta_inc(&a, &b, &x, &y, &wz1, &w, &err);
	
	w1 = wz1 - wx1;
}


double logfishermult(double p, unsigned n) {
	unsigned i;
	double p1, p2, f;
	
	if (p == -numeric_limits<double>::infinity()) return p;
	
	p1 = 0;
	p2 = 1;	
	f = 1;
	for(i = 0; i < n; i++) {
		p1 += p2/f;
		p2 *= -p;
		f *= i+1;
	}
	p1 = p + log(p1);
	
	return p1;
}

double entr(double x) {
	if (x==0 || x==1) return 0;
	else  return (-x*log(x) - (1-x)*log(1-x))/log(2);
}

double q2p(double s) {
	return exp(-log(10)*s/10);
}

double fixprob(double x) {
	if (x<=1e-14) return 0;
	else return x;
}


struct probs0 {
	double q;
	vector<double> f;
	vector<double> p;
};

istream &operator>>(istream &is, probs0 &p) {
	unsigned k;

	is.read((char *)&k, sizeof(unsigned));
	p.p.resize(k, 0);
	p.f.resize(k, 0);
	is.read((char *)&p.f[0], k*sizeof(double));
	is.read((char *)&p.p[0], k*sizeof(double));
}

struct stats1 {
	bool r;
	double p, m, n;
};

istream &operator>>(istream &is, stats1 &s) {
	is >> s.p >> s.n >> s.m;
	s.p = q2p(s.p);
}

typedef vector<stats1> stats2;

istream &operator>>(istream &is, stats2 &s) {
	bool r;
	string l;
	stringstream ss(stringstream::in | stringstream::out);
	stats1 s1;

	is >> s1.r;
	s.clear();
	if (is) {
		getline(is, l);
		ss << l;
		while (true) {
			ss >> s1;
			if (!ss) break;
			s.push_back(s1);
		}
	}
	
	return is;
}


typedef vector<double> probs1;

void get_probs1(probs0 &p0, stats1 &s, probs1 &p) {
	unsigned i;
	double w, a, b, q, r;
	vector<double> dr;
	
	if (s.p == 1) q = p0.q;
	else q = s.p;
	
	if (s.r) a = s.n;
	else a = s.m - s.n;
	b = s.m - a;
	
	a = a + 1;
	b = b + 1;
	
	p.resize(p0.p.size(), 0);
	
	//cerr << "probs1_1.1: " << log(q) << " " << p0.q << " " << s.p << " " << s.n << " " << s.m << endl;
	
	r = beta_log(&a, &b);
	w = -numeric_limits<double>::infinity();
	for(i = 0; i < p0.p.size(); i++) {
		if (p0.f[i] < 0.5) {
			intbetafunc2(p0.f[i], p0.f[i] + q*(1-2*p0.f[i]), a, b, p[i]);
			p[i] = log(p[i]) - log(q*(1-2*p0.f[i]));
		} else if (p0.f[i] > 0.5) {
			intbetafunc2(p0.f[i] - q*(2*p0.f[i]-1), p0.f[i], a, b, p[i]);
			p[i] = log(p[i]) - log(q*(2*p0.f[i]-1));
		} else p[i] = -log(2)*s.m - r;
		p[i] += p0.p[i];
		w = sum_log(w, p[i]);
	}

	for(i = 0; i < p.size(); i++) p[i] -= w;

//	cerr << "probs1_1.2:";
//	for(i = 0; i < p.size(); i++) cerr << "\n" << p[i];
//	cerr << endl;
}

struct probs3 {
	double p1, p2;
};

void set_probs3_p1(double p, probs3 &r) {
	r.p1 = min(p, 0.0);
	r.p2 = diff_log(0, r.p1);
}

void set_probs3_p2(double p, probs3 &r) {
	r.p2 = min(p, 0.0);
	r.p1 = diff_log(0, r.p2);
}


struct scores1 {
	double d, i;
};

ostream &operator<<(ostream &os, scores1 &s) {
	os << (unsigned)s.i << "\t" << fixprob(1-exp(s.d));
	return os;
}

void get_scores1(probs3 &p, scores1 &s) {
	vector<pair<double, unsigned> > v;

	v.resize(2, pair<double, unsigned>(0, 0));
	
	v[0].first = p.p1; v[0].second = 0;
	v[1].first = p.p2; v[1].second = 1;
	
	sort(v.begin(), v.end());

	s.d = v.back().first;
	s.i = v.back().second;
}

struct probs4 {
	double fr;
	vector<double> p;
};

struct opts_probs4 {
	bool f, f2d;
	unsigned i;
};

void get_probs4(vector<probs0> &p0, opts_probs4 &o, stats2 &s, probs4 &t1) {
	unsigned i, j, k, n, m;
	probs1 t1t1, t1t2;
	double d;

	if (o.f) 
		if (o.f2d) {
			t1.p.resize(p0[o.i].p.size(), 0);	
			fill(t1.p.begin(), t1.p.end(), 0);
			get_probs1(p0[o.i], s[o.i], t1t1);
			for(j = 0; j < t1t1.size(); j++) t1.p[j] = t1t1[j];
			
			if (s[o.i].r) n = s[o.i].n;
			else n = s[o.i].m - s[o.i].n;
			
			t1.fr = (double)n/m;
		} else {
			// assumes evenly spaced frequencies starting w 0. todo: extend this to any sequencies
			get_probs1(p0[0], s[0], t1t1);
			get_probs1(p0[1], s[1], t1t2);
			

			t1.p.resize(2*p0[0].p.size()-1, 0);	
			fill(t1.p.begin(), t1.p.end(), -numeric_limits<double>::infinity());
			for(i = 0; i < t1t1.size(); i++) 
				for(j = 0; j < t1t2.size(); j++) {
					k = p0[0].p.size()-1 + i - j;
					t1.p[k] = sum_log(t1.p[k], t1t1[i] + t1t2[j]);
				}
			
			if (s[0].r) n = s[0].n;
			else n = s[0].m - s[0].n;
			m = s[0].m;
			
			t1.fr = (double)n/m;
			
			if (s[1].r) n = s[1].n;
			else n = s[1].m - s[1].n;
			m = s[1].m;
			
			t1.fr = t1.fr - (double)n/m;
		}
	else {
		n = m = 0;
		t1.p.resize(p0[0].p.size(), 0);	
		fill(t1.p.begin(), t1.p.end(), 0);
		for(i = 0; i < s.size(); i++) {
			get_probs1(p0[0], s[i], t1t1);
			for(j = 0; j < t1t1.size(); j++) t1.p[j] += t1t1[j];
			
			if (s[i].r) n += s[i].n;
			else n += s[i].m - s[i].n;
			
			m += s[i].m;
		}
		
		t1.fr = (double)n/m;		
	}
	
	d  = -numeric_limits<double>::infinity();
	for(i = 0; i < t1.p.size(); i++) d = sum_log(d, t1.p[i]);	
	for(i = 0; i < t1.p.size(); i++) t1.p[i] -= d;
}

void get_probs2_1(vector<double> &q, vector<unsigned> &t3, vector<double> &r) {
	bool f1;
	unsigned i;

	r.clear();
	f1 = false;
	for(i = 0; i < t3.size(); ) {
		if (!f1) { r.push_back(q[t3[i]]); f1 = true; i++; }
		else
			if (t3[i] != t3[i-1]+1) { r.push_back(q[t3[i-1]]); f1 = false; }
			else i++;
	}
	if (t3.size() > 0) r.push_back(q[t3[i-1]]);
}

struct probs2 {
	double h, hd, hm, f1, f2, f3;
	
	vector<double> r1;
	double p1;
	
	vector<double> r2;
	double p2;
	
	scores1 s, d;
};


struct probs0_1 {
	unsigned i;
	double p;
};

bool comp_probs0_1(const probs0_1 &a, const probs0_1 &b) {
	return a.p>b.p;
}

void get_probs2(vector<double> &q, probs4 &t1, vector<bool> &f, vector <double> &c, stats2 &s, probs2 &p) {
	unsigned i, z;
	vector<probs0_1> t2;
	vector<unsigned> t3;
	probs3 t4;
	double d;
	
	p.h = 0;
	p.hd = 0;
	p.hm = numeric_limits<double>::infinity();
	p.f1 = t1.fr;
	p.f2 = 0;
	t2.resize(t1.p.size());
	for(i = 0; i < t1.p.size(); i++)  {
		if (t1.p[i] != -numeric_limits<double>::infinity()) {
			p.h -= t1.p[i]*exp(t1.p[i]);
			p.hd += (t1.p[i]*t1.p[i])*exp(t1.p[i]);
			p.hm = min(p.hm, -t1.p[i]);
		}
		t2[i].i = i;
		t2[i].p = t1.p[i];
		p.f2 += q[i] * exp(t1.p[i]);
	}
	p.h = p.h/log(2);		
	p.hd = sqrt(p.hd/(log(2)*log(2)) - p.h*p.h);
	p.hm = p.hm/log(2);
	sort(t2.begin(), t2.end(), comp_probs0_1);
	p.f3 = q[t2[0].i];
	
	if (f[0] || f[4]) {
		if (f[3]) i = (q.size() + 1)/2 - 1;
		else i = 0;
	
		set_probs3_p1(t1.p[i], t4);	
		get_scores1(t4, p.s);
	}
		
	if (f[1]) {
		t3.clear();
		p.p1 = 0;
		for(i = 0; i < t2.size() && 1-p.p1>c[1]; i++) {
			t3.push_back(t2[i].i);
			p.p1 += exp(t2[i].p);
		}
		p.p1 = 1-p.p1;
		sort(t3.begin(), t3.end());
		
		get_probs2_1(q, t3, p.r1);
	}
	
	if (f[5]) {
		t3.clear();
		p.p1 = 0;
		for(i = 0; i < t2.size() && 1-p.p1>c[5]; i++) {
			t3.push_back(t2[i].i);
			p.p1 += exp(t2[i].p);
		}
		p.p1 = 1-p.p1;
		sort(t3.begin(), t3.end());
		
		get_probs2_1(q, t3, p.r1);
	}
	
	if (f[2]) {
		t3.clear();
		p.p2 = 0;
		for(i = 0; i < t2.size() && -t2[i].p/log(2)-(p.h + c[2]*p.hd)<=1e-14; i++) {
			t3.push_back(t2[i].i);
			p.p2 += exp(t2[i].p);
		}
		p.p2 = 1-p.p2;
		sort(t3.begin(), t3.end());
		
		get_probs2_1(q, t3, p.r2);
	}
}


void get_freqs(probs0 &p0, bool f, vector<double> &q) {
	unsigned i;

	if (f) {
		q.resize(2*p0.f.size() - 1, 0);
		for(i = 0; i < q.size(); i++) 
			if (i < p0.f.size()-1) 
				q[i] = -p0.f[p0.f.size()-1 - i];
			else 
				q[i] = p0.f[i - p0.f.size()+1];
		
	} else {
		q.resize(p0.f.size(), 0);
		for(i = 0; i < q.size(); i++) q[i] = p0.f[i];
	}
}
 

void main_savi_poster(int argc, char **argv) {
	unsigned k, i;
	stats2 s;
	vector<probs0> p0, p1;
	probs4 t1;
	opts_probs4 o;
	vector<double> q;
	double xq;
	
	o.f = false;
	o.f2d = false;
	xq = 0.01;
		
	for(i = 1; i < argc;)
		if (strcmp(argv[i], "-pd")==0) {
			ifstream is1(argv[i+1]), is2(argv[i+2]);	
			
			p0.resize(2);

			is1 >> p0[0];
			is2 >> p0[1];
			o.f = true;			
			
			i += 3;
		} else if (strcmp(argv[i], "-p")==0) { 
			ifstream is1(argv[2]);	

			p0.resize(1);
			is1 >> p0[0];
			
			i += 1;
		} else if (strcmp(argv[i], "-q")==0) {
			xq = atof(argv[i+1]);
			
			i += 2;
		} else if (strcmp(argv[i], "-2d")==0) { 
			o.f2d = true;
			
			i += 1;
		} else i += 1;

	if (o.f2d && !o.f) { cerr << "2d needs two priors." << endl; return; }
		
	if (o.f) {
		p0[0].q = xq;
		p0[1].q = xq;
	} else {
		p0[0].q = xq;
	}	
	
	if (o.f2d) {
		get_freqs(p0[0], false, q);
		k = q.size();
		cout.write((char *)&k, sizeof(unsigned));
		cout.write((char *)&q[0], q.size()*sizeof(double));
			
		get_freqs(p0[1], false, q);
		k = q.size();
		cout.write((char *)&k, sizeof(unsigned));
		cout.write((char *)&q[0], q.size()*sizeof(double));
	} else {
		get_freqs(p0[0], o.f, q);
		k = q.size();
		cout.write((char *)&k, sizeof(unsigned));
		cout.write((char *)&q[0], q.size()*sizeof(double));
	}
	
	line = 1;
	while (true) {
		cin >> s;
		if (!cin) return;
		if (o.f && s.size()<2) { cerr << "line " << k << ": not enough data." << endl; continue; }
		
		if (o.f2d) {
			o.i = 0;
			get_probs4(p0, o, s, t1);
			cout.write((char *)&t1.fr, sizeof(double));
			cout.write((char *)&t1.p[0], t1.p.size()*sizeof(double));
			
			o.i = 1;
			get_probs4(p0, o, s, t1);
			cout.write((char *)&t1.fr, sizeof(double));
			cout.write((char *)&t1.p[0], t1.p.size()*sizeof(double));			
		} else {
			get_probs4(p0, o, s, t1);
		
			cout.write((char *)&t1.fr, sizeof(double));
			cout.write((char *)&t1.p[0], t1.p.size()*sizeof(double));
		}
		line++;
	}
}

void main_savi_conf(int argc, char **argv) {
	unsigned k;
	stats2 s;
	probs2 p1;
	probs4 t1;
	unsigned i;
	vector<bool> f;
	vector<double> c, q;
	
	f.resize(5, false);
	c.resize(6, 0);
	for(i = 1; i < argc;) 
		if (strcmp(argv[i], "-z") == 0) {
			f[0] = true;
			
			i += 1;
		} else if (strcmp(argv[i], "-fc") == 0) {
			f[1] = true;
			c[1] = atof(argv[i+1]);
			
			i += 2;
		} else if (strcmp(argv[i], "-fh") == 0) {
			f[2] = true;
			c[2] = atof(argv[i+1]);
			
			i += 2;
		} else if (strcmp(argv[i], "-s") == 0) {
		        f[4] = true;
			c[4] =  atof(argv[i+1]);

			i += 2;
		} else if (strcmp(argv[i], "-fs") == 0) {
			f[5] = true;
			c[5] = atof(argv[i+1]);
			c[6] = atof(argv[i+2]);
			
			i += 3;
		} else i += 1;
	
	cin.read((char *)&k, sizeof(unsigned));
	q.resize(k, 0);
	cin.read((char *)&q[0], k*sizeof(double));
	
	t1.p.resize(k, 0);
	while (true) {
		cin.read((char *)&t1.fr, sizeof(double));
		cin.read((char *)&t1.p[0], k*sizeof(double));
		
		if (!cin) return;
		
		get_probs2(q, t1, f, c, s, p1);
		
		if (f[4]) {
		        cout << fixed << setprecision(c[4]) << 100*p1.f1 << "\t" << 100*p1.f2 << "\t" << 100*p1.f3;
		} else if (f[5]) {
		        cout << fixed << setprecision(c[6]) << 100*p1.f1 << "\t" << 100*p1.f2 << "\t" << 100*p1.f3;
		} else {
		     cout << fixed << setprecision(0) << 100*p1.f1 << "\t" << 100*p1.f2 << "\t" << 100*p1.f3;
		}		

		if (f[0] || f[4]) cout << scientific << setprecision(1) << "\t" << p1.s;
		
		if (f[1]) {
		        cout << fixed << setprecision(0) << "\t" << p1.r1.size();
			for(i = 0; i < p1.r1.size(); i++) cout << "\t" << 100*p1.r1[i];		

			cout << scientific << setprecision(1) << "\t" << fixprob(p1.p1);
		}
		if (f[5]) {
		        cout << fixed << setprecision(c[6]) << "\t" << p1.r1.size();
			for(i = 0; i < p1.r1.size(); i++) cout << "\t" << 100*p1.r1[i];		

			cout << scientific << setprecision(1) << "\t" << fixprob(p1.p1);
		}
		
		if (f[2]) {		
			cout << fixed << setprecision(2) << "\t" << p1.hm << "\t" << p1.h << "\t" <<  p1.hd;
		
			cout << setprecision(0) << "\t" << p1.r2.size();
			for(i = 0; i < p1.r2.size(); i++) cout << "\t" << 100*p1.r2[i];

			cout << scientific << setprecision(1) << "\t" << fixprob(p1.p2);
		}		
		
		cout << endl;
	}
}

void main_savi_comp(int argc, char **argv) {
	unsigned i, k;
	probs2 p1;
	probs4 t1;
	probs3 t2;
	scores1 s;
	unsigned f;
	double c, p;
	vector<double> q;
	
	f = atoi(argv[1]);
	c = atof(argv[2]);
		
	cin.read((char *)&k, sizeof(unsigned));
	q.resize(k, 0);
	cin.read((char *)&q[0], k*sizeof(double));
	
	t1.p.resize(k, 0);
	while (true) {
		cin.read((char *)&t1.fr, sizeof(double));
		cin.read((char *)&t1.p[0], k*sizeof(double));
		
		if (!cin) return;

		p = -numeric_limits<double>::infinity();
		switch (f) {
			case 0:
				for(i = 0; i < q.size(); i++)
					if (q[i] <= c) p = sum_log(p, t1.p[i]);
				break;
				
			case 1:
				for(i = 0; i < q.size(); i++)
					if (q[i] != c) p = sum_log(p, t1.p[i]);
				break;				
				
			case 2:
				for(i = 0; i < q.size(); i++)
					if (q[i] >= c) p = sum_log(p, t1.p[i]);
				break;		
		}
		set_probs3_p2(p, t2);	
		get_scores1(t2, s);
		
		cout << scientific << setprecision(1) << s << endl;
	}
}

void main_savi_poster_accum(int argc, char **argv) {
	unsigned i, j, k, k1, k2;
	probs4 t1, t2;
	double c;
	bool f2d;
	vector<double> r, q1, q2;
		
	f2d = false;
	for(i = 1; i < argc;)
		if (strcmp(argv[i], "-2d")==0) { 
			f2d = true;
			
			i += 1;
		} else i += 1;
	
	cin.read((char *)&k1, sizeof(unsigned));
	q1.resize(k1, 0);
	cin.read((char *)&q1[0], k1*sizeof(double));
	
	if (f2d) {
		cin.read((char *)&k2, sizeof(unsigned));
		q2.resize(k2, 0);
		cin.read((char *)&q2[0], k2*sizeof(double));		
		
		t1.p.resize(k1, 0);	
		t2.p.resize(k2, 0);	
	
		r.resize(k1*k2, -numeric_limits<double>::infinity());		
	} else {
		t1.p.resize(k1, 0);	
	
		r.resize(k1, -numeric_limits<double>::infinity());
	}
	
	line = 1;
	while (true) {
		cin.read((char *)&t1.fr, sizeof(double));
		if (!cin) break;
		cin.read((char *)&t1.p[0], k1*sizeof(double));

		if (f2d) {
			cin.read((char *)&t2.fr, sizeof(double));
			cin.read((char *)&t2.p[0], k2*sizeof(double));
			
			for(i = 0; i < k1; i++)
				for(j = 0; j < k2; j++) {
					k = k1*i+j;
					r[k] = sum_log(r[k], t1.p[i] + t2.p[j]);
				}
		} else	
			for(i = 0; i < k1; i++)
				r[i] = sum_log(r[i], t1.p[i]);
		line++;
	}

	cout.write((char *)&k1, sizeof(unsigned));
	cout.write((char *)&q1[0], q1.size()*sizeof(double));
	if (f2d) {
		cout.write((char *)&k2, sizeof(unsigned));
		cout.write((char *)&q2[0], q2.size()*sizeof(double));
	}
	cout.write((char *)&r[0], r.size()*sizeof(double));
}

void main_savi_poster_merge(int argc, char **argv) {
	unsigned i, j, k, m1;
	probs4 t1;
	double c;
	vector<double> r1, r, q;
	bool fl;
	
	fl = false;
	for(j = 1; j < argc; j++) {
		ifstream is(argv[j]);
		
		//todo: error check
		is.read((char *)&k, sizeof(unsigned));
		q.resize(k, 0);
		is.read((char *)&q[0], k*sizeof(double));
	
		if (j == 1) {
			r.resize(k, -numeric_limits<double>::infinity());	
			r1.resize(k, 0);
			fl = true;
		}

		is.read((char *)&r1[0], k*sizeof(double));		
		for(i = 0; i < r1.size(); i++)
			r[i] = sum_log(r[i], r1[i]);
	}

	if (fl) {
		cout.write((char *)&k, sizeof(unsigned));
		cout.write((char *)&q[0], q.size()*sizeof(double));		
		cout.write((char *)&r[0], r.size()*sizeof(double));
	}
}

void main_savi_poster_prt(int argc, char **argv) {
	unsigned i, j, k, k1, k2;
	double c1, m;
	vector<double> r, q1, q2;
	bool fs, f2d, fm;
		
	fs = false;
	f2d = false;
	fm = false;
	c1 = 1; m = 1;
	
	for(i = 1; i < argc;) 
		if (strcmp(argv[i], "-s") == 0) {
			fs = true;
			
			i += 1;
		} else if (strcmp(argv[i], "-2d")==0) { 
			f2d = true;
			
			i += 1;
		} else if (strcmp(argv[i], "-m")==0) {
			fm = true;
			c1 = atof(argv[i+1]);

			i += 2;
		} else i += 1;
		
		
	cin.read((char *)&k1, sizeof(unsigned));
	q1.resize(k1, 0);
	cin.read((char *)&q1[0], k1*sizeof(double));		
	
	if (f2d) {
		cin.read((char *)&k2, sizeof(unsigned));
		q2.resize(k2, 0);
		cin.read((char *)&q2[0], k2*sizeof(double));	
		
		r.resize(k1*k2, 0);
		cin.read((char *)&r[0], k1*k2*sizeof(double));
	} else {
		r.resize(k1, 0);
		cin.read((char *)&r[0], k1*sizeof(double));
	}
	
	if (fm) {
		m = 0;
		for(i = 0; i < r.size(); i++)  m += exp(r[i]);
	}

	
	if (fs) cout << scientific;	
	if (f2d)
		for(i = 0; i < k1; i++)
			for(j = 0; j < k2; j++) {
				cout << q1[i] << "\t" << q2[j] << "\t";
				k = k1*i + j;
				if (fs) cout << (c1*exp(r[k])/m) << endl;
				else cout << (unsigned int)(c1*exp(r[k])/m) << endl;
			}
	else
		for(i = 0; i < k1; i++) {
			cout << q1[i] << "\t";
			if (fs) cout << (c1*exp(r[i])/m) << endl;
			else cout << (unsigned int)(c1*exp(r[i])/m) << endl;
		}
}

void main_savi_unif_prior(int argc, char **argv) {
	unsigned i, k;
	double d;
	vector<double> r, q;

	d = atof(argv[1]);
	k = 1 + floor(1/d);
	q.resize(k, 0);
	r.resize(k, 0);
	for(i = 0; i < k; i++) q[i] = i*d;
	cout.write((char *)&k, sizeof(unsigned));
	cout.write((char *)&q[0], k*sizeof(double));		
	cout.write((char *)&r[0], k*sizeof(double));
}

void main_savi_txt2prior(int argc, char **argv) {
	unsigned k;
	double f, p;
	vector<double> q, r;

	while (true) {
		cin >> f;
		if (!cin) break;
		cin >> p;
		
		q.push_back(f);
		r.push_back(log(p));
	}
	
	k = q.size();
	cout.write((char *)&k, sizeof(unsigned));
	cout.write((char *)&q[0], q.size()*sizeof(double));		
	cout.write((char *)&r[0], r.size()*sizeof(double));
}

void sighandle(int sig) {
	cerr << "hi: " << line << endl;
}

int main(int argc, char **argv) 
{ 
	signal(SIGUSR1, &sighandle);
	
	main_savi_conf(argc, argv);	
  	return 0;
}
