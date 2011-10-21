#include <boost/serialization/vector.hpp>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <utility>

class pnt {/*{{{*/
	private:
		friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
				ar & x;
				ar & y;
				ar & z;
				ar & dens;
				ar & isBdry;
				ar & idx;
			}
	public:
		double x, y, z;
		double dens;
		int isBdry;
		int idx;


		pnt(double x_, double y_, double z_, double dens_, int isBdry_, int idx_)
			:  x(x_), y(y_), z(z_), dens(dens_), isBdry(isBdry_), idx(idx_) {	}

		pnt(double x_, double y_, double z_, int isBdry_, int idx_)
			:  x(x_), y(y_), z(z_), dens(1.0), isBdry(isBdry_), idx(idx_) {	}

		pnt(double x_, double y_, double z_)
			: x(x_), y(y_), z(z_), dens(1.0), isBdry(0), idx(0) { }

		pnt(double x_, double y_, double z_, int isBdry_)
			: x(x_), y(y_), z(z_), dens(1.0), isBdry(isBdry_), idx(0) { }

		pnt()
			: x(0.0), y(0.0), z(0.0), dens(1.0), isBdry(0), idx(0) { }

		friend pnt operator*(const double d, const pnt &p);
		friend std::ostream & operator<<(std::ostream &os, const pnt &p);
		friend std::istream & operator>>(std::istream &is, pnt &p);

		pnt& operator=(const pnt &p){/*{{{*/
			x = p.x;
			y = p.y;
			z = p.z;
			isBdry = p.isBdry;
			idx = p.idx;
			dens = p.dens;
			return *this;
		}/*}}}*/
		bool operator==(const pnt &p) const {/*{{{*/
			return (x == p.x) & (y == p.y) & (z == p.z);
		}/*}}}*/
		pnt operator-(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x-p.x;
			y_ = y-p.y;
			z_ = z-p.z;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt operator+(const pnt &p) const {/*{{{*/
			double x_, y_, z_;
			double dens_;

			x_ = x+p.x;
			y_ = y+p.y;
			z_ = z+p.z;
			dens_ = dens+p.dens;

			return pnt(x_,y_,z_,dens_,0,0);
		}/*}}}*/
		pnt operator*(double d) const {/*{{{*/
			double x_, y_, z_;
			double dens_;

			x_ = x*d;
			y_ = y*d;
			z_ = z*d;
			dens_ = dens*d;
			return pnt(x_,y_,z_,dens_,0,0);
		}/*}}}*/
		pnt operator/(double d) const {/*{{{*/
			double x_, y_, z_;
			double dens_;

			if(d == 0.0){
				std::cout << "pnt: operator/" << std::endl;
				std::cout << (*this) << std::endl;
			}

			assert(d != 0.0);
			x_ = x/d;
			y_ = y/d;
			z_ = z/d;
			dens_ = dens/d;
			return pnt(x_,y_,z_,dens_,0,0);
		}/*}}}*/
		pnt& operator/=(double d){/*{{{*/
			if(d == 0.0){
				std::cout << "pnt: operator /=" << std::endl << (*this) << std::endl;
			}
			assert(d != 0.0);
			x = x/d;
			y = y/d;
			z = z/d;
			dens = dens/d;
			return *this;
		}/*}}}*/
		pnt& operator+=(const pnt &p){/*{{{*/
			x += p.x;
			y += p.y;
			z += p.z;
			dens += p.dens;
			return *this;
		}/*}}}*/
		double operator[](int i) const {/*{{{*/
			if(i == 0){
				return x;
			} else if(i == 1){
				return y;
			} else {
				return z;
			}
		}/*}}}*/
		void normalize(){/*{{{*/
			double norm;

			norm = x*x + y*y + z*z;
			if(norm == 0){
				std::cout << "pnt: normalize" << std::endl;
				std::cout << x << " " << y << " " << z << " " << isBdry << " " << idx << std::endl;

				assert(norm != 0);
			}	
			norm = sqrt(norm);

			x = x/norm;
			y = y/norm;
			z = z/norm;
		}/*}}}*/
		double dot(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			return junk;
		}/*}}}*/
		double dotForDistance(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			junk = junk - 1.0;

			return fabs(junk);
		}/*}}}*/
		double dotForAngle(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;
			if(junk > 1.0){
				junk = 1.0;
			}

			if(junk < -1.0){
				junk = -1.0;
			}
			return acos(junk);
		}/*}}}*/
		pnt cross(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = y*p.z - p.y*z;
			y_ = z*p.x - p.z*x;
			z_ = x*p.y - p.x*y;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		double magnitude() const{/*{{{*/
			return sqrt(x*x + y*y + z*z);
		}/*}}}*/
		double getLat() const {/*{{{*/
			return asin(z);			
		}/*}}}*/
		double getLon() const {/*{{{*/
			double lon;

			lon = atan2(y,x);

			if(lon < 0){
				return 2.0 * M_PI + lon;
			} else {
				return lon;
			}
		}/*}}}*/
		double magnitude2() const {/*{{{*/
			return x*x + y*y + z*z;
		}/*}}}*/
	struct hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			uint32_t hash; 
			size_t i, key[3] = { (size_t)p.x, (size_t)p.y, (size_t)p.z };
			for(hash = i = 0; i < sizeof(key); ++i) {
				hash += ((uint8_t *)key)[i];
				hash += (hash << 10);
				hash ^= (hash >> 6);
			}
			hash += (hash << 3);
			hash ^= (hash >> 11);
			hash += (hash << 15);
			return hash;
		}
	};/*}}}*/
	struct idx_hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			return (size_t)p.idx;
		}
	};/*}}}*/
};/*}}}*/
class bdry_pnt {/*{{{*/
	private:
		friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
				ar & x;
				ar & y;
				ar & z;
				ar & lat;
				ar & lon;
				ar & idx;
			}
	public:
		double x, y, z;
		double lat, lon;
		int idx;


		bdry_pnt(double x_, double y_, double z_, double lat_, double lon_, int idx_)
			:  x(x_), y(y_), z(z_), lat(lat_), lon(lon_), idx(idx_) {	}

		bdry_pnt(double x_, double y_, double z_, int idx_)
			: x(x_), y(y_), z(z_), idx(idx_) { }

		bdry_pnt(double lat_, double lon_, int idx_)
			: lat(lat_), lon(lon_), idx(idx_) { }

		bdry_pnt(double x_, double y_, double z_)
			: x(x_), y(y_), z(z_), idx(0) { }

		bdry_pnt(double lat_, double lon_)
			: lat(lat_), lon(lon_), idx(0) { }

		bdry_pnt()
			: x(0.0), y(0.0), z(0.0), lat(0.0), lon(0.0), idx(0) { }

		friend std::ostream & operator<<(std::ostream &os, const bdry_pnt &p);
		friend std::istream & operator>>(std::istream &is, bdry_pnt &p);

		bdry_pnt& operator=(const bdry_pnt &p){/*{{{*/
			x = p.x;
			y = p.y;
			z = p.z;
			lat = p.lat;
			lon = p.lon;
			idx = p.idx;
			return *this;
		}/*}}}*/
		bool operator==(const bdry_pnt &p) const {/*{{{*/
			return (x == p.x) & (y == p.y) & (z == p.z);
		}/*}}}*/
		bdry_pnt operator-(const bdry_pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x-p.x;
			y_ = y-p.y;
			z_ = z-p.z;

			return bdry_pnt(x_,y_,z_);
		}/*}}}*/
		bdry_pnt operator+(const bdry_pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x+p.x;
			y_ = y+p.y;
			z_ = z+p.z;

			return bdry_pnt(x_,y_,z_);
		}/*}}}*/
		bdry_pnt& operator/=(double d){/*{{{*/
			if(d == 0.0){
				std::cout << "pnt: operator /=" << std::endl << (*this) << std::endl;
			}
			assert(d != 0.0);
			x = x/d;
			y = y/d;
			z = z/d;
			return *this;
		}/*}}}*/
		bdry_pnt& operator+=(const bdry_pnt &p){/*{{{*/
			x += p.x;
			y += p.y;
			z += p.z;

			return *this;
		}/*}}}*/
		double operator[](int i) const {/*{{{*/
			if(i == 0){
				return x;
			} else if(i == 1){
				return y;
			} else {
				return z;
			}
		}/*}}}*/
		void normalize(){/*{{{*/
			double norm;

			norm = x*x + y*y + z*z;
			if(norm == 0){
				std::cout << "bdry_pnt: normalize" << std::endl;
				std::cout << x << " " << y << " " << z << " " << idx << std::endl;

				assert(norm != 0);
			}	
			norm = sqrt(norm);

			x = x/norm;
			y = y/norm;
			z = z/norm;
		}/*}}}*/
		double dot(const bdry_pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			return junk;
		}/*}}}*/
		double dotForDistance(const bdry_pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			junk = junk - 1.0;

			return fabs(junk);
		}/*}}}*/
		double dotForAngle(const bdry_pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;
			if(junk > 1.0){
				junk = 1.0;
			}

			if(junk < -1.0){
				junk = -1.0;
			}
			return acos(junk);
		}/*}}}*/
		bdry_pnt cross(const bdry_pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = y*p.z - p.y*z;
			y_ = z*p.x - p.z*x;
			z_ = x*p.y - p.x*y;

			return bdry_pnt(x_,y_,z_);
		}/*}}}*/
		double magnitude() const{/*{{{*/
			return sqrt(x*x + y*y + z*z);
		}/*}}}*/
		inline double getLat() const {/*{{{*/
			return lat;
		}/*}}}*/
		inline double getLon() const {/*{{{*/
			return lon;
		}/*}}}*/
		inline double buildLat(){/*{{{*/
			lat = asin(z);

			fixLat();
		}/*}}}*/
		inline double buildLon(){/*{{{*/
			lon = atan2(y,x);

			fixLon();
		}/*}}}*/
		inline double fixLat(){/*{{{*/
			while(lat < -M_PI){
				lat = lat + M_PI;
			}

			while(lat > M_PI){
				lat = lat - M_PI;
			}
		}/*}}}*/
		inline double fixLon(){/*{{{*/
			while(lon < 0.0){
				lon = lon + 2.0*M_PI;
			}

			while(lon > 2.0*M_PI){
				lon = lon - 2.0*M_PI;
			}
		}/*}}}*/
		void convertToLatLon(){/*{{{*/
			double temp;

			temp = getLat();
			temp = getLon();
		}/*}}}*/
		void convertToCart(){/*{{{*/
			x = cos(lon) * cos(lat);
			y = sin(lon) * cos(lat);
			z = sin(lat);
		}/*}}}*/
		double magnitude2() const {/*{{{*/
			return x*x + y*y + z*z;
		}/*}}}*/
};/*}}}*/
class tri {/*{{{*/
	private:
		friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
				ar & vi1;
				ar & vi2;
				ar & vi3;
				ar & idx;
			}
	public:
	int vi1, vi2, vi3;
	int idx;

	tri() : vi1(0), vi2(0), vi3(0), idx(0) { }

	tri(int vi1_, int vi2_, int vi3_)
		: vi1(vi1_), vi2(vi2_), vi3(vi3_), idx(0) { }

	tri(int vi1_, int vi2_, int vi3_, int idx_)
		: vi1(vi1_), vi2(vi2_), vi3(vi3_), idx(idx_) { }


	friend std::ostream & operator<<(std::ostream &os, const tri &t);
	friend std::istream & operator>>(std::istream &is, tri &t);

	tri sortedTri(){/*{{{*/
		int v1, v2, v3, swp_v;
		v1 = vi1;
		v2 = vi2;
		v3 = vi3;

		//Bubble sort on 3 integers
		for(int i = 0; i < 2; i++){
			if(v1 > v2){
				swp_v = v1;
				v1 = v2;
				v2 = swp_v;
			}
			if(v2 > v3){
				swp_v = v2;
				v2 = v3;
				v3 = swp_v;
			}
		}

		return tri(v1,v2,v3);
	}/*}}}*/
	bool operator==(const tri &t) const {/*{{{*/
		return (vi1 == t.vi1) & (vi2 == t.vi2) & (vi3 == t.vi3);
	}/*}}}*/
	struct hasher {/*{{{*/
		size_t operator()(const tri &t) const {
			uint32_t hash; 
			size_t i, key[3] = { t.vi1, t.vi2, t.vi3 };
			for(hash = i = 0; i < sizeof(key); ++i) {
				hash += ((uint8_t *)key)[i];
				hash += (hash << 10);
				hash ^= (hash >> 6);
			}
			hash += (hash << 3);
			hash ^= (hash >> 11);
			hash += (hash << 15);
			return hash;
		}
	};/*}}}*/
};/*}}}*/
class bdry_line {/*{{{*/
	public:
		std::pair<bdry_pnt, bdry_pnt> end_pts;
		int idx;
		bool SCproj;

		double distanceToPoint(const pnt &Q){/*{{{*/
			pnt A, B;
			pnt C, T;
			pnt Q_p, P;
			double r;

			A = pnt(end_pts.first.x, end_pts.first.y, end_pts.first.z);
			B = pnt(end_pts.second.x, end_pts.second.y, end_pts.second.z);

			if(SCproj){

				C = pnt(0, 0, A.z);
				T = A-C;
				r = T.magnitude();

				T = C;

				Q_p = Q - (Q.dot(C))*C - T;
				P = Q_p;

				P.normalize();
				P = P*r + T;

				return Q.dotForAngle(P);
			} else {
				C = A.cross(B);

				return (M_PI/2.0) - C.dotForAngle(Q);
			}
		}/*}}}*/
		pnt projectedPoint(const pnt &Q){/*{{{*/
			pnt A, B;
			pnt C, T;
			pnt Q_p, P;
			double r;

			A = pnt(end_pts.first.x, end_pts.first.y, end_pts.first.z);
			B = pnt(end_pts.second.x, end_pts.second.y, end_pts.second.z);

			if(SCproj){
				C = pnt(0, 0, A.z);
				T = A-C;
				r = T.magnitude();
				C.normalize();

				T = C;
			} else {
				C = A.cross(B);
				C.normalize();
				T = pnt(0, 0, 0);
				r = 1.0;
			}

			Q_p = Q - (Q.dot(C))*C - T;
			P = Q_p;

			P.normalize();
			P = P*r + T;

			return P;
		}/*}}}*/
		bool triangleOnBoundary(const pnt &a, const pnt &b, const pnt &c){/*{{{*/
			pnt begin, end;
			pnt a_vec, b_vec, c_vec, l_vec;
			double a_sign, b_sign, c_sign;

/*			begin = pnt(end_pts.first.x, end_pts.first.y, end_pts.first.z);
			end = pnt(end_pts.second.x, end_pts.second.y, end_pts.second.z);

			l_vec = end-begin;
			a_vec = a-begin;
			b_vec = b-begin;
			c_vec = c-begin;

			a_vec = a_vec.cross(l_vec);
			b_vec = b_vec.cross(l_vec);
			c_vec = c_vec.cross(l_vec);

			a_sign = a_vec.dot(begin);
			b_sign = b_vec.dot(begin);
			c_sign = c_vec.dot(begin);

			a_sign = a_sign/fabs(a_sign);
			b_sign = b_sign/fabs(b_sign);
			c_sign = c_sign/fabs(c_sign);

			if(a_sign == b_sign == c_sign){
				return false;
			} else {
				return true;
			}*/

			pnt proj_a, proj_b, proj_c;
			double a_dist_begin, b_dist_begin, c_dist_begin;
			double a_dist_end, b_dist_end, c_dist_end;
			double ref_dist;
			bool overlapping;

			begin = pnt(end_pts.first.x, end_pts.first.y, end_pts.first.z);
			end = pnt(end_pts.second.x, end_pts.second.y, end_pts.second.z);

			begin.normalize();
			end.normalize();

			proj_a = projectedPoint(a);
			proj_b = projectedPoint(b);
			proj_c = projectedPoint(c);

			a_dist_begin = proj_a.dotForAngle(begin);
			b_dist_begin = proj_b.dotForAngle(begin);
			c_dist_begin = proj_c.dotForAngle(begin);

			a_dist_end = proj_a.dotForAngle(end);
			b_dist_end = proj_b.dotForAngle(end);
			c_dist_end = proj_c.dotForAngle(end);

			ref_dist = begin.dotForAngle(end);

			overlapping = (a_dist_begin <= ref_dist && a_dist_end <= ref_dist) &&
						  (b_dist_begin <= ref_dist && b_dist_end <= ref_dist) &&
						  (c_dist_begin <= ref_dist && c_dist_end <= ref_dist);

			overlapping = overlapping && 
				          (a_dist_begin >= 0 && a_dist_end >= 0) &&
				          (b_dist_begin >= 0 && b_dist_end >= 0) &&
				          (c_dist_begin >= 0 && c_dist_end >= 0);

			return overlapping;



		}/*}}}*/
};/*}}}*/

inline pnt operator*(const double d, const pnt &p){/*{{{*/
	return pnt(d*p.x, d*p.y, d*p.z, 0, 0);
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const pnt &p){/*{{{*/
	//os << '(' << p.x << ", " << p.y << ", " << p.z << ") bdry=" << p.isBdry << " idx=" << p.idx << ' ';
	os  << std::setiosflags(std::ios::fixed) << std::setprecision(16) << p.x << " " << p.y << " " << p.z;
	return os;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, pnt &p){/*{{{*/
	//is >> p.x >> p.y >> p.z >> p.isBdry >> p.idx;
	is >> p.x >> p.y >> p.z;
	return is;
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const tri &t){/*{{{*/
	//return os << t.vi1 << " " << t.vi2 << " " << t.vi3;
	return os << t.vi1 << " " << t.vi2 << " " << t.vi3;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, tri &t){/*{{{*/
	//return is >> t.vi1 >> t.vi2 >> t.vi3;
	return is >> t.vi1 >> t.vi2 >> t.vi3;
}/*}}}*/

void circumcenter(const pnt &A,const pnt &B,const pnt &C, pnt &cent){/*{{{*/
	double a, b, c;
	double pbc, apc, abp;

	a = (B-C).magnitude2();
	b = (C-A).magnitude2();
	c = (A-B).magnitude2();

	pbc = a*(-a + b + c);
	apc = b*( a - b + c);
	abp = c*( a + b - c);

	cent = (pbc*A + apc*B + abp*C)/(pbc + apc + abp);
	cent.dens = (pbc*A.dens + apc*B.dens + abp*C.dens)/(pbc + apc + abp);

//	cent.dens = (A.dens+B.dens+C.dens)/3.0;
}/*}}}*/
double circumradius(const pnt &A, const pnt &B, const pnt &C){/*{{{*/

	pnt ccenter;  
	pnt ac;

	circumcenter(A,B,C,ccenter);
	ac = A - ccenter;

	return ac.magnitude();

/*
	pnt a(0.0,0.0,0.0,0);
	pnt b(0.0,0.0,0.0,0);
	pnt sub(0.0,0.0,0.0,0);
	pnt cross(0.0,0.0,0.0,0);
	double top, bot;

	a = A-C;
	b = B-C;

	sub = a-b;
	cross = a.cross(b);

	top = a.magnitude() * b.magnitude() * sub.magnitude();
	bot = 2.0 * cross.magnitude();

	return top/bot; // */

/*	dx = a.x - b.x;
	dy = a.y - b.y;
	dz = a.z - b.z;

	sa = sqrt(dx*dx + dy*dy + dz*dz);

	dx = b.x - c.x;
	dy = b.y - c.y;
	dz = b.z - c.z;

	sb = sqrt(dx*dx + dy*dy + dz*dz);

	dx = c.x - a.x;
	dy = c.y - a.y;
	dz = c.z - a.z;

	sc = sqrt(dx*dx + dy*dy + dz*dz);

	dx = (sa+sb+sc)/2.0; // Semiperimeter
	dy = 4.0*sqrt(dx*(sa+sb-dx)*(sa+sc-dx)*(sb+sc-dx)); // Bottom of circumradius computation

	return (sa*sb*sc)/dy;*/
}/*}}}*/
double triArea(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	/**************************************************************************
	 * - This function calculates the area of the triangle A,B,C on the
	 *   surface of a sphere.
	 *
	 *   Input: A, B, C
	 *        A: vertex 1 of triangle
	 *        B: vertex 2 of triangle
	 *        C: vertex 3 of triangle
	 *   Output: (returned value area)
	 *   	area: surface area of triangle on sphere.
	 **************************************************************************/
	pnt u12, u23, u31;
	double a, b, c, s, tanqe, area;	
	double sign;

	area = 0.0;

	//Compute Surface normal for triangle to get "sign" of area
	u12 = B - A;
	u23 = C - A;

	u31 = u12.cross(u23);

	u31.normalize();
	sign = u31.magnitude();
	assert(sign != 0.0);

	//dot the surface norm with one of the points to get sign.
	sign = u31.dot(A);
	sign = sign/fabs(sign);

	a = A.dotForAngle(B);
	b = B.dotForAngle(C);
	c = C.dotForAngle(A);

	s = 0.5*(a+b+c);

	tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)));

	area = sign*4.0*atan(tanqe);
//	area = 4.0*atan(tanqe);

	if(isnan(area))
		area = 0.0;

	return area;
}/*}}}*/
int isCcw(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	double sign;
	pnt ab, ac, cross;

	ab = B-A;
	ac = C-A;
	cross = ab.cross(ac);

	sign = cross.dot(A);
	if(sign > 0){
		return 1;
	} else {
		return 0;
	}
}/*}}}*/
pnt pntFromLatLon(const double &lat, const double &lon){/*{{{*/
	pnt temp;
	temp.x = cos(lon) * cos(lat);
	temp.y = sin(lon) * cos(lat);
	temp.z = sin(lat);
	return temp;
}/*}}}*/
