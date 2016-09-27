
/*
   pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

   This file is part of pbrt.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
	met:

	- Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.

	- Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
	IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
	TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
	PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


// shapes/heightfield.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "paramset.h"

#define eps 1e-6
#define VERT(x,y) ((x)+(y)*nx)

// [sx, ex] and [sy, ey]
int Heightfield2::fillLR(int sx, int ex, int sy, int ey){
	int id = num_of_node;
	num_of_node ++;

	Tsx[id] = sx;
	Tex[id] = ex;
	Tsy[id] = sy;
	Tey[id] = ey;

	if((ex - sx) * (ey - sy) < 10 || (ex - sx) * (ey - sy) < nx * ny / 30){
		isleaf[id] = 1;

		float zmax = -INFINITY;
		float zmin = INFINITY;
		L[id] = -1;
		R[id] = -1;
		
		for(int i = sx; i <= ex; i++){
			for(int j = sy; j <= ey; j++){
				zmax = max(zmax, z[VERT(i, j)]);
				zmin = min(zmin, z[VERT(i, j)]);
			}
		}

		float invnxm1 = 1. / (nx - 1);
		float invnym1 = 1. / (ny - 1);

		mybox[id].pMin = Point(sx * invnxm1, sy * invnym1, zmin - eps);
		mybox[id].pMax = Point(ex * invnxm1, ey * invnym1, zmax + eps);
	}
	else{
		isleaf[id] = 0;
		
		int Lsx, Lex, Lsy, Ley;
		int Rsx, Rex, Rsy, Rey;

		if(ex - sx > ey - sy){
			int mx = (sx + ex) / 2;
			Lsx = sx; Lex = mx; Lsy = sy; Ley = ey;
			Rsx = mx; Rex = ex; Rsy = sy; Rey = ey;
		}
		else{
			int my = (sy + ey) / 2;
			Lsx = sx; Lex = ex; Lsy = sy; Ley = my;
			Rsx = sx; Rex = ex; Rsy = my; Rey = ey;
		}

		L[id] = fillLR(Lsx, Lex, Lsy, Ley);
		R[id] = fillLR(Rsx, Rex, Rsy, Rey);
		mybox[id] = Union(mybox[L[id]], mybox[R[id]]);
	}

	return id;
}

// Heightfield Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
		bool ro, int x, int y, const float *zs)
	: Shape(o2w, w2o, ro) {
		nx = x;
		ny = y;
		z = new float[nx*ny];
		memcpy(z, zs, nx*ny*sizeof(float));

		boxminz = (float*)malloc(sizeof(float) * nx * ny);
		boxmaxz = (float*)malloc(sizeof(float) * nx * ny);
		
		for(int i = 0; i < nx - 1; i++){
			for(int j = 0; j < ny - 1; j++){
				float &minz = boxminz[VERT(i, j)];
				float &maxz = boxmaxz[VERT(i, j)];

				minz = min(z[VERT(i, j)], z[VERT(i+1, j)]);
				if(minz > z[VERT(i, j+1)]) minz = z[VERT(i, j+1)];
				if(minz > z[VERT(i+1, j+1)]) minz = z[VERT(i+1, j+1)];
				minz -= 1e-6;

				maxz = max(z[VERT(i, j)], z[VERT(i+1, j)]);
				if(maxz < z[VERT(i, j+1)]) maxz = z[VERT(i, j+1)];
				if(maxz < z[VERT(i+1, j+1)]) maxz = z[VERT(i+1, j+1)];
				maxz += 1e-6;
			}
		}

		L = (int*)malloc(sizeof(int) * (100));
		R = (int*)malloc(sizeof(int) * (100));

		Tsx = (int*)malloc(sizeof(int) * (100));
		Tsy = (int*)malloc(sizeof(int) * (100));
		Tex = (int*)malloc(sizeof(int) * (100));
		Tey = (int*)malloc(sizeof(int) * (100));
		
		isleaf = (bool*)malloc(sizeof(bool) * (100));

		mybox = (BBox*)malloc(sizeof(BBox) * (100));
	
		num_of_node = 0;
		fillLR(0, nx-1, 0, ny-1);

		printf("%d\n", num_of_node);

		printf("%f %f %f %f %f %f\n", mybox[0].pMin.x, mybox[0].pMin.y, mybox[0].pMin.z
									, mybox[0].pMax.x, mybox[0].pMax.y, mybox[0].pMax.z);
	}


Heightfield2::~Heightfield2() {
	delete[] z;
}


BBox Heightfield2::ObjectBound() const {
	float minz = z[0], maxz = z[0];
	for (int i = 1; i < nx*ny; ++i) {
		if (z[i] < minz) minz = z[i];
		if (z[i] > maxz) maxz = z[i];
	}
	return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::CanIntersect() const {
	return true;
}


bool Heightfield2::triangleIntersect(Point p1, Point p2, Point p3, const Ray &ray, float *tHit) const {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector s = ray.o - p1;
	float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	Vector s2 = Cross(s, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
		return false;

	*tHit = t;
	return true;
}


void Heightfield2::gridIntersect(Ray ray, int i, int j, bool &hitsomething, float *tHit, DifferentialGeometry *dg) const{
	float t;
	DifferentialGeometry tdg;
	const Transform *o2w = ObjectToWorld;
	Point p1, p2, p3;
	float invnxm1 = 1. / (nx - 1);
	float invnym1 = 1. / (ny - 1);

	p1 = Point(i * invnxm1, j * invnym1, z[VERT(i, j)]);
	p2 = Point((i+1) * invnxm1, j * invnym1, z[VERT(i+1, j)]);
	p3 = Point((i+1) * invnxm1, (j+1) * invnym1, z[VERT(i+1, j+1)]);

	if(triangleIntersect(p1, p2, p3, ray, &t)){
		hitsomething = true;
		if(*tHit > t){
			*tHit = t;
			Point phit = ray(t);
			*dg = DifferentialGeometry((*o2w)(phit),
					(*o2w)((p2-p1)*(nx-1)), (*o2w)((p3-p2)*(ny-1)),
					Normal(0,0,0), Normal(0,0,0),
					phit.x, phit.y, this);
		}
	}

	// Based on previous p1, p2, p3
	p2 = p3;
	p3 = Point(i * invnxm1, (j+1) * invnym1, z[VERT(i, j+1)]);

	if(triangleIntersect(p1, p2, p3, ray, &t)){
		hitsomething = true;
		if(*tHit > t){
			*tHit = t;
			Point phit = ray(t);
			*dg = DifferentialGeometry((*o2w)(phit),
					(*o2w)((p2-p3)*(nx-1)), (*o2w)((p3-p1)*(ny-1)),
					Normal(0,0,0), Normal(0,0,0),
					phit.x, phit.y, this);
		}
	}
}


// [sx, ex] and [sy, ey]
bool Heightfield2::KDTreeIntersect(Ray Gray, bool &hitsomething, float *tHit, DifferentialGeometry *dg) const{
	
	float invdx = 1. / Gray.d.x, invdy = 1. / Gray.d.y, invdz = 1. / Gray.d.z;
	float invnxm1 = 1. / (nx - 1), invnym1 = 1. / (ny - 1);	
	
	int stack[100], top = 0;
	float Stmin[100], Stmax[100];

	int id = 0;
	float tmin, tmax;
	if(mybox[0].IntersectP(Gray, &tmin, &tmax) == false)
		return false;

	while(1){
		int sx = Tsx[id], sy = Tsy[id], ex = Tex[id], ey = Tey[id];

		if(isleaf[id]){
			float hitt0, hitt1;
			if(mybox[id].IntersectP(Gray, &hitt0, &hitt1) == true){
				Ray ray = Gray;
				ray.mint = hitt0;
				ray.maxt = hitt1;

				int x_going = 1, y_going = 1, outi = ex, outj = ey;
				if(ray.d.x < 0) x_going = -1, outi = sx - 1;
				if(ray.d.y < 0) y_going = -1, outj = sy - 1;
	
				Point startP = ray(ray.mint);
				int i = (int)((startP.x + x_going * eps) * (nx - 1));
				int j = (int)((startP.y + y_going * eps) * (ny - 1));
	
				float txCur = ((i + (x_going+1)/2) * invnxm1 - ray.o.x) * invdx;
				float tyCur = ((j + (y_going+1)/2) * invnym1 - ray.o.y) * invdy;
	
				float txAdd = invnxm1 * invdx * x_going;
				float tyAdd = invnym1 * invdy * y_going;
		
				txCur -= txAdd;
				tyCur -= tyAdd;
	
				// Grid DDA
				while(1){
					float txNxt = txCur + txAdd;
					float tyNxt = tyCur + tyAdd;
	
					float tNear = (boxmaxz[VERT(i, j)] - ray.o[2]) * invdz;
					float  tFar = (boxminz[VERT(i, j)] - ray.o[2]) * invdz;
					if(tNear > tFar) swap(tNear, tFar);
	
					if(txNxt < tyNxt){
						if(max(max(txCur, tyCur), tNear) < min(txNxt, tFar))
							gridIntersect(ray, i, j, hitsomething, tHit, dg);
						if(hitsomething) break;
		
						if(txNxt > ray.maxt) break;
						i += x_going;
						txCur = txNxt;
		
						if(i == outi) break;
					}
					else{
						if(max(max(txCur, tyCur), tNear) < min(tyNxt, tFar))
							gridIntersect(ray, i, j, hitsomething, tHit, dg);
						if(hitsomething) break;
			
						if(tyNxt > ray.maxt) break;
						j += y_going;
						tyCur = tyNxt;
		
						if(j == outj) break;
					}
				}
			
				if(hitsomething) return true;
			}

			if(top > 0){
				--top;
				id = stack[top];
				tmin = Stmin[top];
				tmax = Stmax[top];
			}
			else break;
		}
		else{
			// Split x or y
			float tmid;
			bool belowfirst;
			
			if(ex - sx > ey - sy){
				float mpos = ((sx + ex) / 2) * invnxm1;
				tmid = (mpos - Gray.o.x) * invdx;

				belowfirst = (Gray.o.x <  mpos) ||
							 (Gray.o.x == mpos && Gray.d.x >= 0);
			}
			else{
				float mpos = ((sy + ey) / 2) * invnym1;
				tmid = (mpos - Gray.o.y) * invdy;

				belowfirst = (Gray.o.y <  mpos) ||
							 (Gray.o.y == mpos && Gray.d.y >= 0);
			}
			
			int firstC, secondC;
			if(belowfirst)
				firstC = L[id], secondC = R[id];
			else
				firstC = R[id], secondC = L[id];

			if(tmid > tmax || tmid <= 0){
				id = firstC;
			}
			else if(tmid < tmin)
				id = secondC;
			else{
				stack[top] = secondC;
				Stmin[top] = tmid;
				Stmax[top] = tmax;
				top++;
				
				id = firstC;
				tmax = tmid;
			}
		}
	}
	
	return false;
}


bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg) const {

	bool hitsomething = false;

	Ray ray = (*WorldToObject)(r);
	*tHit = ray.maxt + 1;

	KDTreeIntersect(ray, hitsomething, tHit, dg);
	if(!hitsomething) return false;

	*rayEpsilon = 1e-3f * *tHit;
	return true;
}


bool Heightfield2::IntersectP(const Ray &r) const {
	bool hitsomething = false;
	float tHit;
	DifferentialGeometry dg;
	
	return KDTreeIntersect((*WorldToObject)(r), hitsomething, &tHit, &dg);
}


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
		bool reverseOrientation, const ParamSet &params) {
	int nu = params.FindOneInt("nu", -1);
	int nv = params.FindOneInt("nv", -1);
	int nitems;
	const float *Pz = params.FindFloat("Pz", &nitems);
	Assert(nitems == nu*nv);
	Assert(nu != -1 && nv != -1 && Pz != NULL);
	return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}
