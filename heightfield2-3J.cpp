
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

// Heightfield Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
		bool ro, int x, int y, const float *zs)
	: Shape(o2w, w2o, ro) {
		nx = x;
		ny = y;
		z = new float[nx*ny];
		memcpy(z, zs, nx*ny*sizeof(float));	

		float GLminz = z[0], GLmaxz = z[0];
		for (int i = 1; i < nx*ny; ++i) {
			if (z[i] < GLminz) GLminz = z[i];
			if (z[i] > GLmaxz) GLmaxz = z[i];
		}
		mybox = BBox(Point(0, 0, GLminz-1e-6), Point(1, 1, GLmaxz+1e-6));

		boxminz = (float*)malloc(sizeof(float) * nx * ny);
		boxmaxz = (float*)malloc(sizeof(float) * nx * ny);

		float lenx = 1.f / (nx - 1), leny = 1.f / (ny - 1);
		
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

		bigbox = (BBox*)malloc(sizeof(BBox) * nx * ny);

		for(int i = 0; i < nx - 3; i++){
			for(int j = 0; j < ny - 3; j++){
				float minz = min(boxminz[VERT(i, j)], boxminz[VERT(i+2, j+2)]);
				minz = min(minz, boxminz[VERT(i+2, j)]);
				minz = min(minz, boxminz[VERT(i, j+2)]);

				float maxz = max(boxmaxz[VERT(i, j)], boxmaxz[VERT(i+2, j+2)]);
				maxz = max(maxz, boxmaxz[VERT(i+2, j)]);
				maxz = max(maxz, boxmaxz[VERT(i, j+2)]);

				bigbox[VERT(i, j)] = BBox(Point(    i * lenx,     j * leny, minz),
										  Point((i+3) * lenx, (j+3) * leny, maxz));
			}
		}
	}


Heightfield2::~Heightfield2() {
	delete[] z;
	free(boxminz);
	free(boxmaxz);
	free(bigbox);
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


int out_cnt = 0;
int in_cnt = 0;
int inin_cnt = 0;

bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg) const {
	bool hitsomething = false;

	Ray ray = (*WorldToObject)(r);
	*tHit = ray.maxt + 1;

	float hitt0, hitt1;
	if(mybox.IntersectP(ray, &hitt0, &hitt1) == false)
		return false;

	ray.mint = hitt0;
	ray.maxt = hitt1;

	int x_going = 1, y_going = 1, outi = nx-1, outj = ny-1;
	if(ray.d.x < 0) x_going = -1, outi = -1;
	if(ray.d.y < 0) y_going = -1, outj = -1;

	Point startP = ray(ray.mint);
	int i = (int)((startP.x + x_going * eps) * (nx - 1));
	int j = (int)((startP.y + y_going * eps) * (ny - 1));

	float invdx = 1. / ray.d.x, invdy = 1. / ray.d.y, invdz = 1. / ray.d.z;
	float invnxm1 = 1. / (nx - 1), invnym1 = 1. / (ny - 1);

	float txCur = ((i + (x_going+1)/2) * invnxm1 - ray.o.x) * invdx;
	float tyCur = ((j + (y_going+1)/2) * invnym1 - ray.o.y) * invdy;

	float txAdd = invnxm1 * invdx * x_going;
	float tyAdd = invnym1 * invdy * y_going;

	float inv_txAdd = 1.f / txAdd;
	float inv_tyAdd = 1.f / tyAdd;

	int shiftx = 0;
	if(x_going < 0) shiftx = -2;

	int shifty = 0;
	if(y_going < 0) shifty = -2;

	txCur -= txAdd;
	tyCur -= tyAdd;

	out_cnt ++;

	while(1){
		in_cnt ++;

		if(i + shiftx >= 0 && i + shiftx < nx - 3 &&
		   j + shifty >= 0 && j + shifty < ny - 3 &&
		   !bigbox[VERT(i + shiftx, j + shifty)].IntersectP(ray)){
			float txNxt = txCur + 3 * txAdd;
			float tyNxt = tyCur + 3 * tyAdd;

			if(txNxt < tyNxt){
				if(txNxt > ray.maxt) break;
				int icrJ = (int)((txNxt - tyCur) * inv_tyAdd);
				tyCur += icrJ * tyAdd;
				j += icrJ * y_going;

				txCur = txNxt;
				i += 3 * x_going;

				if(i == outi) break;
			}
			else{
				if(tyNxt > ray.maxt) break;
				int icrI = (int)((tyNxt - txCur) * inv_txAdd);
				txCur += icrI * txAdd;
				i += icrI * x_going;

				tyCur = tyNxt;
				j += 3 * y_going;
			
				if(j == outj) break;
			}

			continue;
		}
	
		float txNxt = txCur + txAdd;
		float tyNxt = tyCur + tyAdd;

		float tNear = (boxmaxz[VERT(i, j)] - ray.o[2]) * invdz;
		float  tFar = (boxminz[VERT(i, j)] - ray.o[2]) * invdz;
		if(tNear > tFar) swap(tNear, tFar);

		if(txNxt < tyNxt){
			if(max(max(txCur, tyCur), tNear) < min(txNxt, tFar))
			{
				inin_cnt ++;
				gridIntersect(ray, i, j, hitsomething, tHit, dg);
			}
			if(hitsomething) break;

			if(txNxt > ray.maxt) break;
			i += x_going;
			txCur = txNxt;

			if(i == outi) break;
		}
		else{
			if(max(max(txCur, tyCur), tNear) < min(tyNxt, tFar))
			{
				inin_cnt ++;
				gridIntersect(ray, i, j, hitsomething, tHit, dg);
			}
			if(hitsomething) break;

			if(tyNxt > ray.maxt) break;
			j += y_going;
			tyCur = tyNxt;

			if(j == outj) break;
		}
	}
	
	//printf("%d / %d = %f, %d / %d = %f\n", in_cnt, out_cnt, 1. * in_cnt / out_cnt, inin_cnt, out_cnt, 1. * inin_cnt / out_cnt);

	if(!hitsomething) return false;

	*rayEpsilon = 1e-3f * *tHit;
	return true;
}


bool Heightfield2::IntersectP(const Ray &r) const {
	DifferentialGeometry dg;
	float tHit, rayEpsilon;

	return Intersect(r, &tHit, &rayEpsilon, &dg);
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
