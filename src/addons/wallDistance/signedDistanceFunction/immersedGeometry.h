// Numerical helpers for masking immersed distance. No STK mutation.
#pragma once
#include <ssdf.hpp>
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace accel::immersedGeometry
{
using Point = std::array<double, 3>;
struct Surface
{
    int dimension;
    std::vector<Point> points;
    std::vector<std::array<int, 3>> faces;
};

inline void collectiveCheck(bool okay, MPI_Comm comm, const std::string& message)
{
    int bad = okay ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, comm);
    if (bad) throw std::runtime_error(message);
}

struct Result
{
    std::vector<double> distance;
    std::vector<uint8_t> inside;
    double distanceSeconds = 0.0;
    double signSeconds = 0.0;
    double combined2DSeconds = 0.0;
};
inline Result evaluate(const Surface& s,const std::vector<Point>& queries)
{
    if (s.points.empty() || s.faces.empty())
        throw std::runtime_error("Empty immersed surface");
    for (const auto& p : s.points)
        for (double v : p)
            if (!std::isfinite(v))
                throw std::runtime_error("Nonfinite immersed surface coordinates");
    for (const auto& f : s.faces)
    {
        for (int k = 0; k < (s.dimension == 3 ? 3 : 2); ++k)
            if (f[k] < 0 || static_cast<std::size_t>(f[k]) >= s.points.size())
                throw std::runtime_error("Invalid immersed surface connectivity");
        if (s.dimension == 2 && s.points[f[0]] == s.points[f[1]])
            throw std::runtime_error("Zero length immersed edge");
    }
    if(queries.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Local query count exceeds SSDF classifier capacity");
    for(const auto& p:queries) for(double v:p)
        if(!std::isfinite(v)) throw std::runtime_error("Nonfinite immersed query coordinates");
    Result r{std::vector<double>(queries.size(),std::numeric_limits<double>::max()),std::vector<uint8_t>(queries.size(),0)};
    if(queries.empty()) return r; // A rank may own no background nodes.
    if(s.dimension==2)
    {
        const double combinedStart = MPI_Wtime();
        // Exact segment projection and half-open ray crossings in the XY plane.
        // Each iteration writes its own output; MPI and STK are outside this loop.
#if defined(_OPENMP) && defined(SSDF_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(queries.size());++i)
        {
            const auto p=queries[i];double best=std::numeric_limits<double>::max();bool inside=false;
            for(const auto& f:s.faces)
            {
                const auto a=s.points[f[0]],b=s.points[f[1]];
                double dx=b[0]-a[0],dy=b[1]-a[1];
                double t=std::clamp(((p[0]-a[0])*dx+(p[1]-a[1])*dy)/(dx*dx+dy*dy),0.0,1.0);
                double x=p[0]-a[0]-t*dx,y=p[1]-a[1]-t*dy;
                best=std::min(best,x*x+y*y);
                if((a[1]>p[1])!=(b[1]>p[1]))
                    if(p[0]<a[0]+(p[1]-a[1])*dx/dy) inside=!inside;
            }
            r.distance[i]=std::sqrt(best);r.inside[i]=inside;
        }
        r.combined2DSeconds = MPI_Wtime() - combinedStart;
        return r;
    }
#ifndef SSDF_ENABLE_CUBIQL
    throw std::runtime_error("3D immersed classification requires SSDF_ENABLE_CUBIQL=ON");
#else
    std::vector<double> x,y,z,sx,sy,sz;std::vector<int>a,b,c;
    for(auto p:queries) { x.push_back(p[0]);y.push_back(p[1]);z.push_back(p[2]); }
    for(auto p:s.points) { sx.push_back(p[0]);sy.push_back(p[1]);sz.push_back(p[2]); }
    for(auto f:s.faces) { a.push_back(f[0]);b.push_back(f[1]);c.push_back(f[2]); }
    const double distanceStart = MPI_Wtime();
    int rc=ssdf::edf_select<double,double,int>(queries.size(),x.data(),y.data(),z.data(),a.size(),a.data(),b.data(),c.data(),sx.size(),sx.data(),sy.data(),sz.data(),r.distance.data());
    r.distanceSeconds = MPI_Wtime() - distanceStart;
    if(rc) throw std::runtime_error("SSDF immersed distance failed");
    std::vector<float> qx(x.begin(),x.end()),qy(y.begin(),y.end()),qz(z.begin(),z.end());
    std::vector<float> vx(sx.begin(),sx.end()),vy(sy.begin(),sy.end()),vz(sz.begin(),sz.end());
    const double signStart = MPI_Wtime();
    rc=ssdf::points_inside_bvh<float,int>(queries.size(),qx.data(),qy.data(),qz.data(),a.size(),a.data(),b.data(),c.data(),sx.size(),vx.data(),vy.data(),vz.data(),r.inside.data());
    r.signSeconds = MPI_Wtime() - signStart;
    if(rc) throw std::runtime_error("SSDF immersed classification failed");
    return r;
#endif
}

// Nearest individual surface + inside-any sign. Exact union SDF only for disjoint solids.
struct Combined
{
    std::vector<double> distance;
    std::vector<uint8_t> inside;
    explicit Combined(std::size_t n):distance(n,std::numeric_limits<double>::max()),inside(n,0){}
    double signedValue(std::size_t i,bool signedOutput) const
    { return distance[i]==0?0.0:((signedOutput && inside[i])?-distance[i]:distance[i]); }
};
}
