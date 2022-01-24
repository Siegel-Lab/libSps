#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <algorithm>
#include <vector>
#include <pybind11/stl.h>

namespace kdtree{

class CacheVector
{
    public:
        typedef std::vector<int64_t> Point;
        typedef std::pair<Point, std::string> DataPoint;
        typedef std::vector<DataPoint> Points;
    private:
        class Node
        {
            public:
                Points vPoints = Points();
                int64_t uiIntegral = 0;
                std::shared_ptr<CacheVector> pCache = nullptr;
        };

        std::vector<Point> vBins;
        std::vector<Node> vData;

        size_t getIdx(const Point& xPoint)
        {
            size_t uiRet = 0;
            for(size_t i = 0; i < vBins.size(); i++)
            {
                if(i > 0)
                    uiRet *= vBins[i-1].size();
                uiRet += (std::upper_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]) - vBins[i].begin()) - 1;
            }
            return uiRet;
        }

        size_t maxIdx()
        {
            size_t uiRet = 1;
            for(size_t i = 0; i < vBins.size(); i++)
                    uiRet *= vBins[i].size();
            return uiRet;
        }

        void integrate(size_t iFixedDim, size_t iDim, Point& vPoint){
            if(iDim == iFixedDim)
                iDim += 1;
            if(iDim >= vBins.size())
            {
                size_t uiIntegral = 0;
                for(size_t uiI : vBins[iFixedDim])
                {
                    vPoint[iFixedDim] = uiI;
                    size_t uiIdx = getIdx(vPoint);
                    uiIntegral += vData[uiIdx].vPoints.size();
                    vData[uiIdx].uiIntegral = uiIntegral;
                }
            }
            else
                for(size_t uiI : vBins[iDim])
                {
                    vPoint[iDim] = uiI;
                    integrate(iFixedDim, iDim + 1, vPoint);
                }
        }

        void integrate()
        {
            Point vPoint(0, vBins.size());
            for(size_t iFixedDim = 0; iFixedDim < vBins.size(); iFixedDim++)
                integrate(iFixedDim, 0, vPoint);
        }

        int64_t count_help(Point& vP, Point vLeft, Point vRight, size_t uiI, size_t uiEdgesFromRight)
        {
            if(uiI >= vBins.size())
            {
                if(uiEdgesFromRight % 2 == 0)
                    return (int64_t)vData[getIdx(vP)].uiIntegral;
                return -(int64_t)vData[getIdx(vP)].uiIntegral;
            }
            int64_t iRet = 0;
            vP[uiI] = vLeft[uiI];
            iRet += count_help(vP, vLeft, vRight, uiI+1, uiEdgesFromRight+1);
            vP[uiI] = vRight[uiI];
            iRet += count_help(vP, vLeft, vRight, uiI+1, uiEdgesFromRight);
            return iRet;
        }

    public:
        CacheVector(Points vvfPoints, std::vector<Point> vBins)
                :
            vBins(vBins),
            vData(maxIdx())
        {
            for(const DataPoint& p: vvfPoints)
                vData[getIdx(p.first)].vPoints.push_back(p);
            integrate();
        }

        size_t count(Point vLeft, Point vRight)
        {
            Point vP(0, vBins.size());
            return (size_t)count_help(vP, vLeft, vRight, 0, 0);
        }
};

}

