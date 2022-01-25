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

                std::string print(){
                    std::string sRet = "";
                    sRet += std::to_string(uiIntegral);
                    return sRet;
                }
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

        std::string poitCoords(Point& vP)
        {
            std::string sRet = "";
            sRet += "(";
            for(size_t uiI = 0; uiI < vP.size(); uiI++)
            {
                sRet += std::to_string(vP[uiI]) ;
                if(uiI + 1 < vP.size())
                    sRet += ", ";
            }
            return sRet + ")";
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
                    //std::cout << "integrate " << poitCoords(vPoint) << ": " << uiIntegral << std::endl;
                    uiIntegral += vData[uiIdx].uiIntegral;
                    vData[uiIdx].uiIntegral = uiIntegral;
                }
                //std::cout << std::endl;
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
            Point vPoint(vBins.size(), 0);
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
            {
                vData[getIdx(p.first)].vPoints.push_back(p);
                ++vData[getIdx(p.first)].uiIntegral;
            }
            integrate();
        }

        size_t count(Point vLeft, Point vRight)
        {
            Point vP(vBins.size(), 0);
            return (size_t)count_help(vP, vLeft, vRight, 0, 0);
        }

        std::string print_help(Point& vP, size_t uiDim)
        {
            std::string sRet = "";
            if(uiDim >= vBins.size())
                sRet = poitCoords(vP) + ": " + vData[getIdx(vP)].print() + "\n";
            else
                for(size_t uiI : vBins[uiDim])
                {
                    vP[uiDim] = uiI;
                    sRet += print_help(vP, uiDim + 1);
                }
            return sRet;
        }
        std::string print()
        {
            Point vP(vBins.size(), 0);
            return print_help(vP, 0);
        }
};

}

