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


class KDTree
{
    public:
        size_t uiThreads;
        size_t uiTRunning;
        size_t uiDims;
        typedef std::pair<std::vector<float>, std::string> Point;
        typedef std::vector<Point> Points;
    private:
        std::mutex xMutex;
        class Node
        {
            public:
                Point vPoint;
                std::shared_ptr<Node> pLeft = nullptr;
                std::shared_ptr<Node> pRight = nullptr;
                size_t uiCount;

                Node() : vPoint(), uiCount(0) {}
                Node(Point vPoint) : vPoint(vPoint), uiCount(1) {}
                Node(Point vPoint, std::shared_ptr<Node> pLeft, std::shared_ptr<Node> pRight) : 
                        vPoint(vPoint), 
                        pLeft(pLeft),
                        pRight(pRight),
                        uiCount(1 + (pLeft == nullptr ? 0 : pLeft->uiCount) + (pRight == nullptr ? 0 : pRight->uiCount))
                {}
        };

        template <class F, class... Args>
        auto enqueue( F&& f, Args&&... args ) -> std::future<typename std::result_of<F( Args... )>::type>
        {
            typedef typename std::result_of<F( Args... )>::type return_type;

            auto task = std::make_shared<std::packaged_task<return_type( )>>(
                std::bind( std::forward<F>( f ), std::forward<Args>( args )... ) );

            std::future<return_type> xFuture = task->get_future( );

            bool binThread = false;
            {
                std::unique_lock<std::mutex> lock( xMutex );
                binThread = uiTRunning < uiThreads;
                if(binThread)
                    uiTRunning += 1;
            }
            if(!binThread)
                ( *task )( );
            else
                std::thread([&](){(*task)(); std::unique_lock<std::mutex> lock( xMutex ); uiTRunning-=1;});
            return xFuture;
        } // method enqueue

        std::shared_ptr<Node> construct(Points::iterator itBegin, Points::iterator itEnd, size_t uiI)
        {
            if(itEnd - itBegin > 1)
            {
                std::sort(itBegin, itEnd, [&] (const Point& rA, const Point& rB){
                    return rA.first[uiI] < rB.first[uiI];
                });
                auto itCenter = itBegin + ((itEnd - itBegin) / 2);
                auto fLeft = enqueue( [&](Points::iterator itA, Points::iterator itB, size_t uiI) 
                { 
                    return construct(itA, itB, uiI); 
                }, itBegin, itCenter, (uiI+1)%uiDims );
                auto fRight = enqueue( [&](Points::iterator itA, Points::iterator itB, size_t uiI) 
                { 
                    return construct(itA, itB, uiI); 
                }, itCenter+1, itEnd, (uiI+1)%uiDims );
 
                return std::make_shared<Node>(*itCenter, fLeft.get(), fRight.get() );
            }
            else if (itEnd - itBegin == 1)
                return std::make_shared<Node>(*itBegin);
            else
                return nullptr;
        }

        size_t count_help(std::shared_ptr<Node> pNode, std::vector<float> vLeft, std::vector<float> vRight,
                     std::vector<bool> vbInLeft, std::vector<bool> vbInRight, size_t uiI)
        {
            if(pNode == nullptr)
                return 0;
            if (
                std::all_of(vbInLeft.begin(), vbInLeft.end(), [](bool v) { return v; }) &&
                std::all_of(vbInRight.begin(), vbInRight.end(), [](bool v) { return v; })
            )
                return pNode->uiCount;
            size_t uiRet = 1;
            for(size_t uiJ = 0; uiJ < uiDims; uiJ++)
                if(vLeft[uiJ] > pNode->vPoint.first[uiJ] || vRight[uiJ] < pNode->vPoint.first[uiJ])
                {
                    uiRet = 0;
                    break;
                }

            std::vector<std::future<size_t>> vTasks;
            if(vRight[uiI] >= pNode->vPoint.first[uiI])
            {
                std::vector<bool> vbInLeft2(vbInLeft);
                if(vLeft[uiI] < pNode->vPoint.first[uiI])
                    vbInLeft2[uiI] = true;
                vTasks.push_back(enqueue([&](std::vector<bool> vbInLeft2){
                    return count_help(pNode->pRight, vLeft, vRight, vbInLeft2, vbInRight, (uiI + 1) % uiDims);
                }, vbInLeft2));
            }
            if(vLeft[uiI] <= pNode->vPoint.first[uiI])
            {
                std::vector<bool> vbInRight2(vbInLeft);
                if(vRight[uiI] > pNode->vPoint.first[uiI])
                    vbInRight2[uiI] = true;
                vTasks.push_back(enqueue([&](std::vector<bool> vbInLeft2){
                    return count_help(pNode->pLeft, vLeft, vRight, vbInLeft, vbInRight2, (uiI + 1) % uiDims);
                }, vbInRight2));
            }

            for(auto& rX : vTasks)
                uiRet += rX.get();

            return uiRet;
        }

        std::shared_ptr<Node> pRoot;
    public:

        KDTree(size_t uiThreads, size_t uiDims, Points vvfPoints) : uiThreads(uiThreads), uiTRunning(0), uiDims(uiDims) 
        {
            pRoot = construct(vvfPoints.begin(), vvfPoints.end(), 0);
        }

        size_t count(std::vector<float> vLeft, std::vector<float> vRight)
        {
            return count_help(pRoot, vLeft, vRight, std::vector<bool>(uiDims, false), std::vector<bool>(uiDims, false),
                         0);
        }
};

class CacheVector
{
    public:
        typedef std::vector<size_t> Point;
        typedef std::pair<Point, std::string> DataPoint;
        typedef std::vector<Point> Points;
        class Node
        {
            public:
                Points vPoints;
                size_t uiIntegral = 0;
                std::shared_ptr<CacheVector> pCache = nullptr;
        };

        std::vector<std::vector<size_t>> vBins;
        std::vector<std::shared_ptr<Node>> vData;

        size_t getIdx(Point& xPoint)
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
            Point vPoint(0, vDims.size());
            for(size_t iFixedDim = 0; iFixedDim < vBins.size(); iFixedDim++)
                integrate(iFixedDim, 0, vPoint);
        }

        CacheVector(Points vvfPoints, std::vector<Point> vBins)
                :
            uiBinSize(uiBinSize),
            vBins(vBins),
            vData(maxIdx())
        {
            for(DataPoint& p: vvfPoints)
                vData[getIdx(p.first)].append(p);
            integrate();
        }
}

}

