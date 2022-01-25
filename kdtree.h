#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <deque>
#include <queue>
#include <stdexcept>
#include <thread>
#include <algorithm>
#include <vector>
#include <pybind11/stl.h>

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif // DEBUG_LEVEL

#if DEBUG_LEVEL >= 1
#define DEBUG( x ) x
#define DEBUG_PARAM( x ) , x
#else // DEBUG_LEVEL
#define DEBUG_PARAM( x )
#define DEBUG( x )
#endif // DEBUG_LEVEL

namespace kdtree
{

/**
 * @brief A Threadpool.
 * @details
 * // create thread pool with 4 worker threads\n
 * ThreadPool pool(4);\n
 *
 * // enqueue and store future\n
 * auto result = pool.enqueue( [](int answer) { return answer; }, 42 );\n
 *
 * // get result from future\n
 * std::cout << result.get() << std::endl;\n
 *
 *
 * if the threadpool is set to have 0 threads, we will execute every task in the main thread
 * immediately, thus we dont need to setup threads at all.
 * this can be used to disable multithreading simply setting the pool to have 0 threads
 */
class ThreadPool
{
  private:
    /* need to keep track of threads so we can join them
     */
    std::vector<std::thread> workers;

    /* the task queue
     */
    std::queue<std::function<void( size_t )>> tasks;

    /* Synchronization
     */
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool bStop;
    const size_t threads;

  public:
    /* External definition
     */
    ThreadPool( size_t );

    /* External definition
     */
    ~ThreadPool( );

    /* add new work item to the pool
     */
    template <class F, class... Args>
    auto enqueue( F&& f, Args&&... args ) -> std::future<typename std::result_of<F( size_t, Args... )>::type>
    {
        typedef typename std::result_of<F( size_t, Args... )>::type return_type;

        /* Don't allow enqueueing after stopping the pool
         */
        if( bStop )
        {
            throw std::runtime_error( "enqueue on stopped ThreadPool" );
        } // if

        /* std::placeholders::_1 and represents some future value (this value will be only known
         * during function definition.
         */
        auto task = std::make_shared<std::packaged_task<return_type( size_t )>>(
            std::bind( std::forward<F>( f ), std::placeholders::_1, std::forward<Args>( args )... ) );

        /* The future outcome of the task. The caller will be later blocked until this result will
         * be available.
         */
        std::future<return_type> xFuture = task->get_future( );

        // if the threadpool is set to have 0 threads then execute the task in the main thread
        if( threads == 0 )
        {
            ( *task )( 0 );
            return xFuture;
        } // if

        {
            /* Mutual access to the task queue has to be synchronized.
             */
            std::unique_lock<std::mutex> lock( queue_mutex );

            /* The task will get delivered its task_id later during its execution.
             */
            tasks.push( [task]( size_t task_id ) { ( *task )( task_id ); } );
        } // end of scope for the lock, so we will release the lock.

        /* Inform some waiting consumer (worker )that we have a fresh task.
         */
        condition.notify_one( );

        return xFuture;
    } // method enqueue
};

/* Constructor just launches some amount of workers
 */
inline ThreadPool::ThreadPool( size_t threads )
    : bStop( false ), // stop must be false in the beginning
      threads( threads )
{
    // if the threadpool is set to have 0 threads, we will execute every task in the main thread
    // immediately, thus we dont need to setup threads at all.
    // this can be used to disable multithreading simply setting the pool to have 0 threads
    if( threads == 0 )
        return;
    for( size_t i = 0; i < threads; ++i )
        workers.emplace_back( [this, i] {
            /* So long we have some data the thread processes these data
             */
            for( ;; )
            {
                /* Synchronization of mutual access to the task queue.
                 */
                std::unique_lock<std::mutex> lock( this->queue_mutex );

                while( !this->bStop && this->tasks.empty( ) )
                {
                    /* We release the lock, so that some producer can push some task into the queue.
                     * The produser will call notify_one(), in order to release
                     */
                    this->condition.wait( lock );
                } // while

                if( this->bStop && this->tasks.empty( ) )
                {
                    /* All work done and signal for termination received.
                     */
                    return;
                } // if

                /* Initialize the function variable task, so that it refers the task at the top of
                 * the queue.
                 */
                std::function<void( size_t )> task( this->tasks.front( ) );
                this->tasks.pop( );

                lock.unlock( );

                try
                {
                    /* Execute the task (that we received as top of the queue)
                     */
                    task( i );
                }
                catch( std::exception e )
                {
                    std::cerr << "exception when executing task:" << std::endl;
                    std::cerr << e.what( ) << std::endl;
                }
                catch( ... )
                {
                    std::cerr << "unknown exception when executing task" << std::endl;
                }
            } // for
        } // lambda
        ); // function call
} // method

/* the destructor joins all threads
 */
inline ThreadPool::~ThreadPool( )
{
    // if the threadpool is set to have 0 threads, we will execute every task in the main thread
    // immediately, thus we dont need to setup threads at all.
    // this can be used to disable multithreading simply setting the pool to have 0 threads
    if( threads == 0 )
        return;

    {
        std::unique_lock<std::mutex> lock( queue_mutex );
        bStop = true;
    } // end of scope for queue_mutex, so we unlock the mutex
    condition.notify_all( );

    /* We wait until all workers finished their job.
     * (Destruction thread is blocked until all workers finished their job.)
     */
    for( size_t i = 0; i < workers.size( ); ++i )
    {
        workers[ i ].join( );
    } // for
}



class CacheVector
{
    public:
        typedef std::vector<int64_t> Point;
        typedef std::pair<Point, std::string> DataPoint;
        typedef std::vector<DataPoint> Points;
        
        CacheVector(Points vvfPoints, std::vector<Point> vBins, size_t uiSizeCache, size_t uiThreads);
    private:
        class Node
        {
            public:
                Points vPoints = Points();
                int64_t uiIntegral = 0;
                std::shared_ptr<CacheVector> pCache = nullptr;

                std::string print(){
                    std::string sRet = "";
                    sRet += std::to_string(vPoints.size());
                    sRet += " -> ";
                    sRet += std::to_string(uiIntegral);
                    return sRet;
                }
        };

        std::vector<Point> vBins;
        std::vector<Node> vData;
        size_t uiSizeCache;
        size_t uiThreads;
        std::deque<Node*> vpCache;

        size_t getIdx(const Point& xPoint, std::vector<bool> vSubstrOne)
        {
            size_t uiRet = 0;
            for(size_t i = 0; i < vBins.size(); i++)
            {
                if(i > 0)
                    uiRet *= vBins[i].size();
                auto itI = std::upper_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]);
                if(vSubstrOne[i] && itI != vBins[i].begin())
                    --itI;
                if(itI == vBins[i].begin())
                    return std::numeric_limits<size_t>::max();
                //bOnSpot = true;
                uiRet += (itI - vBins[i].begin()) - 1;
            }
            return uiRet;
        }

        bool onSpot(const Point& xPoint)
        {
            for(size_t i = 0; i < vBins.size(); i++)
                if(*std::lower_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]) != xPoint[i])
                    return false;
            return true;
        }

        std::vector<Point> makeBins(const Point& xPoint)
        {
            std::vector<Point> vRet;
            for(size_t i = 0; i < vBins.size(); i++)
            {
                auto itI = std::upper_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]);
                Point xP;
                for(int64_t uiI = *(itI-1); uiI <= *itI; uiI++)
                    xP.push_back(uiI);
                DEBUG(std::cerr << "makeBins " << pointCoords(xP) << std::endl;)
                vRet.push_back(xP);
            }
            return vRet;
        }

        size_t maxIdx()
        {
            size_t uiRet = 1;
            for(size_t i = 0; i < vBins.size(); i++)
                    uiRet *= vBins[i].size();
            return uiRet;
        }

        std::string pointCoords(Point& vP)
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

        void integrate(size_t iFixedDim, size_t iDim, Point& vPoint, size_t* uiCnt, std::mutex& rLock){
            if(iDim == iFixedDim)
                iDim += 1;
            if(iDim >= vBins.size())
            {
                size_t uiIntegral = 0;
                for(size_t uiI : vBins[iFixedDim])
                {
                    vPoint[iFixedDim] = uiI;
                    size_t uiIdx = getIdx(vPoint, std::vector<bool>(vBins.size(), false));
                    DEBUG(std::cout << "integrate " << pointCoords(vPoint) << " (" << uiIdx << "/" << vData.size()
                          << ")" << ": " << 
                          vData[uiIdx].uiIntegral << " -> " << uiIntegral << std::endl<< std::endl;)
                    uiIntegral += vData[uiIdx].uiIntegral;
                    vData[uiIdx].uiIntegral = uiIntegral;
                    if(uiCnt != nullptr)
                    {
                        size_t x = (*uiCnt)++;
                        if(x % 1000000 == 0)
                        {
                            std::unique_lock<std::mutex> _(rLock);
                            std::cerr << "\rintegrating " << x << " / " <<  maxIdx() << " = " << 
                                         (100*x)/maxIdx() << "%\033[K";
                        }
                    }
                }
                DEBUG(std::cout << std::endl;)
            }
            else
            {
                ThreadPool xPool((iDim == 0 || (iFixedDim == 0 && iDim == 1)) ? uiThreads : 0);
                for(size_t uiI : vBins[iDim])
                    xPool.enqueue([&](size_t, size_t uiI){
                        Point vP(vPoint);
                        vP[iDim] = uiI;
                        integrate(iFixedDim, iDim + 1, vP, uiCnt, rLock);
                    }, uiI);
            }
        }

        void integrate()
        {
            std::cerr << "\r\033[Kintegrating";
            Point vPoint(vBins.size(), 0);
            size_t uiCnt = 0;
            std::mutex xLock;
            for(size_t iFixedDim = 0; iFixedDim < vBins.size(); iFixedDim++)
                integrate(iFixedDim, 0, vPoint, &uiCnt, xLock);
            std::cerr << "\r\033[K";
        }

        void makeCache(Node& rN, Point& vP)
        {
            if(rN.pCache != nullptr)
                return;
            while(vpCache.size() >= uiSizeCache)
            {
                vpCache.front()->pCache.reset();
                vpCache.pop_front();
            }
            DEBUG(std::cerr << "computing cache " << vpCache.size() << " due to " << pointCoords(vP) << " out of " << rN.vPoints.size() << " elements." << std::endl;)
            rN.pCache = std::make_shared<CacheVector>(rN.vPoints, makeBins(vP), 0, uiThreads);
            DEBUG(std::cerr << rN.pCache->print() << std::endl;)
            vpCache.push_back(&rN);
        }

        int64_t pointVal(Point& vP, std::vector<bool>& vSub)
        {
            auto xIdx = getIdx(vP, vSub);
            if(xIdx == std::numeric_limits<size_t>::max())
                return 0;
            int64_t iRet = 0;
            return (int64_t)vData[xIdx].uiIntegral;
        }

        int64_t count_help(Point& vP, std::vector<bool>& vSub, Point vLeft, Point vRight, 
                           size_t uiI, size_t uiEdgesFromRight)
        {
            if(uiI >= vBins.size())
            {
                int64_t x = 0;
                if(uiEdgesFromRight % 2 == 0)
                    x = pointVal(vP, vSub);
                else
                    x = -pointVal(vP, vSub);
                DEBUG(std::cerr << "count " << pointCoords(vP) << ": " << x << std::endl;)
                return x;
            }
            int64_t iRet = 0;
            vP[uiI] = vLeft[uiI];
            vSub[uiI] = true;
            iRet += count_help(vP, vSub, vLeft, vRight, uiI+1, uiEdgesFromRight+1);
            vP[uiI] = vRight[uiI];
            vSub[uiI] = true;
            iRet += count_help(vP, vSub, vLeft, vRight, uiI+1, uiEdgesFromRight);
            return iRet;
        }

        int64_t count_cache(Point& vP, Point vLeft, Point vRight, size_t i)
        {
            if(i >= vBins.size())
            {
                size_t xIdx = getIdx(vP, std::vector<bool>(vBins.size(), false));
                DEBUG(std::cerr << "pointVal->cache " << xIdx << std::endl;)
                makeCache(vData[xIdx], vP);
                return vData[xIdx].pCache->count(vLeft, vRight);
            }
            else
            {
                int64_t uiRet = 0;
                for(auto itI = std::lower_bound(vBins[i].begin(), vBins[i].end(), vLeft[i]);
                    itI != std::upper_bound(vBins[i].begin(), vBins[i].end(), vRight[i]);
                    ++itI
                    )
                {
                    vP[i] = *itI;
                    uiRet += count_cache(vP, vLeft, vRight, i+1);
                }
                return uiRet;
            }
        }

    public:

        size_t count(Point vLeft, Point vRight)
        {
            if(onSpot(vLeft) && onSpot(vRight))
            {
                Point vP(vBins.size(), 0);
                std::vector<bool> vSub(vBins.size(), false);
                return (size_t)count_help(vP, vSub, vLeft, vRight, 0, 0);
            }
            else 
            {
                return 0;
                DEBUG(std::cerr << "cache " << std::endl;)
                Point vP(vBins.size(), 0);
                return count_cache(vP, vLeft, vRight, 0);
            }
        }

        std::string print_help(Point& vP, size_t uiDim)
        {
            std::string sRet = "";
            if(uiDim >= vBins.size())
                sRet = pointCoords(vP) + ": " + 
                       vData[getIdx(vP, std::vector<bool>(vBins.size(), false))].print() + "\n";
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


CacheVector::CacheVector(Points vvfPoints, std::vector<Point> vBins, size_t uiSizeCache, size_t uiThreads)
        :
    vBins(vBins),
    vData(maxIdx()),
    uiSizeCache(uiSizeCache),
    uiThreads(uiThreads)
{
    {
        std::cerr << "loading";
        size_t uiT = 0;//uiThreads;
        ThreadPool xPool(uiT);
        if(uiT < 1)
            uiT = 1;
        std::vector<std::mutex> vArr(uiT*100);
        for(size_t uiI = 0; uiI < uiT; uiI++)
            xPool.enqueue([&](size_t, size_t uiI){
                for(size_t uiX = uiI; uiX < vvfPoints.size(); uiX+=uiT){
                    if(uiX % 1000000 == 0)
                        std::cerr << "\rloading " << uiX << " / " << vvfPoints.size() << "\033[K";
                    auto p = vvfPoints[uiX];
                    
                    DEBUG(std::cout << "loading " << pointCoords(p.first) << std::endl;)
                    size_t uiIdx = getIdx(p.first, std::vector<bool>(vBins.size(), false));
                    std::unique_lock<std::mutex> xLock(vArr[uiIdx % vArr.size()]);
                    vData[uiIdx].vPoints.push_back(p);
                    vData[uiIdx].uiIntegral = vData[uiIdx].vPoints.size();
                }
            }, uiI);
    }
    std::cerr << "\r\033[K";
    integrate();
}

}

