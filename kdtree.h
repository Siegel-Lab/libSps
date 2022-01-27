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
#include <fstream>
#include <iterator>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

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


template<size_t N>
class CacheVector
{
    public:
        typedef std::array<int64_t, N> Point;
        typedef std::pair<Point, char> DataPoint;
        typedef std::vector<DataPoint> Points;
        
        CacheVector(std::string sFilename, Points vvfPoints, std::array<std::vector<int64_t>, N> vBins,
                    size_t uiSizeCache, size_t uiThreads,
                    size_t uiUntouchedDimensions, size_t uiLayer, size_t uiEleNoCache);
        CacheVector(std::string sFilename, size_t uiSizeCache, size_t uiThreads,
                    size_t uiUntouchedDimensions, size_t uiLayer, size_t uiEleNoCache);
    private:

        std::string sFilename;
        std::array<std::vector<int64_t>, N> vBins;
        size_t uiUntouchedDimensions;
        std::vector<Points> vDataPoints;
        std::vector<std::shared_ptr<CacheVector>> vDataCache;
        std::vector<int64_t> vData;
        size_t uiSizeCache;
        size_t uiEleNoCache;
        size_t uiThreads;
        size_t uiLayer;
        static std::vector<std::deque<std::shared_ptr<CacheVector>>> vpCache;
        static std::mutex xCacheLock;

        template<typename VEC_T>
        static void vectorToFile(std::string sFilename, VEC_T& vector)
        {
            std::ofstream myfile(sFilename, std::ios::binary | std::ios::trunc);
            myfile.write((char*)&vector[0], vector.size() * sizeof(typename VEC_T::value_type));
            myfile.close();
        }

        template<typename VEC_T>
        static void fileToVector(std::string sFilename, VEC_T& vector)
        {
            // open the file:
            std::ifstream file(sFilename, std::ios::binary);

            // Stop eating new lines in binary mode!!!
            file.unsetf(std::ios::skipws);

            file.seekg(0, std::ios::end);
            // reserve capacity
            vector.clear();
            vector.resize(file.tellg() / sizeof(typename VEC_T::value_type));
            file.seekg(0, std::ios::beg);

            // read the data:
            file.read((char*)&vector[0], vector.size() * sizeof(typename VEC_T::value_type));
            file.close();
        }

        void serialize(std::string sFilename, bool bOut)
        {
            for(size_t uiI = 0; uiI < N; uiI++)
            {
                if(bOut)
                    vectorToFile(sFilename + ".axis." + std::to_string(uiI) + ".bin", vBins[uiI]);
                else
                    fileToVector(sFilename + ".axis." + std::to_string(uiI) + ".bin", vBins[uiI]);
            }

            if(bOut)
                vectorToFile(sFilename + ".count.bin", vData);
            else
                fileToVector(sFilename + ".count.bin", vData);

            //std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << (bOut ? "true" : "false") << std::endl;
            if(bOut)
            {
                std::ofstream myfile(sFilename + ".points.bin", std::ios::binary | std::ios::trunc);
                std::vector<size_t> vSizes;
                for(auto& v : vDataPoints)
                    vSizes.push_back(v.size());
                size_t x = vDataPoints.size();
                myfile.write((char*)&x, sizeof(size_t));
                myfile.write((char*)&vSizes[0], vSizes.size() * sizeof(size_t));
                //std::cerr << "yyyyy " << vDataPoints.size() << std::endl;
                for(size_t uiI = 0; uiI < vDataPoints.size(); uiI++)
                {
                    if(vDataPoints[uiI].size() > 0)
                        myfile.write((char*)&vDataPoints[uiI][0], vDataPoints[uiI].size() * sizeof(DataPoint));
                
                    //std::cerr << "yy " << vDataPoints[uiI].size() << std::endl;
                    //for(size_t uiK = 0; uiK < vDataPoints[uiI].size(); uiK++)
                    //    std::cerr << pointCoords(vDataPoints[uiI][uiK].first) << " ";
                    //std::cerr << std::endl;
                }
                std::vector<char> vSuff(8);
                myfile.write(&vSuff[0], vSuff.size());

                myfile.flush();
                myfile.close();
            }
            else
            {
                // open the file:
                std::ifstream file(sFilename + ".points.bin", std::ios::binary);

                // Stop eating new lines in binary mode!!!
                file.unsetf(std::ios_base::skipws);

                size_t uiI = 0;
                file.read((char*)&uiI, sizeof(size_t));
                vDataPoints.resize(uiI);
                std::vector<size_t> vSizes;
                vSizes.resize(uiI);
                file.read((char*)&vSizes[0], uiI * sizeof(size_t));
                for(size_t uiI = 0; uiI < vDataPoints.size(); uiI++)
                {
                    if(vSizes[uiI] > 0)
                    {
                        vDataPoints[uiI].resize(vSizes[uiI]);
                        file.read((char*)&vDataPoints[uiI][0], vSizes[uiI] * sizeof(DataPoint));
                    }
                    //std::cerr << "xx " << vSizes[uiI] << std::endl;
                    //for(size_t uiK = 0; uiK < vSizes[uiI]; uiK++)
                    //    std::cerr << pointCoords(vDataPoints[uiI][uiK].first) << " ";
                    //std::cerr << std::endl;
                }

                // read the data:
                file.close();
            }
        }

        void save_helper(std::string sFlieName, size_t uiMinSizeCacheSave)
        {
            //@todo template this
            serialize(sFlieName, true);
            Point x;
            for(size_t uiI = 0; uiI < vDataCache.size(); uiI++)
                if(vDataCache[uiI] != nullptr || vDataPoints[uiI].size() > uiMinSizeCacheSave)
                {
                    makeCache(uiI, x);
                    vDataCache[uiI]->save_helper(sFlieName + "." + std::to_string(uiI), uiMinSizeCacheSave);
                    vDataCache[uiI].reset();
                }
        }

        size_t getIdx(const Point& xPoint, bool bSubstrOne)
        {
            size_t uiRet = 0;
            for(size_t i = 0; i < vBins.size(); i++)
            {
                uiRet *= vBins[i].size();
                auto itI = std::upper_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]);
                if(bSubstrOne && itI != vBins[i].begin())
                    --itI;
                if(itI == vBins[i].begin())
                    return std::numeric_limits<size_t>::max();
                //bOnSpot = true;
                uiRet += (itI - vBins[i].begin()) - 1;
            }
            return uiRet;
        }

        size_t getIdxPoints(const Point& xPoint, bool bSubstrOne)
        {
            size_t uiRet = 0;
            for(size_t i = 0; i < vBins.size()-uiUntouchedDimensions; i++)
            {
                uiRet *= vBins[i].size();
                auto itI = std::upper_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]);
                if(bSubstrOne && itI != vBins[i].begin())
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
            {
                if(xPoint[i] < vBins[i].front() || xPoint[i] > vBins[i].back())
                    return true;
                if(*std::lower_bound(vBins[i].begin(), vBins[i].end(), xPoint[i]) != xPoint[i])
                {
                    //std::cerr << pointCoords(xPoint) << " not on spot in dimension " << i << std::endl;
                    return false;
                }
            }
            return true;
        }

        std::array<std::vector<int64_t>, N> makeBins(size_t uiPointIdx)
        {
            std::array<std::vector<int64_t>, N> vRet;
            for(size_t i = 0; i < vBins.size(); i++)
            {
                size_t ii = vBins.size() - i - 1;
                std::vector<int64_t> xP;
                if(ii < vBins.size() - uiUntouchedDimensions)
                {
                    size_t j = uiPointIdx % vBins[ii].size();
                    size_t uiAdd = (vBins[i][j+1] - vBins[i][j]) / 500;
                    if(uiAdd < 1)
                        uiAdd = 1;
                    for(int64_t uiI = vBins[i][j]; uiI <= vBins[i][j+1]; uiI+=uiAdd)
                        xP.push_back(uiI);
                }
                else
                    for(int64_t uiI : vBins[ii])
                        xP.push_back(uiI);
                //std::cerr << "makeBins " << vecCoords(xP) << std::endl;
                vRet[ii] = xP;
                if(ii < vBins.size() - uiUntouchedDimensions)
                    uiPointIdx /= vBins[ii].size();
            };
            return vRet;
        }

        size_t maxIdx()
        {
            size_t uiRet = 1;
            for(size_t i = 0; i < vBins.size(); i++)
                    uiRet *= vBins[i].size();
            return uiRet;
        }

        size_t maxIdxPoints()
        {
            size_t uiRet = 1;
            for(size_t i = 0; i < vBins.size()-uiUntouchedDimensions; i++)
                    uiRet *= vBins[i].size();
            return uiRet;
        }

        std::string pointCoords(const Point& vP)
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

        std::string vecCoords(const std::vector<int64_t>& vP)
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

        void integrate(size_t uiTid, size_t iFixedDim, size_t iDim, Point& vPoint, std::vector<size_t>& vCnt, 
                      std::mutex& xLock, size_t& uiCnt){
            if(iDim == iFixedDim)
                iDim += 1;
            if(iDim >= vBins.size())
            {
                size_t uiIntegral = 0;
                for(size_t uiI : vBins[iFixedDim])
                {
                    vPoint[iFixedDim] = uiI;
                    size_t uiIdx = getIdx(vPoint, false);
                    DEBUG(std::cout << "integrate " << pointCoords(vPoint) << " (" << uiIdx << "/" << vData.size()
                          << ")" << ": " << 
                          vData[uiIdx] << " -> " << uiIntegral << std::endl<< std::endl;)
                    uiIntegral += vData[uiIdx];
                    vData[uiIdx] = uiIntegral;
                    vCnt[uiTid]++;
                    if(vCnt[uiTid] >= 1000000 )
                    {
                        std::unique_lock<std::mutex> _(xLock);
                        uiCnt += vCnt[uiTid];
                        vCnt[uiTid] = 0;
                        if(uiTid == 0)
                            std::cerr << "\rintegrating " << uiCnt << " / " << maxIdx()*vBins.size() << " = " << 
                                        (100*uiCnt)/(maxIdx()*vBins.size()) << "%\033[K";
                    }
                }
                DEBUG(std::cout << std::endl;)
            }
            else
            {
                std::vector<size_t> vCnt(uiThreads+1);
                {
                    ThreadPool xPool((iDim == 0 || (iFixedDim == 0 && iDim == 1)) ? uiThreads : 0);
                    for(size_t uiI : vBins[iDim])
                        xPool.enqueue([&](size_t uiTid, size_t uiI){
                            Point vP(vPoint);
                            vP[iDim] = uiI;
                            integrate(uiTid, iFixedDim, iDim + 1, vP, vCnt, xLock, uiCnt);
                        }, uiI);
                }
            }
        }

        void integrate()
        {
            std::cerr << "\r\033[Kintegrating";
            Point vPoint;
            std::vector<size_t> vCnt(1);
            std::mutex xLock;
            size_t uiCnt = 0;
            for(size_t iFixedDim = 0; iFixedDim < vBins.size(); iFixedDim++)
                integrate(0, iFixedDim, 0, vPoint, vCnt, xLock, uiCnt);
            std::cerr << "\r\033[K";
        }

        void makeCache(size_t uiI, Point& vP)
        {
            if(vDataCache[uiI] != nullptr)
                return;
            if(vDataPoints[uiI].size() < uiEleNoCache)
                return;
            {
                std::unique_lock<std::mutex> _(xCacheLock);
                while(vpCache[uiLayer].size() > uiSizeCache)
                {
                    std::cerr << "cache too full" << std::endl;
                    vpCache[uiLayer].front().reset();
                    vpCache[uiLayer].pop_front();
                }
            }
            std::cerr << "computing cache " << uiLayer << ":" << vpCache[uiLayer].size() << " due to " 
                      << pointCoords(vP) << " out of " << vDataPoints[uiI].size() << " elements.\033[K\r";
            //for(auto& rPoint : vDataPoints[uiI])
            //    std::cerr << pointCoords(rPoint.first) << std::endl;
            std::string sNewFileName = sFilename + "." + std::to_string(uiI);
            std::ifstream x((sNewFileName + ".count.bin").c_str());
            bool b = x.good();
            x.close();
            //std::cerr << vDataPoints.size() << " " << uiI << std::endl;
            if(b)
                vDataCache[uiI] = std::make_shared<CacheVector>(sNewFileName, uiSizeCache, uiThreads, 
                                                                uiUntouchedDimensions, uiLayer+1, uiEleNoCache);
            else
                vDataCache[uiI] = std::make_shared<CacheVector>(sNewFileName, 
                                                                vDataPoints[uiI], makeBins(uiI), uiSizeCache, 
                                                                uiThreads, uiUntouchedDimensions, uiLayer+1, 
                                                                uiEleNoCache);  // 
            DEBUG(std::cerr << vDataCache[uiI]->print() << std::endl;)
            
            std::unique_lock<std::mutex> _(xCacheLock);
            vpCache[uiLayer].push_back(vDataCache[uiI]);
        }

        int64_t pointVal(Point& vP)
        {
            auto xIdx = getIdx(vP, true);
            if(xIdx == std::numeric_limits<size_t>::max())
                return 0;
            int64_t iRet = 0;
            return (int64_t)vData[xIdx];
        }

        int64_t count_help(Point& vP, Point vLeft, Point vRight, 
                           size_t uiI, size_t uiEdgesFromRight)
        {
            if(uiI >= vBins.size())
            {
                int64_t x = 0;
                if(uiEdgesFromRight % 2 == 0)
                    x = pointVal(vP);
                else
                    x = -pointVal(vP);
                DEBUG(std::cerr << "count " << pointCoords(vP) << ": " << x << std::endl;)
                return x;
            }
            int64_t iRet = 0;
            vP[uiI] = vLeft[uiI];
            iRet += count_help(vP, vLeft, vRight, uiI+1, uiEdgesFromRight+1);
            vP[uiI] = vRight[uiI];
            iRet += count_help(vP, vLeft, vRight, uiI+1, uiEdgesFromRight);
            return iRet;
        }

        int64_t count_cache(Point& vP, Point vLeft, Point vRight, size_t i)
        {
            if(i >= vBins.size()-uiUntouchedDimensions)
            {
                size_t xIdx = getIdxPoints(vP, false);
                DEBUG(std::cerr << "pointVal->cache " << xIdx << std::endl;)
                makeCache(xIdx, vP);
                if(vDataCache[xIdx] == nullptr)
                {
                    int64_t iRet = 0;
                    for(DataPoint p: vDataPoints[xIdx])
                    {
                        bool bAdd = true;
                        for(size_t i = 0; i < vBins.size() && bAdd; ++i)
                            if(vLeft[i] > p.first[i] || vRight[i] <= p.first[i] )
                                bAdd = false;
                        if(bAdd)
                            ++iRet;
                    }
                    return iRet;
                }
                else
                    return vDataCache[xIdx]->count(vLeft, vRight);
            }
            else
            {
                int64_t uiRet = 0;
                auto xStart = std::upper_bound(vBins[i].begin(), vBins[i].end(), vLeft[i]);
                if(xStart != vBins[i].begin())
                    xStart--;
                for(auto itI = xStart;
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
        void save(size_t uiMinSizeCacheSave)
        {
            save_helper(sFilename, uiMinSizeCacheSave);
        }

        size_t count(Point vLeft, Point vRight)
        {
            if(onSpot(vLeft) && onSpot(vRight))
            {
                Point vP;
                return (size_t)count_help(vP, vLeft, vRight, 0, 0);
            }
            else 
            {
                DEBUG(std::cerr << "cache " << std::endl;)
                Point vP;
                return count_cache(vP, vLeft, vRight, 0);
            }
        }

        std::string print_help(Point& vP, size_t uiDim)
        {
            std::string sRet = "";
            if(uiDim >= vBins.size())
                sRet = pointCoords(vP) + ": " + 
                        std::to_string(vDataPoints[getIdxPoints(vP, false)].size()) + " -> " +
                        std::to_string(vData[getIdx(vP, false)]) + "\n";
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
            Point vP;
            std::string s = "";
            for(size_t uiI = 0; uiI < vDataCache.size(); uiI++)
            {
                if(vDataCache[uiI] != nullptr)
                    s += vDataCache[uiI]->print();
                else
                    s += "nullptr\n";
            }
            return print_help(vP, 0) + s;
        }
};


template <size_t N>
CacheVector<N>::CacheVector(std::string sFileName, Points vvfPoints, std::array<std::vector<int64_t>, N> vBins, 
                        size_t uiSizeCache, size_t uiThreads, size_t uiUntouchedDimensions, size_t uiLayer,
                        size_t uiEleNoCache)
        :
    sFilename(sFileName),
    vBins(vBins),
    uiUntouchedDimensions(uiUntouchedDimensions),
    vData(maxIdx(), 0),
    vDataPoints(maxIdxPoints()),
    vDataCache(maxIdxPoints(), nullptr),
    uiSizeCache(uiSizeCache),
    uiEleNoCache(uiEleNoCache),
    uiThreads(uiThreads),
    uiLayer(uiLayer)
{
    {
        std::unique_lock<std::mutex> _(xCacheLock);
        while(vpCache.size() <= uiLayer)
            vpCache.push_back(std::deque<std::shared_ptr<CacheVector>>());
    }
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
                    
                    DEBUG(std::cout << "loading " << pointCoords(p.first) << " " << getIdxPoints(p.first, false) 
                                    << " " << vDataPoints.size() << std::endl;)
                    size_t uiIdx = getIdx(p.first, false);
                    size_t uiIdx2 = getIdxPoints(p.first, false);
                    if(uiIdx2 == std::numeric_limits<size_t>::max())
                    {
                        std::cerr << "point does not fit in vector " << pointCoords(p.first) << std::endl;
                        for(auto& rBin : vBins)
                            std::cerr << vecCoords(rBin) << std::endl;
                        continue;
                    }

                    std::unique_lock<std::mutex> xLock(vArr[uiIdx % vArr.size()]);
                    vDataPoints[uiIdx2].push_back(p);
                    ++vData[uiIdx];
                }
            }, uiI);
    }
    std::cerr << "\r\033[K";
    integrate();
}


template <size_t N>
CacheVector<N>::CacheVector(std::string sFileName, size_t uiSizeCache, size_t uiThreads, 
                            size_t uiUntouchedDimensions, size_t uiLayer, size_t uiEleNoCache)
        :
    sFilename(sFileName),
    vBins(),
    uiUntouchedDimensions(uiUntouchedDimensions),
    vData(),
    vDataPoints(),
    vDataCache(),
    uiSizeCache(uiSizeCache),
    uiEleNoCache(uiEleNoCache),
    uiThreads(uiThreads),
    uiLayer(uiLayer)
{
    {
        std::unique_lock<std::mutex> _(xCacheLock);
        while(vpCache.size() <= uiLayer)
            vpCache.push_back(std::deque<std::shared_ptr<CacheVector>>());
    }
    serialize(sFileName, false);
    assert(vDataPoints.size() == maxIdxPoints());
    vDataCache.resize(maxIdxPoints(), nullptr);
}

template <size_t N>
typename std::vector<std::deque<std::shared_ptr<CacheVector<N>>>> CacheVector<N>::vpCache;
template <size_t N>
typename std::mutex CacheVector<N>::xCacheLock;

}

template<size_t N>
void export_(pybind11::module& m, std::string sName){
    pybind11::class_<kdtree::CacheVector<N>>( m, sName.c_str() )
        .def( pybind11::init<std::string, typename kdtree::CacheVector<N>::Points, 
                             std::array<std::vector<int64_t>, N>, 
                             size_t, size_t, size_t, size_t, size_t>( ) ) // constructor
        .def( pybind11::init<std::string, size_t, size_t, size_t, size_t, size_t>( ) ) // constructor
        .def( "count", &kdtree::CacheVector<N>::count )
        .def( "save", &kdtree::CacheVector<N>::save )
        .def( "__str__", &kdtree::CacheVector<N>::print );
}

