#pragma once

#include <algorithm>
#include <condition_variable>
#include <deque>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

namespace sps
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
    ThreadPool( );

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
        if( threads <= 1 )
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
            tasks.push( [ task ]( size_t task_id ) { ( *task )( task_id ); } );
        } // end of scope for the lock, so we will release the lock.

        /* Inform some waiting consumer (worker )that we have a fresh task.
         */
        condition.notify_one( );

        return xFuture;
    } // method enqueue

    size_t numThreads( ) const
    {
        return workers.size( );
    }
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
    if( threads <= 1 )
        return;
    for( size_t i = 0; i < threads; ++i )
        workers.emplace_back( [ this, i ] {
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
                catch( std::exception& e )
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

/* Constructor just launches some amount of workers
 */
inline ThreadPool::ThreadPool( ) : ThreadPool( std::thread::hardware_concurrency( ) )
{}

/* the destructor joins all threads
 */
inline ThreadPool::~ThreadPool( )
{
    // if the threadpool is set to have 0 threads, we will execute every task in the main thread
    // immediately, thus we dont need to setup threads at all.
    // this can be used to disable multithreading simply setting the pool to have 0 threads
    if( threads <= 1 )
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

} // namespace sps