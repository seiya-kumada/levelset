//
//  parallel_tools.h
//  LevelSetMethod2
//
//  Created by kumada on 2012/12/19.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#ifndef LevelSetMethod2_parallel_tools_h
#define LevelSetMethod2_parallel_tools_h

#include <thread>
#include <future>
#include <iterator>

namespace lsm
{
        class join_threads
        {
                std::vector<std::thread>& threads;
        public:
                explicit join_threads(std::vector<std::thread>& threads)
                        : threads(threads) {}
                
                ~join_threads()
                {
                        std::for_each(std::begin(threads), std::end(threads),
                                [](std::thread& thread){
                                        if ( thread.joinable() ) {
                                                thread.join();
                                        } else {
                                                std::exit(0);
                                        }
                                }
                        );
                }
        }; // class join_threads

        template<typename Iterator, typename Func>
        inline void parallel_for_each(unsigned long num_threads, Iterator first, Iterator last, Func fun)
        {
                unsigned long const length = std::distance(first, last);
                if ( !length ) {
                        return;
                }
                unsigned long const block_size = length / num_threads;
                std::vector<std::thread> threads {num_threads - 1};
                join_threads joiner {threads};
                Iterator block_start = first;
                for ( unsigned long i = 0; i < num_threads - 1; ++i ) {
                        Iterator block_end = block_start;
                        std::advance(block_end, block_size);
                        std::packaged_task<void(void)> task { [=]{ std::for_each(block_start, block_end, fun); } };
                        threads[i] = std::thread { std::move(task) };
                        block_start = block_end;
                }
                std::for_each(block_start, last, fun);
        }

        template<typename Range, typename Func>
        inline void parallel_for_each(unsigned long num_threads, const Range& range, Func f)
        {
             parallel_for_each(num_threads, std::begin(range), std::end(range), f);
        }
} // namespace lsm

#endif // LevelSetMethod2_parallel_tools_h
