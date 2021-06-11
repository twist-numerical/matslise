//
// Created by toon on 8/13/20.
//

#ifndef MATSLISE_SCOPED_TIMER_H
#define MATSLISE_SCOPED_TIMER_H

#ifdef MATSLISE_PROFILING
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <vector>

namespace matslise {
    class Profiler {
    private:
        Profiler() = default;
        Profiler(Profiler& other) = delete;
        Profiler& operator=(Profiler& other) = delete;
    public:
        static Profiler* profiler;

        std::map<std::string, std::pair<long long, std::chrono::nanoseconds>> executionTimes;
        std::vector<std::string> stack;

        static Profiler& instance() {
            if(!profiler) {
                profiler = new Profiler();
                std::atexit([]() {
                    Profiler::profiler->report();
                    delete Profiler::profiler;
                });
            }
            return *profiler;
        }

        void push(std::string name) {
            stack.push_back(name);
        }

        void pop() {
            stack.pop_back();
        }

        std::string getIdentifier() {
            if (stack.empty())
                return ">";

            std::string id = "";
            for (auto &call : stack)
                id += " > " + call;
            return id;
        }

        void addTime(std::chrono::nanoseconds time) {
            std::string id = getIdentifier();
            if (executionTimes.find(id) == executionTimes.end())
                executionTimes[id] = {1, time};
            else {
                auto &p = executionTimes[id];
                p.first += 1;
                p.second += time;
            }
        }

        void report() {
            long long count;
            std::chrono::nanoseconds time;
            std::cerr << "\n\n**Matslise profiler report:**" << std::endl;
            for (auto &f : executionTimes) {
                tie(count, time) = f.second;
                std::cerr << f.first << "\t" << count << "\t"
                          << std::fixed << std::setprecision(2)
                          << time.count() / 10000 / 100. << "ms" << std::endl;
            }
        }
    };

    inline Profiler* Profiler::profiler = nullptr;

    class ScopedTimer {
    public:
        std::chrono::time_point<std::chrono::high_resolution_clock> start;

        ScopedTimer(const std::string &name) {
            Profiler::instance().push(name);
            start = std::chrono::high_resolution_clock::now();
        }

        ~ScopedTimer() {
            auto end = std::chrono::high_resolution_clock::now();
            Profiler::instance().addTime(end - start);
            Profiler::instance().pop();
        }
    };
}

#define MATSLISE_SCOPED_TIMER(name) ScopedTimer __st(name)
#else
#define MATSLISE_SCOPED_TIMER(name)
#endif //MATSLISE_PROFILING

#endif //MATSLISE_SCOPED_TIMER_H
