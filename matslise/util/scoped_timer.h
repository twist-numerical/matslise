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
        long long insertion_counter = 0;

        Profiler() = default;

        Profiler(Profiler &other) = delete;

        Profiler &operator=(Profiler &other) = delete;

        struct ProfileEntry {
            long long count;
            long long insertion;
            std::chrono::nanoseconds time;
        };

    public:
        static Profiler *profiler;

        std::unordered_map<std::string, ProfileEntry> executionTimes;
        std::vector<std::pair<std::string, long long>> stack;

        static Profiler &instance() {
            if (!profiler) {
                profiler = new Profiler();
                std::atexit([]() {
                    Profiler::profiler->report();
                    delete Profiler::profiler;
                });
            }
            return *profiler;
        }

        void push(std::string name) {
            stack.emplace_back(std::move(name), insertion_counter++);
        }

        void pop() {
            stack.pop_back();
        }

        std::string getIdentifier() {
            if (stack.empty())
                return ">";

            std::string id;
            for (auto &call: stack)
                id += " > " + call.first;
            return id;
        }

        void addTime(std::chrono::nanoseconds time) {
            std::string id = getIdentifier();

            if (executionTimes.find(id) == executionTimes.end())
                executionTimes[id] = ProfileEntry{.count=1, .insertion=stack.back().second, .time= time};
            else {
                auto &p = executionTimes[id];
                p.count += 1;
                p.time += time;
            }
        }

        void report() {
            std::cerr << "\n\n**Matslise profiler report:**" << std::endl;
            std::vector<std::pair<std::string, ProfileEntry>> entries(
                    executionTimes.begin(), executionTimes.end());
            std::sort(entries.begin(), entries.end(),
                      [](const auto &a, const auto &b) { return a.second.insertion < b.second.insertion; });
            for (auto &f: entries) {
                const auto &entry = f.second;
                std::cerr << f.first << "\t" << entry.count << "\t"
                          << std::fixed << std::setprecision(2)
                          << double(entry.time.count() / 10000) / 100. << "ms" << std::endl;
            }
        }
    };

    inline Profiler *Profiler::profiler = nullptr;

    class ScopedTimer {
    public:
        std::chrono::time_point<std::chrono::high_resolution_clock> start;

        explicit ScopedTimer(const std::string &name) {
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

#define MATSLISE_SCOPED_TIMER(name) matslise::ScopedTimer __st(name)
#else
#define MATSLISE_SCOPED_TIMER(name)
#endif //MATSLISE_PROFILING

#endif //MATSLISE_SCOPED_TIMER_H
