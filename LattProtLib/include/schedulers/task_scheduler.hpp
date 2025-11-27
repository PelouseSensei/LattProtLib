#ifndef TASKSCHEDULER_H
#define TASKSCHEDULER_H

#include <vector>
#include <functional>
#include <stdexcept>

template <typename ConfigType>
class TaskScheduler {
public:
    using Task = std::function<void(ConfigType&, int)>;

    void addTask(Task t) {
        tasks.push_back(std::move(t));
    }

    void run(ConfigType& config, int step) {
        for (auto& t : tasks) {
            t(config, step);
        }
    }

    void run_nSteps(ConfigType& config, int nSteps) {
        for (int i = 0; i < nSteps; i++) {
            run(config, i);
        }
    }

    // ----- New methods -----
    void clearTasks() {
        tasks.clear();
    }

    void clearTask(size_t index) {
        if (index >= tasks.size()) {
            throw std::out_of_range("TaskScheduler::clearTask: index out of range");
        }
        tasks.erase(tasks.begin() + index);
    }

private:
    std::vector<Task> tasks;
};

#endif // TASKSCHEDULER_H
