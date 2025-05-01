#ifndef __TYPES_H__
#define __TYPES_H__


#include <cstddef>
#define MAX_CLUSTER_MEMBERS 300
#define MAX_SEARCH_QUEUE_SIZE 300
#define MAX_CLUSTERS 300
#define MAX_NUM_NEIGHBORS 1000
enum class PointLabel
{
    UNVISITED,
    NOISE,
    CLUSTERED
};

struct Point
{
    double x, y;
    int cluster_id = -1;
    PointLabel label = PointLabel::UNVISITED;
};

struct cluster
{
    std::size_t size = 0;
    Point members [MAX_CLUSTER_MEMBERS];
    float centroid_x = 0;
    float centroid_y = 0;
    int cluster_id = -1;
    int bbox_min_x = -1;
    int bbox_max_x = -1;
    int bbox_min_y = -1;
    int bbox_max_y = -1;
};

template<typename T, size_t Capacity>
class StaticQueue {
public:
    StaticQueue() : head_(0), tail_(0), size_(0) {}

    bool push(const T& value) {
        if (size_ >= Capacity) return false; // Queue full
        data_[tail_] = value;
        tail_ = (tail_ + 1) % Capacity;
        ++size_;
        return true;
    }

    bool pop() {
        if (size_ == 0) return false; // Queue empty
        head_ = (head_ + 1) % Capacity;
        --size_;
        return true;
    }

    T& front() {
        return data_[head_];
    }

    const T& front() const {
        return data_[head_];
    }

    bool empty() const {
        return size_ == 0;
    }

    bool full() const {
        return size_ == Capacity;
    }

    size_t size() const {
        return size_;
    }

private:
    T data_[Capacity];
    size_t head_;
    size_t tail_;
    size_t size_;
};
#endif // __TYPES_H__