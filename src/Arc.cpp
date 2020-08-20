//
// Created by carlos on 06/03/19.
//

#include "../headers/Arc.h"

Arc::Arc(int o, int d, int delay, int jitter, int bandwidth, int estimateLinkDuration) {
    this->o = o;
    this->d = d;
    this->delay = delay;
    this->jitter = jitter;
    this->bandwidth = bandwidth;
    this->estimateLinkDuration = estimateLinkDuration;
}

int Arc::getO() const {
    return o;
}

int Arc::getD() const {
    return d;
}

int Arc::getDelay() const {
    return delay;
}

int Arc::getJitter() const {
    return jitter;
}

int Arc::getBandwidth() const {
    return bandwidth;
}

int Arc::getEstimateLinkDuration() const {
    return estimateLinkDuration;
}
